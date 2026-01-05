import pandas as pd
import numpy as np
from math import sqrt

schimm_F8_deltas = [-231.2, -231.2, -166.8, -211, -206.2, -214.2, -166.7, -195.5]
schimm_F8_peaks = [x for x in range(7, 15)]
schimm_F8_lookup = dict(zip(schimm_F8_peaks, schimm_F8_deltas))


def process_2H_vendor_data(exported_data: pd.DataFrame) -> pd.DataFrame:
    """
    Standardise the columns of a dataframe obtained from processing Delta V hydrogen isotope ratio
    `.dxf` files. For information on how to produce this dataframe, see the repository README.

    :param exported_data: A dataframe containing data as exported from hydrogen isotope isodat files
    :return: The dataframe with added/renamed columns
    """

    df = exported_data.rename({"Nr.": "peak_no"}, axis=1)
    df["file_datetime"] = pd.to_datetime(df["file_datetime"])
    df["run_id"] = df["Identifier 1"]
    df["run_id_no"] = df["run_id"].str.extract(r"(\d+)").astype(int)
    df["sample_id"] = df["Identifier 2"]
    df["delta_meas"] = df["d 2H/1H"]

    return df


def get_bracketing_reference_runs(sample_summary_data: pd.DataFrame,
                                  reference_summary_data: pd.DataFrame) -> pd.DataFrame:
    """
    Take a dataframe of sample summary data (one row per sample run containing columns `run_id`,
    `run_id_no`, and `file_datetime`) and a dataframe of reference material summary data (from
    running `reference_scale_normalisation` on the reference material peak data) and add information about
    the bracketing pair of reference runs for each sample run to the sample dataframe.

    :param sample_summary_data: Summary dataframe for samples
    :param reference_summary_data: Summary dataframe for reference materials (including regression metrics)

    :return: Sample summary dataframe with added columns for preceding and following reference material runs
    """

    df_sam_sum = sample_summary_data.copy()
    df_ref_sum = reference_summary_data.copy()

    # Check that the data has been properly processed
    required_cols = ["run_id_no", "file_datetime"]
    for df in (df_sam_sum, df_ref_sum):
        if any(col not in df.columns for col in required_cols):
            raise ValueError("Required columns not found. Please apply a processing function.")

    # For each run, find the immediately preceding and following reference material runs
    for direction in ("backward", "forward"):

        # Rename the reference data columns to reflect how the data relate to the samples
        ref_data = df_ref_sum.copy().rename(
            {x: f"{direction}_ref_{x}" for x in reference_summary_data.columns}, axis=1
            )

        # Join the reference data to the sample data on run ID number in the specified direction
        df_sam_sum = pd.merge_asof(df_sam_sum,
                                   ref_data,
                                   left_on="run_id_no",
                                   right_on=f"{direction}_ref_run_id_no",
                                   direction=direction)

        # Calculate the time delta in seconds between each sample run and each of its bracketing reference runs
        df_sam_sum[f"{direction}_ref_time_delta"] = (
            pd.to_datetime(df_sam_sum[f"{direction}_ref_file_datetime"]) -
            pd.to_datetime(df_sam_sum["file_datetime"])
            ) / np.timedelta64(1, "s")

    return df_sam_sum


def regression_analysis(group: pd.DataFrame) -> pd.Series:
    """
    An aggregation function to convert reference material peak data into linear regression metrics.

    This function is not normally called explicitly; it is most usefully used in a split-apply-combine
    operation, as in `reference_scale_normalisation`.

    :param group: Dataframe containing data for a single reference material run
    :return: A Pandas series with named columns `gradient`, `intercept`, and `rms_error`
    """

    meas_vals = group["delta_meas"]
    known_vals = group["delta_known"]
    grad, incpt = np.polyfit(known_vals, meas_vals, 1)  # Do 1st order linear regression
    norm_meas_vals = [(d - incpt) / grad for d in meas_vals]  # Normalise the measured data with the regression metrics

    delta_diffs = np.subtract(known_vals, norm_meas_vals)
    rms_error = sqrt(sum(diff**2 for diff in delta_diffs) / len(known_vals))  # Calculate mean RMS error

    return pd.Series({"gradient": grad,
                      "intercept": incpt,
                      "rms_error": rms_error})


def reference_scale_normalisation(ref_meas: pd.DataFrame,
                                  ref_info: dict = schimm_F8_lookup
                                  ) -> pd.DataFrame:
    """
    Calculate scale normalisation parameters (gradient, intercept, RMS error) for each reference material run
    in a dataframe. The parameters are of the form

        deltaMeas = *m* · deltaTrue + *c*

    Where deltaMeas is the measured delta value of a given peak, deltaTrue is the "true" (scale corrected)
    value and *m* and *c* are the gradient and intercept of the least squares regression line respectively.

    :param ref_meas: A long-format dataframe of reference material peak data
    :param ref_info: A dictionary of peak_number: known_delta_value pairs for the reference material
    :return: A per-run dataframe containing IDs and scale normalisation parameters for each reference material
    injection
    """

    # Check that the data has been processed
    if "peak_no" not in ref_meas.columns:
        raise ValueError("Required column(s) not found. Please apply a processing function.")

    # Select only the reference peaks, add a column with the known peak values
    df_ref = ref_meas\
        .loc[ref_meas["peak_no"].isin(ref_info.keys())]\
        .assign(delta_known=lambda d: d["peak_no"].map(ref_info))\
        .copy()

    # Perform regression analysis on each run
    group_cols = ["run_id", "run_id_no", "file_datetime"]
    df_ref_reg = df_ref.groupby(group_cols, as_index=False)\
                       .apply(regression_analysis, include_groups=False)

    return df_ref_reg


def normalise_sample_peak_data(sample_peak_data: pd.DataFrame,
                               sample_summary_data: pd.DataFrame,
                               peak_ret_ts: list,
                               ret_t_tolerance: int = 5):
    """
    Select peaks of interest and normalise them according to the correct bracketing standard runs.

    :param sample_peak_data: Long-format multi-run peak data in the format of "vendor data" produced
    by the `isoreader` R package
    :param sample_summary_data: Long-format multi-run sample summary data as produced by
    `get_bracketing_reference_runs`
    :param peak_ret_ts: Expected retention times of peaks of interest (e.g. C16, C18) in seconds
    :param ret_t_tolerance: Tolerance for variation in peak retention times in seconds
    :return: A dataframe of selected peaks, with columns for their raw and normalised delta values
    """

    # Check that both dataframes have been processed correctly
    for df in (sample_peak_data, sample_summary_data):
        if "run_id" not in df.columns:
            raise ValueError("Required column(s) not found. Please apply a processing function.")

    # Find all shared columns in the two dataframes so that we can join on them
    join_cols = list(set(sample_peak_data.columns) & set(sample_summary_data.columns))
    merged_data = pd.merge(sample_peak_data,
                           sample_summary_data,
                           on=join_cols)

    # Select only rows with retention times within margin of the selected RTs
    merged_data = merged_data.loc[(
        np.abs(
            merged_data["Rt"].to_numpy()[:, None] - np.array(peak_ret_ts)
            ) <= ret_t_tolerance).any(axis=1)
    ]

    # Calculate the mean correction coefficients from the bracketing reference material runs
    for coeff in ("gradient", "intercept"):
        merged_data[f"mean_{coeff}"] = merged_data[[f"forward_ref_{coeff}",
                                                    f"backward_ref_{coeff}"]].mean(axis=1)

    # Calculate the normalised delta for each peak of interest
    merged_data["delta_norm"] = (merged_data["delta_meas"] -
                                 merged_data["mean_intercept"]) / merged_data["mean_gradient"]

    return merged_data


def triplicate_analysis(normalised_data: pd.DataFrame,
                        peak_ret_ts: list,
                        ret_t_tolerance: int = 5) -> pd.DataFrame:
    """
    Take long-format normalised peak data and extract the mean and StdDev for each set of triplicates

    :param normalised_data: Output from `normalise_sample_peak_data`
    :param peak_ret_ts: Expected retention times of peaks of interest (e.g. C16, C18) in seconds
    :param ret_t_tolerance: Tolerance for variation in peak retention times in seconds
    :return: Dataframe containing summary delta value data for each sample
    """

    df = normalised_data.copy()
    df["sample_id_short"] = df["sample_id"].str[:-1]  # Remove repeat suffix to allow grouping by sample
    df["approx_RT"] = df["Rt"].apply(  # Round retention times to expected values if within tolerance
        lambda x: int([t for t in peak_ret_ts if abs(x - t) <= ret_t_tolerance][0])
    )

    index_cols = ["sample_id_short",
                  "approx_RT",
                  "forward_ref_run_id",
                  "backward_ref_run_id",
                  "mean_intercept",
                  "mean_gradient"
                  ]

    # Group runs by triplicate grouping and compound (retention time); calculate mean and StdDev, record run count
    df_grouped = df.groupby(index_cols)["delta_norm"]\
                   .agg(["count", "mean", "std"]).reset_index()\
                   .rename(columns={x: "delta_norm_" + x for x in ["mean", "std"]})

    return df_grouped
