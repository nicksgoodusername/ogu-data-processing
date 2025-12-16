import pandas as pd
import numpy as np
from math import sqrt
from re import search


schimm_F8_deltas = [-231.2, -231.2, -166.8, -211, -206.2, -214.2, -166.7, -195.5]
schimm_F8_peaks = [x for x in range(7, 15)]


def get_summary_info(exported_data: pd.DataFrame):

    run_id = exported_data["Identifier 1"].tolist()[0]
    run_id_no = int(search(r"\d+", run_id)[0])
    sample_id = exported_data["Identifier 2"].tolist()[0]
    run_datetime = exported_data["file_datetime"].tolist()[0]

    return {"run_id": run_id,
            "run_id_no": run_id_no,
            "sample_id": sample_id,
            "run_datetime": run_datetime}


def process_vendor_data(exported_data: pd.DataFrame) -> pd.DataFrame:

    exported_data_mod = exported_data.rename({"Nr.": "peak_no"}, axis=1)
    return exported_data_mod


def get_bracketing_reference_runs(sample_summary_data: pd.DataFrame,
                                  reference_summary_data: pd.DataFrame) -> pd.DataFrame:
    """
    Take a dataframe of sample summary data (from `get_summary_info` over several samples)
    and a dataframe of reference material summary data (from `reference_regression` over several
    standards) and append the bracketing pair of reference runs for each sample run to the sample
    dataframe.

    :param sample_summary_data: Summary data for samples
    :param reference_summary_data: Summary data for reference materials (including regression metrics)
    :return: Sample dataframe with columns for preceding and following reference material runs
    """

    for df in (sample_summary_data, reference_summary_data):
        df["run_id_no"] = pd.to_numeric(df["run_id_no"])

    merged_data = sample_summary_data.copy()
    for direction in ("backward", "forward"):
        ref_data = reference_summary_data.copy().rename(
            {x: f"{direction}_ref_{x}" for x in reference_summary_data.columns}, axis=1)
        merged_data = pd.merge_asof(merged_data,
                                    ref_data,
                                    left_on="run_id_no",
                                    right_on=f"{direction}_ref_run_id_no",
                                    direction=direction)

        merged_data[f"{direction}_ref_time_delta"] = (
            pd.to_datetime(merged_data[f"{direction}_ref_run_datetime"]) -
            pd.to_datetime(merged_data["run_datetime"])
            ) / np.timedelta64(1, "s")

    return merged_data


def reference_regression(exported_data: pd.DataFrame,
                         ref_peaks: list = schimm_F8_peaks,
                         ref_known_deltas: list = schimm_F8_deltas
                         ) -> dict:
    """
    Perform a linear regression on deltaD data from a Schimmelmann F8 reference material
    taken on a Delta V instrument. The resultant data will be used to normalise
    sample deltaD values.

    :param peak_data: A Pandas dataframe containing the peak data as exported by my R function
    :param ref_peaks: A list containing the peak numbers of the reference material peaks in order
    :param ref_known_deltas: A list containing the known deltaD values of each reference peak in order
    :return: A dict containing the run ID and results of the linear regression
    """

    out_dict = get_summary_info(exported_data)

    if not len(ref_known_deltas) == len(ref_peaks):
        raise ValueError('The number of reference peaks must be the same as the number of known deltaD values')

    df = process_vendor_data(exported_data)

    df = df.loc[df["peak_no"].isin(ref_peaks)].copy().sort_values(by=["peak_no"])

    x_vals = ref_known_deltas
    y_vals = df["d 2H/1H"].tolist()

    grad, incpt = np.polyfit(x_vals, y_vals, 1).tolist()

    ref_peak_deltas_norm = [(d - incpt) / grad for d in y_vals]
    ref_delta_diffs = np.subtract(ref_known_deltas, ref_peak_deltas_norm)
    rms_error = sqrt(sum([diff**2 for diff in ref_delta_diffs]) / len(ref_peaks))

    out_dict.update({"gradient": grad,
                     "intercept": incpt,
                     "rms_error": rms_error})

    return out_dict


def sel_norm_peak_data(sample_peak_data: pd.DataFrame,
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

    peak_data = process_vendor_data(sample_peak_data)

    merged_data = pd.merge(peak_data,
                           sample_summary_data,
                           left_on="Identifier 1",
                           right_on="run_id")

    # Select only rows with retention times within margin of the selected RTs
    merged_data = merged_data.loc[(
        np.abs(
            merged_data["Rt"].to_numpy()[:, None] - np.array(peak_ret_ts)
            ) <= ret_t_tolerance).any(axis=1)
    ]

    for coeff in ("gradient", "intercept"):
        merged_data[f"mean_{coeff}"] = merged_data[[f"forward_ref_{coeff}",
                                                    f"backward_ref_{coeff}"]].mean(axis=1)

    merged_data["norm d 2H/1H"] = (merged_data["d 2H/1H"] -
                                   merged_data["mean_intercept"]) / merged_data["mean_gradient"]

    return merged_data


def triplicate_analysis(normalised_data: pd.DataFrame,
                        peak_ret_ts: list,
                        ret_t_tolerance: int = 5) -> pd.DataFrame:
    """
    Take long-format normalised peak data and extract the mean and StdDev for each set of triplicates

    :param normalised_data: Output from `sel_norm_peak_data`
    :param peak_ret_ts: Expected retention times of peaks of interest (e.g. C16, C18) in seconds
    :param ret_t_tolerance: Tolerance for variation in peak retention times in seconds
    :return: Dataframe containing summary delta value data for each sample
    """

    df = normalised_data.copy()
    df["sample_id_short"] = df["sample_id"].str[:-1]
    df["approx_RT"] = df["Rt"].apply(
        lambda x: int([t for t in peak_ret_ts if abs(x - t) <= ret_t_tolerance][0])
    )

    index_cols = ["sample_id_short",
                  "approx_RT",
                  "forward_ref_run_id",
                  "backward_ref_run_id",
                  "mean_intercept",
                  "mean_gradient"
                  ]

    df_grouped = df.groupby(["sample_id_short", "approx_RT"])\
                   .describe()["norm d 2H/1H"][["count", "mean", "std"]].reset_index()
    df_grouped.rename(columns={x: x + "_delta_2H" for x in ["mean", "std"]}, inplace=True)

    df_out = df[index_cols].drop_duplicates().merge(df_grouped, on=["sample_id_short", "approx_RT"])

    return df_out
