# IRMS data processing

## Overview

These are Python tools developed by the Organic Geochemistry Unit to help
process isotope ratio mass spectrometry (IRMS) data. Currently, this
process has only been tested on IsoDat `.dxf` files from the Delta V 1,
though we hope to expand the utility further.

The IRMS data must first be prepared into `.csv` files using the
[isoreader](https://isoreader.isoverse.org/index.html) R package. No similar package for
loading IRMS binary files
yet exists in Python.

## Installation

This package can be installed into a Python environment with `pip`, using
the following command:

```bash
pip install git+https://github.com/nicksgoodusername/ogu-data-processing.git
```

To submit modifications/additions, please branch and create a merge request. For help
with this, contact Nick Hall.

## Pre-processing

Currently, this package cannot directly handle output files from IRMS software
(e.g. `.dxf` files).
Handily, the [isoreader](https://isoreader.isoverse.org/index.html) R package has been
developed to handle these files, and can export useful data as easy-to-read `.csv`
files, which this package can then process.

For a directory (folder) containing `.dxf` files, we can create single
`.csv` file containing every run's data as follows:

```R
library(isoreader)

DIR <- "path/to/directory/data"
OUTFILE <- "path/to/directory/all_run_data.csv"

OUTFILE <- paste(DIR, OUTFILE, sep = "")

dxf_files <- list.files(path=DIR, pattern="*.dxf", full.names=TRUE, recursive=FALSE)
file_data <- iso_read_continuous_flow(dxf_files)
output_data <- file_data |> iso_get_vendor_data_table(
    include_file_info=c(file_datetime, starts_with("Identifier"))
    )

write.table(output_data, OUTFILE, sep=",", row.names=FALSE)
```

Where the `.csv` contains the "vendor data" (the per-peak isotope ratio data that you'd
ordinarily copy and paste into an Excel sheet) for each run as well as some extra
information like the date/time the run was started.

This `.csv` is now ready to be processed in Python with this package.
