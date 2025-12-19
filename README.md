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

For a directory (folder) containing `.dxf` files, we can create a matching
`.csv` file for each as follows:

```R
library(isoreader)

DIR <- "path/to/my/directory"

dxf_files <- list.files(path=DIR, pattern="*.dxf", full.names=TRUE, recursive=FALSE)
lapply(dxf_files, function(x) {
    f <- iso_read_continuous_flow(x)
    t <- iso_get_vendor_data_table(f,
                         include_file_info=c(file_datetime, starts_with("Identifier"))
                         )
    write.table(t, gsub(".dxf", ".csv", x), sep=",", row.names=FALSE)
})
```

Where each `.csv` contains the "vendor data" (the per-peak isotope ratio data that you'd
ordinarily copy and paste into an Excel sheet) as well as some extra information like
the date/time the run was started.

These `.csv`s are now ready to be processed in Python with this package.
