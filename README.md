# Processing Fluorcam data into Tcrit, T50 and Tmax

Based on processing steps outlined in [Arnold et al. 2021](https://doi.org/10.1071/FP20344) and refined by Madeline Moran, M.S., this repository provides a workflow for automating detection of thermal limits in high-throughput chlorophyll imaging fluorescence.

## How to use this repository

1.  Clone or download to the desired location on your local computer

2.  Click on `fluorcam-processing.Rproj`, which will open the RStudio project in a new window

3.  Run `renv::restore()` in the console to load package library (see [here](https://rstudio.github.io/renv/articles/collaborating.html) for details about collaborating with 'renv')

4.  From the file pane, open `workflow.R`

5.  Run the first 25 lines of the script, which will load the requisite packages and create the data folders

6.  Move a Fluorcam file (.txt) into `data_raw/`; rename if needed. I recommend keeping the YYYYMMDD as the first part of the file name so files will be chronologically ordered.

7.  Run the remainder of `workflow.R`. When complete, `data_processed/` should contain one folder and one .csv file of the same name as the `data_raw/` input:

    -   The folder will contain a .png image for each well location (A1,A2, ... H4,) plotting the raw and fitted data with thermal limits (Tcrit, T50, Tmax).
    -   The file will contain the thermal limits for each well location in tabular form

### Development notes

This processing workflow is in development. Current features include:

-   detecting outliers in the raw fluorescence data (removes points exceeding the median + 3\*IQR)
-   rescaling fluorescence between minimum and maximum, but only designating minimum if occurred prior to maximum
-   calculating T50 only if it occurred prior to Tmax

Desired features still include:

-   adding conditional to breakpoint for loop if model does not converge

-   joining with dataframe of sample identifiers for each well (not yet supplied)

-   adding `date_run` and `date_processed` columns to final output
