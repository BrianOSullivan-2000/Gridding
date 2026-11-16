
## Using MLE Kriging, I've produced monthly grids dating
## all the way from 1865 to 2023. This script just has
## some rough work pertaining to that (data formatting
## and very small checks/diagnostics)

# %%

library(dplyr)

## First read in our rainfall data
rainfall_data <- read.csv(
    "Data/Monthly_Rainfall/monthly_rain_ROI_NI_1878-2023.csv"
)
## Remove any missing values
rainfall_data <- rainfall_data[!is.na(rainfall_data$rain), ]

## From this we're getting geodata and LTAs
load("Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")

# %%

## Load in all the grids
file_names <- list.files(
    "Results/Monthly_Rainfall/Historic_Grids",
    full.names = TRUE
)

## Combine them together
file_data <- sapply(file_names[2:4], read.csv)
coords <- file_data[[1]][c("east", "north")]

file_data <- lapply(
    file_data, {
        function (df) df |> dplyr::select(-east, -north)
    }
)

## Save each year individually

## Check for number of zero-valued grid points per year-month

## Plot annual trend

## Plot monthly/seasonal trends