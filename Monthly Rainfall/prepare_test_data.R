
## R Script to prepare a test dataset for monthly rain across Ireland
## This file puts together 10 years (2016-2025) and puts together
## train and test folds for each month

## Data is collected from Met Éireann's network (dba.monthly_rain_clean)
## and NI data from the UK Met Office (CEDA Archive)

# %%

library(dplyr)

## Need the read_from_sol function to pull in data
source("utilities/read_from_sol.R")

sql_query <- paste0(
    "select * ",
    "from dba.monthly_rain_clean ",
    "where year >= 2016 and year <= 2025"
)
ROI_rainfall <- read_from_sol(
    sql_query
)

# %%

## Read monthly NI data already collected
## (check in Monthly Rainfall prepare_test_data.R for more info)
NI_rainfall <- read.csv("Data/Monthly_Rainfall/monthly_rain_NI_2016-2025.csv")

## Collate into monthly totals
library(lubridate)

NI_rainfall <- NI_rainfall |>
    group_by(stno, year, month) |>
    filter(n_distinct(day) == days_in_month(make_date(year, month, 1))) |>
    summarise(rain = sum(rain), .groups = "drop")

# %%

## Save these monthly dfs
ROI_rainfall <- ROI_rainfall |>
    dplyr::rename(rain = tot) |>
    dplyr::select(stno, year, month, rain)
# write.csv(
#     ROI_rainfall,
#     "Data/Monthly_Rainfall/monthly_rain_ROI_2016-2025.csv",
#     row.names = FALSE
# )

NI_rainfall <- NI_rainfall |>
    dplyr::select(stno, year, month, rain)
# write.csv(
#     NI_rainfall,
#     "Data/Monthly_Rainfall/monthly_rain_NI_2016-2025.csv",
#     row.names = FALSE
# )

# write.csv(
#     bind_rows(ROI_rainfall, NI_rainfall),
#     "Data/Monthly_Rainfall/monthly_rain_ROI_NI_2016-2025.csv",
#     row.names = FALSE
# )

# %%

## Geodata for stations
station_geodata <- read.csv("Data/General/station_geodata.csv")

## Get station type (for prioritising stations)
sql_query <- paste0(
    "select stno, stationtype ",
    "from dba.ymd_stations"
)
station_types <- read_from_sol(sql_query)
station_types <- station_types |>
    filter(stno %in% c(ROI_rainfall$stno, NI_rainfall$stno))

## Function for dropping pairs of stations that are too close
remove_close_stations <- function(
    data, threshold = 200, return_all = FALSE
) {

    stnos <- unique(data$stno)
    coords <- station_geodata |>
        filter(stno %in% stnos) |>
        dplyr::select(stno, east, north) |>
        distinct(stno, .keep_all = TRUE)

    ## Prioriy order: Synop -> TUSCON -> Climate -> CAMP -> Manual
    coords <- coords |>
        left_join(
            data |>
                select(stno, stationtype) |>
                distinct(),
            by = "stno"
        ) |>
        mutate(
            priority = case_when(
                stationtype == "tucson" ~ 1,
                stationtype == "climate_manual" ~ 2,
                stationtype == "me_automatic_climate_aws" ~ 3,
                stationtype == "aviation_aws" ~ 4,
                stationtype == "rainfall_manual" ~ 5,
                TRUE ~ 6
            )
        ) |>
        arrange(priority, stno)

    ## Get distance matrix and find close pairs
    dist_matrix <- sqrt(
        outer(coords$east, coords$east, "-")^2 +
            outer(coords$north, coords$north, "-")^2
    )
    close_pairs <- which(
        upper.tri(dist_matrix) & dist_matrix < threshold,
        arr.ind = TRUE
    )

    if (return_all) {
        ## Get stnos of each close pair
        close_stnos <-
            matrix(NA, nrow = nrow(close_pairs), ncol = 2)
        for (k in seq_len(nrow(close_pairs))) {
            i <- close_pairs[k, 1]
            j <- close_pairs[k, 2]
            stno1 <- coords$stno[i]
            stno2 <- coords$stno[j]

            close_stnos[k, ] <-
                c(stno1, stno2)
        }

        close_stnos
    } else {
        ## Drop worse station from each pair
        removed_stnos <- c()
        for (k in seq_len(nrow(close_pairs))) {
            i <- close_pairs[k, 1]
            j <- close_pairs[k, 2]
            stno1 <- coords$stno[i]
            stno2 <- coords$stno[j]
            p1 <- coords$priority[coords$stno == stno1]
            p2 <- coords$priority[coords$stno == stno2]

            if (p1 < p2) {
                removed_stnos <- c(removed_stnos, stno2)
            } else if (p2 > p1) {
                removed_stnos <- c(removed_stnos, stno1)
            } else {
                removed_stnos <- c(removed_stnos, stno2)
            }
        }
        removed_stnos
    }
}

# %%

## We'll look at each month
sample_data <-
    bind_rows(
        ROI_rainfall,
        NI_rainfall
    ) |>
    group_by(year, month) |>
    group_split()

## Create a 80/20 train/test split for each month
train_df <- list()
test_df  <- list()

## Loop through each month
for (i in seq_along(sample_data)) {

    ## Drop stations that are too close together
    removed_stations <- remove_close_stations(
        sample_data[[i]] |>
            left_join(station_types, by = "stno"),
        threshold = 200
    )
    sample_data[[i]] <- sample_data[[i]] |>
        filter(!stno %in% removed_stations)

    ## Available stations on current month
    stations <- sample_data[[i]] |>
        pull(stno)

    ## Randomly select 20% stations and assign to train/test
    test_stations <- sort(
        sample(
            stations,
            size = ceiling(length(stations) * 0.20),
            replace = FALSE
        )
    )
    train_df[[i]] <- sample_data[[i]] |>
        filter(!stno %in% test_stations)
    test_df[[i]] <- sample_data[[i]] |>
        filter(stno %in% test_stations)
}
test_df <- bind_rows(test_df)
train_df <- bind_rows(train_df)

# write.csv(
#     train_df,
#     "Data/Monthly_Rainfall/folds/train_split_2016-2025.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Monthly_Rainfall/folds/test_split_2016-2025.csv",
#     row.names = FALSE
# )

# %%

## The 80-20 split was just a preliminary validation method
## for screening methods. I have identified the best ones
## and need to run a more thorough 10-fold CV for them

rainfall <- read.csv(
    "Data/Monthly_Rainfall/monthly_rain_ROI_NI_2016-2025.csv"
)

close_stno_pairs <-
    remove_close_stations(
        rainfall |>
            group_by(stno) |>
            slice(1) |>
            left_join(
                station_geodata |>
                    dplyr::select(east, north, stno),
                by = "stno"
            ) |>
            left_join(station_types, by = "stno"),
        threshold = 200,
        return_all = TRUE
    )

## Folds need to be defined according to station
set.seed(222)
nfolds <- 10
station_folds <- sample(unique(rainfall$stno))
station_folds <- data.frame(
    stno = station_folds,
    fold = cut(
        seq_along(station_folds),
        breaks = nfolds,
        labels = FALSE
    )
)


## For close stations, manually set them to be in the same fold
close_stno_pairs <-
    remove_close_stations(
        rainfall |>
            group_by(stno) |>
            slice(1) |>
            left_join(
                station_geodata |>
                    dplyr::select(east, north, stno),
                by = "stno"
            ) |>
            left_join(station_types, by = "stno"),
        threshold = 200,
        return_all = TRUE
    )


for (i in seq_len(nrow(close_stno_pairs))) {
    stno_pair <- close_stno_pairs[i, ]
    first_stno_fold <- station_folds[
        station_folds$stno == stno_pair[1],
    ]$fold
    station_folds <- station_folds |>
        mutate(
            fold = if_else(
                stno == stno_pair[2],
                first_stno_fold,
                fold
            )
        )
}

rainfall <- rainfall |>
    dplyr::left_join(station_folds, by = "stno") |>
    left_join(
        station_geodata |>
            dplyr::select(stno, east, north),
        by = "stno"
    )  |>
    dplyr::distinct(year, month, east, north, .keep_all = TRUE) |>
    dplyr::select(stno, year, month, rain, fold)

# write.csv(
#     rainfall,
#     paste0(
#         "Data/Monthly_Rainfall/folds/",
#         nfolds,
#         "fold_2016-2025.csv"
#     ),
#     row.names = FALSE
# )
