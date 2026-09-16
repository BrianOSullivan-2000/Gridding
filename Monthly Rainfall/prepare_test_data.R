
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

## Read daily NI data already collected
## (check in Daily Rainfall prepare_test_data.R for more info)
NI_rainfall <- read.csv("Data/Daily_Rainfall/daily_rain_NI_2016-2025.csv")

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

## Get station type (for prioritising stations)
sql_query <- paste0(
    "select stno, stationtype ",
    "from dba.ymd_stations"
)
station_types <- read_from_sol(sql_query)
station_types <- station_types |>
    filter(stno %in% c(ROI_rainfall$stno, NI_rainfall$stno))

## Function for dropping pairs of stations that are too close
remove_close_stations <- function(data, threshold = 200) {
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
