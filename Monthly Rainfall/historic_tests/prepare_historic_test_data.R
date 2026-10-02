
## R Script to prepare a test dataset for monthly rain across Ireland
## This file puts together 57 years (1883-1940), samples 300 random
## and creates train and test folds for each of those months

## Data is collected from Met Éireann's historic data
## used to prepare monthly grids

# %%

library(dplyr)

ROI_rainfall <- read.csv(
    paste0(
        "Data/Monthly_rainfall/historic_data/",
        "roi_monthly_rain_1792-2024.csv"
    )
)
NI_rainfall <- read.csv(
    paste0(
        "Data/Monthly_rainfall/historic_data/",
        "ni_monthly_rain_1836-2023.csv"
    )
)

## Save these monthly dfs
ROI_rainfall <- ROI_rainfall |>
    dplyr::rename(rain = tot) |>
    dplyr::select(stno, year, month, rain, ind) |>
    filter(year >= 1883, year <= 1940)

NI_rainfall <- NI_rainfall |>
    dplyr::rename(rain = tot, stno = nistno) |>
    dplyr::select(stno, year, month, rain, ind) |>
    filter(year >= 1883, year <= 1940)

# write.csv(
#     bind_rows(ROI_rainfall, NI_rainfall),
#     "Data/Monthly_Rainfall/monthly_rain_ROI_NI_1883-1940.csv",
#     row.names = FALSE
# )

# %%

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
                select(stno) |>
                distinct(),
            by = "stno"
        )

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
        removed_stnos <- c(removed_stnos, stno2)
    }
    removed_stnos
}

# %%

## Randomly sample 300 months
set.seed(222)
sample_300_months <- sample(
    seq(as.Date("1883-01-01"), as.Date("1940-12-31"), by = "month"),
    300,
    replace = FALSE
)
sample_300_months <- sort(sample_300_months)
sample_dates <- data.frame(
    year = as.integer(format(sample_300_months, "%Y")),
    month = as.integer(format(sample_300_months, "%m"))
)

## Pull those days from the data
sample_data <-
    bind_rows(
        ROI_rainfall,
        NI_rainfall
    ) |>
    dplyr::semi_join(sample_dates, by = c("year", "month")) |>
    group_by(year, month) |>
    group_split()

## We'll create a 80/20 train/test split for each of the 300 months
train_df <- list()
test_df <- list()

## Loop through each month
for (i in seq_along(sample_data)) {

    ## Drop stations that are too close together
    removed_stations <- remove_close_stations(
        sample_data[[i]],
        threshold = 300
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
#     "Data/Monthly_Rainfall/folds/train_split_1883-1940.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Monthly_Rainfall/folds/test_split_1883-1940.csv",
#     row.names = FALSE
# )
