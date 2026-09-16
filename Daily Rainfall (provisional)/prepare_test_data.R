
## R Script to prepare a test dataset for daily rain across Ireland
## This file puts together all provisional daily rain data
## from bcoonan.daily_rain_provisional_clean as of 16-09-2026

# %%

library(dplyr)

## Need the read_from_sol function to pull in data
source("utilities/read_from_sol.R")

sql_query <- paste0(
    "select * ",
    "from bcoonan.daily_rain_provisional_clean"
)
rainfall <- read_from_sol(
    sql_query
)


# %%

## Save these monthly dfs
rainfall <- rainfall |>
    filter(!ind %in% c("qc_ok", "qc_ok_man", "bootstrapping", "regression")) |>
    dplyr::rename(rain = rainfall) |>
    dplyr::select(stno, year, month, day, rainfall)

write.csv(
    rainfall,
    "Data/Daily_Rainfall (provisional)/daily_rain_provisional.csv",
    row.names = FALSE
)



# %%

## Randomly sample 200 days
set.seed(222)
sample_200_days <- sample(
    seq(as.Date("2024-10-13"), as.Date("2026-09-15"), by = "day"),
    200,
    replace = FALSE
)
sample_200_days <- sort(sample_200_days)
sample_dates <- data.frame(
    year = as.integer(format(sample_200_days, "%Y")),
    month = as.integer(format(sample_200_days, "%m")),
    day = as.integer(format(sample_200_days, "%d"))
)

## Pull those days from the data
sample_data <- rainfall |>
    dplyr::semi_join(sample_dates, by = c("year", "month", "day")) |>
    group_by(year, month, day) |>
    group_split()

## We'll create a 80/20 train/test split for each of the 200 days
train_df <- list()
test_df  <- list()

## Loop through each day
for (i in seq_along(sample_data)) {

    ## Available stations on current day
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
#     "Data/Daily_Rainfall (provisional)/folds/train_split.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Daily_Rainfall (provisional)/folds/test_split.csv",
#     row.names = FALSE
# )
