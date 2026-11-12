
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
    as.data.frame() |>
    filter(ind %in% c("qc_ok", "qc_ok_man", "bootstrapping", "regression")) |>
    dplyr::rename(rain = rainfall) |>
    distinct(east, north, year, month, day, .keep_all = TRUE) |>
    dplyr::select(stno, year, month, day, rain)

# write.csv(
#     rainfall,
#     "Data/Daily_Rainfall/daily_rain_provisional.csv",
#     row.names = FALSE
# )

# %%

## Some stations don't have geodata or LTAs
## interpolate them from the grid

## Do geodata first

library(sf)
library(gstat)

station_geodata <- read.csv("Data/General/station_geodata.csv")
grid <- read.csv("Data/General/grid_geodata.csv")

missing_geodata_stnos <-
    unique(rainfall$stno)[
        !(unique(rainfall$stno) %in% unique(station_geodata$stno))
    ]

## Read in grid
grid_sf <- st_as_sf(grid, coords = c("east", "north"), crs = NA)
stno_sf <- rainfall |>
    group_by(stno) |>
    slice(1) |>
    filter(stno %in% missing_geodata_stnos) |>
    st_as_sf(
        coords = c("east", "north"), crs = NA
    )

## Add empty columns
geodata_vars <- names(grid)[4:ncol(grid)]
missing_stno_coords <- rainfall |>
    group_by(stno) |>
    slice(1) |>
    filter(stno %in% missing_geodata_stnos) |>
    dplyr::select(stno, east, north)

## Interpolate from grid
for (v in geodata_vars) {
    result <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = grid_sf,
        newdata = stno_sf,
        idp = 2,
        nmax = 8
    )
    missing_stno_coords[[v]] <- result$var1.pred
}

## If you have the stno name, add it here, e.g.
# missing_stno_coords$name <- c("BALLYCANEW", "KILKENNY ACS")
missing_stno_coords <- missing_stno_coords |>
    dplyr::select(all_of(names(station_geodata)))

station_geodata <- bind_rows(
    station_geodata,
    missing_stno_coords
)

# write.csv(
#     station_geodata,
#     "Data/General/station_geodata.csv",
#     row.names = FALSE
# )

# %%

## Do the same for LTAs

station_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)

LTAs_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Rainfall/",
        "LTAs_Rainfall_grid_1961-1990.csv"
    )
)
LTAs_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Rainfall/",
        "LTAs_Rainfall_grid_1991-2020.csv"
    )
)

## NOTE - the code assumes anything missing in 6190 LTAs
## are also going to be missing in the 9120 LTAs
missing_LTA_stnos <-
    unique(rainfall$stno)[
        !(unique(rainfall$stno) %in% unique(station_LTAs_6190$stno))
    ]

LTA_vars <- c(
    paste0("m_", 1:12),
    "Ann",
    "spring", "summer", "autumn", "winter"
)

missing_stno_coords_6190 <- rainfall |>
    group_by(stno) |>
    slice(1) |>
    filter(stno %in% missing_LTA_stnos) |>
    dplyr::select(stno, east, north)
missing_stno_coords_9120 <- missing_stno_coords_6190
missing_stno_coords_6190[LTA_vars] <- NA
missing_stno_coords_9120[LTA_vars] <- NA

LTAs_6190_sf <- st_as_sf(LTAs_6190, coords = c("east", "north"), crs = NA)
LTAs_9120_sf <- st_as_sf(LTAs_9120, coords = c("east", "north"), crs = NA)
stno_6190_sf <- st_as_sf(
    missing_stno_coords_6190, coords = c("east", "north"), crs = NA
)
stno_9120_sf <- st_as_sf(
    missing_stno_coords_9120, coords = c("east", "north"), crs = NA
)

for (v in LTA_vars) {
    result_6190 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_6190_sf,
        newdata = stno_6190_sf,
        idp = 2,
        nmax = 8
    )
    result_9120 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_9120_sf,
        newdata = stno_9120_sf,
        idp = 2,
        nmax = 8
    )
    missing_stno_coords_6190[[v]] <- result_6190$var1.pred
    missing_stno_coords_9120[[v]] <- result_9120$var1.pred
}

station_LTAs_6190 <- bind_rows(
    station_LTAs_6190,
    missing_stno_coords_6190
)
station_LTAs_9120 <- bind_rows(
    station_LTAs_9120,
    missing_stno_coords_9120
)

# write.csv(
#     station_LTAs_6190,
#     "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv",
#     row.names = FALSE
# )

# write.csv(
#     station_LTAs_9120,
#     "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv",
#     row.names = FALSE
# )

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
#     "Data/Daily_Rainfall/folds/train_split_80_20_provisional.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Daily_Rainfall/folds/test_split_80_20_provisional.csv",
#     row.names = FALSE
# )
