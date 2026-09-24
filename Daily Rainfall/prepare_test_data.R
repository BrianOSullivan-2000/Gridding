
## R Script to prepare a test dataset for daily rain across Ireland
## This file puts together 10 year of daily rain (2016-2025), randomly
## samples 200 days, and puts together train and test folds

## Data is collected from Met Éireann's network (dba.daily_rain_clean)
## and NI data from the UK Met Office (CEDA Archive)

# %%

library(dplyr)

## Need the read_from_sol function to pull in data
source("utilities/read_from_sol.R")

sql_query <- paste0(
    "select * ",
    "from dba.daily_rain_clean ",
    "where year >= 2016 and year <= 2025"
)
ROI_rainfall <- read_from_sol(
    sql_query
)


# %%

## Next we'll get NI data from the CEDA Archive
## Get metadata
NI_metadata <- read.csv(
    paste0(
        "Data/Daily_Rainfall/NI_Data/",
        "midas-open_uk-daily-rain-obs_dv-202607_station-metadata.csv"
    ),
    skip = 48,
    header = TRUE
)
counties <- c("antrim", "armagh", "down", "derry", "tyrone", "fermanagh")

urls <- c()


for (county in counties) {

    stnos <- unique(NI_metadata[NI_metadata$historic_county == county, ]$src_id)
    stations <- unique(
        NI_metadata[NI_metadata$historic_county == county, ]$station_file_name
    )

    for (i in seq_along(stnos)) {

        stno <- stnos[i]
        station <- stations[i]

        for (year in 2016:2025) {
            url <- paste0(
                "https://dap.ceda.ac.uk/badc/ukmo-midas-open/",
                "data/uk-daily-rain-obs/dataset-version-202607/",
                county, "/", stno, "_", station,
                "/qc-version-1/midas-open_uk-daily-rain-obs_dv-202607_",
                county, "_", stno, "_", station, "_qcv-1_", year, ".csv"
            )

            urls <- c(urls, url)
        }
    }
}

# writeLines(urls, "Data/Daily_Rainfall/NI_Data/NI-2016-2025.txt")

# %%

## Get geodata and LTAs
## when getting elevations get more accurate ones directly from UK Met files

## Get names of all files
NI_files <- list.files(
    "Data/Daily_Rainfall/NI_Data/",
    pattern = "\\.csv$", full.names = TRUE
)

## Open each file, convert to similar format as ROI data
NI_dfs <- lapply(NI_files, read.csv, skip = 64)
NI_dfs <- lapply(
    NI_dfs,
    function(df) {
        df$ob_date <- as.Date(df$ob_date)
        df$t <- as.numeric(as.Date(df$ob_date) - as.Date("2016-01-01")) + 1
        df <- df |>
            dplyr::rename(stno = src_id, rain = prcp_amt) |>
            filter(ob_day_cnt == 1) |>
            dplyr::select(stno, rain, t)
        df
    }
)

## Get longitude, latitude, and elevation of each station
NI_coords <- lapply(
    NI_files,
    function(filepath) {
        coords <- read.csv(filepath, header = FALSE)

        stno <- as.numeric(coords[70, 4])
        name <- coords[11, 3]
        longitude <- as.numeric(coords[17, 4])
        latitude <- as.numeric(coords[17, 3])
        elevation <- as.numeric(coords[18, 3])

        data.frame(
            stno = stno, name = name,
            lon = longitude, lat = latitude, elev = elevation
        )
    }
)
NI_coords <- bind_rows(NI_coords) |>
    group_by(stno) |>
    slice(1)

# %%

## Now we need to add a few more geodata values
## First get easting and northing by changing the projection
library(sf)
library(gstat)
NI_coords <- NI_coords |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(29903) |>
    mutate(
        east  = st_coordinates(geometry)[, 1],
        north = st_coordinates(geometry)[, 2]
    ) |>
    st_drop_geometry()

## For other geodata, assign each station geodata values by
## interpolating from the underlying grid

## Read in grid
grid <- read.csv("Data/General/grid_geodata.csv")
grid_sf <- st_as_sf(grid, coords = c("east", "north"), crs = NA)
NI_sf <- st_as_sf(NI_coords, coords = c("east", "north"), crs = NA)

# Add empty columns
geodata_vars <- names(grid)[5:ncol(grid)]
NI_coords[geodata_vars] <- NA

for (v in geodata_vars) {
    result <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = grid_sf,
        newdata = NI_sf,
        idp = 2,
        nmax = 8
    )
    NI_coords[[v]] <- result$var1.pred
}

# %%

## Finally we can add the geodata to observation data
NI_df <- bind_rows(NI_dfs)
NI_df <- NI_df |> left_join(NI_coords, by = "stno")
NI_df <- NI_df[complete.cases(NI_df), ]

## Get specific month and day back
NI_df$year <- as.integer(
    format(as.Date(NI_df$t - 2, origin = "2016-01-01"), "%Y")
)
NI_df$month <- as.integer(
    format(as.Date(NI_df$t - 2, origin = "2016-01-01"), "%m")
)
NI_df$day <- as.integer(
    format(as.Date(NI_df$t - 2, origin = "2016-01-01"), "%d")
)

## Add 99 to these station IDs
NI_df$stno <- NI_df$stno + 9900000
NI_coords$stno <- NI_coords$stno + 9900000
NI_df$niflag <- TRUE

# %%

## Similarly, we need to get LTAs for each station
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

NI_coords_6190 <- NI_coords |>
    dplyr::select(stno, east, north) |>
    mutate(stno = stno - 9900000)
NI_coords_9120 <- NI_coords_6190
LTA_vars <- c(
    paste0("m_", 1:12),
    "Ann",
    "spring", "summer", "autumn", "winter"
)
NI_coords_6190[LTA_vars] <- NA
NI_coords_9120[LTA_vars] <- NA

LTAs_6190_sf <- st_as_sf(LTAs_6190, coords = c("east", "north"), crs = NA)
LTAs_9120_sf <- st_as_sf(LTAs_9120, coords = c("east", "north"), crs = NA)
NI_6190_sf <- st_as_sf(NI_coords_6190, coords = c("east", "north"), crs = NA)
NI_9120_sf <- st_as_sf(NI_coords_9120, coords = c("east", "north"), crs = NA)

for (v in LTA_vars) {
    result_6190 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_6190_sf,
        newdata = NI_6190_sf,
        idp = 2,
        nmax = 8
    )
    result_9120 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_9120_sf,
        newdata = NI_9120_sf,
        idp = 2,
        nmax = 8
    )
    NI_coords_6190[[v]] <- result_6190$var1.pred
    NI_coords_9120[[v]] <- result_9120$var1.pred
}

NI_coords_6190 <- NI_coords_6190 |>
    dplyr::select(-east, -north) |>
    mutate(stno = stno + 9900000)
NI_coords_9120 <- NI_coords_9120 |>
    dplyr::select(-east, -north) |>
    mutate(stno = stno + 9900000)

# %%

## Update station geodata to include NI stations
station_geodata <- read.csv("Data/General/station_geodata.csv")

station_geodata <- bind_rows(
    station_geodata,
    NI_coords |>
        anti_join(station_geodata, by = "stno")
)

# write.csv(
#     station_geodata,
#     "Data/General/station_geodata.csv",
#     row.names = FALSE
# )

## Do the same for LTAs
station_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)

station_LTAs_6190 <- bind_rows(
    station_LTAs_6190,
    NI_coords_6190 |>
        anti_join(station_LTAs_6190, by = "stno")
)

station_LTAs_9120 <- bind_rows(
    station_LTAs_9120,
    NI_coords_9120 |>
        anti_join(station_LTAs_9120, by = "stno")
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

## Save rainfall observations

ROI_rainfall <- ROI_rainfall |>
    dplyr::select(stno, year, month, day, rain, ind)
# write.csv(
#     ROI_rainfall,
#     "Data/Daily_Rainfall/daily_rain_ROI_2016-2025.csv",
#     row.names = FALSE
# )

NI_rainfall <- NI_df |>
    mutate(ind = NA) |>
    dplyr::select(stno, year, month, day, rain, ind)
# write.csv(
#     NI_rainfall,
#     "Data/Daily_Rainfall/daily_rain_NI_2016-2025.csv",
#     row.names = FALSE
# )

# write.csv(
#     bind_rows(ROI_rainfall, NI_rainfall),
#     "Data/Daily_Rainfall/daily_rain_ROI_NI_2016-2025.csv",
#     row.names = FALSE
# )

# %%

## We're going to drop stations that are too close (< 200m)
## Given we're testing gridding methods and most grid points
## are not right next to other stations, this should make
## our results more realistic

## First we'll also need data on the station type
## so we know which station to prioritise if there is a conflict
sql_query <- paste0(
    "select stno, stationtype ",
    "from dba.ymd_stations"
)
station_types <- read_from_sol(sql_query)
station_types <- station_types |>
    filter(stno %in% c(ROI_rainfall$stno, NI_rainfall$stno))

## Function for checking nearby stations
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

## Randomly sample 200 days
set.seed(222)
sample_200_days <- sample(
    seq(as.Date("2016-01-01"), as.Date("2025-12-31"), by = "day"),
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
sample_data <-
    bind_rows(
        ROI_rainfall,
        NI_rainfall
    ) |>
    dplyr::semi_join(sample_dates, by = c("year", "month", "day")) |>
    group_by(year, month, day) |>
    group_split()

## We'll create a 80/20 train/test split for each of the 200 days
train_df <- list()
test_df  <- list()

## Loop through each day
for (i in seq_along(sample_data)) {

    ## Drop stations that are too close together
    removed_stations <- remove_close_stations(
        sample_data[[i]] |>
            left_join(station_types, by = "stno"),
        threshold = 200
    )
    sample_data[[i]] <- sample_data[[i]] |>
        filter(!stno %in% removed_stations)

    ## Available stations on current day
    stations <- sample_data[[i]] |>
        filter(ind == 0 | is.na(ind)) |>
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
#     "Data/Daily_Rainfall/folds/train_split_2016-2025.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Daily_Rainfall/folds/test_split_2016-2025.csv",
#     row.names = FALSE
# )
