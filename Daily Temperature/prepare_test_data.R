
## R Script to prepare a test dataset for daily max/min temperature
## This file puts together maximum and minimum temperatures
## for 2025 and puts together train and test folds

## Data is collected from Met Éireann's network
## This includes QC data from climate stations (dba.allclimat)
## QC data from CAMP (bcoonan.daily_temperature_provisional_clean)
## and NI QC data from the UK Met Office (CEDA Archive)

# %%

library(dplyr)
library(tidyr)

## Need the read_from_sol function to pull in data
source("utilities/read_from_sol.R")

## Read stations from allclimat
sql_query <- paste0(
    "select  stno, year, month, day, maxt, kmaxt, mint, kmint ",
    "from dba.allclimat ",
    "where year = 2025"
)
allclimat_temperature <- read_from_sol(sql_query)

# %%

sql_query <- paste0(
    "select stno, stn_type, year, month, day, element, temperature, ind ",
    "from bcoonan.daily_temperature_provisional_clean ",
    "where year = 2025"
)
camp_temperature <- read_from_sol(sql_query)

# %%

## Next we'll get NI data from the CEDA Archive
## Get metadata
NI_metadata <- read.csv(
    paste0(
        "Data/Daily_Temperature/NI_Data/",
        "midas-open_uk-daily-temperature-obs_dv-202607_station-metadata.csv"
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

        for (year in 2025) {
            url <- paste0(
                "https://dap.ceda.ac.uk/badc/ukmo-midas-open/",
                "data/uk-daily-temperature-obs/dataset-version-202607/",
                county, "/", stno, "_", station,
                "/qc-version-1/midas-open_uk-daily-temperature-obs_dv-202607_",
                county, "_", stno, "_", station, "_qcv-1_", year, ".csv"
            )

            urls <- c(urls, url)
        }
    }
}

# writeLines(urls, "Data/Daily_Temperature/NI_Data/NI-2025.txt")

# %%

library(lubridate)

## Get geodata and LTAs
## when getting elevations get more accurate ones directly from UK Met files

## Get names of all files
NI_files <- list.files(
    "Data/Daily_Temperature/NI_Data/",
    pattern = "\\.csv$", full.names = TRUE
)

## Open each file, convert to similar format as ROI data
NI_dfs <- lapply(NI_files, read.csv, skip = 93)

NI_dfs <- lapply(
    NI_dfs,
    function(df) {
        main_interval <-
            as.integer(
                names(sort(table(df$ob_hour_count), decreasing = TRUE))[1]
            )

        ## Combine 12 hour periods into 24 hour
        if (main_interval == 12) {
            df <- df |>
                filter(ob_hour_count == 12) |>
                mutate(
                    ob_end_time = ymd_hms(ob_end_time),
                    day = as.Date(ob_end_time - hours(9))
                ) |>
                group_by(day) |>
                mutate(
                    max_air_temp = max(max_air_temp, na.rm = TRUE),
                    min_air_temp = min(min_air_temp, na.rm = TRUE),
                    .groups = "drop"
                )
        } else if (main_interval == 24) {
            df <- df |>
                filter(ob_hour_count == 24)
        } else {
            return(NULL)
        }

        df$t <- as.numeric(as.Date(df$ob_end_time) - as.Date("2025-01-01")) + 1
        df <- df |>
            dplyr::rename(
                stno = src_id,
                maxt = max_air_temp,
                mint = min_air_temp
            ) |>
            group_by(t) |>
            summarise(
                stno = first(stno),
                maxt = first(maxt),
                mint = first(mint),
                .groups = "drop"
            )
        df
    }
)

## Get longitude, latitude, and elevation of each station
NI_coords <- lapply(
    NI_files,
    function(filepath) {
        coords <- read.csv(filepath, header = FALSE)

        stno <- as.numeric(coords[101, 3])
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

library(gstat)

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
    format(as.Date(NI_df$t - 1, origin = "2025-01-01"), "%Y")
)
NI_df$month <- as.integer(
    format(as.Date(NI_df$t - 1, origin = "2025-01-01"), "%m")
)
NI_df$day <- as.integer(
    format(as.Date(NI_df$t - 1, origin = "2025-01-01"), "%d")
)

## Add 99 to these station IDs
NI_df$stno <- NI_df$stno + 9900000
NI_coords$stno <- NI_coords$stno + 9900000
NI_df$niflag <- TRUE

# %%

## Similarly, we need to get LTAs for each station
LTAs_max_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_grid_1961-1990.csv"
    )
)
LTAs_max_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_grid_1991-2020.csv"
    )
)
LTAs_min_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1961-1990.csv"
    )
)
LTAs_min_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1991-2020.csv"
    )
)

NI_coordinates <- NI_coords |>
    dplyr::select(stno, east, north) |>
    mutate(stno = stno - 9900000)
LTA_vars <- c(
    paste0("m_", 1:12),
    "ann"
)
NI_coordinates[LTA_vars] <- NA
NI_coords_max_6190 <- NI_coordinates
NI_coords_max_9120 <- NI_coordinates
NI_coords_min_6190 <- NI_coordinates
NI_coords_min_9120 <- NI_coordinates

LTAs_max_6190_sf <-
    st_as_sf(LTAs_max_6190, coords = c("east", "north"), crs = NA)
LTAs_max_9120_sf <-
    st_as_sf(LTAs_max_9120, coords = c("east", "north"), crs = NA)
LTAs_min_6190_sf <-
    st_as_sf(LTAs_min_6190, coords = c("east", "north"), crs = NA)
LTAs_min_9120_sf <-
    st_as_sf(LTAs_min_9120, coords = c("east", "north"), crs = NA)

NI_max_6190_sf <- st_as_sf(
    NI_coords_max_6190,
    coords = c("east", "north"), crs = NA
)
NI_max_9120_sf <- st_as_sf(
    NI_coords_max_9120,
    coords = c("east", "north"), crs = NA
)
NI_min_6190_sf <- st_as_sf(
    NI_coords_min_6190,
    coords = c("east", "north"), crs = NA
)
NI_min_9120_sf <- st_as_sf(
    NI_coords_min_9120,
    coords = c("east", "north"), crs = NA
)


for (v in LTA_vars) {
    result_max_6190 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_max_6190_sf,
        newdata = NI_max_6190_sf,
        idp = 2,
        nmax = 8
    )
    result_max_9120 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_max_9120_sf,
        newdata = NI_max_9120_sf,
        idp = 2,
        nmax = 8
    )
    result_min_6190 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_min_6190_sf,
        newdata = NI_min_6190_sf,
        idp = 2,
        nmax = 8
    )
    result_min_9120 <- idw(
        formula = as.formula(paste(v, "~ 1")),
        locations = LTAs_min_9120_sf,
        newdata = NI_min_9120_sf,
        idp = 2,
        nmax = 8
    )
    NI_coords_max_6190[[v]] <- result_max_6190$var1.pred
    NI_coords_max_9120[[v]] <- result_max_9120$var1.pred
    NI_coords_min_6190[[v]] <- result_min_6190$var1.pred
    NI_coords_min_9120[[v]] <- result_min_9120$var1.pred
}

NI_coords_max_6190 <- NI_coords_max_6190 |>
    dplyr::select(-east, -north) |>
    mutate(stno = stno + 9900000)
NI_coords_max_9120 <- NI_coords_max_9120 |>
    dplyr::select(-east, -north) |>
    mutate(stno = stno + 9900000)
NI_coords_min_6190 <- NI_coords_min_6190 |>
    dplyr::select(-east, -north) |>
    mutate(stno = stno + 9900000)
NI_coords_min_9120 <- NI_coords_min_9120 |>
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
station_LTAs_max_6190 <- read.csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1961-1990.csv"
)
station_LTAs_max_9120 <- read.csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1991-2020.csv"
)
station_LTAs_min_6190 <- read.csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1961-1990.csv"
)
station_LTAs_min_9120 <- read.csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1991-2020.csv"
)


station_LTAs_max_6190 <- bind_rows(
    station_LTAs_max_6190,
    NI_coords_max_6190 |>
        anti_join(station_LTAs_max_6190, by = "stno")
)
station_LTAs_max_9120 <- bind_rows(
    station_LTAs_max_9120,
    NI_coords_max_9120 |>
        anti_join(station_LTAs_max_9120, by = "stno")
)
station_LTAs_min_6190 <- bind_rows(
    station_LTAs_min_6190,
    NI_coords_min_6190 |>
        anti_join(station_LTAs_min_6190, by = "stno")
)
station_LTAs_min_9120 <- bind_rows(
    station_LTAs_min_9120,
    NI_coords_min_9120 |>
        anti_join(station_LTAs_min_9120, by = "stno")
)

# write.csv(
#     station_LTAs_max_6190,
#     paste0(
#         "Data/General/LTAs/Temperature/",
#         "LTAs_Temperature_max_stations_1961-1990.csv"
#     ),
#     row.names = FALSE
# )

# write.csv(
#     station_LTAs_max_9120,
#     paste0(
#         "Data/General/LTAs/Temperature/",
#         "LTAs_Temperature_max_stations_1991-2020.csv"
#     ),
#     row.names = FALSE
# )

# write.csv(
#     station_LTAs_min_6190,
#     paste0(
#         "Data/General/LTAs/Temperature/",
#         "LTAs_Temperature_min_stations_1961-1990.csv"
#     ),
#     row.names = FALSE
# )

# write.csv(
#     station_LTAs_min_9120,
#     paste0(
#         "Data/General/LTAs/Temperature/",
#         "LTAs_Temperature_min_stations_1991-2020.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Combine temperature observations

allclimat_temperature <- allclimat_temperature |>
    filter(kmint == 0, kmaxt == 0) |>
    dplyr::select(stno, year, month, day, mint, maxt)

camp_temperature <- camp_temperature |>
    filter(!stno %in% unique(allclimat_temperature$stno)) |>
    filter(ind %in% c("qc_ok", "bootstrapping", "regression")) |>
    pivot_wider(
        names_from = element,
        values_from = temperature
    ) |>
    dplyr::rename(maxt = maxt9, mint = mint9) |>
    dplyr::select(stno, year, month, day, mint, maxt)

## Save temperature observations

ROI_temperature <- bind_rows(
    allclimat_temperature,
    camp_temperature
) |>
    filter(!is.na(maxt), !is.na(mint))

# write.csv(
#     ROI_temperature,
#     "Data/Daily_Temperature/daily_temperature_ROI_2025.csv",
#     row.names = FALSE
# )

NI_temperature <- NI_df |>
    dplyr::select(stno, year, month, day, mint, maxt) |>
    filter(!is.na(maxt), !is.na(mint))

# write.csv(
#     NI_temperature,
#     "Data/Daily_Temperature/daily_temperature_NI_2025.csv",
#     row.names = FALSE
# )

# write.csv(
#     bind_rows(ROI_temperature, NI_temperature),
#     "Data/Daily_Temperature/daily_temperature_ROI_NI_2025.csv",
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
    filter(stno %in% c(ROI_temperature$stno, NI_temperature$stno))

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
    seq(as.Date("2025-01-01"), as.Date("2025-12-31"), by = "day"),
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
        ROI_temperature,
        NI_temperature
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
#     "Data/Daily_Temperature/folds/train_split_2025.csv",
#     row.names = FALSE
# )
# write.csv(
#     test_df,
#     "Data/Daily_Temperature/folds/test_split_2025.csv",
#     row.names = FALSE
# )
