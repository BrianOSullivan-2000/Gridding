
# %%

## I need to figure this out

## QUESTIONS ##

# 1. Are all the stations that are in my test data present within both the geodata and the LTAs?
# Two separate checks here - for ROI stations and for NI stations


# 2. Is the geodata for each station correct? NI elevs should be from metadata, everything
# else is sampled from grid.


# 3. Are the LTAs from each station correct? These should all be sampled from grid
# there are a lot of LTAs (rainfall + all the different temperature types), so
# you need to check to make sure it was done properly (no duplicates!)


# %%

# Question 1

library(dplyr)

rain_stnos <- bind_rows(
    read.csv("Data/Daily_Rainfall/folds/test_split_2016-2025.csv"),
    read.csv("Data/Daily_Rainfall/folds/train_split_2016-2025.csv")
) |>
    group_by(stno) |>
    slice(1)

rain_prov_stnos <- bind_rows(
    read.csv("Data/Daily_Rainfall (provisional)/folds/test_split.csv"),
    read.csv("Data/Daily_Rainfall (provisional)/folds/train_split.csv")
) |>
    group_by(stno) |>
    slice(1)

temp_stnos <- bind_rows(
    read.csv("Data/Daily_Temperature/folds/test_split_2025.csv"),
    read.csv("Data/Daily_Temperature/folds/train_split_2025.csv")
) |>
    group_by(stno) |>
    slice(1)

geodata_stnos <- read.csv("Data/General/station_geodata.csv") |>
    group_by(stno) |>
    slice(1)

rain_LTA_stnos <- bind_rows(
    read.csv("Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"),
    read.csv("Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv")
) |>
    group_by(stno) |>
    slice(1)

temp_max_LTA_stnos <- bind_rows(
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1961-1990.csv"),
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1991-2020.csv")
) |>
    group_by(stno) |>
    slice(1)

temp_max_LTA_stnos <- bind_rows(
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1961-1990.csv"),
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1991-2020.csv")
) |>
    group_by(stno) |>
    slice(1)

temp_mean_LTA_stnos <- bind_rows(
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_mean_stations_1961-1990.csv"),
    read.csv("Data/General/LTAs/Temperature/LTAs_Temperature_mean_stations_1991-2020.csv")
) |>
    group_by(stno) |>
    slice(1)



rain_stnos |>
    filter(!stno %in% rain_LTA_stnos$stno)


## WHEN YOU GET BACk
## you have added stnos 2016, 2203, 2309, 2425, 4236, 6131, and 6427 to the station rainfall LTAs