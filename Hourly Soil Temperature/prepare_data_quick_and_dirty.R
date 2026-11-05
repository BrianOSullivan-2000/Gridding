
## I want to try and do a very rough CV for hourly soil temperature
## Will come back to this sooner to do something a bit more sophisitcated

## The core idea is we will get Rain and Temperature data too, plug it into
## a regression-kriging or spatio-temporal XGBoost and see how it does

# %%

## First load in some soil temperature data

source("utilities/read_from_sol.R")

soil_temperature_query <- paste("")