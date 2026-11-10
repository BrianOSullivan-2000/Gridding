
## Validation tests for monthly rainfall grids ##

## Spatially embedded XGBoost (seXGBoost, need a better name) ##
## Predicts Monthly Rainfall using a standard XGBoost model   ##
## except spatial basis functions are also included to add    ##
## a spatial effect                                           ##

# %%

## Setup and data

## Setup and data
from pathlib import Path
from time import perf_counter

import dill
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.spatial import cKDTree

from utilities.add_basis_functions import add_basis_functions
from utilities.make_mesh import make_mesh

dill.load_session("Data/Monthly_Rainfall/train_test_80_20_2016-2025.dill")

session = globals()
rain_data = session["rain_data"]
dates = session["dates"]
island_outline = session["island_outline"]
get_monthly_rain_data = session["get_monthly_rain_data"]
update_monthly_predictions = session["update_monthly_predictions"]
collect_metrics = session["collect_metrics"]

# %%

## Spatial basis functions
nodelist = [
    make_mesh(island_outline, a=0.15, q=30, buffer=0, simplify=10000, plot=False),
    make_mesh(island_outline, a=0.02, q=30, buffer=2000, simplify=10000, plot=False),
]

threshold = 4000
for index in range(1, len(nodelist)):
    previous_nodes = np.vstack(nodelist[:index])
    distances, _ = cKDTree(previous_nodes).query(nodelist[index], k=1)
    nodelist[index] = nodelist[index][distances > threshold]

station_locations = (
    pd.concat([rain_data["train"], rain_data["test"]], ignore_index=True)
    .drop_duplicates(subset="stno")
    .loc[:, ["stno", "east", "north"]]
)
basis_columns = [
    f"B{index + 1}"
    for index in range(sum(len(nodes) for nodes in nodelist))
]
station_basis = pd.DataFrame(
    add_basis_functions(
        station_locations[["east", "north"]].to_numpy(),
        nodelist,
        supportlist=[20, 10],
        bfs_type="Mat",
        parlist=[
            (1, 60000, 0, 0.5),
            (1, 25000, 0, 0.5),
        ],
    ),
    columns=basis_columns
)
station_basis.index = station_locations["stno"].to_numpy()

for split in ("train", "test"):
    split_data = rain_data[split].assign(
        eastone=rain_data[split]["east"] / 1000,
        northone=rain_data[split]["north"] / 1000,
        aspectsin=np.sin(rain_data[split]["aspect_10km"]),
        aspectcos=np.cos(rain_data[split]["aspect_10km"]),
    )
    rain_data[split] = split_data

# %% 
 
## XGBoost model setup
geodata_features = [
    "eastone", "northone",
    "points5", "dist2c",
    "n5", "e5", "w5", "s5",
    "exp25k",
    "aspectsin", "aspectcos",
]
model_features = geodata_features + basis_columns
xgb_params = {
    "max_depth": 3,
    "eta": 0.05,
    "min_child_weight": 1,
    "objective": "reg:squarederror",
    "subsample": 0.2,
    "colsample_bytree": 0.1,
    "lambda": 50,
}
num_boost_round = 1600

## Prediction function
def predict_values(train_data, target_data):
    train_data = train_data.dropna(subset=["y"])

    Xtrain = xgb.DMatrix(
        train_data[model_features],
        label=train_data["y"].to_numpy(dtype=float),
    )
    model = xgb.train(
        xgb_params,
        Xtrain,
        num_boost_round=num_boost_round,
    )
    Xtest = xgb.DMatrix(target_data[model_features])
    return model.predict(Xtest)

def predict_monthly_values(monthly_data):
    return predict_values(monthly_data["train"], monthly_data["test"])


# %% 

## Validation test

for date_index in range(len(dates)):

    start_time = perf_counter()

    monthly_data = get_monthly_rain_data(rain_data, dates, date_index)

    for split in ("train", "test"):
        split_data = monthly_data[split]
        monthly_data[split] = split_data.join(station_basis, on="stno")
    
    predicted_y = predict_monthly_values(monthly_data)
    rain_data = update_monthly_predictions(
        rain_data,
        predicted_y,
        dates,
        date_index,
    )
    print(date_index)
    print(perf_counter() - start_time)

metrics_table = collect_metrics(
    None,
    "Spatial XGBoost",
    rain_data["test"]["rain"],
    rain_data["test"]["predicted_rain"],
)
print(metrics_table.to_string(index=False))

metrics_table.to_csv(
    "Results/Monthly_Rainfall/train_test_80_20_2016-2025/spatial_XGBoost.csv",
    index=False,
)

# %%

## Actual gridding, first need to set up basis functions
prediction_grid = session["grid_geodata"].copy()
grid_LTAs_9120 = session["grid_LTAs_9120"]
monthly_rain_plot = session["monthly_rain_plot"]

grid_basis = pd.DataFrame(
    add_basis_functions(
        prediction_grid[["east", "north"]].to_numpy(),
        nodelist,
        supportlist=[20, 10],
        bfs_type="Mat",
        parlist=[
            (1, 60000, 0, 0.5),
            (1, 25000, 0, 0.5),
        ],
    ),
    columns=basis_columns
)

## Prepare other predictors
prediction_grid = prediction_grid.assign(
    eastone=prediction_grid["east"] / 1000,
    northone=prediction_grid["north"] / 1000,
    aspectsin=np.sin(prediction_grid["aspect_10km"]),
    aspectcos=np.cos(prediction_grid["aspect_10km"]),
)

prediction_grid = pd.concat((prediction_grid, grid_basis), axis=1)

# %% 

## Make grids and check computation

grid_times = np.full(len(dates), np.nan, dtype=float)
figure_directory = Path("Figures/Monthly_Rainfall/spatial_XGBoost")
figure_directory.mkdir(parents=True, exist_ok=True)

# for date_index in range(len(dates)):
for date_index in range(10):
    start_time = perf_counter()
    current_date = dates.iloc[date_index]
    current_day = int(current_date["day"])
    current_month = int(current_date["month"])
    current_year = int(current_date["year"])

    monthly_rain_data = get_monthly_rain_data(
        rain_data,
        dates,
        date_index,
        experiment_type="all_data",
    )
    monthly_rain_data = monthly_rain_data.join(station_basis, on="stno")

    predicted_y = predict_values(monthly_rain_data, prediction_grid)
    prediction_grid["rain"] = (
        np.expm1(predicted_y)
        * grid_LTAs_9120[f"m_{current_month}"].to_numpy()
    )

    plot_destination = figure_directory / (
        f"spatial_XGBoost_{current_year}_{current_month:02d}_{current_day:02d}.jpg"
    )
    monthly_rain_plot(
        prediction_grid,
        monthly_rain_data,
        plot_destination=str(plot_destination),
    )

    grid_times[date_index] = perf_counter() - start_time
    print(date_index)

print("Mean Times")
print(np.nanmean(grid_times))

# %%

## Things to check

## Different basis functions (number, support, and type)
## Higher max_depths, higher lambdas
## Maybe subsample at 0.3, or colsubsample at 0.3?

## Use LTA as predictor instead of in normalization
## (So at the moment, it's worse. But, I want to run
## it through the hyperparameter optimization scheme
## that I have set up in the rough work)
## See it's overall performance is worse, but I
## also think it's better at not overfitting
