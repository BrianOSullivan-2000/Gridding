
# %%

## Code used to set up starting environment for each gridding method ##

try:
    from IPython import get_ipython
except ImportError:
    get_ipython = None

if get_ipython is not None:
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic("reset", "-f")

import warnings

import numpy as np
import pandas as pd

from utilities.metric_functions import (
    JSD,
    ME10,
    R2,
    mae,
    me,
    rmse,
)
from visualisation.base_grid_plot import base_grid_plot

## Station geodata and LTAs
station_geodata = pd.read_csv("Data/General/station_geodata.csv")
station_LTAs_6190 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)

## Grid geodata and LTAs
grid_geodata = pd.read_csv("Data/General/grid_geodata.csv")
grid_LTAs_6190 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1961-1990.csv"
)
grid_LTAs_9120 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1991-2020.csv"
)


## Island Outline for plotting
island_outline = pd.read_csv("Data/General/Irlcoast.csv")

## 80/20 Test (2016-2025)
train_data = pd.read_csv("Data/Monthly_Rainfall/folds/train_split_2016-2025.csv")
test_data = pd.read_csv("Data/Monthly_Rainfall/folds/test_split_2016-2025.csv")

rain_data = {
    "train": train_data.merge(station_geodata, on="stno", how="left")
    .merge(station_LTAs_9120, on="stno", how="left")
    .assign(predicted_rain=float("nan")),
    "test": test_data.merge(station_geodata, on="stno", how="left")
    .merge(station_LTAs_9120, on="stno", how="left")
    .assign(predicted_rain=float("nan")),
}

dates = (
    rain_data["train"][["year", "month"]]
    .drop_duplicates()
    .sort_values(["year", "month"])
    .reset_index(drop=True)
)

## Template for metrics
metrics_table = pd.DataFrame(
    columns=[
        "Name",
        "RMSE",
        "ME",
        "Computation_Time",
        "R2",
        "ME10",
    ]
)


## Function for getting metrics
def collect_metrics(metrics_table, name, y, y_hat):

    ## If any missing data in y_hat, drop it but give a warning
    y = np.asarray(y)
    y_hat = np.asarray(y_hat)
    missing_predictions = pd.isna(y_hat)
    missing_count = missing_predictions.sum()
    if missing_count:
        warnings.warn(
            f"Found {missing_count} NA values in y_hat; excluding them from metric calculations.",
            UserWarning,
        )
        y = y[~missing_predictions]
        y_hat = y_hat[~missing_predictions]

    new_row = pd.DataFrame(
        [
            {
                "Name": name,
                "RMSE": rmse(y_hat, y),
                "ME": me(y_hat, y),
                "Computation_Time": float("nan"),
                "MAE": mae(y_hat, y),
                "R2": R2(y_hat, y),
                "ME10": ME10(y_hat, y),
            }
        ]
    )

    if metrics_table is None:
        return new_row

    return pd.concat([metrics_table, new_row], ignore_index=True)


## Function that prepares data for a given day
def get_monthly_rain_data(
    rain_data,
    dates,
    date_index,
    experiment_type="train_test",
):
    if experiment_type not in ("train_test", "all_data"):
        raise ValueError("experiment_type must be 'train_test' or 'all_data'")

    current_date = dates.iloc[date_index]
    current_month = current_date["month"]
    current_year = current_date["year"]

    def prepare_monthly_data(data):
        monthly_data = data.loc[
            (data["month"] == current_month)
            & (data["year"] == current_year)
        ].copy()
        monthly_data["LTA"] = monthly_data[f"m_{int(current_month)}"]
        monthly_data["y"] = monthly_data["rain"] / monthly_data["LTA"]
        return monthly_data

    if experiment_type == "train_test":
        return {
            "train": prepare_monthly_data(rain_data["train"]),
            "test": prepare_monthly_data(rain_data["test"]),
        }

    all_data = pd.concat(
        [rain_data["train"], rain_data["test"]],
        ignore_index=True,
    )
    return prepare_monthly_data(all_data)


## Function that takes y_hat (predicted rain) and adds it to test set
def update_monthly_predictions(rain_data, y_hat, dates, date_index):
    monthly_rain_data = get_monthly_rain_data(
        rain_data,
        dates,
        date_index,
    )
    monthly_test_data = monthly_rain_data["test"]
    predicted_rain = (
        np.asarray(y_hat) * monthly_test_data["LTA"].to_numpy()
    )

    current_date = dates.iloc[date_index]
    prediction_mask = (
        (rain_data["test"]["month"] == current_date["month"])
        & (rain_data["test"]["year"] == current_date["year"])
    )
    updated_rain_data = rain_data.copy()
    updated_rain_data["test"] = rain_data["test"].copy()
    updated_rain_data["test"].loc[prediction_mask, "predicted_rain"] = predicted_rain

    return updated_rain_data


def monthly_rain_plot(grid, monthly_rain_data, plot_destination):
    return base_grid_plot(
        grid,
        stations=monthly_rain_data,
        response="rain",
        breaks=[
            0,
            20,
            50,
            100,
            150,
            200,
            250,
            300,
            400,
            500,
            600,
            800,
        ],
        bias=0.5,
        save_file=plot_destination,

        station_size=35
    )

# import dill

# dill.dump_session("Data/Monthly_Rainfall/train_test_80_20_2016-2025.dill")
