
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
    # Heavy_Wet_Day_Bias,
    # Rain_Day_Bias,
    Wet_Day_Bias,
    mae,
    me,
    rmse,
)
from visualisation.base_grid_plot import base_grid_plot

# from visualisation.prediction_plot import prediction_plot

## Station geodata and LTAs
station_geodata = pd.read_csv("Data/General/station_geodata.csv").assign(
    aspect_10km_cos=lambda data: np.cos(data["aspect_10km"]),
    aspect_10km_sin=lambda data: np.sin(data["aspect_10km"]),
)
station_LTAs_6190 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)

## Grid geodata and LTAs
grid_geodata = pd.read_csv("Data/General/grid_geodata.csv").assign(
    aspect_10km_cos=lambda data: np.cos(data["aspect_10km"]),
    aspect_10km_sin=lambda data: np.sin(data["aspect_10km"]),
)
grid_LTAs_6190 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1961-1990.csv"
)
grid_LTAs_9120 = pd.read_csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1991-2020.csv"
)


## Island Outline for plotting
island_outline = pd.read_csv("Data/General/Irlcoast.csv")

## 80/20 Test
train_data = pd.read_csv("Data/Daily_Rainfall/folds/train_split_80_20_provisional.csv")
test_data = pd.read_csv("Data/Daily_Rainfall/folds/test_split_80_20_provisional.csv")

rain_data = {
    "train": train_data.merge(station_geodata, on="stno", how="left")
    .merge(station_LTAs_9120, on="stno", how="left")
    .assign(predicted_rain=float("nan")),
    "test": test_data.merge(station_geodata, on="stno", how="left")
    .merge(station_LTAs_9120, on="stno", how="left")
    .assign(predicted_rain=float("nan")),
}

dates = (
    rain_data["train"][["year", "month", "day"]]
    .drop_duplicates()
    .sort_values(["year", "month", "day"])
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
        "JSD",
        "ME10",
        "Wet_Day_Bias",
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
                "JSD": JSD(y, y_hat),
                "ME10": ME10(y_hat, y),
                "Wet_Day_Bias": Wet_Day_Bias(y_hat, y),
            }
        ]
    )

    if metrics_table is None:
        return new_row

    return pd.concat([metrics_table, new_row], ignore_index=True)


## Function that prepares data for a given day
def get_daily_rain_data(
    rain_data,
    dates,
    date_index,
    experiment_type="train_test",
):
    if experiment_type not in ("train_test", "all_data"):
        raise ValueError("experiment_type must be 'train_test' or 'all_data'")

    current_date = dates.iloc[date_index]
    current_day = current_date["day"]
    current_month = current_date["month"]
    current_year = current_date["year"]

    def prepare_daily_data(data):
        daily_data = data.loc[
            (data["day"] == current_day)
            & (data["month"] == current_month)
            & (data["year"] == current_year)
        ].copy()
        daily_data["LTA"] = daily_data[f"m_{int(current_month)}"]
        daily_data["normalized_rain"] = daily_data["rain"] / daily_data["LTA"]
        daily_data["y"] = np.log1p(daily_data["normalized_rain"])
        return daily_data

    if experiment_type == "train_test":
        return {
            "train": prepare_daily_data(rain_data["train"]),
            "test": prepare_daily_data(rain_data["test"]),
        }

    all_data = pd.concat(
        [rain_data["train"], rain_data["test"]],
        ignore_index=True,
    )
    return prepare_daily_data(all_data)


## Function that takes y_hat (predicted rain) and adds it to test set
def update_daily_predictions(rain_data, y_hat, dates, date_index):
    daily_rain_data = get_daily_rain_data(
        rain_data,
        dates,
        date_index,
    )
    daily_test_data = daily_rain_data["test"]
    predicted_rain = (
        np.expm1(np.asarray(y_hat)) * daily_test_data["LTA"].to_numpy()
    )

    current_date = dates.iloc[date_index]
    prediction_mask = (
        (rain_data["test"]["day"] == current_date["day"])
        & (rain_data["test"]["month"] == current_date["month"])
        & (rain_data["test"]["year"] == current_date["year"])
    )
    updated_rain_data = rain_data.copy()
    updated_rain_data["test"] = rain_data["test"].copy()
    updated_rain_data["test"].loc[prediction_mask, "predicted_rain"] = predicted_rain

    return updated_rain_data


def daily_rain_plot(grid, daily_rain_data, plot_destination):
    return base_grid_plot(
        grid,
        stations=daily_rain_data,
        response="rain",
        breaks=[
            0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.75,
            1,
            2,
            3,
            4,
            5,
            7.5,
            10,
            20,
            30,
            40,
            50,
            75,
            100,
        ],
        bias=0.5,
        save_file=plot_destination,
    )

# import dill

# dill.dump_session("Data/Daily_Rainfall/train_test_80_20_provisional.dill")
