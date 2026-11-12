
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
    ME10,
    R2,
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
station_LTAs_max_6190 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1961-1990.csv"
)
station_LTAs_max_9120 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_stations_1991-2020.csv"
)
station_LTAs_min_6190 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1961-1990.csv"
)
station_LTAs_min_9120 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_stations_1991-2020.csv"
)

## Grid geodata and LTAs
grid_geodata = pd.read_csv("Data/General/grid_geodata.csv").assign(
    aspect_10km_cos=lambda data: np.cos(data["aspect_10km"]),
    aspect_10km_sin=lambda data: np.sin(data["aspect_10km"]),
)
LTAs_max_6190 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_grid_1961-1990.csv"
).merge(grid_geodata[["east", "north"]], on=["east", "north"], how="inner")

LTAs_max_9120 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_max_grid_1991-2020.csv"
).merge(grid_geodata[["east", "north"]], on=["east", "north"], how="inner")

LTAs_min_6190 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_grid_1961-1990.csv"
).merge(grid_geodata[["east", "north"]], on=["east", "north"], how="inner")

LTAs_min_9120 = pd.read_csv(
    "Data/General/LTAs/Temperature/LTAs_Temperature_min_grid_1991-2020.csv"
).merge(grid_geodata[["east", "north"]], on=["east", "north"], how="inner")


## Island Outline for plotting
island_outline = pd.read_csv("Data/General/Irlcoast.csv")

## 80/20 Test (2016-2025)
train_data = pd.read_csv("Data/Daily_Temperature/folds/train_split_2025.csv")
test_data = pd.read_csv("Data/Daily_Temperature/folds/test_split_2025.csv")

temperature_data = {}
for split_name, data in {"train": train_data, "test": test_data}.items():
    temperature_data[split_name] = (
        data.merge(station_geodata, on="stno", how="left")
        .merge(station_LTAs_max_9120, on="stno", how="left")
        .drop(columns="ann")
        .rename(
            columns={
                column: f"maxt_{column}"
                for column in station_LTAs_max_9120.columns
                if column.startswith("m_")
            }
        )
        .merge(station_LTAs_min_9120, on="stno", how="left")
        .drop(columns="ann")
        .rename(
            columns={
                column: f"mint_{column}"
                for column in station_LTAs_min_9120.columns
                if column.startswith("m_")
            }
        )
        .assign(predicted_maxt=float("nan"), predicted_mint=float("nan"))
    )

dates = (
    temperature_data["train"][["year", "month", "day"]]
    .drop_duplicates()
    .sort_values(["year", "month", "day"])
    .reset_index(drop=True)
)

## Template for metrics
metrics_table = pd.DataFrame(
    columns=[
        "Element",
        "Name",
        "RMSE",
        "ME",
        "Computation_Time",
        "R2",
        "ME10",
    ]
)


## Function for getting metrics
def collect_metrics(metrics_table, name, y_maxt, y_hat_maxt, y_mint, y_hat_mint):

    ## If any missing data in predictions, drop it but give a warning
    y_maxt = np.asarray(y_maxt)
    y_hat_maxt = np.asarray(y_hat_maxt)
    missing_predictions = pd.isna(y_hat_maxt)
    missing_count = missing_predictions.sum()
    if missing_count:
        warnings.warn(
            f"Found {missing_count} NA values in y_hat_maxt; excluding them from metric calculations.",
            UserWarning,
        )
        y_maxt = y_maxt[~missing_predictions]
        y_hat_maxt = y_hat_maxt[~missing_predictions]

    y_mint = np.asarray(y_mint)
    y_hat_mint = np.asarray(y_hat_mint)
    missing_predictions = pd.isna(y_hat_mint)
    missing_count = missing_predictions.sum()
    if missing_count:
        warnings.warn(
            f"Found {missing_count} NA values in y_hat_mint; excluding them from metric calculations.",
            UserWarning,
        )
        y_mint = y_mint[~missing_predictions]
        y_hat_mint = y_hat_mint[~missing_predictions]

    new_row_mint = pd.DataFrame(
        [
            {
                "Element": "Min Temperature",
                "Name": name,
                "RMSE": rmse(y_hat_mint, y_mint),
                "ME": me(y_hat_mint, y_mint),
                "Computation_Time": float("nan"),
                "MAE": mae(y_hat_mint, y_mint),
                "R2": R2(y_hat_mint, y_mint),
                "ME10": ME10(y_hat_mint, y_mint),
            }
        ]
    )
    new_row_maxt = pd.DataFrame(
        [
            {
                "Element": "Max Temperature",
                "Name": name,
                "RMSE": rmse(y_hat_maxt, y_maxt),
                "ME": me(y_hat_maxt, y_maxt),
                "Computation_Time": float("nan"),
                "MAE": mae(y_hat_maxt, y_maxt),
                "R2": R2(y_hat_maxt, y_maxt),
                "ME10": ME10(y_hat_maxt, y_maxt),
            }
        ]
    )
    new_rows = pd.concat([new_row_mint, new_row_maxt], ignore_index=True)

    if metrics_table is None:
        return new_rows

    return pd.concat([metrics_table, new_rows], ignore_index=True)


## Function that prepares data for a given day
def get_daily_temperature_data(
    temperature_data,
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
        daily_data["LTA_maxt"] = daily_data[f"maxt_m_{int(current_month)}"]
        daily_data["LTA_mint"] = daily_data[f"mint_m_{int(current_month)}"]
        daily_data["y_maxt"] = daily_data["maxt"] - daily_data["LTA_maxt"]
        daily_data["y_mint"] = daily_data["mint"] - daily_data["LTA_mint"]
        return daily_data

    if experiment_type == "train_test":
        return {
            "train": prepare_daily_data(temperature_data["train"]),
            "test": prepare_daily_data(temperature_data["test"]),
        }

    all_data = pd.concat(
        [temperature_data["train"], temperature_data["test"]],
        ignore_index=True,
    )
    return prepare_daily_data(all_data)


## Function that takes y_hat (predicted temperature) and adds it to test set
def update_daily_predictions(temperature_data, y_hat_maxt, y_hat_mint, dates, date_index):
    daily_temperature_data = get_daily_temperature_data(
        temperature_data,
        dates,
        date_index,
    )
    daily_test_data = daily_temperature_data["test"]
    predicted_maxt = (
        np.asarray(y_hat_maxt) + daily_test_data["LTA_maxt"].to_numpy()
    )
    predicted_mint = (
        np.asarray(y_hat_mint) + daily_test_data["LTA_mint"].to_numpy()
    )

    current_date = dates.iloc[date_index]
    prediction_mask = (
        (temperature_data["test"]["day"] == current_date["day"])
        & (temperature_data["test"]["month"] == current_date["month"])
        & (temperature_data["test"]["year"] == current_date["year"])
    )
    updated_temperature_data = temperature_data.copy()
    updated_temperature_data["test"] = temperature_data["test"].copy()
    updated_temperature_data["test"].loc[prediction_mask, "predicted_maxt"] = predicted_maxt
    updated_temperature_data["test"].loc[prediction_mask, "predicted_mint"] = predicted_mint

    return updated_temperature_data


def daily_temperature_plot(grid, daily_temperature_data, plot_destination):
    import matplotlib.pyplot as plt

    temperature_palette = [
        "#502020", "#611838", "#721050", "#830868", "#940080",
        "#000080", "#00009F", "#0000BF", "#0000DF", "#0000FF",
        "#566CF2", "#6781F0", "#7997ED", "#8AACEB", "#9BC2E8",
        "#91B7CE", "#87ACB5", "#7DA19B", "#739682", "#698B69",
        "#76AC76", "#7CBC7C", "#83CD83", "#89DD89", "#90EE90",
        "#A8D56A", "#B4C958", "#C1BD45", "#CDB132", "#DAA520",
        "#B1700D", "#C5740A", "#D87806", "#EB7C03", "#FF8000",
        "#FF5000", "#DF3C00", "#BF2800", "#9F1400", "#800000",
        "#993300", "#B26600", "#CC9900", "#E5CC00", "#FFFF00",
    ]

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for response, ax in zip(("maxt", "mint"), axes):
        lower_break = max(-10, int(np.trunc(grid[response].min())) - 1)
        upper_break = min(34, int(np.trunc(grid[response].max())) + 1)
        breaks = np.arange(lower_break, upper_break + 1)
        color_palette = temperature_palette[
            lower_break + 10 : upper_break + 11
        ]

        base_grid_plot(
            grid,
            stations=daily_temperature_data,
            response=response,
            breaks=breaks,
            color_palette=color_palette,
            station_size=15,
            ax=ax,
            include_axis=False,
        )
        fig.axes[-1].tick_params(labelsize=8)

    fig.savefig(plot_destination, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)

# import dill

# dill.dump_session("Data/Daily_Temperature/train_test_80_20_2025.dill")
