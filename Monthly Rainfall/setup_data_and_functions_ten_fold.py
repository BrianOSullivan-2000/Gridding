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

import dill
import numpy as np
import pandas as pd

from utilities.metric_functions import (
    ME10,
    R2,
    mae,
    me,
    rmse,
)

# Restore shared objects created by setup_data_and_functions.py.
dill.load_session("Data/Monthly_Rainfall/train_test_80_20_2016-2025.dill")
session_objects = globals()
station_geodata = session_objects["station_geodata"]
station_LTAs_6190 = session_objects["station_LTAs_6190"]
station_LTAs_9120 = session_objects["station_LTAs_9120"]
grid_geodata = session_objects["grid_geodata"]
grid_LTAs_6190 = session_objects["grid_LTAs_6190"]
grid_LTAs_9120 = session_objects["grid_LTAs_9120"]
island_outline = session_objects["island_outline"]

## Ten-fold cross-validation (2016-2025)
rain_data = (
    pd.read_csv("Data/Monthly_Rainfall/folds/10fold_2016-2025.csv")
    .merge(station_geodata, on="stno", how="left")
    .merge(station_LTAs_9120, on="stno", how="left")
    .assign(predicted_rain=float("nan"))
)

dates = (
    rain_data[["year", "month"]]
    .drop_duplicates()
    .sort_values(["year", "month"])
    .reset_index(drop=True)
)

## Template for metrics
metrics_table = pd.DataFrame(
    columns=["Name", "RMSE", "ME", "MAE", "R2", "ME10"]
)


## Function that prepares data for a given month and held-out fold
def get_monthly_rain_data(
    rain_data,
    dates,
    date_index=None,
    current_fold=None,
    experiment_type="train_test",
    date_indices=None,
):
    if experiment_type not in ("train_test", "all_data"):
        raise ValueError("experiment_type must be 'train_test' or 'all_data'")

    if date_indices is None:
        date_indices = date_index
    if date_indices is None:
        raise ValueError("Provide one or more date_indices")

    if np.isscalar(date_indices):
        date_indices = [date_indices]
    date_indices = list(date_indices)
    if not date_indices:
        raise ValueError("Provide one or more date_indices")

    selected_dates = dates.iloc[date_indices][["year", "month"]].drop_duplicates()
    selected_months = rain_data.merge(
        selected_dates,
        on=["year", "month"],
        how="inner",
        sort=False,
    ).copy()
    selected_months["LTA"] = selected_months.apply(
        lambda row: row[f"m_{int(row['month'])}"],
        axis=1,
    )
    selected_months["y"] = selected_months["rain"] / selected_months["LTA"]

    if experiment_type == "all_data":
        return selected_months

    if current_fold is None:
        raise ValueError("Need to give fold for train-test split")

    return {
        "train": selected_months.loc[
            selected_months["fold"] != current_fold
        ].copy(),
        "test": selected_months.loc[
            selected_months["fold"] == current_fold
        ].copy(),
    }


## Function that adds predictions to the held-out fold
def update_monthly_predictions(
    rain_data,
    y_hat,
    dates,
    date_index=None,
    current_fold=None,
    date_indices=None,
):
    if current_fold is None:
        raise ValueError("Need to give fold for train-test split")

    if date_indices is None:
        date_indices = date_index

    monthly_rain_data = get_monthly_rain_data(
        rain_data,
        dates,
        current_fold=current_fold,
        date_indices=date_indices,
    )
    test_data = monthly_rain_data["test"]
    predicted_rain = np.asarray(y_hat, dtype=float) * test_data["LTA"].to_numpy()
    if len(predicted_rain) != len(test_data):
        raise ValueError(
            f"Expected {len(test_data)} predictions, got {len(predicted_rain)}"
        )

    selected_dates = dates.iloc[
        [date_indices] if np.isscalar(date_indices) else list(date_indices)
    ][["year", "month"]].drop_duplicates()
    target_rows = rain_data.merge(
        selected_dates,
        on=["year", "month"],
        how="left",
        indicator=True,
        sort=False,
        )["_merge"].eq("both").to_numpy(copy=True)
    target_rows = target_rows & (rain_data["fold"].to_numpy() == current_fold)

    updated_rain_data = rain_data.copy()
    updated_rain_data.loc[target_rows, "predicted_rain"] = predicted_rain
    return updated_rain_data


## Function that calculates the mean metric across folds
def collect_metrics(metrics_table, name, rain_data):
    fold_RMSEs = []
    fold_MEs = []
    fold_MAEs = []
    fold_R2s = []
    fold_ME10s = []

    for current_fold in rain_data["fold"].drop_duplicates():
        fold_data = rain_data.loc[rain_data["fold"] == current_fold]
        y = fold_data["rain"].to_numpy(dtype=float)
        y_hat = fold_data["predicted_rain"].to_numpy(dtype=float)
        missing_predictions = pd.isna(y_hat)
        missing_count = missing_predictions.sum()
        if missing_count:
            warnings.warn(
                f"Found {missing_count} missing predictions in fold "
                f"{current_fold}; excluding them from metric calculations.",
                UserWarning,
            )
            y = y[~missing_predictions]
            y_hat = y_hat[~missing_predictions]

        fold_RMSEs.append(rmse(y_hat, y))
        fold_MEs.append(me(y_hat, y))
        fold_MAEs.append(mae(y_hat, y))
        fold_R2s.append(R2(y_hat, y))
        fold_ME10s.extend(
            float(error)
            for error in ME10(y_hat, y).split(", ")
            if error != "NA"
        )

    new_row = pd.DataFrame(
        [
            {
                "Name": name,
                "RMSE": np.mean(fold_RMSEs),
                "ME": np.mean(fold_MEs),
                "MAE": np.mean(fold_MAEs),
                "R2": np.mean(fold_R2s),
                "ME10": ", ".join(
                    f"{error:g}"
                    for error in sorted(fold_ME10s, reverse=True)[:10]
                ),
            }
        ]
    )

    if metrics_table is None:
        return new_row

    return pd.concat([metrics_table, new_row], ignore_index=True)

dill.dump_session("Data/Monthly_Rainfall/train_test_ten_fold_2016-2025.dill")