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

from utilities.metric_functions import ME10, R2, mae, me, rmse

# Restore shared objects created by setup_data_and_functions.py.
dill.load_session("Data/Daily_Temperature/train_test_80_20_2025.dill")
session_objects = globals()
station_geodata = session_objects["station_geodata"]
station_LTAs_max_9120 = session_objects["station_LTAs_max_9120"]
station_LTAs_min_9120 = session_objects["station_LTAs_min_9120"]
grid_geodata = session_objects["grid_geodata"]
LTAs_max_9120 = session_objects["LTAs_max_9120"]
LTAs_min_9120 = session_objects["LTAs_min_9120"]
island_outline = session_objects["island_outline"]

## Setup for the five-fold data
temperature_data = (
    pd.read_csv("Data/Daily_Temperature/folds/ST_5fold_2025.csv")
    .merge(station_geodata, on="stno", how="left")
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
    temperature_data[["year", "month", "day"]]
    .drop_duplicates()
    .sort_values(["year", "month", "day"])
    .reset_index(drop=True)
)

## Template for metrics
metrics_table = pd.DataFrame(
    columns=["Element", "Name", "RMSE", "ME", "MAE", "R2", "ME10"]
)


def get_daily_temperature_data(
    temperature_data,
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

    if any(index < 0 or index >= len(dates) for index in date_indices):
        raise ValueError("date_indices contains an invalid index")

    selected_dates = dates.iloc[date_indices][["year", "month", "day"]].drop_duplicates()
    daily_temperature_data = temperature_data.merge(
        selected_dates,
        on=["year", "month", "day"],
        how="inner",
        sort=False,
    ).copy()

    daily_temperature_data["LTA_maxt"] = daily_temperature_data.apply(
        lambda row: row[f"maxt_m_{int(row['month'])}"],
        axis=1,
    )
    daily_temperature_data["LTA_mint"] = daily_temperature_data.apply(
        lambda row: row[f"mint_m_{int(row['month'])}"],
        axis=1,
    )
    daily_temperature_data["y_maxt"] = (
        daily_temperature_data["maxt"] - daily_temperature_data["LTA_maxt"]
    )
    daily_temperature_data["y_mint"] = (
        daily_temperature_data["mint"] - daily_temperature_data["LTA_mint"]
    )

    if experiment_type == "all_data":
        return daily_temperature_data

    if current_fold is None:
        raise ValueError("Need to give fold for train-test split")

    return {
        "train": daily_temperature_data.loc[
            daily_temperature_data["fold"] != current_fold
        ].copy(),
        "test": daily_temperature_data.loc[
            daily_temperature_data["fold"] == current_fold
        ].copy(),
    }


def get_temperature_data(
    temperature_data,
    current_fold=None,
    experiment_type="train_test",
    date_index=None,
    dates=None,
    date_indices=None,
):
    if experiment_type not in ("train_test", "all_data"):
        raise ValueError("experiment_type must be 'train_test' or 'all_data'")

    if date_indices is None:
        date_indices = date_index

    prepared_temperature_data = temperature_data.copy()
    prepared_temperature_data["LTA_maxt"] = prepared_temperature_data.apply(
        lambda row: row[f"maxt_m_{int(row['month'])}"],
        axis=1,
    )
    prepared_temperature_data["LTA_mint"] = prepared_temperature_data.apply(
        lambda row: row[f"mint_m_{int(row['month'])}"],
        axis=1,
    )
    prepared_temperature_data["y_maxt"] = (
        prepared_temperature_data["maxt"] - prepared_temperature_data["LTA_maxt"]
    )
    prepared_temperature_data["y_mint"] = (
        prepared_temperature_data["mint"] - prepared_temperature_data["LTA_mint"]
    )

    selected_dates = None
    if date_indices is not None:
        if dates is None:
            raise ValueError("Provide dates when using date_indices")
        if np.isscalar(date_indices):
            date_indices = [date_indices]
        date_indices = list(date_indices)
        if any(index < 0 or index >= len(dates) for index in date_indices):
            raise ValueError("date_indices contains an invalid index")
        selected_dates = dates.iloc[date_indices][["year", "month", "day"]].drop_duplicates()

    if experiment_type == "train_test":
        if current_fold is None:
            raise ValueError("Need to give fold for train-test split")

        prepared_temperature_data = {
            "train": prepared_temperature_data.loc[
                prepared_temperature_data["fold"] != current_fold
            ].copy(),
            "test": prepared_temperature_data.loc[
                prepared_temperature_data["fold"] == current_fold
            ].copy(),
        }

        if selected_dates is not None:
            for split in ("train", "test"):
                prepared_temperature_data[split] = prepared_temperature_data[split].merge(
                    selected_dates,
                    on=["year", "month", "day"],
                    how="inner",
                    sort=False,
                )
        return prepared_temperature_data

    if selected_dates is not None:
        prepared_temperature_data = prepared_temperature_data.merge(
            selected_dates,
            on=["year", "month", "day"],
            how="inner",
            sort=False,
        )
    return prepared_temperature_data


def update_daily_predictions(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    dates,
    date_index=None,
    current_fold=None,
    date_indices=None,
):
    if current_fold is None:
        raise ValueError("Need to give fold for train-test split")

    if date_indices is None:
        date_indices = date_index
    if date_indices is None:
        date_indices = list(range(len(dates)))

    if np.isscalar(date_indices):
        date_indices = [date_indices]

    daily_temperature_data = get_daily_temperature_data(
        temperature_data,
        dates,
        current_fold=current_fold,
        date_indices=date_indices,
    )
    test_data = daily_temperature_data["test"]
    predicted_maxt = (
        np.asarray(y_hat_maxt, dtype=float) + test_data["LTA_maxt"].to_numpy()
    )
    predicted_mint = (
        np.asarray(y_hat_mint, dtype=float) + test_data["LTA_mint"].to_numpy()
    )
    if len(predicted_maxt) != len(test_data) or len(predicted_mint) != len(test_data):
        raise ValueError("Prediction count does not match the selected test rows")

    selected_dates = dates.iloc[list(date_indices)][["year", "month", "day"]].drop_duplicates()
    target_rows = (
        temperature_data.merge(
            selected_dates,
            on=["year", "month", "day"],
            how="left",
            indicator=True,
            sort=False,
        )["_merge"]
        .eq("both")
        .to_numpy(copy=True)
    )
    target_rows &= temperature_data["fold"].to_numpy() == current_fold

    updated_temperature_data = temperature_data.copy()
    updated_temperature_data.loc[target_rows, "predicted_maxt"] = predicted_maxt
    updated_temperature_data.loc[target_rows, "predicted_mint"] = predicted_mint
    return updated_temperature_data


def update_predictions(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    current_fold,
    dates=None,
    date_indices=None,
):
    prepared_temperature_data = get_temperature_data(
        temperature_data,
        current_fold=current_fold,
        dates=dates,
        date_indices=date_indices,
    )
    test_data = prepared_temperature_data["test"]

    predicted_maxt = (
        np.asarray(y_hat_maxt, dtype=float) + test_data["LTA_maxt"].to_numpy()
    )
    predicted_mint = (
        np.asarray(y_hat_mint, dtype=float) + test_data["LTA_mint"].to_numpy()
    )
    if len(predicted_maxt) != len(test_data) or len(predicted_mint) != len(test_data):
        raise ValueError("Prediction count does not match the selected test rows")

    if date_indices is None:
        target_data = temperature_data.loc[
            temperature_data["fold"] == current_fold
        ].copy()
    else:
        if np.isscalar(date_indices):
            date_indices = [date_indices]
        selected_dates = dates.iloc[list(date_indices)][["year", "month", "day"]].drop_duplicates()
        target_data = temperature_data.loc[
            temperature_data["fold"] == current_fold
        ].merge(
            selected_dates,
            on=["year", "month", "day"],
            how="inner",
            sort=False,
        )

    if len(target_data) != len(test_data):
        raise ValueError("Prediction count does not match the selected test rows")

    if date_indices is None:
        target_mask = temperature_data["fold"].to_numpy() == current_fold
    else:
        target_mask = (
            temperature_data.merge(
                selected_dates,
                on=["year", "month", "day"],
                how="left",
                indicator=True,
                sort=False,
            )["_merge"]
            .eq("both")
            .to_numpy(copy=True)
        )
        target_mask &= temperature_data["fold"].to_numpy() == current_fold

    updated_temperature_data = temperature_data.copy()
    updated_temperature_data.loc[target_mask, "predicted_maxt"] = predicted_maxt
    updated_temperature_data.loc[target_mask, "predicted_mint"] = predicted_mint
    return updated_temperature_data


def collect_metrics(metrics_table, name, temperature_data):
    rows = []
    for element_name, element in (
        ("Min Temperature", "mint"),
        ("Max Temperature", "maxt"),
    ):
        fold_RMSEs = []
        fold_MEs = []
        fold_MAEs = []
        fold_R2s = []
        fold_ME10s = []

        for current_fold in temperature_data["fold"].dropna().drop_duplicates():
            fold_data = temperature_data.loc[
                temperature_data["fold"] == current_fold
            ]
            y = fold_data[element].to_numpy(dtype=float)
            y_hat = fold_data[f"predicted_{element}"].to_numpy(dtype=float)
            missing_predictions = pd.isna(y_hat)
            missing_count = missing_predictions.sum()
            if missing_count:
                warnings.warn(
                    f"Found {missing_count} missing predictions for {element_name} in fold "
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

        rows.append(
            {
                "Element": element_name,
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
        )

    new_rows = pd.DataFrame(rows)
    if metrics_table is None:
        return new_rows
    return pd.concat([metrics_table, new_rows], ignore_index=True)

dill.dump_session("Data/Daily_Temperature/train_test_five_fold_2025.dill")
