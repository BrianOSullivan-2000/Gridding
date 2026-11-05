
# # IDW Regression for Daily Rainfall
#
# Fit a linear trend for each date, then interpolate its residuals using IDW.

# %% 

# ## Setup and data
import dill
import numpy as np
from scipy.spatial import cKDTree

dill.load_session("Data/Daily_Rainfall/train_test_80_20_2016-2025.dill")

session = globals()
rain_data = session["rain_data"]
dates = session["dates"]
metrics_table = session["metrics_table"]
get_daily_rain_data = session["get_daily_rain_data"]
update_daily_predictions = session["update_daily_predictions"]
collect_metrics = session["collect_metrics"]

# %% 

## Regression Model
PREDICTORS = ("east", "north", "points5", "dist2c", "exp25k")
IDW_POWER = 2
IDW_MAX_NEIGHBORS = 15


# %% 

## IDW prediction
def idw_predict(residuals, coordinates, target_coordinates, power=2, nmax=15):
    neighbor_count = min(nmax, len(coordinates))
    distances, neighbors = cKDTree(coordinates).query(
        target_coordinates,
        k=neighbor_count,
    )

    if neighbor_count == 1:
        distances = distances[:, np.newaxis]
        neighbors = neighbors[:, np.newaxis]

    neighbor_residuals = residuals[neighbors]
    predictions = np.empty(len(target_coordinates), dtype=float)
    exact_matches = distances == 0
    has_exact_match = exact_matches.any(axis=1)

    for row_index in np.flatnonzero(has_exact_match):
        predictions[row_index] = neighbor_residuals[
            row_index, exact_matches[row_index]
        ].mean()

    non_exact = ~has_exact_match
    weights = np.zeros_like(distances, dtype=float)
    np.divide(1.0, distances**power, out=weights, where=distances > 0)
    weight_totals = weights[non_exact].sum(axis=1)
    predictions[non_exact] = (
        (weights[non_exact] * neighbor_residuals[non_exact]).sum(axis=1)
        / weight_totals
    )

    return predictions


def predict_daily_values(daily_data):
    train_data = daily_data["train"].dropna(subset=["y", *PREDICTORS])
    test_data = daily_data["test"]
    predictions = np.full(len(test_data), np.nan, dtype=float)
    valid_test = test_data.loc[:, PREDICTORS].notna().all(axis=1).to_numpy()

    if train_data.empty or not valid_test.any():
        return predictions

    training_features = train_data.loc[:, PREDICTORS].to_numpy(dtype=float)
    training_design = np.column_stack(
        (np.ones(len(train_data)), training_features)
    )
    coefficients, _, _, _ = np.linalg.lstsq(
        training_design,
        train_data["y"].to_numpy(dtype=float),
        rcond=None,
    )
    training_residuals = (
        train_data["y"].to_numpy(dtype=float) - training_design @ coefficients
    )
    training_coordinates = train_data.loc[:, ["east", "north"]].to_numpy()
    test_coordinates = test_data.loc[valid_test, ["east", "north"]].to_numpy()
    interpolated_residuals = idw_predict(
        training_residuals,
        training_coordinates,
        test_coordinates,
        power=IDW_POWER,
        nmax=IDW_MAX_NEIGHBORS,
    )

    test_features = test_data.loc[valid_test, PREDICTORS].to_numpy(dtype=float)
    test_design = np.column_stack((np.ones(len(test_features)), test_features))
    regression_predictions = test_design @ coefficients
    predictions[valid_test] = regression_predictions + interpolated_residuals
    return predictions

# %% Run all dates

# Reinitialize from setup data so any one-date experimentation above is not
# applied twice.
for date_index in range(len(dates)):
    daily_data = get_daily_rain_data(rain_data, dates, date_index)
    predicted_y = predict_daily_values(daily_data)
    rain_data = update_daily_predictions(
        rain_data,
        predicted_y,
        dates,
        date_index,
    )
    print(date_index)

# %% Metrics
metrics_table = collect_metrics(
    None,
    "IDW MLR",
    rain_data["test"]["rain"],
    rain_data["test"]["predicted_rain"],
)
print(metrics_table.to_string(index=False))

# metrics_table.to_csv(
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_MLR_python.csv",
#     index=False,
# )

# %%
