
## Validation tests for daily rainfall grids ##

## Spatially embedded XGBoost (seXGBoost, need a better name) ##
## Predicts Daily Rainfall using a standard XGBoost model     ##
## except spatial basis functions are also included to add    ##
## a spatial effect                                           ##

# %%

## Setup and data

## Setup and data
from itertools import product
from pathlib import Path

import dill
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.spatial import cKDTree

from utilities.add_basis_functions import add_basis_functions
from utilities.make_mesh import make_mesh

dill.load_session("Data/Daily_Rainfall/train_test_80_20_provisional.dill")

session = globals()
rain_data = session["rain_data"]
dates = session["dates"]
island_outline = session["island_outline"]
get_daily_rain_data = session["get_daily_rain_data"]
update_daily_predictions = session["update_daily_predictions"]
collect_metrics = session["collect_metrics"]

# %%

## Spatial basis functions
nodelist = [
    make_mesh(
        island_outline, a = 0.1, q = 30, 
        buffer = 0, simplify = 10000, plot = False
    ),
    make_mesh(
        island_outline, a = 0.04, q = 30,
        buffer = 500, simplify = 10000, plot = False
    ),
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
        supportlist=[12, 12],
        bfs_type="Mat",
        parlist=[
            (1, 15000, 0, 2),
            (1, 7000, 0, 2),
        ],
    ),
    columns=basis_columns
)

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
    "eta": 0.1, 
    "objective": "reg:squarederror",
    "subsample": 0.2,
    "colsample_bytree": 0.2,
    "lambda": 5,
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

def predict_daily_values(daily_data):
    return predict_values(daily_data["train"], daily_data["test"])


# %%

## Actual gridding, first need to set up basis functions
prediction_grid = session["grid_geodata"].copy()
grid_LTAs_9120 = session["grid_LTAs_9120"]
daily_rain_plot = session["daily_rain_plot"]

grid_basis = pd.DataFrame(
    add_basis_functions(
        prediction_grid[["east", "north"]].to_numpy(),
        nodelist,
        supportlist=[12, 12],
        bfs_type="Mat",
        parlist=[
            (1, 15000, 0, 2),
            (1, 7000, 0, 2),
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

## First-date hyperparameter grid search

hyp_figure_directory = Path("Figures/Daily_Rainfall_provisional/spatial_XGBoost_hyp_2")
hyp_figure_directory.mkdir(parents=True, exist_ok=True)

first_date = dates.iloc[0]
first_date_index = 0
first_date_data = get_daily_rain_data(
    rain_data,
    dates,
    first_date_index,
    experiment_type="all_data",
)

station_basis_by_station = station_basis.copy()
station_basis_by_station.index = station_locations["stno"].to_numpy()
first_date_data = first_date_data.join(station_basis_by_station, on="stno")
first_date_data = first_date_data.dropna(subset=["y"]).reset_index(drop=True)

split_rng = np.random.default_rng(42)
shuffled_indices = split_rng.permutation(len(first_date_data))
validation_size = max(1, round(len(first_date_data) * 0.15))
validation_size = min(validation_size, len(first_date_data) - 1)
validation_indices = shuffled_indices[:validation_size]
training_indices = shuffled_indices[validation_size:]
hyperparameter_training_data = first_date_data.iloc[training_indices]
hyperparameter_validation_data = first_date_data.iloc[validation_indices]

training_matrix = xgb.DMatrix(
    hyperparameter_training_data[model_features],
    label=hyperparameter_training_data["y"].to_numpy(dtype=float),
)
validation_matrix = xgb.DMatrix(
    hyperparameter_validation_data[model_features]
)
validation_lta = hyperparameter_validation_data["LTA"].to_numpy(dtype=float)
validation_observed_rain = hyperparameter_validation_data["rain"].to_numpy(
    dtype=float
)

hyperparameter_axes = {
    "lambda": [1, 5, 10, 30, 50, 70],
    "max_depth": [3, 4, 5],
    "eta": [0.05, 0.1, 0.2],
    "min_child_weight": [1],
    "subsample": [0.2, 0.3, 0.4],
    "colsample_bytree": [0.1, 0.2],
    "num_boost_round": [800, 1600, 2000],
}
hyperparameter_names = list(hyperparameter_axes)
hyperparameter_grid = [
    dict(zip(hyperparameter_names, values))
    for values in product(*hyperparameter_axes.values())
]


first_date_label = (
    f"{int(first_date['year'])}_{int(first_date['month']):02d}_"
    f"{int(first_date['day']):02d}"
)
print(
    f"Hyperparameter search for {first_date_label}: "
    f"{len(hyperparameter_grid)} combinations, "
    f"{len(training_indices)} training stations, "
    f"{len(validation_indices)} validation stations.",
    flush=True,
)

hyperparameter_results = []
for combination_index, combination in enumerate(hyperparameter_grid, start=1):
    model_params = {
        "objective": "reg:squarederror",
        "max_depth": combination["max_depth"],
        "eta": combination["eta"],
        "min_child_weight": combination["min_child_weight"],
        "subsample": combination["subsample"],
        "colsample_bytree": combination["colsample_bytree"],
        "lambda": combination["lambda"],
        "seed": 42,
    }
    model = xgb.train(
        model_params,
        training_matrix,
        num_boost_round=combination["num_boost_round"],
    )

    predicted_validation_rain = (
        np.expm1(model.predict(validation_matrix)) * validation_lta
    )
    valid_validation_predictions = np.isfinite(predicted_validation_rain) & np.isfinite(
        validation_observed_rain
    )
    validation_rmse = np.sqrt(
        np.mean(
            (
                predicted_validation_rain[valid_validation_predictions]
                - validation_observed_rain[valid_validation_predictions]
            )
            ** 2
        )
    )

    print(
        f"{combination_index}/{len(hyperparameter_grid)} "
        f"{combination} validation RMSE={validation_rmse:.4f}",
        flush=True,
    )

    hyperparameter_results.append(
        {
            **combination,
            "validation_rmse": validation_rmse,
        }
    )

hyperparameter_results = pd.DataFrame(hyperparameter_results).sort_values(
    "validation_rmse"
).reset_index(drop=True)
results_path = hyp_figure_directory / f"hyperparameter_results_{first_date_label}.csv"
hyperparameter_results.to_csv(results_path, index=False)
print(f"Saved sorted validation results to {results_path}", flush=True)
print("Best hyperparameter combinations:")
print(hyperparameter_results.head(10).to_string(index=False))

# %%

hyp2_figure_directory = Path("Figures/Daily_Rainfall_provisional/spatial_XGBoost_hyp")
hyp2_figure_directory.mkdir(parents=True, exist_ok=True)

first_day_results = pd.read_csv(results_path)
selected_hyperparameters = first_day_results.loc[
    first_day_results["validation_rmse"] < 0.8
].copy()

validation_date_count = 20
validation_sets_by_date = []
print(
    f"Preparing 15% validation sets for the first {validation_date_count} dates "
    f"and {len(selected_hyperparameters)} selected parameter combinations.",
    flush=True,
)

for date_index in range(validation_date_count):
    current_date = dates.iloc[date_index]
    current_date_label = (
        f"{int(current_date['year'])}_{int(current_date['month']):02d}_"
        f"{int(current_date['day']):02d}"
    )
    daily_data = get_daily_rain_data(
        rain_data,
        dates,
        date_index,
        experiment_type="all_data",
    )
    daily_data = daily_data.join(station_basis_by_station, on="stno")
    daily_data = daily_data.dropna(subset=["y"]).reset_index(drop=True)

    date_rng = np.random.default_rng(42 + date_index)
    date_indices = date_rng.permutation(len(daily_data))
    date_validation_size = max(1, round(len(daily_data) * 0.15))
    date_validation_size = min(date_validation_size, len(daily_data) - 1)
    date_validation_indices = date_indices[:date_validation_size]
    date_training_indices = date_indices[date_validation_size:]
    training_data = daily_data.iloc[date_training_indices]
    validation_data = daily_data.iloc[date_validation_indices]

    validation_sets_by_date.append(
        {
            "date_label": current_date_label,
            "daily_data": daily_data,
            "training_matrix": xgb.DMatrix(
                training_data[model_features],
                label=training_data["y"].to_numpy(dtype=float),
            ),
            "validation_matrix": xgb.DMatrix(validation_data[model_features]),
            "validation_lta": validation_data["LTA"].to_numpy(dtype=float),
            "validation_rain": validation_data["rain"].to_numpy(dtype=float),
            "grid_lta": grid_LTAs_9120[
                f"m_{int(current_date['month'])}"
            ].to_numpy(dtype=float),
        }
    )
    print(
        f"  {current_date_label}: {len(date_training_indices)} training, "
        f"{len(date_validation_indices)} validation stations",
        flush=True,
    )

hyp2_run_id = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
daily_rmse_path = hyp2_figure_directory / (
    f"validation_rmse_by_day_{hyp2_run_id}.csv"
)
mean_rmse_path = hyp2_figure_directory / (
    f"validation_rmse_mean_{hyp2_run_id}.csv"
)
parameter_columns = [
    "max_depth",
    "eta",
    "min_child_weight",
    "subsample",
    "colsample_bytree",
    "lambda",
    "num_boost_round",
]

selected_combinations = []
for combination_index, (_, selected_row) in enumerate(
    selected_hyperparameters.iterrows(), start=1
):
    combination = {
        name: selected_row[name]
        for name in parameter_columns
    }
    combination["max_depth"] = int(combination["max_depth"])
    combination["num_boost_round"] = int(combination["num_boost_round"])
    first_day_rmse = float(selected_row["validation_rmse"])
    model_params = {
        "objective": "reg:squarederror",
        **combination,
        "seed": 42,
    }
    selected_combinations.append(
        {
            "combination_index": combination_index,
            "combination_id": combination_index,
            "combination": combination,
            "first_day_rmse": first_day_rmse,
            "model_params": model_params,
            "daily_rmse": [],
        }
    )

print(
    f"Starting date-major validation for {len(validation_sets_by_date)} dates "
    f"and {len(selected_combinations)} combinations.",
    flush=True,
)
for date_index, date_validation in enumerate(validation_sets_by_date, start=1):
    date_rows = []
    print(
        f"Date {date_index}/{len(validation_sets_by_date)} "
        f"{date_validation['date_label']}: fitting all combinations.",
        flush=True,
    )
    for selected_combination in selected_combinations:
        combination = selected_combination["combination"]
        model = xgb.train(
            selected_combination["model_params"],
            date_validation["training_matrix"],
            num_boost_round=combination["num_boost_round"],
        )
        predicted_validation_rain = (
            np.expm1(model.predict(date_validation["validation_matrix"]))
            * date_validation["validation_lta"]
        )
        observed_validation_rain = date_validation["validation_rain"]
        valid_predictions = np.isfinite(predicted_validation_rain) & np.isfinite(
            observed_validation_rain
        )
        if valid_predictions.any():
            daily_rmse = np.sqrt(
                np.mean(
                    (
                        predicted_validation_rain[valid_predictions]
                        - observed_validation_rain[valid_predictions]
                    )
                    ** 2
                )
            )
        else:
            daily_rmse = float("nan")
        selected_combination["daily_rmse"].append(daily_rmse)

        date_rows.append(
            {
                "combination_id": selected_combination["combination_id"],
                "date": date_validation["date_label"],
                **combination,
                "first_day_validation_rmse": selected_combination[
                    "first_day_rmse"
                ],
                "validation_rmse": daily_rmse,
                "validation_station_count": int(valid_predictions.sum()),
            }
        )
        print(
            f"  Combination {selected_combination['combination_index']}/"
            f"{len(selected_combinations)}: RMSE={daily_rmse:.4f}",
            flush=True,
        )

    pd.DataFrame(date_rows).to_csv(
        daily_rmse_path,
        mode="a",
        header=not daily_rmse_path.exists(),
        index=False,
    )

combination_summaries = []
for selected_combination in selected_combinations:
    combination_rmse = selected_combination["daily_rmse"]
    mean_rmse = float(np.nanmean(combination_rmse))
    combination = selected_combination["combination"]
    combination_summaries.append(
        {
            "combination_id": selected_combination["combination_id"],
            **combination,
            "first_day_validation_rmse": selected_combination["first_day_rmse"],
            "mean_validation_rmse": mean_rmse,
            "validation_day_count": len(combination_rmse),
        }
    )
    print(
        f"Combination {selected_combination['combination_index']}/"
        f"{len(selected_combinations)} mean validation RMSE across "
        f"{len(combination_rmse)} dates: {mean_rmse:.4f}",
        flush=True,
    )

pd.DataFrame(combination_summaries).to_csv(mean_rmse_path, index=False)

print(f"Saved per-day RMSEs to {daily_rmse_path}", flush=True)
print(f"Saved mean RMSEs to {mean_rmse_path}", flush=True)

# %%

station_basis_by_station = station_basis.copy()
station_basis_by_station.index = station_locations["stno"].to_numpy()

## Plot grids for the six best configurations across the first ten dates
validation_rmse_mean_file = Path(
    r"Figures/Daily_Rainfall_provisional/spatial_XGBoost_hyp/validation_rmse_mean_20261112_220622.csv"
)
figure_output_directory = Path(
    r"Figures/Daily_Rainfall_provisional/spatial_XGBoost_hyp"
)

if not validation_rmse_mean_file.is_file():
    raise FileNotFoundError(
        f"Set validation_rmse_mean_file to an existing CSV: "
        f"{validation_rmse_mean_file}"
    )

top_configurations = pd.read_csv(validation_rmse_mean_file)
top_configurations = top_configurations.loc[
    np.isfinite(top_configurations["mean_validation_rmse"])
].nsmallest(6, "mean_validation_rmse")
if top_configurations.empty:
    raise ValueError("No configurations have a finite mean validation RMSE")

figure_output_directory.mkdir(parents=True, exist_ok=True)
prediction_grid_matrix = xgb.DMatrix(prediction_grid[model_features])
grid_date_count = min(10, len(dates))

print(
    f"Plotting {len(top_configurations)} configurations across "
    f"{grid_date_count} dates to {figure_output_directory}",
    flush=True,
)

for date_index in range(grid_date_count):
    current_date = dates.iloc[date_index]
    current_day = int(current_date["day"])
    current_month = int(current_date["month"])
    current_year = int(current_date["year"])
    current_date_label = (
        f"{current_year}_{current_month:02d}_{current_day:02d}"
    )

    daily_rain_data = get_daily_rain_data(
        rain_data,
        dates,
        date_index,
        experiment_type="all_data",
    )
    daily_rain_data = daily_rain_data.join(station_basis_by_station, on="stno")
    daily_rain_data = daily_rain_data.dropna(subset=["y"])
    training_matrix = xgb.DMatrix(
        daily_rain_data[model_features],
        label=daily_rain_data["y"].to_numpy(dtype=float),
    )

    for configuration_rank, (_, configuration) in enumerate(
        top_configurations.iterrows(),
        start=1,
    ):
        model_params = {
            "objective": "reg:squarederror",
            "max_depth": int(configuration["max_depth"]),
            "eta": float(configuration["eta"]),
            "min_child_weight": float(configuration["min_child_weight"]),
            "subsample": float(configuration["subsample"]),
            "colsample_bytree": float(configuration["colsample_bytree"]),
            "lambda": float(configuration["lambda"]),
            "seed": 42,
        }
        model = xgb.train(
            model_params,
            training_matrix,
            num_boost_round=int(configuration["num_boost_round"]),
        )
        predicted_y = model.predict(prediction_grid_matrix)
        prediction_grid["rain"] = (
            np.expm1(predicted_y)
            * grid_LTAs_9120[f"m_{current_month}"].to_numpy()
        )

        configuration_id = int(configuration["combination_id"])
        plot_destination = figure_output_directory / (
            f"spatial_XGBoost_top{configuration_rank:02d}_"
            f"config{configuration_id}_{current_date_label}.jpg"
        )
        daily_rain_plot(
            prediction_grid,
            daily_rain_data,
            plot_destination=str(plot_destination),
        )
        print(
            f"Date {current_date_label}: plotted configuration "
            f"{configuration_rank}/{len(top_configurations)} "
            f"(ID {configuration_id})",
            flush=True,
        )
