
## Validation tests for daily temperature grids ##

## Spatially embedded XGBoost (seXGBoost, need a better name) ##
## Predicts Daily Temperature using a standard XGBoost model  ##
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
from utilities.metric_functions import rmse

# Load in general objects (for CV experiment)
dill.load_session("Data/Daily_Temperature/train_test_80_20_2025.dill")
session = globals()
temperature_data = session["temperature_data"]
dates = session["dates"]
island_outline = session["island_outline"]
get_daily_temperature_data = session["get_daily_temperature_data"]
update_daily_predictions = session["update_daily_predictions"]
collect_metrics = session["collect_metrics"]

# %%

## Spatial basis functions

## Make mesh (node locations)
nodelist = [
    make_mesh(island_outline, a=0.15, q=30, buffer=0, simplify=10000, plot=False),
    make_mesh(island_outline, a=0.02, q=30, buffer=2000, simplify=10000, plot=False),
]

## Remove nodes that are too close to each other
threshold = 4000
for index in range(1, len(nodelist)):
    previous_nodes = np.vstack(nodelist[:index])
    distances, _ = cKDTree(previous_nodes).query(nodelist[index], k=1)
    nodelist[index] = nodelist[index][distances > threshold]

## Calculate basis functions between stations and nodes
station_locations = (
    pd.concat([temperature_data["train"], temperature_data["test"]], ignore_index=True)
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

## Prepare some predictors
for split in ("train", "test"):
    split_data = temperature_data[split].assign(
        eastone=temperature_data[split]["east"] / 1000,
        northone=temperature_data[split]["north"] / 1000,
        aspectsin=np.sin(temperature_data[split]["aspect_10km"]),
        aspectcos=np.cos(temperature_data[split]["aspect_10km"]),
    )
    temperature_data[split] = split_data

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

## Hyperparameters (best defaults)
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
def predict_values(train_data, target_data, response):
    train_data = train_data.dropna(subset=[response])

    ## Fit model
    Xtrain = xgb.DMatrix(
        train_data[model_features],
        label=train_data[response].to_numpy(dtype=float),
    )
    model = xgb.train(
        xgb_params,
        Xtrain,
        num_boost_round=num_boost_round,
    )

    ## Predict values
    Xtest = xgb.DMatrix(target_data[model_features])
    return model.predict(Xtest)

## Predict values for maxt and mint
def predict_daily_values(daily_data):
    predicted_y_maxt = predict_values(
        daily_data["train"], daily_data["test"], "y_maxt"
    )
    predicted_y_mint = predict_values(
        daily_data["train"], daily_data["test"], "y_mint"
    )
    return predicted_y_maxt, predicted_y_mint


# %%

## Actual gridding, first need to set up basis functions
prediction_grid = session["grid_geodata"].copy()
grid_coordinates = prediction_grid[["east", "north"]]
grid_LTAs_max_9120 = grid_coordinates.merge(
    session["LTAs_max_9120"], on=["east", "north"], how="left"
)
grid_LTAs_min_9120 = grid_coordinates.merge(
    session["LTAs_min_9120"], on=["east", "north"], how="left"
)
daily_temperature_plot = session["daily_temperature_plot"]

## Get basis functions
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

## All predictors
prediction_grid = pd.concat((prediction_grid, grid_basis), axis=1)

# %%

## Hyperparameter grid search

## Make a directory to save results
hyp_figure_directory = Path("Figures/Daily_Temperature/spatial_XGBoost_hyp")
hyp_figure_directory.mkdir(parents=True, exist_ok=True)

## Hyperparameters to search over
hyperparameter_axes = {
    "max_depth": [2, 3, 4, 6],
    "eta": [0.02, 0.05, 0.1],
    "min_child_weight": [1],
    "subsample": [0.1, 0.2, 0.4],
    "colsample_bytree": [0.1, 0.2, 0.3, 0.5],
    "lambda": [1, 10, 30, 50],
    "num_boost_round": [1600],
}
hyperparameter_names = list(hyperparameter_axes)
hyperparameter_grid = [
    dict(zip(hyperparameter_names, values))
    for values in product(*hyperparameter_axes.values())
]

## We'll look at 5 days
NO_DATES_TO_CHECK = 5
date_rng = np.random.default_rng(42)
date_indices = np.sort(
    date_rng.choice(len(dates), size=NO_DATES_TO_CHECK, replace=False)
)

## Write the daily data to here (as train and validation splits)
daily_data_by_date = {}
print(
    f"Setting up validation tests for {len(date_indices)} days.",
    flush=True,
)

## For each date, add basis functions and do 75/25 split
for date_index in date_indices:
    current_date = dates.iloc[date_index]
    date_label = (
        f"{int(current_date['year'])}_{int(current_date['month']):02d}_"
        f"{int(current_date['day']):02d}"
    )
    daily_data_by_date[date_label] = get_daily_temperature_data(
        temperature_data,
        dates,
        date_index,
        experiment_type="all_data",
    ).join(station_basis, on="stno")

    ## Make 75/25 validation split
    val_rng = np.random.default_rng(42)
    shuffled_indices = val_rng.permutation(
        len(daily_data_by_date[date_label])
    )
    val_size = round(len(daily_data_by_date[date_label]) * 0.15)
    val_idx = shuffled_indices[:val_size]
    train_idx = shuffled_indices[val_size:]

    daily_data_by_date[date_label] = {
        "train": daily_data_by_date[date_label].iloc[train_idx],
        "val": daily_data_by_date[date_label].iloc[val_idx],
    }

print("Daily data prepared", flush=True)
print(
    f"Hyperparameter search over {len(hyperparameter_grid)} combinations",
    flush=True,
)

## If the RMSE falls below a threshold, it's good enough to plot and look at
MAXT_RMSE_GRID_THRESHOLD = None
MINT_RMSE_GRID_THRESHOLD = None

## There's a little extra code here so you can do checkpoints
## as we go, hyperparameter results will continually get saved, but this way
## we won't accidentally overwrite anything
checkpoint_paths = sorted(
    hyp_figure_directory.glob("hyperparameter_results_checkpoint_*.csv")
)
previous_results = (
    pd.read_csv(checkpoint_paths[-1]) if checkpoint_paths else pd.DataFrame()
)
hyperparameter_results = previous_results.to_dict("records")
start_combination_index = (
    int(previous_results["combination_index"].max()) + 1
    if not previous_results.empty
    else 1
)
run_timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S_%f")
results_path = hyp_figure_directory / (
    f"hyperparameter_results_checkpoint_{run_timestamp}.csv"
)
if hyperparameter_results:
    pd.DataFrame(hyperparameter_results).to_csv(results_path, index=False)


## Start the tests
print(f"Starting at combination {start_combination_index}", flush=True)

## Loop over each hyperparameter combination
for combination_index, combination in enumerate(
    hyperparameter_grid[start_combination_index - 1 :],
    start=start_combination_index,
):

    print(
        f"Start combination {combination_index}/"
        f"{len(hyperparameter_grid)} {combination}",
        flush=True,
    )

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

    ## Write RMSEs to here for each date
    daily_maxt_rmses = np.empty(len(date_indices))
    daily_mint_rmses = np.empty(len(date_indices))
    models_by_date = []

    ## Loop over each date
    for daily_index, date_index in enumerate(date_indices):

        current_date = dates.iloc[date_index]
        date_label = (
            f"{int(current_date['year'])}_{int(current_date['month']):02d}_"
            f"{int(current_date['day']):02d}"
        )

        ## XGBoost setup
        train_data = daily_data_by_date[date_label]["train"]
        val_data = daily_data_by_date[date_label]["val"]
        X_train_maxt = xgb.DMatrix(
            train_data[model_features],
            label=train_data["y_maxt"].to_numpy(dtype=float),
        )
        X_train_mint = xgb.DMatrix(
            train_data[model_features],
            label=train_data["y_mint"].to_numpy(dtype=float),
        )
        X_val = xgb.DMatrix(
            val_data[model_features]
        )

        ## Train!
        model_maxt = xgb.train(
            model_params,
            X_train_maxt,
            num_boost_round=combination["num_boost_round"],
        )
        model_mint = xgb.train(
            model_params,
            X_train_mint,
            num_boost_round=combination["num_boost_round"],
        )

        ## Predict!
        pred_maxt = model_maxt.predict(X_val) + val_data["LTA_maxt"]
        pred_mint = model_mint.predict(X_val) + val_data["LTA_mint"]

        ## Get and print RMSEs
        daily_maxt_rmses[daily_index] = rmse(pred_maxt, val_data["maxt"])
        daily_mint_rmses[daily_index] = rmse(pred_mint, val_data["mint"])
        models_by_date.append(
            (date_label, model_maxt, model_mint)
        )

        print(
            f"  {date_label}: max RMSE={daily_maxt_rmses[daily_index]:.4f}, "
            f"min RMSE={daily_mint_rmses[daily_index]:.4f}",
            flush=True,
        )

    ## Final RMSEs for the hyperparameter combination
    rmse_maxt = np.mean(daily_maxt_rmses)
    rmse_mint = np.mean(daily_mint_rmses)
    print(
        f"Max Temp RMSE={rmse_maxt:.4f},"
        f" Min Temp RMSE={rmse_mint:.4f}",
        flush=True,
    )

    ## If RMSEs are low enough, plot a figure
    if (
        MAXT_RMSE_GRID_THRESHOLD is not None
        and MINT_RMSE_GRID_THRESHOLD is not None
        and rmse_maxt < MAXT_RMSE_GRID_THRESHOLD
        and rmse_mint < MINT_RMSE_GRID_THRESHOLD
    ):
        grid_matrix = xgb.DMatrix(prediction_grid[model_features])
        for date_label, month, model_maxt, model_mint in models_by_date:
            grid_prediction = prediction_grid.copy()
            grid_prediction["maxt"] = model_maxt.predict(grid_matrix) + (
                grid_LTAs_max_9120[f"m_{month}"].to_numpy(dtype=float)
            )
            grid_prediction["mint"] = model_mint.predict(grid_matrix) + (
                grid_LTAs_min_9120[f"m_{month}"].to_numpy(dtype=float)
            )
            daily_temperature_data = pd.concat(
                (
                    daily_data_by_date[date_label]["train"],
                    daily_data_by_date[date_label]["val"],
                )
            ).sort_index()
            figure_path = hyp_figure_directory / (
                f"spatial_XGBoost_hyp_{date_label}_{combination_index:03d}.jpg"
            )
            daily_temperature_plot(
                grid_prediction,
                daily_temperature_data,
                plot_destination=str(figure_path),
            )

    ## Save results and update
    hyperparameter_results.append(
        {
            "combination_index": combination_index,
            **combination,
            "rmse_maxt": rmse_maxt,
            "rmse_mint": rmse_mint,
        }
    )
    pd.DataFrame(hyperparameter_results).to_csv(results_path, index=False)

print(f"Saved hyperparameter results to {results_path}", flush=True)
