## Validation tests for daily temperature grids ##

## Spatially embedded XGBoost (seXGBoost, need a better name) ##
## Predicts Daily Temperature using a standard XGBoost model     ##
## except spatial basis functions are also included to add    ##
## a spatial effect                                           ##

# %%

## Load in data and functions

from pathlib import Path
from time import perf_counter

import dill
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.spatial import cKDTree

from utilities.add_basis_functions import add_basis_functions
from utilities.make_mesh import make_mesh

# Run setup_data_and_functions_five_fold.py first to create this session.
dill.load_session("Data/Daily_Temperature/train_test_five_fold_2025.dill")

session = globals()
temperature_data = session["temperature_data"]
dates = session["dates"]
island_outline = session["island_outline"]
get_daily_temperature_data = session["get_daily_temperature_data"]
update_daily_predictions = session["update_daily_predictions"]
collect_metrics = session["collect_metrics"]

# %%

## Save final results here

final_results = temperature_data.loc[:, [
    "stno", "year", "month", "day", "mint", "maxt", "fold"
]].copy()

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

station_locations = temperature_data.drop_duplicates(subset="stno").loc[
    :, ["stno", "east", "north"]
]
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
    columns=basis_columns,
)
station_basis.index = station_locations["stno"].to_numpy()

temperature_data = temperature_data.assign(
    eastone=temperature_data["east"] / 1000,
    northone=temperature_data["north"] / 1000,
    aspectsin=np.sin(temperature_data["aspect_10km"]),
    aspectcos=np.cos(temperature_data["aspect_10km"]),
)

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
def predict_values(train_data, target_data, response):
    train_data = train_data.dropna(subset=[response])

    Xtrain = xgb.DMatrix(
        train_data[model_features],
        label=train_data[response].to_numpy(dtype=float),
    )
    model = xgb.train(
        xgb_params,
        Xtrain,
        num_boost_round=num_boost_round,
    )
    Xtest = xgb.DMatrix(target_data[model_features])
    return model.predict(Xtest)


def predict_daily_values(daily_data):
    predicted_y_maxt = predict_values(
        daily_data["train"], daily_data["test"], "y_maxt"
    )
    predicted_y_mint = predict_values(
        daily_data["train"], daily_data["test"], "y_mint"
    )
    return predicted_y_maxt, predicted_y_mint

# %%

## Spatial XGBoost five-fold CV
for fold_number, current_fold in enumerate(
    sorted(temperature_data["fold"].dropna().unique()),
    start=1,
):
    print(f"Starting fold {fold_number} (fold value {current_fold})")
    fold_start = perf_counter()

    for date_index in range(len(dates)):
        daily_data = get_daily_temperature_data(
            temperature_data,
            dates,
            date_index=date_index,
            current_fold=current_fold,
        )

        for split in ("train", "test"):
            daily_data[split] = daily_data[split].join(
                station_basis,
                on="stno",
            )

        predicted_y_maxt, predicted_y_mint = predict_daily_values(daily_data)
        temperature_data = update_daily_predictions(
            temperature_data,
            predicted_y_maxt,
            predicted_y_mint,
            dates,
            date_index=date_index,
            current_fold=current_fold,
        )
        print(f"Date index {date_index} completed")

    print(f"Fold {fold_number} completed in {perf_counter() - fold_start:.1f} seconds")

metrics_table = collect_metrics(None, "Spatial XGBoost", temperature_data)
print(metrics_table.to_string(index=False))

final_results["Spatial XGBoost Min"] = temperature_data["predicted_mint"].to_numpy()
final_results["Spatial XGBoost Max"] = temperature_data["predicted_maxt"].to_numpy()

results_file = Path(
    "Results/Daily_Temperature/five_fold_CV_2025/spatial_XGBoost.csv"
)
results_file.parent.mkdir(parents=True, exist_ok=True)
metrics_table.to_csv(results_file, index=False)

final_results_file = results_file.with_name("final_results_python.csv")
final_results.to_csv(final_results_file, index=False)
