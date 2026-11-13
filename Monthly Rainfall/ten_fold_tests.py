## Validation tests for monthly rainfall grids
## Ten-fold cross-validation

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

dill.load_session("Data/Monthly_Rainfall/train_test_ten_fold_2016-2025.dill")

session = globals()
rain_data = session["rain_data"]
dates = session["dates"]
island_outline = session["island_outline"]
get_monthly_rain_data = session["get_monthly_rain_data"]
update_monthly_predictions = session["update_monthly_predictions"]
collect_metrics = session["collect_metrics"]

# %% 

## Save final results here

final_results = rain_data.loc[:, ["stno", "year", "month", "rain", "fold"]].copy()

# %%

## For XGBoost, setup spatial basis functions ahead of time
nodelist = [
    make_mesh(island_outline, a=0.15, q=30, buffer=0, simplify=10000, plot=False),
    make_mesh(island_outline, a=0.02, q=30, buffer=2000, simplify=10000, plot=False),
]

threshold = 4000
for index in range(1, len(nodelist)):
    previous_nodes = np.vstack(nodelist[:index])
    distances, _ = cKDTree(previous_nodes).query(nodelist[index], k=1)
    nodelist[index] = nodelist[index][distances > threshold]

station_locations = rain_data.drop_duplicates(subset="stno").loc[
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

rain_data = rain_data.assign(
    eastone=rain_data["east"] / 1000,
    northone=rain_data["north"] / 1000,
    aspectsin=np.sin(rain_data["aspect_10km"]),
    aspectcos=np.cos(rain_data["aspect_10km"]),
)

# %% 

## spatial XGBoost ten-fold CV
geodata_features = [
    "eastone",
    "northone",
    "points5",
    "dist2c",
    "n5",
    "e5",
    "w5",
    "s5",
    "exp25k",
    "aspectsin",
    "aspectcos",
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


def predict_values(train_data, target_data):
    train_data = train_data.dropna(subset=["y"])
    training_matrix = xgb.DMatrix(
        train_data[model_features],
        label=train_data["y"].to_numpy(dtype=float),
    )
    model = xgb.train(
        xgb_params,
        training_matrix,
        num_boost_round=num_boost_round,
    )
    target_matrix = xgb.DMatrix(target_data[model_features])
    return model.predict(target_matrix)


for fold_number, current_fold in enumerate(
    sorted(rain_data["fold"].dropna().unique()),
    start=1,
):
    print(f"Starting fold {fold_number} (fold value {current_fold})")
    fold_start = perf_counter()

    for date_index in range(len(dates)):
        monthly_data = get_monthly_rain_data(
            rain_data,
            dates,
            date_index=date_index,
            current_fold=current_fold,
        )
        for split in ("train", "test"):
            monthly_data[split] = monthly_data[split].join(
                station_basis,
                on="stno",
            )

        predicted_y = predict_values(monthly_data["train"], monthly_data["test"])
        rain_data = update_monthly_predictions(
            rain_data,
            predicted_y,
            dates,
            date_index=date_index,
            current_fold=current_fold,
        )
        print(f"Date index {date_index} completed")

    print(f"Fold {fold_number} completed in {perf_counter() - fold_start:.1f} seconds")

metrics_table = collect_metrics(None, "Spatial XGBoost", rain_data)
print(metrics_table.to_string(index=False))

final_results["Spatial XGBoost"] = rain_data["predicted_rain"].to_numpy()

results_file = Path(
    "Results/Monthly_Rainfall/ten_fold_2016-2025/spatial_XGBoost.csv"
)
results_file.parent.mkdir(parents=True, exist_ok=True)
metrics_table.to_csv(results_file, index=False)

final_results_file = results_file.with_name("final_results_python.csv")
final_results.to_csv(final_results_file, index=False)
