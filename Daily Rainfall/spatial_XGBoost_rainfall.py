
## Spatially embedded XGBoost modelling (seXGBoost, I need a better name)

## XGBoost is a standard, but powerful ML model, by including spatial
## basis functions as predictors, you can incorporate spatial effects
## into the model in a simple, effective way

# %% 

# ## Setup and data
from pathlib import Path
from time import perf_counter

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
grid_geodata = session["grid_geodata"]
grid_LTAs_9120 = session["grid_LTAs_9120"]
daily_rain_plot = session["daily_rain_plot"]
