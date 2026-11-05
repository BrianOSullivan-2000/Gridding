import matplotlib.pyplot as plt
import numpy as np


def prediction_plot(observed, predicted, metrics=False, unit=""):
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)

    limits = np.concatenate((observed.ravel(), predicted.ravel()))
    min_limit = np.nanmin(limits)
    max_limit = np.nanmax(limits)

    _, ax = plt.subplots()
    ax.scatter(observed, predicted)
    ax.plot(
        [min_limit, max_limit],
        [min_limit, max_limit],
        color="red",
        linestyle="--",
    )
    ax.set_xlim(min_limit, max_limit)
    ax.set_ylim(min_limit, max_limit)
    ax.set_xlabel("Observed")
    ax.set_ylabel("Predicted")

    if metrics:
        difference = predicted - observed
        rmse_value = np.round(np.sqrt(np.mean(difference**2)), 3)
        me_value = np.round(np.mean(difference), 3)
        r2_value = np.round(
            1
            - np.sum(difference**2)
            / np.sum((observed - np.mean(predicted)) ** 2),
            3,
        )
        mae_value = np.round(np.mean(np.abs(difference)), 3)

        ax.set_title(
            f"RMSE = {rmse_value}{unit}     R2 = {r2_value}"
            f"     Mean Error = {me_value}     MAE = {mae_value}"
        )

    return ax