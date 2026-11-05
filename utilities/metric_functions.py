import numpy as np


## Root Mean Squared Error
def rmse(pred, obs):
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    return np.sqrt(np.mean((pred - obs) ** 2))


## Mean Absolute Error
def mae(pred, obs):
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    return np.mean(np.abs(pred - obs))


## Coefficient of Determination
def R2(pred, obs):
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    return 1 - np.sum((obs - pred) ** 2) / np.sum((obs - np.mean(obs)) ** 2)


## Mean Error (Bias)
def me(pred, obs):
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    return np.mean(pred - obs)


## Jensen-Shannon Divergence (compares empirical distributions)
def JSD(y, y_hat):
    y = np.asarray(y, dtype=float)
    y_hat = np.asarray(y_hat, dtype=float)

    p = y / np.sum(y)
    q = y_hat / np.sum(y_hat)
    middle = 0.5 * (p + q)

    p_positive = p > 0
    q_positive = q > 0
    kl_p_middle = np.sum(
        p[p_positive] * np.log2(p[p_positive] / middle[p_positive])
    )
    kl_q_middle = np.sum(
        q[q_positive] * np.log2(q[q_positive] / middle[q_positive])
    )

    return 0.5 * kl_p_middle + 0.5 * kl_q_middle


## Ten largest errors
def ME10(y_hat, y):
    y_hat = np.asarray(y_hat, dtype=float)
    y = np.asarray(y, dtype=float)
    errors = np.sort(np.round(y_hat - y, 2).ravel())[::-1][:10]
    largest_errors = [f"{error:g}" for error in errors]
    largest_errors.extend(["NA"] * (10 - len(largest_errors)))
    return ", ".join(largest_errors)


## Check bias for Wet Days at different thresholds (Daily Rainfall)
def _wet_day_bias(y_hat, y, threshold):
    y_hat = np.asarray(y_hat, dtype=float)
    y = np.asarray(y, dtype=float)
    predicted_wet_days = np.sum(y_hat > threshold)
    observed_wet_days = np.sum(y > threshold)
    return (predicted_wet_days - observed_wet_days) * 100 / observed_wet_days


def Rain_Day_Bias(y_hat, y):
    return _wet_day_bias(y_hat, y, 0.2)

def Wet_Day_Bias(y_hat, y):
    return _wet_day_bias(y_hat, y, 1)

def Heavy_Wet_Day_Bias(y_hat, y):
    return _wet_day_bias(y_hat, y, 10)