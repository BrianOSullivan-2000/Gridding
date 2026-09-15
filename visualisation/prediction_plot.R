
################################################################################

# Plot one value against another (observed vs. predicted)
# This is just a super simple function

################################################################################

prediction_plot <- function(
    observed,
    predicted,
    metrics = FALSE,
    unit = ""
) {

    # observed          -observed values
    # predicted         -predicted values
    # metrics           -if TRUE, show RMSE and R2 values
    # unit              -character if RMSE has a unit

    require(ggplot2)

    min_limit <- min(observed, predicted, na.rm = TRUE)
    max_limit <- max(observed, predicted, na.rm = TRUE)
    plot_df <- data.frame(observed, predicted)

    pred_plot <- ggplot(data = plot_df) +
        geom_point(x = observed, y = predicted) +
        geom_abline(
            slope = 1,
            intercept = 0,
            color = "red",
            linetype = "dashed"
        ) +
        xlim(c(min_limit, max_limit)) +
        ylim(c(min_limit, max_limit)) +
        xlab("Observed") + ylab("Predicted")

    if (metrics) {
        require(Metrics)
        R2 <- function(y1, y2) {
            1 - sum((y2 - y1)^2) / sum((y2 - mean(y1))^2)
        }
        me <- function(y1, y2) {
            mean(y1 - y2)
        }

        rmse_value <- round(rmse(predicted, observed), 3)
        me_value <- round(me(predicted, observed), 3)
        R2_value <- round(R2(predicted, observed), 3)
        mae_value <- round(mae(predicted, observed), 3)

        pred_plot <-
            pred_plot +
            ggtitle(
                paste0(
                    "RMSE = ", rmse_value, unit,
                    "     ", "R2 = ", R2_value,
                    "     ", "Mean Error = ", me_value,
                    "     ", "MAE = ", mae_value
                )
            )
    }

    pred_plot
}