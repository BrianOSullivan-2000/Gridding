
#-----------------------------------------------------------#
## AUTHORED BY BARRY COONAN - 18th November 2026 ##
#-----------------------------------------------------------#

## Brian O'Sullivan - just changed formatting to match this project's
## coding style, and added some print statments for timing

# calculate the 1991-2020 daily lta using month lta for each station
# uses a spline fitted to adjusted monthly normals
# adjusted normals ensure average of daily splines = monthly normal

# paramnorm <- rainnorm
# parameter <- "rain"

month_lta_to_day_lta <- function(paramnorm, parameter, yearis = 2020) {

    require(lubridate)
    require(stringr)

    function_start <- proc.time()[["elapsed"]]
    station_count <- nrow(paramnorm)

    print_phase_timing <- function(phase_name, phase_start) {
        elapsed <- proc.time()[["elapsed"]] - phase_start
        cat(sprintf(
            "TIMING %s: %.3f sec (station count: %d)\n",
            phase_name, elapsed, station_count
        ))
        invisible(elapsed)
    }

    #################################################################
    # spline interpolation of 'daily' values using month normal means
    #################################################################

    # model of days of year rain/temp values using month norm of
    # rain/temp as middle of month value and fit a spline

    # create a dataframe with months, days of year
    for (monthis in 1:12) {
        aaa <- cbind(
            rep(monthis, days_in_month(monthis)),
            seq(1:days_in_month(monthis))
        )

        if (monthis == 1) {
            year_frame <- aaa
        } else {
            year_frame <- rbind(year_frame, aaa)
        }
    }

    year_frame <- data.frame(year_frame)
    # add day number
    year_frame$day_num <- seq(1, nrow(year_frame))
    colnames(year_frame) <- c("month", "day", "day_num")

    # to do a fit of rain/temp per day create a dataset with the norm
    # value in the middle of each month
    # year_frame$value <- 0
    year_frame <- cbind(
        year_frame,
        matrix(
            rep(0, nrow(year_frame) * nrow(paramnorm)),
            nrow = nrow(year_frame)
        )
    )
    year_day <- 0

    # for each month place the monthly average at the middle of each month
    for (monthis in 1:12) {
        day_of_year <- year_day + floor(days_in_month(monthis) / 2)
        if (parameter == "temp") {
            year_frame[day_of_year, 4:ncol(year_frame)] <-
                paramnorm[[monthis + 3]]
        } else {
            year_frame[day_of_year, 4:ncol(year_frame)] <-
                paramnorm[[monthis + 3]] / days_in_month(monthis)
        }
        year_day <- days_in_month(monthis) + year_day
    }
    colnames(year_frame)[4:ncol(year_frame)] <- paste0("stn_", paramnorm$stno)

    # create a 3 year long dataset with normal value at the middle of
    # each month and exclude all other days
    # this is so spline calc will be continious dec/jan
    # also n-day average can include days in the previous/future year
    data_days <- year_frame[year_frame[, 4] > 0, ]

    # # save this -> testing
    # write.csv(
    #     data_days,
    #     paste0(
    #         homedir,
    #         "OUT/monthly_for_spline.csv"
    #     ),
    #     row.names = FALSE
    # )

    data_days2 <- data_days
    data_days3 <- data_days
    data_days2$day_num <- data_days2$day_num + max(data_days$day_num) + 15
    data_days3$day_num <- data_days3$day_num + max(data_days2$day_num) + 15
    # put these together
    data_days <- rbind(data_days, data_days2, data_days3)

    # reset the day number so dataset includes negative days
    # (days before year start) and days into new year
    # need these extra days to calculate a rolling average
    data_days$day_num <- data_days$day_num - 364

    # spline fit using normals at middle of months
    # the spline is basically a moving monthly average at this stage
    for (stn in seq_len(nrow(paramnorm))) {
        days_spline <- spline(
            data_days$day_num, data_days[, stn + 3],
            method = "natural", xmin = -360, xmax = 360 * 2,
            n = (nrow(year_frame) * 3) - 14
        )

        if (stn == 1) {
            spline_calc <- days_spline
        } else {
            spline_calc <- merge(spline_calc, days_spline, by = c("x"))
        }
    }

    spline_calc <- data.frame(spline_calc)
    colnames(spline_calc) <- c("day_num", paste0("stn_", paramnorm$stno))

    if (exists("initial_spline_start")) {
        print_phase_timing("initial spline", initial_spline_start)
    }

    # # just to see it makes sense
    # plot(spline_calc$day_num, spline_calc$stn_1275, pwd = 0.5, type = "l",
    #      main = "using month lta", ylim = c(0, 5))
    # points(
    #     data_days$day_num,
    #     data_days$stn_1275,
    #     col = "red", pch = 19, type = "p"
    # )

    # append a date column
    spline_calc$date <- as.Date(paste0(yearis, "-01-01"), format("%Y-%m-%d"))
    spline_calc$date <- spline_calc$date + spline_calc$day_num
    spline_calc$month <- month(spline_calc$date)
    spline_calc$year <- year(spline_calc$date)

    # re arrange
    spline_calc <- spline_calc[
        ,
        c((ncol(spline_calc) - 2):ncol(spline_calc),
          1:(ncol(spline_calc) - 3))
    ]
    # this is a rolling average of monthly long term average

    # to convert to a daily estimate need to ensure that the mean of the
    #  predicted ('spline') daily values for each month = the monthly mean
    # so calculate the average of each month as determined using the spline
    # of the monthly average
    spline_avg <- spline_calc[spline_calc$year == yearis, ]
    spline_avg <- aggregate(spline_avg, by = list(spline_avg$month), mean)
    spline_avg <- spline_avg[, c(3, 6:ncol(spline_avg))]

    # get the monthly long term averages
    norm_avg <- data_days[data_days$day_num > 0 & data_days$day_num < 366, ]
    norm_avg <- norm_avg[, c(1, 4:ncol(norm_avg))]
    colnames(norm_avg) <- c(
        "month",
        paste0("nor_", str_sub(colnames(norm_avg)[2:ncol(norm_avg)], 5, 10))
    )

    # calculate the difference between the month mean of spline dailys
    # and the normal month mean this needs to be corrected for
    avg_diff <- norm_avg - spline_avg
    avg_diff <- data.frame(avg_diff)
    avg_diff$month <- seq(1, 12)
    # write.csv(
    #     avg_diff,
    #     paste0(
    #         homedir,
    #         "OUT/diff_month-lta_daily-lta_stations_",
    #         parameter, ".csv"
    #     ),
    #     row.names = FALSE
    # )

    # RMSE norm value versus calculated daily values per month
    diff_rmse <- apply(
        avg_diff[, 2:ncol(avg_diff)],
        2,
        function(x) sqrt(mean(x^2))
    )
    diff_rmse_raw <- diff_rmse
    stnos <- as.integer(str_replace(
        colnames(avg_diff)[2:ncol(avg_diff)],
        pattern = "nor_",
        replacement = ""
    ))
    month_norm_avg <- apply(norm_avg[, 2:ncol(norm_avg)], 2, mean)
    diff_rmse_raw <- data.frame(t(rbind(stnos, month_norm_avg, diff_rmse_raw)))

    # this is the rmse of daily lta versus monthly lta with
    # no adjustment to the month lta (0)
    diff_rmse <- c(0, diff_rmse)

    # to get better daily estimates adjust the month normal value using
    # difference between norm monthly values and the spline predicted monthly
    # mean of daily values

    #########################################################################
    # spline interpolation of 'daily' values using adjusted month normal means
    #########################################################################

    initial_spline_start <- proc.time()[["elapsed"]]

    # model of days of year versus rain/temp using month norm as
    # middle of month value  and fit a spline
    # adjust_factor <- 1 # testing

    adjusted_loop_start <- proc.time()[["elapsed"]]
    adjustment_factors <- seq(0.1, 3, 0.1)

    for (adjustment_index in seq_along(adjustment_factors)) {
        adjust_factor <- adjustment_factors[adjustment_index]

        # if(exists("year_frame")) rm(year_frame)
        # create a dataframe with months, days of year
        for (monthis in 1:12) {
            aaa <- cbind(
                rep(monthis, days_in_month(monthis)),
                seq(1:days_in_month(monthis))
            )

            if (monthis == 1) {
                year_frame <- aaa
            } else {
                year_frame <- rbind(year_frame, aaa)
            }
        }
        year_frame <- data.frame(year_frame)
        # add day number
        year_frame$day_num <- seq(1, nrow(year_frame))
        colnames(year_frame) <- c("month", "day", "day_num")

        # to do a fit of rain/temp per day create a dataset with the norm
        # value in the middle of each month
        # year_frame$value <- 0
        year_frame <- cbind(
            year_frame,
            matrix(
                rep(0, nrow(year_frame) * nrow(paramnorm)),
                nrow = nrow(year_frame)
            )
        )
        year_day <- 0
        monthis <- 1 # testing

        # put adjusted normals monthly average at the middle of each month
        for (monthis in 1:12) {
            day_of_year <- year_day + floor(days_in_month(monthis) / 2)

            if (parameter == "temp") {
                year_frame[day_of_year, 4:ncol(year_frame)] <-
                    paramnorm[[monthis + 3]] +
                    adjust_factor *
                    avg_diff[monthis, c(2:ncol(avg_diff))]
            } else {
                year_frame[day_of_year, 4:ncol(year_frame)] <-
                    (paramnorm[[monthis + 3]] / days_in_month(monthis)) +
                    adjust_factor * avg_diff[monthis, c(2:ncol(avg_diff))]
            }

            year_day <- days_in_month(monthis) + year_day
        }
        colnames(year_frame)[4:ncol(year_frame)] <-
            paste0("stn_", paramnorm$stno)

        # create a 3 year long dataset with normal value at the middle of each
        # month
        # this is so start and end of year could use rolling average of
        # calculated 'daily' values
        data_days_adj <- year_frame[year_frame[, 4] > 0, ]
        data_days_adj2 <- data_days_adj
        data_days_adj3 <- data_days_adj
        data_days_adj2$day_num <- data_days_adj2$day_num +
            max(data_days_adj$day_num) + 15
        data_days_adj3$day_num <- data_days_adj3$day_num +
            max(data_days_adj2$day_num) + 15
        # put these together
        data_days_adj <- rbind(data_days_adj, data_days_adj2, data_days_adj3)

        # reset these so dataset includes negative days
        # (days before year start) and days into new year
        # need these extra days to calculate a rolling average
        data_days_adj$day_num <- data_days_adj$day_num - 364

        # spline fit using adjusted normals at middle of months
        for (stn in seq_len(nrow(paramnorm))) {
            days_spline <- spline(
                data_days_adj$day_num, data_days_adj[, stn + 3],
                method = "natural", xmin = -360, xmax = 360 * 2,
                n = (nrow(year_frame) * 3) - 14
            )

            if (stn == 1) {
                spline_calc_adj <- days_spline
            } else {
                spline_calc_adj <- merge(
                    spline_calc_adj,
                    days_spline, by = c("x")
                )
            }
        }

        spline_calc_adj <- data.frame(spline_calc_adj)
        colnames(spline_calc_adj) <-
            c("day_num", paste0("stn_", paramnorm$stno))

        # # just to see it makes sense
        day_num <- days_in_month(1:12)
        day_num <- cumsum(day_num)

        # if (adjust_factor %in% c(0.1, 1, 1.2, 1.6, 2, 3))
        #   {
        #   png(paste0(homedir, "stno_1475_", adjust_factor, "_", parameter,
        #   ".png"))
        #   plot(spline_calc_adj$day_num, spline_calc_adj$stn_1185, lwd = 1,
        #   type = "l", lty = 1, col = "red",
        #        main = "stno 1475 gurteen using adjusted month lta red, month
        #        lta black\ngray lines are months",
        #        ylim = c(1, 4), xlim = c(10, 365),
        #        xlab = "day of year", ylab = "daily rain")
        #   points(spline_calc$day_num, spline_calc$stn_1185, lwd = 1, type =
        #   "l", lty = 1)
        #   points(data_days$day_num, data_days$stn_1185, pch = 1, type = "p")
        #   abline(v = day_num, col = "gray80")
        #   dev.off()
        # }
        # points(data_days_adj$day_num, data_days_adj$stn_485, col = "red", pch
        # = 19, type = "p")

        # append the date column
        spline_calc_adj$date <-
            as.Date(paste0(yearis, "-01-01"), format("%Y-%m-%d"))
        spline_calc_adj$date <- spline_calc_adj$date + spline_calc_adj$day_num
        spline_calc_adj$month <- month(spline_calc_adj$date)
        spline_calc_adj$year <- year(spline_calc_adj$date)

        # re arrange
        spline_calc_adj <- spline_calc_adj[
            ,
            c(
                (ncol(spline_calc_adj) - 2):ncol(spline_calc_adj),
                1:(ncol(spline_calc_adj) - 3)
            )
        ]
        # this is a rolling average of monthly long term average
        # to convert to a daily estimate need to endure the predicted
        # ('spline') mean
        # for each month = the monthly mean
        spline_avg <- spline_calc_adj[spline_calc_adj$year == yearis, ]
        spline_avg <- aggregate(spline_avg, by = list(spline_avg$month), mean)
        spline_avg <- spline_avg[, c(3, 6:ncol(spline_avg))]

        # calculate the difference between the spline derived month mean and
        # the normal month mean
        avg_diff2 <- norm_avg - spline_avg

        # put the mnth lta versus calculated daily lta together
        if (adjust_factor == 0.1) {
            month_diff <- cbind(rep(adjust_factor, 12), seq(1, 12), avg_diff2)
        } else {
            month_diff <- rbind(
                month_diff,
                cbind(rep(adjust_factor, 12), seq(1, 12), avg_diff2)
            )
        }

        if (adjustment_index == 1 || adjustment_index %% 5 == 0 ||
                adjustment_index == length(adjustment_factors)) {

            elapsed_factor_loop <-
                proc.time()[["elapsed"]] - adjusted_loop_start
            estimated_factor_loop <-
                elapsed_factor_loop * length(adjustment_factors) /
                adjustment_index
            cat(sprintf(
                paste0(
                    "TIMING adjusted loop: %.0f/%.0f factors,",
                    "elapsed %.3f sec ",
                    "(station count: %d, full-loop estimate:",
                    "%.1f sec / %.1f min)\n"
                ),
                adjustment_index, length(adjustment_factors),
                elapsed_factor_loop,
                station_count,
                estimated_factor_loop, estimated_factor_loop / 60
            ))
        }
    }

    print_phase_timing("adjusted-factor loop", adjusted_loop_start)


    # difference bwteen the mean of daily month average and the normals monthly
    # lta
    month_diff_abs <- month_diff[, c(1, 2, 4:ncol(month_diff))]

    # make these difference absolute values
    month_diff_abs[, 3:ncol(month_diff_abs)] <-
        sapply(month_diff_abs[, 3:ncol(month_diff_abs)], abs)
    colnames(month_diff_abs)[1:2] <- c("adjust", "month")

    # get the minimum difference for each station month = best aggreement daily
    # predictions and monthly lta
    month_mins <- month_diff_abs
    month_mins <- aggregate(
        month_mins[, c(2:ncol(month_mins))],
        by = list(month_mins$month),
        min
    )[, 2:ncol(month_mins)]

    # identify the adjusment corresponding with these month minimum values
    optimization_start <- proc.time()[["elapsed"]]
    month_diff_mins <- month_diff_abs

    stn <- 2
    monthis <- 1

    for (stn in 2:(ncol(month_diff_abs) - 1)) {
        for (monthis in 1:12) {
            # for the particular station/month extract the difference daily to
            # monthly
            month_stn <-
                month_diff_abs[month_diff_abs$month == monthis, c(1:2, stn + 1)]

            # determine the 'best' adjustment for this month which is the
            # middle adjustment -> rank 15
            month_stn <- month_stn[month_stn[, 3] == min(month_stn[, 3]), ]

            if (stn == 2 && monthis == 1) {
                # add the station number
                month_stn$stno <-
                    as.integer(str_replace(
                        colnames(month_stn)[3],
                        pattern = "nor_", replacement = ""
                    ))
                colnames(month_stn)[3] <- "diff"
                optimised_adjust <- month_stn
            } else {
                # add the station number
                month_stn$stno <-
                    as.integer(str_replace(
                        colnames(month_stn)[3],
                        pattern = "nor_", replacement = ""
                    ))
                colnames(month_stn)[3] <- "diff"
                optimised_adjust <- rbind(optimised_adjust, month_stn)
            }
        }
    }

    ############################################################################
    ##########
    #
    # rerun calculation of daily values based on this optimum adjustment of
    # mean values
    #
    ############################################################################
    ##########

    final_spline_start <- proc.time()[["elapsed"]]

    # create a dataframe with months, days of year
    for (monthis in 1:12) {
        aaa <- cbind(
            rep(monthis, days_in_month(monthis)),
            seq(1:days_in_month(monthis))
        )

        if (monthis == 1) {
            year_frame <- aaa
        } else {
            year_frame <- rbind(year_frame, aaa)
        }
    }
    year_frame <- data.frame(year_frame)
    # add day number
    year_frame$day_num <- seq(1, nrow(year_frame))
    colnames(year_frame) <- c("month", "day", "day_num")

    # to do a fit of rain/temp per day create a dataset with the norm value in
    # the middle of each month
    # year_frame$value <- 0
    year_frame <- cbind(
        year_frame,
        matrix(
            rep(0, nrow(year_frame) * nrow(paramnorm)),
            nrow = nrow(year_frame)
        )
    )
    year_day <- 0

    # put adjusted normals monthly average at the middle of each month
    for (monthis in 1:12) {
        day_of_year <- year_day + floor(days_in_month(monthis) / 2)

        if (parameter == "temp") {
            # select this months adjustment of monthly to better predict daily
            optimum_adjust <-
                optimised_adjust$adjust[optimised_adjust$month == monthis]
            # adjust monthly normal by station
            year_frame[day_of_year, 4:ncol(year_frame)] <-
                paramnorm[[monthis + 3]] +
                optimum_adjust * avg_diff[monthis, c(2:ncol(avg_diff))]
        } else {
            # select this months adjustment of monthly to better predict daily
            optimum_adjust <-
                optimised_adjust$adjust[optimised_adjust$month == monthis]
            # adjust monthly normal by station
            year_frame[day_of_year, 4:ncol(year_frame)] <-
                (paramnorm[[monthis + 3]] / days_in_month(monthis)) +
                optimum_adjust * avg_diff[monthis, c(2:ncol(avg_diff))]
        }

        year_day <- days_in_month(monthis) + year_day
    }
    colnames(year_frame)[4:ncol(year_frame)] <- paste0("stn_", paramnorm$stno)

    # create a 3 year long dataset with normal value at the middle of each month
    # this is so start and end of year coul use rolling average of calculated
    # 'daily' values
    data_days_adj <- year_frame[year_frame[, 4] > 0, ]

    data_days_adj2 <- data_days_adj
    data_days_adj3 <- data_days_adj
    data_days_adj2$day_num <- data_days_adj2$day_num +
        max(data_days_adj$day_num) + 15
    data_days_adj3$day_num <- data_days_adj3$day_num +
        max(data_days_adj2$day_num) + 15
    # put these together
    data_days_adj <- rbind(data_days_adj, data_days_adj2, data_days_adj3)

    # reset these so dataset includes negative days (days before year start)
    # and days into new year
    # need these extra days to calculate a rolling average
    data_days_adj$day_num <- data_days_adj$day_num - 364

    # spline fit using normals at middle of months
    # the spline is basically a moving monthly average at this stage
    for (stn in seq_len(nrow(paramnorm))) {
        days_spline <- spline(
            data_days_adj$day_num, data_days_adj[, stn + 3],
            method = "natural", xmin = -360, xmax = 360 * 2,
            n = (nrow(year_frame) * 3) - 14
        )

        if (stn == 1) {
            spline_calc_adj <- days_spline
        } else {
            spline_calc_adj <- merge(spline_calc_adj, days_spline, by = c("x"))
        }
    }

    # these are the daily predicted values
    spline_calc_adj <- data.frame(spline_calc_adj)
    colnames(spline_calc_adj) <- c("day_num", paste0("stn_", paramnorm$stno))

    # # just to see it makes sense
    day_num <- days_in_month(1:12)
    day_num <- cumsum(day_num)

    # append the date column
    spline_calc_adj$date <-
        as.Date(paste0(yearis, "-01-01"), format("%Y-%m-%d"))
    spline_calc_adj$date <- spline_calc_adj$date + spline_calc_adj$day_num
    spline_calc_adj$month <- month(spline_calc_adj$date)
    spline_calc_adj$year <- year(spline_calc_adj$date)

    # re arrange
    spline_calc_adj <-
        spline_calc_adj[
            ,
            c(
                (ncol(spline_calc_adj) - 2):ncol(spline_calc_adj),
                1:(ncol(spline_calc_adj) - 3)
            )
        ]

    # calculate the daily rainfall value
    spline_calc_adj$monthdays <- days_in_month(spline_calc_adj$month)

    # best fit monthly avgs based on calculated daily vlaues
    spline_avg <- spline_calc_adj[spline_calc_adj$year == yearis, ]
    spline_avg <- aggregate(spline_avg, by = list(spline_avg$month), mean)
    spline_avg <- spline_avg[, c(3, 6:ncol(spline_avg))]

    # # best fit
    # write.csv(spline_avg, paste0(homedir,
    # "OUT/adjusted_daily_mean_monthly_for_spline.csv"), row.names = FALSE)

    # save the stats
    avg_diff2 <- norm_avg - spline_avg[, 1:(ncol(spline_avg) - 1)]
    colnames(avg_diff2) <- colnames(spline_avg)[seq_len(ncol(avg_diff2))]
    avg_diff2$month <- seq(1, 12)
    # write.csv(avg_diff2, paste0(homedir,
    # "OUT/adjusted_diff_month-lta_daily-lta_stations_", parameter, ".csv"),
    # row.names = FALSE)

    # RMSE norm value versus calculated daily values per month
    adj_diff_rmse <-
        apply(avg_diff2[, 2:ncol(avg_diff2)], 2, function(x) sqrt(mean(x^2)))

    # save rmse calculation
    stnos <- as.integer(
        str_replace(
            colnames(avg_diff2)[2:ncol(avg_diff2)],
            pattern = "stn_", replacement = ""
        )
    )
    adj_diff_rmse <- data.frame(t(rbind(stnos, adj_diff_rmse)))
    # add in the result of just using unadjusted monthly normals
    rmse <- merge(diff_rmse_raw, adj_diff_rmse, by = "stnos")
    # write.csv(rmse, paste0(homedir, "/OUT/rmse_normal_vs_spline_daily_",
    # parameter, ".csv"), row.names= FALSE)

    # png(paste0(homedir, "stno_1185_spline_fitted_daily_", parameter, ".png"),
    # width = 900, height = 600)
    # plot(spline_calc_adj$day_num, spline_calc_adj$stn_1185, lwd = 1, type =
    # "l", lty = 1, col = "red",
    #      main = "stno 1185 piltown using adjusted month lta red, month lta
    #      black\ngray lines are months",
    #      ylim = c(2, 4.5), xlim = c(10, 365), cex = 18,
    #      xlab = "day of year", ylab = "daily rain")
    # points(spline_calc$day_num, spline_calc$stn_1185, lwd = 1, type = "l",
    # lty = 1)
    # points(data_days$day_num, data_days$stn_1185, pch = 1, type = "p")
    # abline(v = day_num, col = "gray80")
    # dev.off()

    # write.csv(
    #     spline_calc_adj[, c(1:(ncol(spline_calc_adj) - 1))],
    #     paste0(homedir, "OUT/daily_norms_", parameter, ".csv"),
    #     row.names = FALSE
    # )

    print_phase_timing("final spline and statistics", final_spline_start)

    station_columns <- paste0("stn_", paramnorm$stno)
    daily_normals <- spline_calc_adj[
        spline_calc_adj$year == yearis,
        c("date", "month", "year", "day_num", station_columns)
    ]

    print_phase_timing("total function", function_start)

    daily_normals
}
