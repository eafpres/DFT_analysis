#
  source("dft_analysis.R")
  library(dplyr)
  library(ggplot2)
#
# test only
#
  time_units <- 730
#
  period_1 <- 90
  period_2 <- 10000
  period_3 <- 10000
#
  granularity <- 10
  time <- seq(1, time_units, 1 / granularity)
#
# sinusoid with period "period_1" units
#
  X1 <- sin(2 * pi * time / period_1)
#
# sinusoid with period "period_2" units
#
  X2 <- sin(2 * pi * time / period_2)
#
# sinusoid with period "period_3" units
# and a phase offset of of "offset time units
#
  offset <- 17
  X3 <- sin(offset / period_3 + 2 * pi * time / period_3)
#
# linear growth with time
#  
  X4 <- 0.0 * time
#
# Y is a linear function of the two sinusoids plus the growth term
#
  Y <- 2.5 * X1 + X4
  my_data <- data.frame(time = time,
                        X1 = X1,
                        X2 = X2,
                        X3 = X3,
                        X4 = X4,
                        Y = Y,
                        stringsAsFactors = FALSE)
#
# plot the first 300 time units
#
  my_data %>%
    filter(time < 300) %>%
    ggplot(aes(x = time, y = Y)) +
    geom_line(color = "blue", size = 1)
#
# run the DFT analysis
#
  x_var <- "time"
  y_var <- "Y"
  low_cutoff <- 1000
  dc_threshold <- low_cutoff - 0
  noise_threshold <- 0.05
  significance_threshold <- 1
  peak_threshold_ratio <- 0.01
  cycles = 3
#
  dft_results <-
    dft_analysis(data = my_data, 
                 x_var = x_var,
                 y_var = y_var,
                 low_cutoff = low_cutoff,
                 dc_threshold = dc_threshold,
                 peak_threshold_ratio = peak_threshold_ratio,
                 cycles = cycles)
#
# dft_results has the periods and the model variable names
#
  print(dft_results)