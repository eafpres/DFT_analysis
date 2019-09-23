#
# demo how to use DFT function
#
# clear the environment
#
  rm(list = ls())
#
# get the function
#
  source("dft_analysis.R")
#
# libraries
# 
# for plotting
#
  library(ggplot2)
#
# for pipes and data processing
#
  library(dplyr)
#
  time_units <- 730
  period_1 <- 3.5
  period_2 <- 13
#
# here we generate 10 points per time unit granularity
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
# linear growth with time (to demonstrate DFT still works with trend)
#  
  X4 <- 0.15 * time
#
# Y is a linear function of the two sinusoids plus the growth term
#
  Y <- 2.5 * X1 + 3 * X2 + X4
  my_data <- data.frame(time = time,
                        X1 = X1,
                        X2 = X2,
                        X4 = X4,
                        Y = Y,
                        stringsAsFactors = FALSE)
#
# plot the first 100 time steps
#
  my_data %>%
    filter(time < 100) %>%
    ggplot(aes(x = time, y = Y)) +
    geom_line(color = "blue", size = 1)
#
# run the DFT analysis
#
# define some common parameters for reuse
#
  x_var <- "time"
  y_var <- "Y"
  low_cutoff <- 3000
  dc_threshold <- low_cutoff - 50
  peak_threshold_ratio <- 0.01
  cycles <- 10
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
#
# to demonstrate what happens if we too aggressively cut
# out low frequencies we increase the impact of dc_threshold
#
  dc_threshold <- low_cutoff - 100
  dft_results <-
    dft_analysis(data = my_data, 
                 x_var = x_var,
                 y_var = y_var,
                 low_cutoff = low_cutoff,
                 dc_threshold = dc_threshold,
                 peak_threshold_ratio = peak_threshold_ratio,
                 cycles = cycles)
#
# repeat with periods of even days and lower granularity
#
  time_units <- 730
#
# this makes one periodic component essentially DC
#
  period_1 <- 10000
#  
  period_2 <- 30
#
  granularity <- 1
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
# linear growth with time (to demonstrate DFT still works with trend)
#  
  X4 <- 0.0 * time
#
# Y is a linear function of the two sinusoids plus the growth term
#
  Y <- 2.5 * X1 + 3 * X2 + X4
  my_data <- data.frame(time = time,
                        X1 = X1,
                        X2 = X2,
                        X4 = X4,
                        Y = Y,
                        stringsAsFactors = FALSE)
#
# plot the first 100 time steps
#
  my_data %>%
    filter(time < 100) %>%
    ggplot(aes(x = time, y = Y)) +
    geom_line(color = "blue", size = 1)
#
# run the DFT analysis
#
  low_cutoff <- 365
  dc_threshold <- low_cutoff - 5
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
#
# repeat with mor complex series
#
  time_units <- 730
#
  period_1 <- 90
  period_2 <- 35
  period_3 <- 43
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
  X4 <- 0.07 * time
#
# Y is a linear function of the two sinusoids plus the growth term
#
  Y <- 2.5 * X1 + 3 * X2 + 4.31 * X3 + X4
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
  low_cutoff <- 0
  dc_threshold <- low_cutoff - 0
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
  