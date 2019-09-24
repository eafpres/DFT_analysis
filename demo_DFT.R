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
# test DFT with all defaults
#
  x_var <- "time"
  y_var <- "Y"
#
  dft_results <-
    dft_analysis(data = my_data, 
                 x_var = x_var,
                 y_var = y_var)
#
# dft_results has the periods and the model variable names
#
  print(dft_results)
#
# to demonstrate what happens if we too aggressively cut
# out low frequencies
# low_cutoff = 100 means 100 / 10 = 10 days is the longest
# period we will keep in the DFT results
#
# we also increase cycles from the defualt becuase
# in the DFT function it plots 'cycles' number of periods
# of the longest period, and if we cut of the longer
# period we won't plot very much data at the default
# of 3 cycles
#
  cycles = 10
  low_cutoff = 100
  dft_results <-
    dft_analysis(data = my_data, 
                 x_var = x_var,
                 y_var = y_var,
                 low_cutoff = low_cutoff,
                 dc_threshold = dc_threshold,
                 cycles = cycles)
#
# repeat with periods of even days and lower granularity
# and only use a pure tone  
#
  time_units <- 730
#  
  period_2 <- 30
#
  granularity <- 1
  time <- seq(1, time_units, 1 / granularity)
#
# sinusoid with period "period_2" units
#
  X2 <- sin(2 * pi * time / period_2)
#
# Y is a linear function of the two sinusoids plus the growth term
#
  Y <- 3 * X2
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
                 dc_threshold = dc_threshold)
#
# dft_results has the periods and the model variable names
#
# note that we suffer some accuracy loss due to the 
# lower time granularity; this cna be significant for
# predicting longer time scales
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
  x_var <- "time"
  y_var <- "Y"
  low_cutoff <- 1000
  dc_threshold <- low_cutoff - 0
  peak_threshold_ratio <- 0.01
  noise_threshold <- 0.05
  significance_threshold <- 1
  cycles <- 3
  smooth_n <- 1
#
  dft_results <-
    dft_analysis(data = my_data, 
                 x_var = x_var,
                 y_var = y_var,
                 low_cutoff = low_cutoff,
                 dc_threshold = dc_threshold,
                 peak_threshold_ratio = peak_threshold_ratio,
                 cycles = cycles,
                 smooth_n = smooth_n)
#
# dft_results has the periods and the model variable names
#
  print(dft_results)
  