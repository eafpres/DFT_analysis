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
  time_steps <- 730
  period_1 <- 3.5
  period_2 <- 13
#
  time <- seq(1:time_steps)
#
# sinusoid with period "period_1" units
#
  X1 <- sin(2 * pi * seq(1:time_steps) / period_1)
#
# sinusoid with period "period_2" units
#
  X2 <- sin(2 * pi * seq(1:time_steps) / period_2)
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
  dft_results <-
    dft_analysis(my_data, 
                 "time",
                 "Y",
                 low_cutoff = 100,
                 peak_threshold_ratio = 0.01,
                 cycles = 10)
#
# dft_results has the periods and the model variable names
#
  print(dft_results)
#
# repeat with periods of even days
#
  time_steps <- 730
  period_1 <- 7
  period_2 <- 31
#
  time <- seq(1:time_steps)
#
# sinusoid with period "period_1" units
#
  X1 <- sin(2 * pi * seq(1:time_steps) / period_1)
#
# sinusoid with period "period_2" units
#
  X2 <- sin(2 * pi * seq(1:time_steps) / period_2)
#
# linear growth with time (to demonstrate DFT still works with trend)
#  
  X4 <- 0.05 * time
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
  dft_results <-
    dft_analysis(my_data, 
                 "time",
                 "Y",
                 low_cutoff = 100,
                 high_cutoff = 2,
                 hf_threshold = 2,
                 peak_threshold_ratio = 0,
                 cycles = 10)
#
# dft_results has the periods and the model variable names
#
  print(dft_results)
  