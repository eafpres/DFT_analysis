#
  dft_analysis <- function(data, 
                           x_var, y_var, 
                           low_cutoff = 366,
                           dc_threshold = low_cutoff,
                           peak_threshold_ratio = 0.01,
                           noise_threshold = 0.05,
                           significance_threshold = 1,
                           cycles = 3,
                           smooth_n = 1) {
#
# some needed libraries if not loaded in the calling code
#
    require(dplyr)
    require(zoo)
#
# data is a data.frame or an object coercible to a data.frame
# x_var is a character string naming a column in data that represents time
# y_var is a character string naming a column in data for analysis
#    
# peak_threshold_ratio is a numeric which determines which FFT peaks are kept
# low_cutoff is the longest period for which we will keep a found peak
#    
# dc_threshold can be used to remove more very low frequency values
# to improve finding of legitimate peaks
#
# noise_threshold is used in the peak find section as a crude filter
#    
# significance_threshold is used in a simple linear model to drop
# any parameters with p values less than significance threshold
#
# cycles is the number of cycles of the longest period found
# to be used in the model fit visualization
#
# smooth_n is the number of points over which to perform optional
# moving average smoothing; this can 'sometimes' help with noisy
# data but in general will negatively impact the resulting fit
#    
# NOTE that significance_threshold only affects the model 
# which is a side effect and NOT returned by the function
#    
# similarly, the visualizations are side effects and NOT 
# returned by the function
#
# the function assumes the units of variable x_var are the time units
#
    data <- as.data.frame(data)
    period <- 
      mean(data[2:nrow(data), x_var] - 
             data[1:(nrow(data) - 1), x_var])
    fft_data <- data[, y_var]
    fft_raw <- fft(fft_data)
#
# f.Nyquist is the maximum sampling frequency given the 1/day period
#
    f.Nyquist <- 0.5 * (1 / period)
    freq <- 
      c(0,
        f.Nyquist * c(seq(length(fft_data) / 2), 
                    -rev(seq(length(fft_data) / 2)))/(length(fft_data) / 2))
#
# a property of the discrete FFT is that it reflects
# the data to negative frequencies, so we strip those off
#
    real_points <- which(freq >= 0)
    freq <- freq[real_points]
    fft_raw <- fft_raw[real_points]
#
# now we do smoothing over smooth_n periods, centered, to eliminate
# noise to avoid duplicating peaks in adjacent samples
#
    FFT_processed <- 
      fft_raw %>%
      as.data.frame() %>%
      rename(fft_raw = names(.)[1]) %>%
      mutate(smooth_fft = rollmean(Mod(fft_raw), k = smooth_n, fill = "extend", align = "center"))
#
# find the first peak then get some properties
# of all the remaining data
#
    found_start_first_peak <- FALSE
    index <- 1
    while (!found_start_first_peak & index < length(fft_raw)) {
      index <- index + 1 
      if (FFT_processed$smooth_fft[index] > (1 + noise_threshold) * FFT_processed$smooth_fft[index - 1]) {
        found_start_first_peak <- TRUE
        start_first_peak <- index - 1
      }
    }
    found_end_last_peak <- FALSE
    index <- nrow(FFT_processed)
    while (!found_end_last_peak & index > 1) {
      index <- index - 1
      if (FFT_processed$smooth_fft[index] > (1 + noise_threshold) * FFT_processed$smooth_fft[index + 1]) {
        found_end_last_peak <- TRUE
        end_last_peak <- index + 1
      }
    }
    max_retained_fft <-
      max(FFT_processed$smooth_fft[start_first_peak:nrow(FFT_processed)])
    last_retained_fft_index <-
      min(nrow(FFT_processed), end_last_peak + 3)
    min_retained_fft <-
      min(FFT_processed$smooth_fft[start_first_peak:nrow(FFT_processed)])
    first_retained_fft_index <-
      max(1, start_first_peak - 3)
#    
# visualize the range selected with peaks
#
    FFT_processed %>%
      mutate(X = row_number()) %>%
      filter(X >= first_retained_fft_index) %>%
      filter(X <= last_retained_fft_index) %>%
      ggplot(aes(x = X, y = smooth_fft)) +
      geom_point(color = "red", size = 2) +
      ylim(c(0.99 * min_retained_fft, 1.01 * max_retained_fft)) +
      geom_smooth(color = "black", size = 0.25, linetype = 2, se = FALSE, span = 0.2)
#
# calculate the first derivative on the left and right of each point excluding the ends
#
    derivatives <- 
      data.frame(left_derivative = numeric(length = nrow(FFT_processed)),
                 right_derivative = numeric(length = nrow(FFT_processed)))
    for (index in 1:nrow(FFT_processed)) {
      if (index == 1) {
        derivatives$left_derivative[index] <- NA
        derivatives$right_derivative[index] <-
          (FFT_processed$smooth_fft[index + 1] - FFT_processed$smooth_fft[index]) / period
      } else if (index == length(fft_raw)) {
        derivatives$left_derivative[index] <-
          (FFT_processed$smooth_fft[index] - FFT_processed$smooth_fft[index - 1]) / period
        derivatives$right_derivative[index] <- NA
      } else {
        derivatives$left_derivative[index] <-
          (FFT_processed$smooth_fft[index] - FFT_processed$smooth_fft[index - 1]) / period
        derivatives$right_derivative[index] <-
          (FFT_processed$smooth_fft[index + 1] - FFT_processed$smooth_fft[index]) / period
      }
    }
    FFT_processed <-
      FFT_processed %>%
      mutate(fft = Mod(fft_raw)) %>%
      mutate(left_derivative = derivatives$left_derivative) %>%
      mutate(right_derivative = derivatives$right_derivative) %>%
      mutate(PEAK = case_when(
        left_derivative > 0 & right_derivative <= 0 ~ TRUE,
        left_derivative >= 0 & right_derivative < 0 ~ TRUE,
        TRUE ~ FALSE
      )) %>%
      mutate(frequency = freq)
#
# the peak_threshold_ratio is just a cutoff to find peaks in the
# power spectrum vs. frequency
#
# we look only at points from first_peak or the value determined
# as an offset from low_cutoff using dc_threshold, i.e. if 
# low_cutoff is 365 (longest period is 365 time steps) and 
# dc_threshold is default, then dc_threshold = low_cutoff
# and therefore the first data point we look at is at f > 1 / 365 or first_peak,
# whichever is greater, but we look back 1 step from first_peak in the 
# frequency spectrum to get the start of the peak
#
    label_points <- 
      which(FFT_processed$PEAK)
    if (min(low_cutoff, dc_threshold) > 0) {
      label_points <- 
        label_points[which(FFT_processed$frequency[label_points] > 
                             (1 / min(low_cutoff * period, dc_threshold * period)))]
    }
#
# there can be multiple samples in a peak, so we need to 
# select just the actual peak point; the following logic
# takes care of that by looking for closely spaced samples
# and choosing the one with the maximum power
#
    final_label_points <- integer()
    keep_points <- integer()
    if (length(label_points) > 1) {
      for (index in 1:(length(label_points))) {
        if (FFT_processed$smooth_fft[label_points[index]] > 
            (1 + noise_threshold) * FFT_processed$smooth_fft[label_points[index] - 1] &
            FFT_processed$smooth_fft[label_points[index] + 1] < 
            (1 - noise_threshold) * FFT_processed$smooth_fft[label_points[index]]) {
          keep_points <- c(keep_points, index)
        }
      }
      keep_points <- label_points[keep_points]
    } else {
      keep_points <- label_points
    }
    if (length(keep_points) > 0) {
      final_label_points <- keep_points
    }
#
# if we didn't find anything, return NULL
#
    if (length(final_label_points) == 0) {
      dft_results <- NULL
      return(dft_results)
    }
#
# we use all the frequencies above the threshold
# to fit the seasonal data to sine an cosine series
#
    fit_frequencies <- freq[final_label_points]
    seasonal_data_model <- 
      as.data.frame(matrix(nrow = length(fft_data), 
                           ncol = 1 + 2 * length(fit_frequencies)))
    seasonal_data_model[, 1] <- data[, x_var]
    seasonal_data_model[, 1] <- seasonal_data_model[, 1]
    for (i in 1:length(fit_frequencies)) {
      seasonal_data_model[, 2 * i] <- 
        sin(2 * pi * fit_frequencies[i] * seasonal_data_model[, 1])
      seasonal_data_model[, 2 * i + 1] <-
        cos(2 * pi * fit_frequencies[i] * seasonal_data_model[, 1])
    }
    seasonal_data_model <- cbind(seasonal_data_model, 
                                 data[, y_var])
    sdm_colnames <- character()
    for (i in 1:length(fit_frequencies)) {
      sdm_colnames <- c(sdm_colnames,
                        c(paste0("sin_", ceiling(100 / fit_frequencies[i]) / 100),
                        paste0("cos_", ceiling(100 / fit_frequencies[i]) / 100)))
    }
    sdm_colnames <- c("time", sdm_colnames, "fft_data")
    colnames(seasonal_data_model) <- sdm_colnames
#
# fit over the requested number cycles of the longest period
# note that depending on the data and filtering settings
# this parameter can have a significant impact when using 
# the results over long numbers of time steps
#
    fit_range <- cycles / min(fit_frequencies) / period
    prelim_model <- lm(fft_data ~ ., data = seasonal_data_model[1:fit_range, ])
#
# build a formula to fit using only sigificant factors
#
    threshold <- significance_threshold
    signif_form <- 
      as.formula(paste("fft_data ~ ",
                       paste(names(which((summary(prelim_model)$coefficients[
                         2:(nrow(summary(prelim_model)$coefficients)), 4] < 
                           threshold) == TRUE)), 
                         collapse = "+")))
#
# refit with only the signficant factors
#
    final_model <- lm(signif_form, data = seasonal_data_model[1:fit_range, ])
#
# convert the labels to period in days
#
    label_frequencies <- freq[final_label_points]
    label_values <- as.character(ceiling(100 / label_frequencies) / 100)
    plot(x = freq[max(1, min(final_label_points) - 3):
                    min(length(freq), 2 * max(final_label_points))], 
         y = FFT_processed$smooth_fft[max(1, min(final_label_points) - 3):
                                 min(length(freq), min(length(freq), 2 * max(final_label_points)))], 
         type = "l", 
         xaxt = "n",
         yaxt = "n",
         xlab = "period (time units)",
         ylab = "")
    title(ylab = "relative power in frequency", mgp = c(2, 1, 0))
    par(mgp = c(3, 1, 0))
    text(x = freq[final_label_points], 
         y = FFT_processed$smooth_fft[final_label_points], 
         label_values, 
         pos = 4)
    low_cutoff_x <-  min(min(freq), min(1 / (low_cutoff * period), 1 / (dc_threshold * period)))
    abline(v = low_cutoff_x, col = "red")
    text(x = low_cutoff_x,
         y = 0.95 * max(FFT_processed$fft),
         "low cutoff",
         pos = 4, 
         col = "red")
#
# create the x-axis in period vs. frequency
#
    axis_points <- 
      seq(0, max(freq[1:min(length(freq), 2 * max(final_label_points))]), 
          length = 10)
    axis(1, at = axis_points,
         label = c("", 
                   ceiling(10 / axis_points[2:length(axis_points)]) / 10))
    axis(2, at = 
           seq(min(FFT_processed$smooth_fft[max(1, min(final_label_points) - 3):
                     min(length(freq), min(length(freq), 2 * max(final_label_points)))]),
               max(FFT_processed$smooth_fft[max(1, min(final_label_points) - 3):
                     min(length(freq), min(length(freq), 2 * max(final_label_points)))]), 
               length = 10),
         label = rep("", 10))
#
# plot the prediction
#
    periods <- 1 / label_frequencies / period
    plot_time_units <- ceiling(cycles * max(periods))
#
# generate plot data from the model
#
    plot_data_names <-
      attr(summary(final_model)[["terms"]], "term.labels")
    plot_data <- 
      as.data.frame(matrix(0, 
                           ncol = length(plot_data_names),
                           nrow = nrow(data)),
                    stringsAsFactors = FALSE)
    colnames(plot_data) <- 
      plot_data_names
    plot_data[, x_var] <- data[, x_var]
    for (i in which(!(colnames(plot_data) == x_var))) {
      plot_data[, i] <-
        ifelse(rep(strsplit(plot_data_names[i], "_")[[1]][1] == "sin", nrow(plot_data)),
               sin(2 * pi * plot_data[, x_var] / as.numeric(strsplit(plot_data_names[i], "_")[[1]][2])),
               cos(2 * pi * plot_data[, x_var] / as.numeric(strsplit(plot_data_names[i], "_")[[1]][2])))
    }
    plot(data[1:plot_time_units, x_var], data[1:plot_time_units, y_var], type = "l")
    points(data[1:plot_time_units, x_var], predict(final_model, newdata = plot_data)[1:plot_time_units], col = "red")
    lines(data[1:plot_time_units, x_var], predict(final_model, newdata = plot_data)[1:plot_time_units], col = "red")  
#
    dft_results <- list(periods = 1 / label_frequencies,
                        labels = sdm_colnames)
    return(dft_results)
  }
  