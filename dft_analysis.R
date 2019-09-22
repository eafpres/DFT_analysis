#
  dft_analysis <- function(data, 
                           x_var, y_var, 
                           peak_threshold_ratio,
                           low_cutoff = 366,
                           dc_threshold = low_cutoff - 5,
                           noise_threshold = 0.05,
                           significance_threshold = 1,
                           cycles = 3) {
#
# data is a data.frame or an object coercible to a data.frame
# x_var is a character string naming a column in data that represents time
# y_var is a character string naming a column in data for analysis
#    
# peak_threshold_ration is a numeric which determines which FFT peaks are kept
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
      f.Nyquist * c(seq(length(fft_data) / 2), 
                    -rev(seq(length(fft_data) / 2)))/(length(fft_data) / 2)
#
# a property of the discrete FFT is that it reflects
# the data to negative frequencies, so we strip those off
#
    real_points <- which(freq >= 0)
    freq <- freq[real_points]
    fft_raw <- fft_raw[real_points]
#
# find the first point where the value increases
# this is used to reject the typically huge values at long periods (DC)
#
    first_peak <- 1
    provisional_peak <- 0
    for (i in 2:(length(Mod(fft_raw) - 1))) {
      if (Mod(fft_raw)[i] > (1 + noise_threshold) * Mod(fft_raw)[i - 1]) {
        provisional_peak <- i
        for (j in (provisional_peak + 1):(length(Mod(fft_raw) - 1))) {
          if  (Mod(fft_raw)[j] > (1 - noise_threshold) * Mod(fft_raw)[j - 1]) {
            provisional_peak <- j
          } else {
            first_peak <- provisional_peak
            cat("exiting inner loop\n")
            break()
          }
        }
        if (first_peak > 1) {
          cat("exiting 2rd inner loop\n")
          break()
        }
      }
      if (first_peak > 1) {
        cat("exiting outer loop\n")
        break()
      }
    }
#
# the peak_threshold_ratio is just a cutoff to find peaks in the
# power spectrum vs. frequency
#
# we look only at points from first_peak or the value determined
# as an offset from low_cutoff using dc_threshold, i.e. if 
# low_cutoff is 365 (longest period is 365 time steps) and 
# dc_threshold is default, then dc_threshold = 365 - 5 = 360
# and therefore the first data point we look at is 5 or first_peak,
# whichever is greater, but we look back 3 steps from first_peak in the 
# frequency spectrum to get the start of the peak
#
    labels_start <- max((low_cutoff - dc_threshold), (first_peak - 3))
    label_points <- 
      which(Mod(fft_raw)[labels_start:length(fft_raw)] > 
              peak_threshold_ratio * max(Mod(fft_raw)[labels_start:length(fft_raw)])) + labels_start
#
# there can be multiple samples in a peak, so we need to 
# select just the actual peak point; the following logic
# takes care of that by looking for closely spaced samples
# and choosing the one with the maximum power
#
    final_label_points <- integer()
#
# find start of first real peak
#
    first_label_start <- as.integer(NA)
    for (possible_first_label_start in 1:(length(label_points) - 1)) {
      if (Mod(fft_raw)[label_points[possible_first_label_start + 1]] > 
          (1 + noise_threshold) * Mod(fft_raw)[label_points[possible_first_label_start]]) {
        first_label_start <- possible_first_label_start
        if (!(is.na(first_label_start))) {
          break()
        }
      }
    }
    final_label_points <- label_points[first_label_start:length(label_points)]
    keep_points <- integer()
    for (index in 2:(length(final_label_points) - 1)) {
      if (Mod(fft_raw)[index] > (1 + noise_threshold) * Mod(fft_raw)[index - 1] &
          Mod(fft_raw)[index + 1] < (1 - noise_threshold) * Mod(fft_raw)[index]) {
        keep_points <- c(keep_points, index)
      }
    }
    if (length(keep_points) > 0) {
      final_label_points <- final_label_points[which(final_label_points %in% keep_points)]
    } else {
      final_label_points <-
        final_label_points[c(1, length(final_label_points))]
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
    plot(x = freq[1:length(freq)], 
         y = Mod(fft_raw[1:length(freq)]), 
         type = "l", 
         xaxt = "n",
         yaxt = "n",
         xlab = "period (time units)",
         ylab = "")
    title(ylab = "relative power in frequency", mgp = c(2, 1, 0))
    par(mgp = c(3, 1, 0))
    text(x = freq[final_label_points], 
         y = Mod(fft_raw[final_label_points]), 
         label_values, 
         pos = 4)
    abline(v = freq[low_cutoff - dc_threshold], col = "red")
    text(x = freq[low_cutoff - dc_threshold],
         y = 0.95 * max(Mod(fft_raw[1:length(freq)])),
         "low cutoff",
         pos = 4, 
         col = "red")
#
# create the x-axis in period vs. frequency
#
    axis_points <- seq(0, ceiling(max(freq)), length = 10)
    axis(1, at = axis_points,
         label = c("", 
                   ceiling(10 / axis_points[2:length(axis_points)]) / 10))
    axis(2, at = 
           seq(min(Mod(fft_raw)), 
               max(Mod(fft_raw)), length = 10),
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
  