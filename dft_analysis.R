#
  dft_analysis <- function(data, 
                           x_var, y_var, 
                           peak_threshold_ratio,
                           low_cutoff = 366,
                           high_cutoff = 1,
                           dc_threshold = 5 * high_cutoff,
                           noise_threshold = 0.5,
                           significance_threshold = 1,
                           cycles = 30) {
#
# data is a data.frame or an object coercible to a data.frame
# x_var is a character string naming a column in data that represents time
# y_var is a character string naming a column in data for analysis
#    
# peak_threshold_ration is a numeric which determines which FFT peaks are kept
# low_cutoff is the longest period for which we will keep a found peak
# high_cutoff is the shortest period for which we will keep a found peak
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
    noise_threshold <- 0.5
    for (i in 2:(length(Mod(fft_raw) - 1))) {
      if (Mod(fft_raw)[i] > (1 + noise_threshold) * Mod(fft_raw)[i - 1]) {
        provisional_peak <- i
        for (j in (provisional_peak + 1):(length(Mod(fft_raw) - 1)) {
          if  (Mod(fft_raw)[j] > (1 + noise_threshold) * Mod(fft_raw)[j - 1]) {
            provisional_peak <- j
          } else {
            first_peak <- provisional_peak
          }
          break()
        }
      }
    }
#
# the peak_threshold_ratio is just a cutoff to find peaks in the
# power spectrum vs. frequency
#
    label_points <- 
      which(Mod(fft_raw)[dc_threshold:length(fft_raw)] > 
              peak_threshold_ratio * max(Mod(fft_raw)[dc_threshold:length(fft_raw)]))
#
# there can be multiple samples in a peak, so we need to 
# select just the actual peak point; the following logic
# takes care of that by looking for closely spaced samples
# and choosing the one with the maximum power
#
    final_label_points <- integer()
    keep <- 1
    for (index in 2:length(label_points)) {
      if (abs(label_points[index] - label_points[index - 1]) < 2) {
        if (Mod(fft_raw[label_points[index]]) >
            Mod(fft_raw[label_points[keep]])) {
          keep <- index
        }
        if (index == length(label_points)) {
          final_label_points <- c(final_label_points, label_points[keep])
        }
      } else {
        final_label_points <- c(final_label_points, label_points[keep])
        keep <- index
        if (index == length(label_points)) {
          final_label_points <- c(final_label_points, label_points[keep])
        }
      }
    }
    if (length(final_label_points) == 0) {
      final_label_points <- label_points[keep]
    }
#
# this function assumes we don't want the lowest frequency
#
    final_label_points <- 
      final_label_points[!(final_label_points == 1)]
    final_label_points <-
      final_label_points[!(1/(freq[final_label_points]) > low_cutoff)]
    final_label_points <-
      final_label_points[!(1/(freq[final_label_points]) < high_cutoff)]
    if (length(final_label_points) == 0) {
      dft_results <- NULL
      return(dft_results)
    }
#
# we use all the frequencies above thethreshold
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
# fit over the requested cycles of the longest period
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
    plot(x = freq[dc_threshold:length(freq)], 
         y = Mod(fft_raw[dc_threshold:length(freq)]), 
         type = "l", 
         xaxt = "n",
         yaxt = "n",
         xlab = "period (time steps)",
         ylab = "")
    title(ylab = "relative power in frequency", mgp = c(2, 1, 0))
    par(mgp = c(3, 1, 0))
    text(x = freq[final_label_points], 
         y = Mod(fft_raw[final_label_points]), 
         label_values, pos = 4)
    abline(v = 1 / low_cutoff, col = "red")
    abline(v = 1 / high_cutoff, col = "red")
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
  