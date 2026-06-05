#' 2D Kalman Filter for GPS Coordinates
#'
#' Filters noisy GPS coordinates using a 2D Kalman filter implemented in C.
#' The filter uses position and velocity as state variables, and automatically
#' detects and rejects outliers (points more than 30 meters from the predicted
#' position). Supports automatic conversion from semicircles to decimal degrees.
#'
#' @param time Numeric vector with time in seconds
#' @param lat Numeric vector with latitude in decimal degrees or semicircles
#' @param lon Numeric vector with longitude in decimal degrees or semicircles
#' @param speed Numeric vector with speed in m/s
#' @param sigma_a Process noise parameter - controls how much acceleration
#'   is expected between steps (default 0.8)
#' @param sigma_gps GPS measurement noise in meters - higher values trust
#'   the model more than the GPS signal (default 1.0)
#'
#' @return A data.frame with two columns:
#' \describe{
#'   \item{clean_lat}{Filtered latitude in decimal degrees}
#'   \item{clean_lon}{Filtered longitude in decimal degrees}
#' }
#' @examples
#' time <- c(0, 1, 2, 3)
#' lat <- c(52.2, 52.2001, 52.2002, 52.2003)
#' lon <- c(21.0, 21.0001, 21.0002, 21.0003)
#' speed <- c(13.0,13.0,13.0,13.0)
#' kalman_filter_2d(time, lat, lon, speed)
#' @export
#' @useDynLib gpscleaner, .registration = TRUE
gps_noise_cleaner <- function(time, lat, lon, speed, sigma_a = 0.8, sigma_gps = 1.0) {
  stopifnot(is.numeric(time), is.numeric(lat), is.numeric(lon), is.numeric(speed))
  stopifnot(is.numeric(sigma_a), is.numeric(sigma_gps))

  n <- length(lat)
  if (n == 0) {
    return(data.frame(clean_lat = numeric(0), clean_lon = numeric(0)))
  }

  stopifnot(length(time) == n, length(lon) == n, length(speed) == n)

  if (max(abs(lat), na.rm = TRUE) > 180 || max(abs(lon), na.rm = TRUE) > 180) {
    message("Detected semicircles format - converting to decimal degrees")
    lat <- lat * (180 / 2^31)
    lon <- lon * (180 / 2^31)
  }

  if (any(abs(lat) > 90, na.rm = TRUE))
    stop("Latitude values must be between -90 and 90 degrees")
  if (any(abs(lon) > 180, na.rm = TRUE))
    stop("Longitude values must be between -180 and 180 degrees")
  if (any(speed < 0, na.rm = TRUE)) {
    warning("Negative speed values detected - taking absolute values")
    speed <- abs(speed)
  }

  # valid_idx musi byc zdefiniowane przed uyciem
  valid_idx <- which(!is.na(time) & !is.na(lat) & !is.na(lon) & !is.na(speed))

  if (any(diff(time[valid_idx]) <= 0))
    warning("Time vector is not strictly increasing - results may be unreliable")

  if (length(valid_idx) < 2) {
    return(data.frame(clean_lat = rep(NA_real_, n),
                      clean_lon = rep(NA_real_, n)))
  }

  raw_result <- .Call("kalman",
                      as.numeric(time[valid_idx]),
                      as.numeric(lat[valid_idx]),
                      as.numeric(lon[valid_idx]),
                      as.numeric(speed[valid_idx]),
                      as.numeric(sigma_a),
                      as.numeric(sigma_gps),
                      PACKAGE = "gpscleaner")

  clean_lat <- rep(NA_real_, n)
  clean_lon <- rep(NA_real_, n)
  clean_lat[valid_idx] <- raw_result[[1]]
  clean_lon[valid_idx] <- raw_result[[2]]

  return(data.frame(clean_lat = clean_lat, clean_lon = clean_lon))
}
#' Linear Interpolation of Missing GPS Coordinates with Speed Recalculation
#'
#' Interpolates missing latitude and longitude values using linear interpolation
#' based on timestamps. For points where coordinates were missing, speed is
#' recalculated from the newly interpolated positions. Original speed values
#' are preserved where coordinates were valid.
#'
#' @param time Numeric vector with time in seconds
#' @param lat Numeric vector with latitude (may contain NA)
#' @param lon Numeric vector with longitude (may contain NA)
#' @param speed Numeric vector with speed in m/s (may contain NA)
#'
#' @return A data.frame with three columns:
#' \describe{
#'   \item{interp_lat}{Interpolated latitude in decimal degrees}
#'   \item{interp_lon}{Interpolated longitude in decimal degrees}
#'   \item{interp_speed}{Original speed where coordinates were valid,
#'   recalculated from interpolated positions where coordinates were missing}
#' }
#' @examples
#' time <- c(0, 1, 2, 3, 4)
#' lat <- c(52.2, NA, NA, 52.2003, 52.2004)
#' lon <- c(21.0, NA, NA, 21.0003, 21.0004)
#' speed <- c(13.0,13.0,13.0,13.0,13.0)
#' interpolate_gps(time, lat, lon, speed)
#' @export
interpolate_gps <- function(time, lat, lon, speed) {
  stopifnot(is.numeric(time), is.numeric(lat), 
            is.numeric(lon), is.numeric(speed))
  n <- length(lat)
  stopifnot(length(time) == n, length(lon) == n, length(speed) == n)

  if (n == 0)
    return(data.frame(interp_lat   = numeric(0),
                      interp_lon   = numeric(0),
                      interp_speed = numeric(0)))

  if (sum(!is.na(lat)) < 2 || sum(!is.na(lon)) < 2) {
    warning("Not enough valid points for interpolation")
    return(data.frame(interp_lat   = lat,
                      interp_lon   = lon,
                      interp_speed = speed))
  }
  if (max(abs(lat), na.rm = TRUE) > 180 || max(abs(lon), na.rm = TRUE) > 180) {
    message("Detected semicircles format - converting to decimal degrees")
    lat <- lat * (180 / 2^31)
    lon <- lon * (180 / 2^31)
  }

  raw_result <- .Call("interpolate_gps",
                      as.numeric(time),
                      as.numeric(lat),
                      as.numeric(lon),
                      as.numeric(speed),
                      PACKAGE = "gpscleaner")

  return(data.frame(interp_lat   = raw_result[[1]],
                    interp_lon   = raw_result[[2]],
                    interp_speed = raw_result[[3]]))
}

