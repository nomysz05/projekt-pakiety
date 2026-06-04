#' 2D Kalman Filter for GPS Coordinates
#' @useDynLib mojpakiet, .registration = TRUE
#' @param time Numeric vector with time in seconds
#' @param lat Numeric vector with latitude
#' @param lon Numeric vector with longitude
#' @param speed Numeric vector with speed in km/h
#' @param sigma_a Process noise parameter (default 0.8)
#' @param sigma_gps GPS measurement noise in meters (default 1.0)
#'
#' @return A data.frame with filtered coordinates: clean_lat and clean_lon
#' @examples
#' time <- c(0, 1, 2, 3)
#' lat <- c(52.2, 52.2001, 52.2002, 52.2003)
#' lon <- c(21.0, 21.0001, 21.0002, 21.0003)
#' speed <- c(10, 10, 10, 10)
#' kalman_filter_2d(time, lat, lon, speed)
#' @export
kalman_filter_2d <- function(time, lat, lon, speed, sigma_a = 0.8, sigma_gps = 1.0) {

  stopifnot(is.numeric(time), is.numeric(lat), is.numeric(lon), is.numeric(speed))
  stopifnot(is.numeric(sigma_a), is.numeric(sigma_gps))
  
  n <- length(lat)
  if (n == 0) {
    return(data.frame(clean_lat = numeric(0), clean_lon = numeric(0)))
  }
  
  stopifnot(length(time) == n, length(lon) == n, length(speed) == n)
  
  raw_result <- .Call("kalman", 
                      as.numeric(time), 
                      as.numeric(lat), 
                      as.numeric(lon), 
                      as.numeric(speed), 
                      as.numeric(sigma_a), 
                      as.numeric(sigma_gps),
                      PACKAGE = "mojpakiet")
  
  output <- data.frame(
    clean_lat = raw_result[[1]],
    clean_lon = raw_result[[2]]
  )
  
  return(output)
}

