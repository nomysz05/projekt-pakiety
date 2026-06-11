# 2D Kalman Filter for GPS Coordinates

Filters noisy GPS coordinates using a 2D Kalman filter implemented in C.
The filter uses position and velocity as state variables, and
automatically detects and rejects outliers (points more than 30 meters
from the predicted position). Supports automatic conversion from
semicircles to decimal degrees.

## Usage

``` r
gps_noise_cleaner(time, lat, lon, speed, sigma_a = 0.8, sigma_gps = 1)
```

## Arguments

- time:

  Numeric vector with time in seconds

- lat:

  Numeric vector with latitude in decimal degrees or semicircles

- lon:

  Numeric vector with longitude in decimal degrees or semicircles

- speed:

  Numeric vector with speed in m/s

- sigma_a:

  Process noise parameter - controls how much acceleration is expected
  between steps (default 0.8)

- sigma_gps:

  GPS measurement noise in meters - higher values trust the model more
  than the GPS signal (default 1.0)

## Value

A data.frame with two columns:

- clean_lat:

  Filtered latitude in decimal degrees

- clean_lon:

  Filtered longitude in decimal degrees

## Examples

``` r
time <- c(0, 1, 2, 3)
lat <- c(52.2, 52.2001, 52.2002, 52.2003)
lon <- c(21.0, 21.0001, 21.0002, 21.0003)
speed <- c(13.0,13.0,13.0,13.0)
gps_noise_cleaner(time, lat, lon, speed)
#> --- URUCHAMIAM CZYSTY FILTR KALMANA (WSTĘPNA FILTRACJA W R) --- 
#>   clean_lat clean_lon
#> 1   52.2000   21.0000
#> 2   52.2001   21.0001
#> 3   52.2002   21.0002
#> 4   52.2003   21.0003
```
