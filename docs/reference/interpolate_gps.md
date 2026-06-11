# Linear Interpolation of Missing GPS Coordinates with Speed Recalculation

Interpolates missing latitude and longitude values using linear
interpolation based on timestamps. For points where coordinates were
missing, speed is recalculated from the newly interpolated positions.
Original speed values are preserved where coordinates were valid.

## Usage

``` r
interpolate_gps(time, lat, lon, speed)
```

## Arguments

- time:

  Numeric vector with time in seconds

- lat:

  Numeric vector with latitude (may contain NA)

- lon:

  Numeric vector with longitude (may contain NA)

- speed:

  Numeric vector with speed in m/s (may contain NA)

## Value

A data.frame with three columns:

- interp_lat:

  Interpolated latitude in decimal degrees

- interp_lon:

  Interpolated longitude in decimal degrees

- interp_speed:

  Original speed where coordinates were valid, recalculated from
  interpolated positions where coordinates were missing

## Examples

``` r
time <- c(0, 1, 2, 3, 4)
lat <- c(52.2, NA, NA, 52.2003, 52.2004)
lon <- c(21.0, NA, NA, 21.0003, 21.0004)
speed <- c(13.0,13.0,13.0,13.0,13.0)
interpolate_gps(time, lat, lon, speed)
#>   interp_lat interp_lon interp_speed
#> 1    52.2000    21.0000     13.00000
#> 2    52.2001    21.0001     13.03201
#> 3    52.2002    21.0002     13.03201
#> 4    52.2003    21.0003     13.00000
#> 5    52.2004    21.0004     13.00000
```
