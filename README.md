# mojpakiet

A 2D Kalman Filter for GPS trajectory smoothing, implemented in R with a C backend via `.Call` interface.

## Installation

```r
install.packages("mojpakiet")
```

## Usage

```r
library(mojpakiet)
result <- kalman_filter_2d(time, lat, lon, speed)
```

## Parameters
- `time` - numeric vector of timestamps in seconds
- `lat` - numeric vector of latitudes
- `lon` - numeric vector of longitudes  
- `speed` - numeric vector of speed in km/h
