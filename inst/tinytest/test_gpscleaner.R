library("gpscleaner")
time  <- c(0, 1, 2, 3)
lat   <- c(52.2000, 52.2100, 52.2002, 52.2003)
lon   <- c(21.0000, 21.0001, 21.0002, 21.0003)
speed <- c(13.0, 13.0, 13.0, 13.0)

res <- gps_noise_cleaner(time, lat, lon, speed)

expect_true(is.data.frame(res))
expect_equal(names(res), c("clean_lat", "clean_lon"))
expect_equal(nrow(res), 4)
expect_equal(res$clean_lat[1], 52.2000)
expect_true(abs(res$clean_lat[1]-res$clean_lat[2])<abs(lat[1]-lat[2]))
expect_error(
  gps_noise_cleaner(c(1, 2), c(95, 52), c(21, 21), c(5, 5)), 
  pattern = "Latitude values must be between -90 and 90 degrees"
)




time_na  <- c(0, 1, 2, 3, 4)
lat_na   <- c(52.2, NA, NA, 52.5, 52.6)
lon_na   <- c(21.0, NA, NA, 21.3, 21.4)
speed_na <- c(13.0, 13.0, 13.0, 13.0, 13.0)

res_interp <- interpolate_gps(time_na, lat_na, lon_na, speed_na)

expect_false(any(is.na(res_interp$interp_lat)))
expect_false(any(is.na(res_interp$interp_lon)))
expect_equal(res_interp$interp_lat[1], 52.2)
expect_equal(res_interp$interp_lat[5], 52.6)
expect_equal(res_interp$interp_lat[2], 52.3)
expect_equal(res_interp$interp_lon[2], 21.1)


time_out  <- c(0, 1, 2, 3)
lat_out   <- c(52.2000, 52.2001, 52.9000, 52.2003)
lon_out   <- c(21.0000, 21.0001, 21.0002, 21.0003)
speed_out <- c(5.0, 5.0, 5.0, 5.0)

res_out <- gps_noise_cleaner(time_out, lat_out, lon_out, speed_out)

expect_true(res_out$clean_lat[3] < 52.3)


time_interp  <- c(0, 1, 2)
lat_interp   <- c(52.2, NA, 52.2009)
lon_interp   <- c(21.0, NA, 21.0000)
speed_interp <- c(0, NA, 0)

res_spd <- interpolate_gps(time_interp, lat_interp, lon_interp, speed_interp)
expect_true(res_spd$interp_speed[2] > 0)


time_noise  <- c(0, 1, 2, 3, 4)

lat_noise   <- c(52.200, 52.202, 52.199, 52.202, 52.200)
lon_noise   <- c(21.000, 21.001, 21.002, 21.003, 21.004)
speed_noise <- c(5.0, 5.0, 5.0, 5.0, 5.0)

res_noise <- gps_noise_cleaner(time_noise, lat_noise, lon_noise, speed_noise)

oblicz_calkowity_dystans <- function(lats, lons) {
  n_pts <- length(lats)
  dystans_suma <- 0
  for (i in 1:(n_pts - 1)) {
    dlat <- lats[i+1] - lats[i]
    dlon <- lons[i+1] - lons[i]
    dystans_suma <- dystans_suma + sqrt(dlat^2 + dlon^2)
  }
  return(dystans_suma)
}

dystans_przed <- oblicz_calkowity_dystans(lat_noise, lon_noise)
dystans_po    <- oblicz_calkowity_dystans(res_noise$clean_lat, res_noise$clean_lon)

expect_true(dystans_po < dystans_przed)