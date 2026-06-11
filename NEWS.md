# gpscleaner 1.0.0

* **Package ready:** First official release of the fully functional package.
* **Core functions:** * `gps_noise_cleaner()` – smoothing telemetry noise using a 2D Kalman Filter.
  * `interpolate_gps()` – reconstructing missing points in a GPS track (`NA` gaps) along with automatic speed recalculation.
* **Documentation update:** The project website has been updated with new usage examples and plots verifying the algorithms' performance (including a U-shaped trajectory test case with sharp turns and an injected anomaly).
