# geoconv - convert between different coordinate systems

[crates.io](https://crates.io/crates/geoconv)
[docs.rs](https://docs.rs/geoconv/latest/geoconv/)

`geoconv` implements support for converting between some basic coordinate
systems. This package contains support for WGS84
Latitude/Longitude/Elevation (LLE) geodetic coordinates, Earth-centered,
Earth-fixed (XYZ), and local tangent plane (both East/North/Up (ENU)
and Azimuth/Elevation/Range (AER)) systems.

> [!WARNING]
> **This is also absolutely not ready to be used for navigational purposes**.
> Please do not use this library in any situation that may cause harm to life or property.

This package is particularly useful if you know your location on Earth,
and want to geolocate an object that you can observe, or if you know your
location on Earth, and an object's location on Earth and want to determine
your relative positions.

This also includes Haversine-distance using the circumference of the earth
for approximate great-circle distance between two Lat/Lon points, but only
on the "Earth's surface". In some cases (long distances) this is a way
better approach for distance. If both points are line-of-sight, you likely
care more about using the local tangent plane to work out the Range
to the target in 3-space by converting both lat/lons to Earth-North-Up or
an Azimuth/Elevation/Range.

# `no_std`

This module is `no_std`, using `libm` for the underlying math.
