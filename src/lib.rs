// {{{ Copyright (c) Paul R. Tagliamonte <paultag@gmail.com>, 2023
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE. }}}

#![forbid(unsafe_code)]

//! This crate implements support for converting between geodetic,
//! relative and cartesian coordinates.

mod wgs84;

pub use crate::wgs84::{WGS84};

/// Meters represent the SI unit of measure, Meter
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Meters(pub f64);

/// Degrees is an angular measure that ranges from 0 to 360 (sometimes negative)
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Degrees(pub f64);

impl Degrees {
    // to_radians will convert the Degrees to Radians.
    pub fn to_radians(self) -> Radians {
        return Radians(std::f64::consts::PI / 180.0 * self.0);
    }
}

/// Radians is an angular measure that ranges from 0 to 2π (𝜏).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Radians(pub f64);

impl Radians {
    /// to_radians will convert the Radians to Degrees.
    pub fn to_degrees(self) -> Degrees {
        return Degrees(180.0 / std::f64::consts::PI * self.0)
    }
}

/// LLE or Latitude, Longitude, Elevation, is a location somewhere around Earth.
/// Elevation is the number of Meters above the CoordinateSystem's ellipsoid,
/// *not* the altitude above or below the surface.
///
/// This is a *absolute* and *angular* measure.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LLE {
    pub latitude: Degrees,
    pub longitude: Degrees,
    pub elevation: Meters,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MustNotHaveElevationError;

const EARTH_RADIUS: Meters = Meters(6371000.0);

impl LLE {

    /// haversine_distance will compute the great-circle distance between two
    /// points. Since implementing this in 3-D is fairly tricky, this is the
    /// *2D* distance between two Lat/Lon, not the 3D great-circle distance
    /// between two Lat/Lon.
    ///
    /// As such, passing in any Elevation other than '0.0' will result in an
    /// error ('MustNotHaveElevationError').
    pub fn haversine_distance(self, other: LLE) -> Result<Meters, MustNotHaveElevationError> {
        if self.elevation.0 != 0.0 || other.elevation.0 != 0.0 {
            return Err(MustNotHaveElevationError)
        }

        let self_lat = self.latitude.to_radians().0;
        let self_lon = self.longitude.to_radians().0;
        let other_lat = other.latitude.to_radians().0;
        let other_lon = other.longitude.to_radians().0;

        let delta_lat = self_lat - other_lat;
        let delta_lon = self_lon - other_lon;

        let a = (delta_lat / 2.0).sin().powf(2.0) + self_lat.cos() * other_lat.cos() * (delta_lon / 2.0).sin().powf(2.0);
        let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());

        Ok(Meters(EARTH_RADIUS.0 * c))
    }
}

/// XYZ is the earth-centric XYZ point system. XYZ locations can be turned into
/// points on the Earth's ellipsoid, but plotted using cartesian coordinates
/// relative to Earth, rather than angular LLE measurements.
///
/// This is a *absolute* and *cartesian* measure.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct XYZ {
    pub x: Meters,
    pub y: Meters,
    pub z: Meters,
}

/// AER represents an Azimuth, Elevation, and Range measurement.
///
/// Azimuth/Elevation (or Az/El) is a common way of locating objects measured
/// at a specific location on Earth in the local tangent plane
/// (for instance, from a RADAR)
///
/// This is a *relative* and *angular* measure.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct AER {
    pub azimuth: Degrees,
    pub elevation: Degrees,
    pub range: Meters,
}

impl AER {
    /// to_enu will convert to the Earth-North-Up cartesian local tangent plane
    /// system.
    pub fn to_enu(self) -> ENU {
        let az_rad = self.azimuth.to_radians();
        let el_rad = self.elevation.to_radians();
        let r = Meters(self.range.0 * el_rad.0.cos());
        ENU{
            east: Meters(r.0 * az_rad.0.sin()),
            north: Meters(r.0 * az_rad.0.cos()),
            up: Meters(self.range.0 * el_rad.0.sin()),
        }
    }
}

/// ENU is East, North, Up in Meters. These measures are in the local
/// tangent plane, which is to say, increasing "North" will get further and
/// further away from the Earth's surface (well, unless there's a mountain
/// range in front of you).
///
/// This is a *relative* and *cartesian* measure.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ENU {
    pub east: Meters,
    pub north: Meters,
    pub up: Meters,
}

impl ENU {
    /// to_aer will convert to the Azimuth-Elevation-Range angular local tangent
    /// plane system.
    pub fn to_aer(self) -> AER {
        let r = (self.east.0 * self.east.0 + self.north.0 * self.north.0).sqrt();
        let tau = std::f64::consts::PI * 2.0;
        AER{
            azimuth: Radians(self.east.0.atan2(self.north.0) % tau).to_degrees(),
            elevation: Radians(self.up.0.atan2(r)).to_degrees(),
            range: Meters((r * r + self.up.0 * self.up.0).sqrt()),
        }
    }
}

/// CoordinateSystem is a trait to enable converstion between locations
/// (usually Latitude and Longitude, in the form of LLA objects) to absolute
/// points in space and vice versa.
///
/// Different systems have different measurements of Earth's surface, and
/// a Latitude / Longitude must be understood within its CoordinateSystem,
/// or significant errors can be introduced.
pub trait CoordinateSystem {
    /// lle_to_xyz will convert between Lat/Lon and ECEF style coordiantes.
    fn lle_to_xyz(geo: LLE) -> XYZ;

    /// xyz_to_lle will convert between ECEF and Lat/Lon style coordinates.
    fn xyz_to_lle(c: XYZ) -> LLE;

    /// xyz_to_enu will convert an ECEF coordinate as observed at a Lat/Lon
    /// into a local tangent plane East/North/Up coordinate.
    fn xyz_to_enu(refr: LLE, c: XYZ) -> ENU;

    /// enu_to_xyz will convert an East/North/Up coordinate as observed at a
    /// Lat/Lon into a ECEF coordinate.
    fn enu_to_xyz(refr: LLE, c: ENU) -> XYZ;

    /// lle_to_enu will convert a Lat/Lon as seen by another Lat/Lon into a
    /// local tangent plane East/North/Up coordinate.
    fn lle_to_enu(refr: LLE, geo: LLE) -> ENU;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degrees_to_radians() {
        let result = Degrees(180.0).to_radians();
        assert_eq!(result.0, std::f64::consts::PI);
    }

    #[test]
    fn radians_to_degrees() {
        let result = Radians(std::f64::consts::PI).to_degrees();
        assert_eq!(result.0, 180.0);
    }

    macro_rules! assert_in_eps {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d && $y - $x < $d) { panic!(); }
        }
    }

    #[test]
    fn haversine_known() {
        let pointa = LLE{
            latitude: Degrees(22.55),
            longitude: Degrees(43.12),
            elevation: Meters(0.0),
        };
        let pointb = LLE{
            latitude: Degrees(13.45),
            longitude: Degrees(100.28),
            elevation: Meters(0.0),
        };
        let result = pointa.haversine_distance(pointb).expect("distance");
        assert_in_eps!(6094544.408786774, result.0, 1e-6);

        let pointa = LLE{
            latitude: Degrees(51.510357),
            longitude: Degrees(-0.116773),
            elevation: Meters(0.0),
        };
        let pointb = LLE{
            latitude: Degrees(38.889931),
            longitude: Degrees(-77.009003),
            elevation: Meters(0.0),
        };
        let result = pointa.haversine_distance(pointb).expect("distance");
        assert_in_eps!(5897658.288856054, result.0, 1e-6);
    }

    #[test]
    fn haversine_elevation() {
        let pointa = LLE{
            latitude: Degrees(51.510357),
            longitude: Degrees(-0.116773),
            elevation: Meters(0.0),
        };
        let pointb = LLE{
            latitude: Degrees(38.889931),
            longitude: Degrees(-77.009003),
            elevation: Meters(10.0),
        };

        assert!(pointa.haversine_distance(pointb).is_err());
        assert!(pointb.haversine_distance(pointa).is_err());
    }
}

// vim: foldmethod=marker
