// {{{ Copyright (c) Paul R. Tagliamonte <paultag@gmail.com>, 2023-2024
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

use crate::wgs::Wgs;

/// [Wgs84] ("World Geodetic System 1984") is a standard published and maintained
/// by the United States National Geospatial-Intelligence Agency. This is the
/// CoordinateSystem used by most systems "by default", including GPS.
///
/// If you have a "Latitude and Longitude" without further information about
/// the coordinate system, it's not a bad guess to try it out with [Wgs84] first.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Wgs84 {}

impl Wgs for Wgs84 {
    const A: f64 = 6378137.0;
    const B: f64 = 6356752.314245;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AER, CoordinateSystem, Degrees, ENU, LLE, Meters, Wgs72, XYZ};

    type Wgs84LLE = LLE<Wgs84, Degrees>;
    type Wgs72LLE = LLE<Wgs72, Degrees>;

    macro_rules! assert_in_eps {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d && $y - $x < $d) {
                panic!();
            }
        };
    }

    #[test]
    fn lle_to_xyz() {
        let g = Wgs84LLE::new(
            Degrees::new(34.00000048),
            Degrees::new(-117.3335693),
            Meters::new(251.702),
        );
        let refx = Wgs84::lle_to_xyz(&g);
        assert_in_eps!(-2430601.8, refx.x.as_float(), 0.1);
        assert_in_eps!(-4702442.7, refx.y.as_float(), 0.1);
        assert_in_eps!(3546587.4, refx.z.as_float(), 0.1);
    }

    #[test]
    fn around_the_world() {
        let refr = Wgs84LLE::new(
            Degrees::new(38.897957),
            Degrees::new(-77.036560),
            Meters::new(30.0),
        );

        let position = Wgs84LLE::new(
            Degrees::new(38.8709455),
            Degrees::new(-77.0552551),
            Meters::new(100.0),
        );

        let positionx = Wgs84::lle_to_xyz(&position);
        let positionenu = Wgs84::xyz_to_enu(&refr, &positionx);
        let positionxx1 = Wgs84::enu_to_xyz(&refr, &positionenu);
        let position1: Wgs84LLE = Wgs84::xyz_to_lle(&positionxx1);

        assert_in_eps!(
            position.latitude.as_float(),
            position1.latitude.as_float(),
            1e-7
        );
        assert_in_eps!(
            position.longitude.as_float(),
            position1.longitude.as_float(),
            1e-7
        );
        assert_in_eps!(
            position.elevation.as_float(),
            position1.elevation.as_float(),
            1e-7
        );
    }

    #[test]
    fn lle_to_enu() {
        let refr = Wgs84LLE::new(
            Degrees::new(34.00000048),
            Degrees::new(-117.3335693),
            Meters::new(251.702),
        );
        let refx = Wgs84::lle_to_xyz(&refr);

        let point = XYZ {
            x: Meters::new(refx.x.as_float() + 1.0),
            y: refx.y,
            z: refx.z,
        };
        let pointenu = Wgs84::xyz_to_enu(&refr, &point);

        assert_in_eps!(0.88834836, pointenu.east.as_float(), 0.1);
        assert_in_eps!(0.25676467, pointenu.north.as_float(), 0.1);
        assert_in_eps!(-0.38066927, pointenu.up.as_float(), 0.1);

        let point = XYZ {
            x: refx.x,
            y: Meters::new(refx.y.as_float() + 1.0),
            z: refx.z,
        };
        let pointenu = Wgs84::xyz_to_enu(&refr, &point);
        assert_in_eps!(-0.45917011, pointenu.east.as_float(), 0.1);
        assert_in_eps!(0.49675810, pointenu.north.as_float(), 0.1);
        assert_in_eps!(-0.73647416, pointenu.up.as_float(), 0.1);

        let point = XYZ {
            x: refx.x,
            y: refx.y,
            z: Meters::new(refx.z.as_float() + 1.0),
        };
        let pointenu = Wgs84::xyz_to_enu(&refr, &point);
        assert_eq!(0.0, pointenu.east.as_float());
        assert_in_eps!(0.82903757, pointenu.north.as_float(), 0.1);
        assert_in_eps!(0.55919291, pointenu.up.as_float(), 0.1);
    }

    #[test]
    fn enu_aer_round_trip() {
        let enu = ENU {
            east: Meters::new(10.0),
            north: Meters::new(20.0),
            up: Meters::new(30.0),
        };

        let aer: AER<Degrees> = enu.into();
        let enu1: ENU = aer.into();

        assert_in_eps!(10.0, enu1.east.as_float(), 1e-6);
        assert_in_eps!(20.0, enu1.north.as_float(), 1e-6);
        assert_in_eps!(30.0, enu1.up.as_float(), 1e-6);
    }

    #[test]
    fn aer_enu_round_trip() {
        let refr = Wgs84LLE::new(
            Degrees::new(38.897957),
            Degrees::new(-77.036560),
            Meters::new(30.0),
        );

        let position = Wgs84LLE::new(
            Degrees::new(38.8709455),
            Degrees::new(-77.0552551),
            Meters::new(30.0),
        );

        let positionx = Wgs84::lle_to_xyz(&position);
        let positionenu = Wgs84::xyz_to_enu(&refr, &positionx);

        let aed: AER<Degrees> = positionenu.into();
        let positionenu1: ENU = aed.into();

        let positionx1 = Wgs84::enu_to_xyz(&refr, &positionenu1);
        let position1: Wgs84LLE = Wgs84::xyz_to_lle(&positionx1);

        assert_in_eps!(
            position.latitude.as_float(),
            position1.latitude.as_float(),
            1e-7
        );
        assert_in_eps!(
            position.longitude.as_float(),
            position1.longitude.as_float(),
            1e-7
        );
        assert_in_eps!(
            position.elevation.as_float(),
            position1.elevation.as_float(),
            1e-7
        );
    }

    #[test]
    fn round_trip_wgs84_to_wgs72() {
        let origin = Wgs84LLE::new(Degrees::new(10.0), Degrees::new(-20.0), Meters::new(30.0));
        let translated: Wgs72LLE = origin.translate();
        let round_trip: Wgs84LLE = translated.translate();

        assert_in_eps!(
            origin.latitude.as_float(),
            round_trip.latitude.as_float(),
            1e-7
        );

        assert_in_eps!(
            origin.longitude.as_float(),
            round_trip.longitude.as_float(),
            1e-7
        );

        assert_in_eps!(
            origin.elevation.as_float(),
            round_trip.elevation.as_float(),
            1e-7
        );
    }
}

// vim: foldmethod=marker
