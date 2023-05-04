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

use crate::{
  CoordinateSystem,
  LLE, XYZ, ENU,
  Meters,
  Radians,
};

/// WGS84 ("World Geodetic System 1984") is a standard published and maintained
/// by the United States National Geospatial-Intelligence Agency. This is the
/// coordinate system used by most systems "by default", including GPS.
///
/// If you have a "Latitude and Longitude" without further information about
/// the coordinate system, it's not a bad guess to try it out with WGS84 first.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct WGS84;

const WGS84_A: f64 = 6378137.0;
const WGS84_B: f64 = 6356752.314245;
const WGS84_F: f64 = (WGS84_A - WGS84_B) / WGS84_A;
const WGS84_E_SQ: f64 = WGS84_F * (2.0 - WGS84_F);

impl CoordinateSystem for WGS84 {
    fn lle_to_xyz(self, g: LLE) -> XYZ {
        let lambda = g.latitude.to_radians().0;
        let phi = g.longitude.to_radians().0;
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let n = WGS84_A / (1.0 - WGS84_E_SQ * sin_lambda * sin_lambda).sqrt();

        XYZ{
            x: Meters((g.elevation.0 + n) * cos_lambda * cos_phi),
            y: Meters((g.elevation.0 + n) * cos_lambda * sin_phi),
            z: Meters((g.elevation.0 + (1.0 - WGS84_E_SQ) * n) * sin_lambda),
        }
    }

    fn xyz_to_lle(self, x: XYZ) -> LLE {
        let eps = WGS84_E_SQ / (1.0 - WGS84_E_SQ);
        let p = (x.x.0 * x.x.0 + x.y.0 * x.y.0).sqrt();
        let q = (x.z.0 * WGS84_A).atan2(p * WGS84_B);

        let sin_q = q.sin();
        let cos_q = q.cos();

        let sin_q_3 = sin_q * sin_q * sin_q;
        let cos_q_3 = cos_q * cos_q * cos_q;

        let phi = (x.z.0 + eps * WGS84_B * sin_q_3).atan2(p - WGS84_E_SQ * WGS84_A * cos_q_3);
        let lambda = x.y.0.atan2(x.x.0);
        let v = WGS84_A / (1.0 - WGS84_E_SQ * phi.sin() * phi.sin()).sqrt();
        let h = Meters((p / phi.cos()) - v);

        LLE{
            latitude: Radians(phi).to_degrees(),
            longitude: Radians(lambda).to_degrees(),
            elevation: h,
        }
    }

    fn xyz_to_enu(self, g: LLE, x: XYZ) -> ENU {
        let lambda = g.latitude.to_radians().0;
        let phi = g.longitude.to_radians().0;
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        let xref = self.lle_to_xyz(g);
        let xd = x.x.0 - xref.x.0;
        let yd = x.y.0 - xref.y.0;
        let zd = x.z.0 - xref.z.0;

        ENU{
            east: Meters(-sin_phi*xd + cos_phi*yd),
            north: Meters(-cos_phi*sin_lambda*xd - sin_lambda*sin_phi*yd + cos_lambda*zd),
            up: Meters(cos_lambda*cos_phi*xd + cos_lambda*sin_phi*yd + sin_lambda*zd),
        }
    }

    fn enu_to_xyz(self, g: LLE, lt: ENU) -> XYZ {
        let lambda = g.latitude.to_radians().0;
        let phi = g.longitude.to_radians().0;
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let n = WGS84_A / (1.0 - WGS84_E_SQ * sin_lambda * sin_lambda).sqrt();

        let x0 = (g.elevation.0 + n) * cos_lambda * cos_phi;
        let y0 = (g.elevation.0 + n) * cos_lambda * sin_phi;
        let z0 = (g.elevation.0 + (1.0 - WGS84_E_SQ) * n) * sin_lambda;

        let east = lt.east.0;
        let north = lt.north.0;
        let up = lt.up.0;

        let xd = -sin_phi * east - cos_phi * sin_lambda * north + cos_lambda * cos_phi * up;
        let yd = cos_phi * east - sin_lambda * sin_phi * north + cos_lambda * sin_phi * up;
        let zd = cos_lambda * north + sin_lambda * up;

        let x = xd + x0;
        let y = yd + y0;
        let z = zd + z0;

        XYZ{
            x: Meters(x),
            y: Meters(y),
            z: Meters(z),
        }
    }

    fn lle_to_enu(self, r: LLE, geo: LLE) -> ENU {
        let xyz = self.lle_to_xyz(geo);
        self.xyz_to_enu(r, xyz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        LLE,
        Degrees,
    };

    macro_rules! assert_in_eps {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d && $y - $x < $d) { panic!(); }
        }
    }

    #[test]
    fn lle_to_xyz() {
        let g = LLE{
            latitude: Degrees(34.00000048),
            longitude: Degrees(-117.3335693),
            elevation: Meters(251.702),
        };
        let refx = WGS84.lle_to_xyz(g);
        assert_in_eps!(-2430601.8, refx.x.0, 0.1);
        assert_in_eps!(-4702442.7, refx.y.0, 0.1);
        assert_in_eps!(3546587.4, refx.z.0, 0.1);
    }

    #[test]
    fn around_the_world() {
        let refr = LLE{
            latitude: Degrees(38.897957),
            longitude: Degrees(-77.036560),
            elevation: Meters(30.0),
        };

        let position = LLE{
            latitude: Degrees(38.8709455),
            longitude: Degrees(-77.0552551),
            elevation: Meters(100.0),
        };

        let positionx = WGS84.lle_to_xyz(position);
        let positionenu = WGS84.xyz_to_enu(refr, positionx);
        let positionxx1 = WGS84.enu_to_xyz(refr, positionenu);
        let position1 = WGS84.xyz_to_lle(positionxx1);

        assert_in_eps!(position.latitude.0, position1.latitude.0, 1e-7);
        assert_in_eps!(position.longitude.0, position1.longitude.0, 1e-7);
        assert_in_eps!(position.elevation.0, position1.elevation.0, 1e-7);
    }

    #[test]
    fn lle_to_enu() {
        let refr = LLE{
            latitude: Degrees(34.00000048),
            longitude: Degrees(-117.3335693),
            elevation: Meters(251.702),
        };
        let refx = WGS84.lle_to_xyz(refr);

        let point = XYZ{
            x: Meters(refx.x.0 + 1.0),
            y: refx.y,
            z: refx.z,
        };
        let pointenu = WGS84.xyz_to_enu(refr, point);

        assert_in_eps!(0.88834836, pointenu.east.0, 0.1);
        assert_in_eps!(0.25676467, pointenu.north.0, 0.1);
        assert_in_eps!(-0.38066927, pointenu.up.0, 0.1);

        let point = XYZ{
            x: refx.x,
            y: Meters(refx.y.0 + 1.0),
            z: refx.z,
        };
        let pointenu = WGS84.xyz_to_enu(refr, point);
        assert_in_eps!(-0.45917011, pointenu.east.0, 0.1);
        assert_in_eps!(0.49675810, pointenu.north.0, 0.1);
        assert_in_eps!(-0.73647416, pointenu.up.0, 0.1);

        let point = XYZ{
            x: refx.x,
            y: refx.y,
            z: Meters(refx.z.0 + 1.0),
        };
        let pointenu = WGS84.xyz_to_enu(refr, point);
        assert_eq!(0.0, pointenu.east.0);
        assert_in_eps!(0.82903757, pointenu.north.0, 0.1);
        assert_in_eps!(0.55919291, pointenu.up.0, 0.1);
    }

    #[test]
    fn enu_aer_round_trip() {
        let enu = ENU{
            east:  Meters(10.0),
            north: Meters(20.0),
            up:    Meters(30.0),
        };
        let enu1 = enu.to_aer().to_enu();

        assert_in_eps!(10.0, enu1.east.0, 1e-6);
        assert_in_eps!(20.0, enu1.north.0, 1e-6);
        assert_in_eps!(30.0, enu1.up.0, 1e-6);
    }

    #[test]
    fn aer_enu_round_trip() {
        let refr = LLE{
            latitude: Degrees(38.897957),
            longitude: Degrees(-77.036560),
            elevation: Meters(30.0),
        };

        let position = LLE{
            latitude: Degrees(38.8709455),
            longitude: Degrees(-77.0552551),
            elevation: Meters(30.0),
        };

        let positionx = WGS84.lle_to_xyz(position);
        let positionenu = WGS84.xyz_to_enu(refr, positionx);

        let aed = positionenu.to_aer();
        let positionenu1 = aed.to_enu();

        let positionx1 = WGS84.enu_to_xyz(refr, positionenu1);
        let position1 = WGS84.xyz_to_lle(positionx1);

        assert_in_eps!(position.latitude.0, position1.latitude.0, 1e-7);
        assert_in_eps!(position.longitude.0, position1.longitude.0, 1e-7);
        assert_in_eps!(position.elevation.0, position1.elevation.0, 1e-7);
    }
}

// vim: foldmethod=marker
