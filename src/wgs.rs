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

use crate::{
    CoordinateSystem, Enu, Lle, Meters, Radians, Xyz,
    math::{atan2, cos, sin, sqrt},
};

/// Wgs ("World Geodetic System") is a set of standards published and maintained
/// by the United States National Geospatial-Intelligence Agency.
///
/// If you see Wgs, you likely actually want Wgs84.
///
/// I'm unclear how right this actually is. I have a hunch that it's
/// not as simple as this, but we'll see.
pub(crate) trait Wgs
where
    Self: Copy,
{
    const A: f64;
    const B: f64;
    const F: f64 = (Self::A - Self::B) / Self::A;
    const E_SQ: f64 = Self::F * (2.0 - Self::F);
}

impl<T, AngularMeasure> CoordinateSystem<AngularMeasure> for T
where
    AngularMeasure: From<Radians> + Copy,
    Radians: From<AngularMeasure> + Copy,
    T: Wgs,
{
    fn lle_to_xyz(g: &Lle<Self, AngularMeasure>) -> Xyz {
        let lambda: Radians = g.latitude.into();
        let phi: Radians = g.longitude.into();
        let lambda = lambda.as_float();
        let phi = phi.as_float();

        let sin_lambda = sin(lambda);
        let cos_lambda = cos(lambda);
        let sin_phi = sin(phi);
        let cos_phi = cos(phi);
        let n = Self::A / sqrt(1.0 - Self::E_SQ * sin_lambda * sin_lambda);

        Xyz {
            x: Meters::new((g.elevation.as_float() + n) * cos_lambda * cos_phi),
            y: Meters::new((g.elevation.as_float() + n) * cos_lambda * sin_phi),
            z: Meters::new((g.elevation.as_float() + (1.0 - Self::E_SQ) * n) * sin_lambda),
        }
    }

    fn xyz_to_lle(x: &Xyz) -> Lle<Self, AngularMeasure> {
        let eps = Self::E_SQ / (1.0 - Self::E_SQ);
        let p = sqrt(x.x.as_float() * x.x.as_float() + x.y.as_float() * x.y.as_float());
        let q = atan2(x.z.as_float() * Self::A, p * Self::B);

        let sin_q = sin(q);
        let cos_q = cos(q);

        let sin_q_3 = sin_q * sin_q * sin_q;
        let cos_q_3 = cos_q * cos_q * cos_q;

        let phi = atan2(
            x.z.as_float() + eps * Self::B * sin_q_3,
            p - Self::E_SQ * Self::A * cos_q_3,
        );
        let lambda = atan2(x.y.as_float(), x.x.as_float());
        let v = Self::A / sqrt(1.0 - Self::E_SQ * sin(phi) * sin(phi));
        let h = Meters::new((p / cos(phi)) - v);

        Lle::<Self, AngularMeasure>::new(Radians::new(phi).into(), Radians::new(lambda).into(), h)
    }

    fn xyz_to_enu(g: &Lle<Self, AngularMeasure>, x: &Xyz) -> Enu {
        let lambda: Radians = g.latitude.into();
        let phi: Radians = g.longitude.into();
        let lambda = lambda.as_float();
        let phi = phi.as_float();

        let sin_lambda = sin(lambda);
        let cos_lambda = cos(lambda);
        let sin_phi = sin(phi);
        let cos_phi = cos(phi);

        let xref = Self::lle_to_xyz(g);
        let xd = x.x.as_float() - xref.x.as_float();
        let yd = x.y.as_float() - xref.y.as_float();
        let zd = x.z.as_float() - xref.z.as_float();

        Enu {
            east: Meters::new(-sin_phi * xd + cos_phi * yd),
            north: Meters::new(
                -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd,
            ),
            up: Meters::new(
                cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd,
            ),
        }
    }

    fn enu_to_xyz(g: &Lle<Self, AngularMeasure>, lt: &Enu) -> Xyz {
        let lambda: Radians = g.latitude.into();
        let phi: Radians = g.longitude.into();
        let lambda = lambda.as_float();
        let phi = phi.as_float();

        let sin_lambda = sin(lambda);
        let cos_lambda = cos(lambda);
        let sin_phi = sin(phi);
        let cos_phi = cos(phi);
        let n = Self::A / sqrt(1.0 - Self::E_SQ * sin_lambda * sin_lambda);

        let x0 = (g.elevation.as_float() + n) * cos_lambda * cos_phi;
        let y0 = (g.elevation.as_float() + n) * cos_lambda * sin_phi;
        let z0 = (g.elevation.as_float() + (1.0 - Self::E_SQ) * n) * sin_lambda;

        let east = lt.east.as_float();
        let north = lt.north.as_float();
        let up = lt.up.as_float();

        let xd = -sin_phi * east - cos_phi * sin_lambda * north + cos_lambda * cos_phi * up;
        let yd = cos_phi * east - sin_lambda * sin_phi * north + cos_lambda * sin_phi * up;
        let zd = cos_lambda * north + sin_lambda * up;

        let x = xd + x0;
        let y = yd + y0;
        let z = zd + z0;

        Xyz {
            x: Meters::new(x),
            y: Meters::new(y),
            z: Meters::new(z),
        }
    }
}

// vim: foldmethod=marker
