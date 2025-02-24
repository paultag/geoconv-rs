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

/// [Wgs72] ("World Geodetic System 1972") is a standard published by the United
/// States National Geospatial-Intelligence Agency.
///
/// This CoordinateSystem isn't used very much any more (if at all), but is
/// implemented to assist with the conversion of old geospatial data
/// to the Wgs84 system, what you actually want.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Wgs72 {}

impl Wgs for Wgs72 {
    const A: f64 = 6378135.0;
    const B: f64 = 6356750.52;
}

#[cfg(test)]
mod tests {
    use super::Wgs72;
    use crate::{Degrees, LLE, Meters, Wgs84};

    type Wgs72LLE = LLE<Wgs72, Degrees>;
    type Wgs84LLE = LLE<Wgs84, Degrees>;

    macro_rules! assert_in_eps {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d && $y - $x < $d) {
                panic!();
            }
        };
    }

    #[test]
    fn round_trip_wgs72_to_wgs84() {
        let origin = Wgs72LLE::new(Degrees::new(10.0), Degrees::new(-20.0), Meters::new(30.0));
        let translated: Wgs84LLE = origin.translate();
        let round_trip: Wgs72LLE = translated.translate();

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
