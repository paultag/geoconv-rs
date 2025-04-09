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

//! A Maidenhead Locator is a way of describing a location on earth commonly
//! used by amateur radio operators. These are usually either 4 or 6 bytes
//! long, with increasing precision as the number of bytes increases.

use super::{Degrees, Lle, Meters, Wgs84};
use std::ops::Range;

/// A Maidenhead Grid Locator is the actual location somewhere on earth. They
/// look something like `BL11bh` or `FM18`. The more digits, the more precision
/// on the location.
///
/// Currently, only `Grid<6>` is implemented, but I'd like to sit down and
/// work through `Grid<4>` and maybe `Grid<8>`.
pub struct Grid<const N: usize>([u8; N]);

/// Possible errors that can be returned when parsing a Maidenhead Grid Locator.
#[derive(Copy, Clone, Debug)]
pub enum Error {
    /// One of the ASCII values is invalid according to the how the scheme
    /// encodes Latitude or Longitude Degrees.
    OutOfBoundsGrid,

    /// Returned if the underlying Grid type (for instance, `Grid<4>`)
    /// disagrees with the string to be parsed.
    InvalidGridLength,

    /// Returned if bytes outside the Maidenhead scheme are used in the
    /// Grid string.
    Malformed,
}

impl std::str::FromStr for Grid<6> {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != 6 {
            return Err(Error::InvalidGridLength);
        }
        Ok(Self(s.as_bytes().try_into().unwrap()))
    }
}

impl Grid<6> {
    /// Compute the center of the Maidenhead grid square. All points will
    /// be at 0 meters (on the ellipsoid).
    pub fn center_lat_lon(&self) -> Result<Lle<Wgs84>, Error> {
        let top_left = self.raw_top_left_lat_lon()?;
        Ok(Lle::new(
            Degrees::new(top_left.latitude.as_float() + (1.0 / 48.0) - 90.0),
            Degrees::new(top_left.longitude.as_float() + (1.0 / 24.0) - 180.0),
            Meters::new(0.0),
        ))
    }

    /// Compute the "top left" and "bottom right" of the Maidenhead grid square.
    /// All points will be at 0 meters (on the ellipsoid).
    pub fn corners_lat_lon(&self) -> Result<[Lle<Wgs84>; 2], Error> {
        let top_left = self.raw_top_left_lat_lon()?;

        let next_lat = top_left.latitude.as_float() + (1.0 / 24.0) - 90.0;
        let next_lon = top_left.longitude.as_float() + (1.0 / 12.0) - 180.0;
        let lat = top_left.latitude.as_float() - 90.0;
        let lon = top_left.longitude.as_float() - 180.0;

        Ok([
            Lle::new(Degrees::new(lat), Degrees::new(lon), Meters::new(0.0)),
            Lle::new(
                Degrees::new(next_lat),
                Degrees::new(next_lon),
                Meters::new(0.0),
            ),
        ])
    }

    fn raw_top_left_lat_lon(&self) -> Result<Lle<Wgs84>, Error> {
        let [
            // First is lon/lat; "F" and "M" (A-Z).
            lon,
            lat,
            // Next is our square - "1" and "8" (0-9)
            square_lon,
            square_lat,
            // Finally, our subsquare, "l" and "v" (a-z)
            subsquare_lon,
            subsquare_lat,
        ] = self.0;

        fn ascii_to_int(range: Range<u8>, v: u8) -> Result<u8, Error> {
            if !(range.contains(&v)) {
                return Err(Error::Malformed);
            }
            Ok(v - range.start)
        }

        let lon: f64 = ascii_to_int(65..91, lon)? as f64;
        let lat: f64 = ascii_to_int(65..91, lat)? as f64;

        let square_lon: f64 = ascii_to_int(48..58, square_lon)? as f64;
        let square_lat: f64 = ascii_to_int(48..58, square_lat)? as f64;

        let subsquare_lon: f64 = ascii_to_int(97..123, subsquare_lon)? as f64;
        let subsquare_lat: f64 = ascii_to_int(97..123, subsquare_lat)? as f64;

        // the 1/48 and 1/24 here will put it in the center of the grid.
        // need to implement bounding box style

        Ok(Lle::new(
            Degrees::new((lat * 10.0) + (square_lat) + (subsquare_lat / 24.0)),
            Degrees::new((lon * 20.0) + (square_lon * 2.0) + (subsquare_lon / 12.0)),
            Meters::new(0.0),
        ))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_example_from_internet() {
        let grid: Grid<6> = "IO93ob".parse().unwrap();
        let center = grid.center_lat_lon().unwrap();

        assert_eq!(-0.7916666666666856, center.longitude.as_float());
        assert_eq!(53.0625, center.latitude.as_float());
    }
}

// vim: foldmethod=marker
