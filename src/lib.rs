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

#![forbid(unsafe_code)]
#![deny(missing_docs)]
#![deny(missing_copy_implementations)]
#![deny(trivial_casts)]
#![deny(trivial_numeric_casts)]
#![deny(unused_import_braces)]
#![deny(unused_qualifications)]
#![deny(rustdoc::broken_intra_doc_links)]
#![deny(rustdoc::private_intra_doc_links)]
#![no_std]

//! `geoconv` implements support for converting between some basic coordinate
//! systems. This package contains support for [Wgs84]
//! Latitude/Longitude/Elevation ([Lle]) geodetic coordinates, Earth-centered,
//! Earth-fixed ([Xyz]), and local tangent plane (both East/North/Up ([Enu])
//! and Azimuth/Elevation/Range ([Aer])) systems.
//!
//! <div class="warning">
//! <b>This is also absolutely not ready to be used for navigational purposes</b>.
//! Please do not use this library in any situation that may cause harm to life
//! or property.
//! </div>
//!
//! This package is particularly useful if you know your location on Earth,
//! and want to geolocate an object that you can observe, or if you know your
//! location on Earth, and an object's location on Earth and want to determine
//! your relative positions.
//!
//! This also includes Haversine-distance using the circumference of the earth
//! for approximate great-circle distance between two Lat/Lon points, but only
//! on the "Earth's surface". In some cases (long distances) this is a way
//! better approach for distance. If both points are line-of-sight, you likely
//! care more about using the local tangent plane to work out the Range
//! to the target in 3-space by converting both lat/lons to Earth-North-Up or
//! an Azimuth/Elevation/Range.
//!
//! # Supported Coordinate Systems
//!
//! I only implemented what I was interested in ([Wgs84]), but I threw
//! in another system ([Wgs72]) since it should be pretty straight-forward to
//! support both. I have no reference data, so some additional testing for
//! that coordinate system would be most welcome.
//!
//! | Coordinate System | State       | Note                     |
//! | ----------------- | ----------- | ------------------------ |
//! | [Wgs84]           | Tested      | **You likely want this** |
//! | [Wgs72]           | Implemented |                          |
//!
//! # Types
//!
//! Within `geoconv`, the main types that you'll be interacting with for location
//! data are listed below, along with a short description.
//!
//! | Name  | Cartesian or Angular | Where is `0`, `0`, `0`?               | Description                                                  |
//! | ----- | -------------------- | ------------------------------------- | ------------------------------------------------------------ |
//! | [Lle] | Angular              | Null Island                           | Latitude, Longitude, and Elevation                           |
//! | [Xyz] | Cartesian            | Earth's Core                          | X, Y, Z (sometimes called ECEF, Earth-centered, Earth-fixed) |
//! | [Aer] | Angular              | Some point on the local tangent plane | Azimuth, Elevation and Range                                 |
//! | [Enu] | Cartesian            | Some point on the local tangent plane | East-North-Up                                                |
//!
//! # Determining the Azimuth, Elevation and Rage between two points
//!
//! Using `geoconv`, we can take some Latitude, Longitude and Elevation ([Lle])
//! and figure out where another [Lle] point is in reference to our local
//! tangent plane -- for instance, what the bearing (azimuth), elevation (up
//! and down) and range to the "other" point is from where "we" are.
//!
//! # `no_std`
//!
//! This module is `no_std`, using `libm` for the underlying math.
//!
//! ```rust
//! use geoconv::{Lle, Wgs84, Degrees, Meters, Aer};
//!
//! // Alias for a lat/lon/elevation in the Wgs84 coordinate system,
//! // used by (among many others), GPS -- this is usually what you
//! // want when you're processing lat/lon information.
//! type LleW84 = Lle<Wgs84>;
//!
//! // "My" point on earth.
//! let me = LleW84::new(
//!     Degrees::new(42.352211),
//!     Degrees::new(-71.051315),
//!     Meters::new(0.0),
//! );
//!
//! // "Your" point on earth.
//! let you = LleW84::new(
//!     Degrees::new(42.320239),
//!     Degrees::new(-70.929482),
//!     Meters::new(100.0),
//! );
//!
//! // Compute in what direction I'd need to look to see you.
//! let look: Aer<Degrees> = me.aer_to(&you);
//! ```
//!
//! # Determine the coordinates of something you can range
//!
//! Using `geoconv`, we can take some observation taken from a point,
//! and figure out where that point is. Let's work through
//! taking a reading in Azimuth, Elevation and Range ([Aer]) and
//! turning that back into Latitude, Longitude and Elevation ([Lle])
//! given our location.
//!
//! ```rust
//! use geoconv::{Lle, Wgs84, Degrees, Meters, Aer};
//!
//! type LleW84 = Lle<Wgs84>;
//!
//! // "My" point on earth.
//! let me = LleW84::new(
//!     Degrees::new(42.352211),
//!     Degrees::new(-71.051315),
//!     Meters::new(0.0),
//! );
//!
//! // I see something straight ahead of me, 45 degrees in elevation (up),
//! // and 30 meters away.
//! let observation = Aer {
//!     azimuth: Degrees::new(0.0),
//!     elevation: Degrees::new(45.0),
//!     range: Meters::new(30.0),
//! };
//!
//! // Assuming I have a perfect reading on where that object is, where
//! // is that object as a latitude/longitude?
//! let observed_lle: LleW84 = observation.to_lle(&me);
//! ```
//!
//! # Haversine (Great Circle) distance between two points
//!
//! In addition to operations on the local tangent plane, sometimes
//! the two points being measured are far enough away where the range
//! between two points is *not* the distance you'd travel between them.
//!
//! For instance, when traving, from San Francisco to Paris, we need to
//! travel along the curve of the earth, rather than directly through
//! the crust of Earth.
//!
//! To compute distance beyond line of sight between to latitude and longitude
//! points, we can use the Haversine formula to approximate the distance.
//!
//! Haversine distance **does not** easily take into account elevation,
//! so the [haversine_distance] function does not accept a [Lle], rather,
//! a tuple of AngularMeasures, such as [Degrees] or [Radians].
//!
//! ```rust
//! use geoconv::{Degrees, haversine_distance, Meters};
//!
//! // Latitude, then Longitude (in that order) for the tuples
//! let lat_lon_a = (Degrees::new(22.55), Degrees::new(43.12));
//! let lat_lob_b = (Degrees::new(13.45), Degrees::new(100.28));
//!
//! let result: Meters = haversine_distance(lat_lon_a, lat_lob_b);
//! ```

pub mod maidenhead;
mod wgs;
mod wgs72;
mod wgs84;

pub use wgs72::Wgs72;
pub use wgs84::Wgs84;

pub(crate) mod std {
    // pub use alloc::*;
    pub use core::*;
}

use libm::{atan2, cos, pow, sin, sqrt};
use std::marker::PhantomData;

/// [Meters] represent the SI unit of measure, Meter
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Meters(f64);

impl Meters {
    /// Create a new Meters struct, initialized with the provided
    /// distance, in meters.
    pub const fn new(meters: f64) -> Meters {
        Meters(meters)
    }

    /// Return the number of meters as a floating point number.
    /// This number may be less than 0 (for distances less than a meter), or
    /// very large (for kilometers, etc).
    pub const fn as_float(self) -> f64 {
        self.0
    }
}

/// [Degrees] is an angular measure that ranges from 0 to 360 (sometimes negative)
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Degrees(f64);

impl Degrees {
    /// Create a new [Degrees] struct, initialized with the provided
    /// number of degrees, either between 0 and 360 or -180 to +180.
    pub const fn new(degrees: f64) -> Degrees {
        Degrees(degrees)
    }

    /// Return the [Degrees] as a floating point number.
    pub const fn as_float(self) -> f64 {
        self.0
    }
}

impl From<Degrees> for Radians {
    fn from(deg: Degrees) -> Self {
        Radians::new(std::f64::consts::PI / 180.0 * deg.as_float())
    }
}

/// [Radians] is an angular measure that ranges from 0 to 2π (𝜏).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Radians(f64);

impl Radians {
    /// Create a new [Radians] struct, initialized with the provided
    /// number of radians, ranging between 0 and to 2π (𝜏).
    pub const fn new(radians: f64) -> Radians {
        // TODO(paultag): put this into the range 0 -> 𝜏
        Radians(radians)
    }

    /// Return the [Radians] as a floating point number from 0
    /// to 2π (𝜏).
    pub const fn as_float(self) -> f64 {
        self.0
    }
}

impl From<Radians> for Degrees {
    fn from(rad: Radians) -> Self {
        Degrees::new(180.0 / std::f64::consts::PI * rad.as_float())
    }
}

/// [Lle] or Latitude, Longitude, Elevation, is a location somewhere around
/// Earth, in refernce to some [CoordinateSystem]'s ellipsoid.
///
/// Relatedly, this means that the `elevation` is the number of [Meters]
/// above the [CoordinateSystem]'s ellipsoid, *not* the altitude above or
/// below the true surface of the Earth.
///
/// The provided `CoordinateSystem` should be the [CoordinateSystem]
/// that the Latitude, Longitude and Elevation are in refernce to, usually
/// [Wgs84]. If you don't know otherwise, it's a safe to assume that you
/// want to use `Lle<Wgs84>`. If you know otherwise, you may not have
/// needed this warning.
///
/// Additionally, the `latitude` and `longitude` fields can be specified in
/// either [Degrees] or [Radians]. Unless you know otherwise, you likely
/// want to work with [Degrees], and is the default if nothing is specified.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Lle<CoordinateSystem, AngularMeasure = Degrees>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    _unused: PhantomData<CoordinateSystem>,

    /// Latitude, in degrees, for some geodetic coordinate system.
    pub latitude: AngularMeasure,

    /// Longitude, in degrees, for some geodetic coordinate system.
    pub longitude: AngularMeasure,

    /// Elevation above the ellipsoid for some geodetic coordinate system.
    /// This is *not* the same as altitude above the surface.
    pub elevation: Meters,
}

/// Alias for [LLE]. This struct was renamed to better match Rust style,
/// and this alias will be removed in a future release.
#[deprecated]
pub type LLE<CoordinateSystem> = Lle<CoordinateSystem>;

impl<CoordinateSystem, AngularMeasure> Lle<CoordinateSystem, AngularMeasure>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Create a new Lle (latitude, longitude, elevation) from parts.
    pub const fn new(
        latitude: AngularMeasure,
        longitude: AngularMeasure,
        elevation: Meters,
    ) -> Lle<CoordinateSystem, AngularMeasure> {
        Lle {
            latitude,
            longitude,
            elevation,
            _unused: PhantomData,
        }
    }

    /// Compute the East-North-Up of some provided point (`other`) given
    /// some refernce location `self`.
    pub fn enu_to(&self, other: &Lle<CoordinateSystem, AngularMeasure>) -> Enu {
        CoordinateSystem::lle_to_enu(self, other)
    }

    /// Compute the Az-El-Range of some provided point (`other`) given
    /// some reference location `self`.
    pub fn aer_to<AerAngularMeasure>(
        &self,
        other: &Lle<CoordinateSystem, AngularMeasure>,
    ) -> Aer<AerAngularMeasure>
    where
        Aer<AerAngularMeasure>: From<Enu>,
    {
        self.enu_to(other).into()
    }
}

/// [Xyz] is the earth-centric Xyz point system.
///
/// `{ x: 0.0, y: 0.0, z: 0.0 }` is the center of the CoordinateSystem (which
/// is to say, inside Earth's Core). X/Y/Z coordinates are cartesian, and
/// as they grow larger in absolute value, will get further from the center of
/// Earth.
///
/// Xyz locations can be turned into points in reference to Earth's ellipsoid
/// by using some [CoordinateSystem], such as [Wgs84], either directly
/// or via [Lle].
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Xyz {
    /// Distance from the center of the ECEF coordinate system. The X axis
    /// intersects the equator along the semi-major axis.
    pub x: Meters,

    /// Distance from the center of the ECEF coordinate system. The Y axis
    /// intersects the equator along the semi-minor axis.
    pub y: Meters,

    /// Distance from the center of the ECEF coordinate system. The Z axis
    /// intersects the north and south poles.
    pub z: Meters,
}

/// Alias for [XYZ]. This struct was renamed to better match Rust style,
/// and this alias will be removed in a future release.
#[deprecated]
pub type XYZ = Xyz;

// Implemnt conversions to/from Lle and Xyz coordinates.

impl<CoordinateSystem, AngularMeasure> From<Lle<CoordinateSystem, AngularMeasure>> for Xyz
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Convert some [Lle] into an [Xyz] coordinate using [Lle]'s
    /// [CoordinateSystem].
    fn from(lle: Lle<CoordinateSystem, AngularMeasure>) -> Self {
        CoordinateSystem::lle_to_xyz(&lle)
    }
}

impl<CoordinateSystem, AngularMeasure> From<Xyz> for Lle<CoordinateSystem, AngularMeasure>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Convert some [Xyz] into an [Lle] coordinate using [Lle]'s
    /// [CoordinateSystem].
    fn from(xyz: Xyz) -> Self {
        CoordinateSystem::xyz_to_lle(&xyz)
    }
}

/// [Aer] represents an Azimuth, Elevation, and Range measurement.
///
/// Azimuth/Elevation (or Az/El) is a common way of locating objects measured
/// at a specific location on Earth in the local tangent plane
/// (for instance, from a RADAR). That reference location is not included in
/// the [Aer] struct, so re-contextualizing this requires knowing the point
/// at which the observation was taken from or was in reference to.
///
/// An [Aer] can be passed an "`AngularMeasure`" generic argument, which means,
/// for all practical purposes, either [Degrees] or [Radians]. The default is
/// [Degrees], but can be overriden by passing [Radians] into the [Aer].
/// [Degrees] tend to be a lot easier to use as an API and for humans, but
/// when passing data between third party code and this library, you may find
/// it's less work to skip conversions through [Degrees] when reading or
/// writing [Aer] data. If you don't have some overriding conviction, using
/// [Degrees] is a good idea.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Aer<AngularMeasure = Degrees> {
    /// Angle in azimuth (left and right) in relation to some fixed point and
    /// orientation.
    pub azimuth: AngularMeasure,

    /// Angle in elevation (up and down) in relation to some fixed point and
    /// orientation.
    pub elevation: AngularMeasure,

    /// Range to some point in space along the defined azimuth and elevation,
    /// in Meters.
    pub range: Meters,
}

/// Alias for [AER]. This struct was renamed to better match Rust style,
/// and this alias will be removed in a future release.
#[deprecated]
pub type AER = Aer;

impl<AngularMeasure> Aer<AngularMeasure>
where
    Enu: From<Aer<AngularMeasure>>,
    Self: Copy,
{
    /// Take some observation (`self`), and contextualize it into some
    /// absolute [Lle] given some reference point `obs` that this [Aer] is
    /// in reference to.
    pub fn to_lle<CoordinateSystem, LleAngularMeasure>(
        &self,
        obs: &Lle<CoordinateSystem, LleAngularMeasure>,
    ) -> Lle<CoordinateSystem, LleAngularMeasure>
    where
        Radians: From<LleAngularMeasure> + Copy,
        LleAngularMeasure: From<Radians> + Copy,
        CoordinateSystem: crate::CoordinateSystem<LleAngularMeasure>,
    {
        let enu: Enu = (*self).into();
        enu.to_lle(obs)
    }
}

// Here we need to bite the bullet and just implement manual conversions.
// We can't do a blanket implementation using `From` and generics since
// the stdlib implements From<X> for X, which winds up getting a conflicting
// implementation.
//
// As a result, every new AngularMeasure is going to need to grow stupid
// conversions. If I need to add a new angle type, this should all be
// turned into macros. It's the literal same code.

impl From<Aer<Degrees>> for Aer<Radians> {
    /// Convert an [Aer] in [Degrees] into [Radians]
    fn from(aer: Aer<Degrees>) -> Self {
        Self {
            azimuth: aer.azimuth.into(),
            elevation: aer.elevation.into(),
            range: aer.range,
        }
    }
}

impl From<Aer<Radians>> for Aer<Degrees> {
    /// Convert an [Aer] in [Radians] into [Degrees].
    fn from(aer: Aer<Radians>) -> Self {
        Self {
            azimuth: aer.azimuth.into(),
            elevation: aer.elevation.into(),
            range: aer.range,
        }
    }
}

// Implement conversions from Aer -> Enu; for Radians and Degrees. Again,
// this can't be done easily using generics, so I'm just going to make it
// explicit.

impl From<Aer<Radians>> for Enu {
    /// Convert to [Enu] cartesian local tangent plane coordinates from an
    /// [Aer] local tangent plane angular measurement in [Radians].
    fn from(aer: Aer<Radians>) -> Self {
        let az_rad: Radians = aer.azimuth;
        let el_rad: Radians = aer.elevation;
        let r = Meters::new(aer.range.as_float() * cos(el_rad.as_float()));
        Enu {
            east: Meters::new(r.as_float() * sin(az_rad.as_float())),
            north: Meters::new(r.as_float() * cos(az_rad.as_float())),
            up: Meters::new(aer.range.as_float() * sin(el_rad.as_float())),
        }
    }
}

impl From<Aer<Degrees>> for Enu {
    /// Convert to [Enu] cartesian local tangent plane coordinates from an
    /// [Aer] local tangent plane angular measurement in [Degrees].
    fn from(aer: Aer<Degrees>) -> Self {
        let aer: Aer<Radians> = aer.into();
        aer.into()
    }
}

/// [Enu] is East, North, Up in Meters.
///
/// East-North-Up are cartesian coordinates rooted at some reference
/// point's local tangent plane. That reference location is not included in
/// the [Enu] struct, so re-contextualizing an [Enu] to a [Lle] or [Xyz]
/// requires knowing the point at which the observation was taken from or
/// was in reference to.
///
/// These measures are cartesian in the local tangent plane, which is to say,
/// increasing "North" will continue straight ahead, but not following the
/// curve of the Earth -- so you'll slowly get further and further away from
/// Earth's surface.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Enu {
    /// East is the distance in the easterly direction on some local tangent
    /// plane reference. This may be negative to indicate a westerly
    /// distance.
    pub east: Meters,

    /// North is the distance in the northernly direction on some local tangent
    /// plane reference. This may be negative to indicate a southern
    /// distance.
    pub north: Meters,

    /// Up is the distance above the local tangent plane that intersects with
    /// the coordinate system ellipsoid (this is *not* the same as altitude
    /// above the surface, this is elevation above the ellipsoid).
    pub up: Meters,
}

/// Alias for [Enu]. This struct was renamed to better match Rust style,
/// and this alias will be removed in a future release.
#[deprecated]
pub type ENU = Enu;

impl Enu {
    /// Take some observation (`self`), and contextualize it into some
    /// absolute [Lle] given some reference point `obs` that this [Aer]
    /// is in reference to.
    pub fn to_lle<CoordinateSystem, LleAngularMeasure>(
        &self,
        obs: &Lle<CoordinateSystem, LleAngularMeasure>,
    ) -> Lle<CoordinateSystem, LleAngularMeasure>
    where
        Radians: From<LleAngularMeasure> + Copy,
        LleAngularMeasure: From<Radians> + Copy,
        CoordinateSystem: crate::CoordinateSystem<LleAngularMeasure>,
    {
        CoordinateSystem::enu_to_lle(obs, self)
    }
}

// Implement conversions from Aer -> Enu; for Radians and Degrees. Again,
// this can't be done easily using generics, so I'm just going to make it
// explicit.

impl From<Enu> for Aer<Radians> {
    /// Convert to [Aer] cartesian local tangent plane angular measurement in
    /// [Radians] from [Enu] local tangent plane coordinates.
    fn from(enu: Enu) -> Self {
        let r = sqrt(
            enu.east.as_float() * enu.east.as_float() + enu.north.as_float() * enu.north.as_float(),
        );

        let tau = std::f64::consts::PI * 2.0;
        Self {
            azimuth: Radians::new(atan2(enu.east.as_float(), enu.north.as_float()) % tau),
            elevation: Radians::new(atan2(enu.up.as_float(), r)),
            range: Meters::new(sqrt(r * r + enu.up.as_float() * enu.up.as_float())),
        }
    }
}

impl From<Enu> for Aer<Degrees> {
    /// Convert to [Aer] cartesian local tangent plane angular measurement in
    /// [Degrees] from [Enu] local tangent plane coordinates.
    fn from(enu: Enu) -> Self {
        let aer: Aer<Radians> = enu.into();
        aer.into()
    }
}

/// [CoordinateSystem] is a trait to enable converstion between locations
/// (usually Latitude and Longitude, in the form of [Lle] objects) to absolute
/// points in space and vice versa.
///
/// Different systems have different measurements of Earth's surface, and
/// a Latitude / Longitude must be understood within its [CoordinateSystem],
/// or significant errors can be introduced.
pub trait CoordinateSystem<AngularMeasure>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    Self: Sized,
{
    /// Convert [Xyz] to [Lle] in our [CoordinateSystem].
    fn xyz_to_lle(c: &Xyz) -> Lle<Self, AngularMeasure>;

    /// Convert [Xyz] to [Enu], referenced to `refr` ([Lle]) in our
    /// [CoordinateSystem].
    fn xyz_to_enu(refr: &Lle<Self, AngularMeasure>, c: &Xyz) -> Enu;

    /// Convert [Enu] to [Xyz], referenced to `refr` ([Lle]) in our
    /// [CoordinateSystem].
    fn enu_to_xyz(refr: &Lle<Self, AngularMeasure>, c: &Enu) -> Xyz;

    /// Convert [Lle] to [Xyz] in our [CoordinateSystem].
    fn lle_to_xyz(geo: &Lle<Self, AngularMeasure>) -> Xyz;

    /// Convert [Enu] to [Lle], referenced to `refr` ([Lle]) in our
    /// [CoordinateSystem].
    fn enu_to_lle(refr: &Lle<Self, AngularMeasure>, c: &Enu) -> Lle<Self, AngularMeasure> {
        let xyz = Self::enu_to_xyz(refr, c);
        Self::xyz_to_lle(&xyz)
    }

    /// Convert [Lle] to [Enu], referenced to `refr` ([Lle]) in our
    /// [CoordinateSystem].
    fn lle_to_enu(r: &Lle<Self, AngularMeasure>, geo: &Lle<Self, AngularMeasure>) -> Enu {
        let xyz = Self::lle_to_xyz(geo);
        Self::xyz_to_enu(r, &xyz)
    }
}

impl<FromCoordinateSystem, FromAngularMeasure> Lle<FromCoordinateSystem, FromAngularMeasure>
where
    Radians: From<FromAngularMeasure> + Copy,
    FromAngularMeasure: From<Radians> + Copy,
    FromCoordinateSystem: CoordinateSystem<FromAngularMeasure>,
    Self: Copy,
{
    /// Translate from one [CoordinateSystem] to another [CoordinateSystem].
    ///
    /// This isn't a From trait since Rust doesn't have the ability to implement
    /// conversions between the same type. In the future when Rust can handle
    /// specialization or negating a blanket generic implementation
    /// (such as `where FromCoordinateSystem != ToCoordinateSystem`), this
    /// method will be deprecated for a plain `.into()`.
    ///
    /// <div class="warning">
    /// It's not clear to me that, as things are implemented right now,
    /// that <code>Xyz</code>'s <code>(0.0, 0.0, 0.0)</code> is the same
    /// exact point in space for all the <code>CoordinateSystem</code>s,
    /// although this library <u>will assume they are</u>. As such, I'm not
    /// confident these conversions are implemented properly, and
    /// <b>must not be relied upon for anything approaching important</b>.
    /// </div>
    pub fn translate<ToCoordinateSystem, ToAngularMeasure>(
        &self,
    ) -> Lle<ToCoordinateSystem, ToAngularMeasure>
    where
        Radians: From<ToAngularMeasure> + Copy,
        ToAngularMeasure: From<Radians> + Copy,
        ToCoordinateSystem: CoordinateSystem<ToAngularMeasure>,
    {
        let xyz = FromCoordinateSystem::lle_to_xyz(self);
        ToCoordinateSystem::xyz_to_lle(&xyz)
    }
}

/// For both [haversine_distance] and [bearing] we need to know these values
/// in this format, and writing this every time is a bit tedious.
fn angular_measure_to_radian_delta<AngularMeasure>(
    self_lat_lon: (AngularMeasure, AngularMeasure),
    other_lat_lon: (AngularMeasure, AngularMeasure),
) -> ((f64, f64), (f64, f64), (f64, f64))
where
    Radians: From<AngularMeasure> + Copy,
{
    let (self_lat, self_lon) = self_lat_lon;
    let (self_lat, self_lon): (Radians, Radians) = (self_lat.into(), self_lon.into());
    let (self_lat, self_lon) = (self_lat.as_float(), self_lon.as_float());

    let (other_lat, other_lon) = other_lat_lon;
    let (other_lat, other_lon): (Radians, Radians) = (other_lat.into(), other_lon.into());
    let (other_lat, other_lon) = (other_lat.as_float(), other_lon.as_float());

    let delta_lat = self_lat - other_lat;
    let delta_lon = self_lon - other_lon;

    (
        (delta_lat, delta_lon),
        (self_lat, self_lon),
        (other_lat, other_lon),
    )
}

/// Haversine (Great Circle) distance between two points returns the
/// distance between to Latitude/Longitude points along the great circle
/// path along Earth's "surface", in either [Degrees] or [Radians].
///
/// "Surface" has quotes around it beacuse great circle paths as
/// implemented by the haversine formula (and therefore [haversine_distance])
/// will treat earth as a perfect sphere, sure to annoy both flat-earthers
/// as well as globeheads. This means that, in reality, the point on the
/// sphere's latitude/longitude at elevation `0` may be above or below the
/// true surface of earth, some [CoordinateSystem] ellipsoid (such as the
/// very likely case these are [Wgs84] coordinates) -- and absolutely won't
/// account for changes in the true elevation on Earth's surface.
///
/// But hey, what's [804672](https://www.youtube.com/watch?v=tbNlMtqrYS0)
/// [Meters] between friends?
pub fn haversine_distance<AngularMeasure>(
    self_lat_lon: (AngularMeasure, AngularMeasure),
    other_lat_lon: (AngularMeasure, AngularMeasure),
) -> Meters
where
    Radians: From<AngularMeasure> + Copy,
{
    const EARTH_RADIUS: Meters = Meters(6371000.0);

    let (
        (delta_lat, delta_lon),
        // f64 of the measure in radians.
        (self_lat, _self_lon),
        (other_lat, _other_lon),
    ) = angular_measure_to_radian_delta(self_lat_lon, other_lat_lon);

    let a = pow(sin(delta_lat / 2.0), 2.0)
        + cos(self_lat) * cos(other_lat) * pow(sin(delta_lon / 2.0), 2.0);
    let c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    Meters::new(EARTH_RADIUS.as_float() * c)
}

/// Compute the bearing (compass heading) from one point to another, provided
/// in Latitude/Longitude points along Earth's "surface" in either [Degrees]
/// or [Radians].
pub fn bearing<AngularMeasure>(
    self_lat_lon: (AngularMeasure, AngularMeasure),
    other_lat_lon: (AngularMeasure, AngularMeasure),
) -> AngularMeasure
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
{
    let (
        (_delta_lat, delta_lon),
        // f64 of the measure in radians.
        (self_lat, _self_lon),
        (other_lat, _other_lon),
    ) = angular_measure_to_radian_delta(self_lat_lon, other_lat_lon);

    let x = cos(self_lat) * sin(other_lat) - sin(self_lat) * cos(other_lat) * cos(delta_lon);
    let y = sin(delta_lon) * cos(other_lat);

    Radians::new(atan2(y, x) + std::f64::consts::PI).into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degrees_as_radians() {
        let result: Radians = Degrees::new(180.0).into();
        assert_eq!(result.as_float(), std::f64::consts::PI);
    }

    #[test]
    fn radians_as_degrees() {
        let result: Degrees = Radians::new(std::f64::consts::PI).into();
        assert_eq!(result.as_float(), 180.0);
    }

    macro_rules! assert_in_eps {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d && $y - $x < $d) {
                panic!("{} is not ~= to {}", $x, $y);
            }
        };
    }

    #[test]
    fn haversine_known() {
        let result = haversine_distance(
            (Degrees::new(22.55), Degrees::new(43.12)),
            (Degrees::new(13.45), Degrees::new(100.28)),
        );
        assert_in_eps!(6094544.408786774, result.as_float(), 1e-6);

        let result = haversine_distance(
            (Degrees::new(51.510357), Degrees::new(-0.116773)),
            (Degrees::new(38.889931), Degrees::new(-77.009003)),
        );
        assert_in_eps!(5897658.288856054, result.as_float(), 1e-6);
    }

    #[test]
    fn bearing_known() {
        // somewhere in south america's bearing to null island (strictly
        // due east).
        let result = bearing(
            (Degrees::new(0.0), Degrees::new(-88.0)),
            (Degrees::new(0.0), Degrees::new(0.0)),
        );
        assert_in_eps!(90.0, result.as_float(), 1e-6);

        let result = bearing(
            (Degrees::new(0.0), Degrees::new(0.0)),
            (Degrees::new(0.0), Degrees::new(-88.0)),
        );
        assert_in_eps!(270.0, result.as_float(), 1e-6);
    }
}

// vim: foldmethod=marker
