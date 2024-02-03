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

//! `geoconv` implements support for converting between some basic coordinate
//! systems. This package contains support for [WGS84]
//! Latitude/Longitude/Elevation ([LLE]) geodetic coordinates, Earth-centered,
//! Earth-fixed ([XYZ]), and local tangent plane (both East/North/Up ([ENU])
//! and Azimuth/Elevation/Range ([AER])) systems.
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
//! I only implemented what I was interested in ([WGS84]), but I threw
//! in another system ([WGS72]) since it should be pretty straight-forward to
//! support both. I have no reference data, so some additional testing for
//! that coordinate system would be most welcome.
//!
//! | Coordinate System | State       | Note                     |
//! | ----------------- | ----------- | ------------------------ |
//! | [WGS84]           | Tested      | **You likely want this** |
//! | [WGS72]           | Implemented |                          |
//!
//! # Types
//!
//! Within `geoconv`, the main types that you'll be interacting with for location
//! data are listed below, along with a short description.
//!
//! | Name  | Cartesian or Angular | Where is `0`, `0`, `0`?               | Description                                                  |
//! | ----- | -------------------- | ------------------------------------- | ------------------------------------------------------------ |
//! | [LLE] | Angular              | Null Island                           | Latitude, Longitude, and Elevation                           |
//! | [XYZ] | Cartesian            | Earth's Core                          | X, Y, Z (sometimes called ECEF, Earth-centered, Earth-fixed) |
//! | [AER] | Angular              | Some point on the local tangent plane | Azimuth, Elevation and Range                                 |
//! | [ENU] | Cartesian            | Some point on the local tangent plane | East-North-Up                                                |
//!
//! # Determining the Azimuth, Elevation and Rage between two points
//!
//! Using `geoconv`, we can take some Latitude, Longitude and Elevation ([LLE])
//! and figure out where another [LLE] point is in reference to our local
//! tangent plane -- for instance, what the bearing (azimuth), elevation (up
//! and down) and range to the "other" point is from where "we" are.
//!
//! ```rust
//! use geoconv::{LLE, WGS84, Degrees, Meters, AER};
//!
//! // Alias for a lat/lon/elevation in the WGS84 coordinate system,
//! // used by (among many others), GPS -- this is usually what you
//! // want when you're processing lat/lon information.
//! type LLEWGS84 = LLE<WGS84>;
//!
//! // "My" point on earth.
//! let me = LLEWGS84::new(
//!     Degrees::new(42.352211),
//!     Degrees::new(-71.051315),
//!     Meters::new(0.0),
//! );
//!
//! // "Your" point on earth.
//! let you = LLEWGS84::new(
//!     Degrees::new(42.320239),
//!     Degrees::new(-70.929482),
//!     Meters::new(100.0),
//! );
//!
//! // Compute in what direction I'd need to look to see you.
//! let look: AER<Degrees> = me.aer_to(&you);
//! ```
//!
//! # Determine the coordinates of something you can range
//!
//! Using `geoconv`, we can take some observation taken from a point,
//! and figure out where that point is. Let's work through
//! taking a reading in Azimuth, Elevation and Range ([AER]) and
//! turning that back into Latitude, Longitude and Elevation ([LLE])
//! given our location.
//!
//! ```rust
//! use geoconv::{LLE, WGS84, Degrees, Meters, AER};
//!
//! type LLEWGS84 = LLE<WGS84>;
//!
//! // "My" point on earth.
//! let me = LLEWGS84::new(
//!     Degrees::new(42.352211),
//!     Degrees::new(-71.051315),
//!     Meters::new(0.0),
//! );
//!
//! // I see something straight ahead of me, 45 degrees in elevation (up),
//! // and 30 meters away.
//! let observation = AER {
//!     azimuth: Degrees::new(0.0),
//!     elevation: Degrees::new(45.0),
//!     range: Meters::new(30.0),
//! };
//!
//! // Assuming I have a perfect reading on where that object is, where
//! // is that object as a latitude/longitude?
//! let observed_lle: LLEWGS84 = observation.to_lle(&me);
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
//! so the [haversine_distance] function does not accept a [LLE], rather,
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

mod wgs;
mod wgs72;
mod wgs84;

pub use wgs72::WGS72;
pub use wgs84::WGS84;

use std::marker::PhantomData;

/// [Meters] represent the SI unit of measure, Meter
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Meters(f64);

impl Meters {
    /// Create a new Meters struct, initialized with the provided
    /// distance, in meters.
    pub fn new(meters: f64) -> Meters {
        Meters(meters)
    }

    /// Return the number of meters as a floating point number.
    /// This number may be less than 0 (for distances less than a meter), or
    /// very large (for kilometers, etc).
    pub fn as_float(self) -> f64 {
        self.0
    }
}

/// [Degrees] is an angular measure that ranges from 0 to 360 (sometimes negative)
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Degrees(f64);

impl Degrees {
    /// Create a new [Degrees] struct, initialized with the provided
    /// number of degrees, either between 0 and 360 or -180 to +180.
    pub fn new(degrees: f64) -> Degrees {
        Degrees(degrees)
    }

    /// Return the [Degrees] as a floating point number.
    pub fn as_float(self) -> f64 {
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
    pub fn new(radians: f64) -> Radians {
        // TODO(paultag): put this into the range 0 -> 𝜏
        Radians(radians)
    }

    /// Return the [Radians] as a floating point number from 0
    /// to 2π (𝜏).
    pub fn as_float(self) -> f64 {
        self.0
    }
}

impl From<Radians> for Degrees {
    fn from(rad: Radians) -> Self {
        Degrees::new(180.0 / std::f64::consts::PI * rad.as_float())
    }
}

/// [LLE] or Latitude, Longitude, Elevation, is a location somewhere around
/// Earth, in refernce to some [CoordinateSystem]'s ellipsoid.
///
/// Relatedly, this means that the `elevation` is the number of [Meters]
/// above the [CoordinateSystem]'s ellipsoid, *not* the altitude above or
/// below the true surface of the Earth.
///
/// The provided `CoordinateSystem` should be the [CoordinateSystem]
/// that the Latitude, Longitude and Elevation are in refernce to, usually
/// [WGS84]. If you don't know otherwise, it's a safe to assume that you
/// want to use `LLE<WGS84>`. If you know otherwise, you may not have
/// needed this warning.
///
/// Additionally, the `latitude` and `longitude` fields can be specified in
/// either [Degrees] or [Radians]. Unless you know otherwise, you likely
/// want to work with [Degrees], and is the default if nothing is specified.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LLE<CoordinateSystem, AngularMeasure = Degrees>
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

impl<CoordinateSystem, AngularMeasure> LLE<CoordinateSystem, AngularMeasure>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Create a new LLE (latitude, longitude, elevation) from parts.
    pub fn new(
        latitude: AngularMeasure,
        longitude: AngularMeasure,
        elevation: Meters,
    ) -> LLE<CoordinateSystem, AngularMeasure> {
        LLE {
            latitude,
            longitude,
            elevation,
            _unused: PhantomData,
        }
    }

    /// Compute the East-North-Up of some provided point (`other`) given
    /// some refernce location `self`.
    pub fn enu_to(&self, other: &LLE<CoordinateSystem, AngularMeasure>) -> ENU {
        CoordinateSystem::lle_to_enu(self, other)
    }

    /// Compute the Az-El-Range of some provided point (`other`) given
    /// some reference location `self`.
    pub fn aer_to<AERAngularMeasure>(
        &self,
        other: &LLE<CoordinateSystem, AngularMeasure>,
    ) -> AER<AERAngularMeasure>
    where
        AER<AERAngularMeasure>: From<ENU>,
    {
        self.enu_to(other).into()
    }
}

/// [XYZ] is the earth-centric XYZ point system.
///
/// `{ x: 0.0, y: 0.0, z: 0.0 }` is the center of the CoordinateSystem (which
/// is to say, inside Earth's Core). X/Y/Z coordinates are cartesian, and
/// as they grow larger in absolute value, will get further from the center of
/// Earth.
///
/// XYZ locations can be turned into points in reference to Earth's ellipsoid
/// by using some [CoordinateSystem], such as [WGS84], either directly
/// or via [LLE].
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct XYZ {
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

// Implemnt conversions to/from LLE and XYZ coordinates.

impl<CoordinateSystem, AngularMeasure> From<LLE<CoordinateSystem, AngularMeasure>> for XYZ
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Convert some [LLE] into an [XYZ] coordinate using [LLE]'s
    /// [CoordinateSystem].
    fn from(lle: LLE<CoordinateSystem, AngularMeasure>) -> Self {
        CoordinateSystem::lle_to_xyz(&lle)
    }
}

impl<CoordinateSystem, AngularMeasure> From<XYZ> for LLE<CoordinateSystem, AngularMeasure>
where
    Radians: From<AngularMeasure> + Copy,
    AngularMeasure: From<Radians> + Copy,
    CoordinateSystem: crate::CoordinateSystem<AngularMeasure>,
{
    /// Convert some [XYZ] into an [LLE] coordinate using [LLE]'s
    /// [CoordinateSystem].
    fn from(xyz: XYZ) -> Self {
        CoordinateSystem::xyz_to_lle(&xyz)
    }
}

/// [AER] represents an Azimuth, Elevation, and Range measurement.
///
/// Azimuth/Elevation (or Az/El) is a common way of locating objects measured
/// at a specific location on Earth in the local tangent plane
/// (for instance, from a RADAR). That reference location is not included in
/// the [AER] struct, so re-contextualizing this requires knowing the point
/// at which the observation was taken from or was in reference to.
///
/// An [AER] can be passed an "`AngularMeasure`" generic argument, which means,
/// for all practical purposes, either [Degrees] or [Radians]. The default is
/// [Degrees], but can be overriden by passing [Radians] into the [AER].
/// [Degrees] tend to be a lot easier to use as an API and for humans, but
/// when passing data between third party code and this library, you may find
/// it's less work to skip conversions through [Degrees] when reading or
/// writing [AER] data. If you don't have some overriding conviction, using
/// [Degrees] is a good idea.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct AER<AngularMeasure = Degrees> {
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

impl<AngularMeasure> AER<AngularMeasure>
where
    ENU: From<AER<AngularMeasure>>,
    Self: Copy,
{
    /// Take some observation (`self`), and contextualize it into some
    /// absolute [LLE] given some reference point `obs` that this [AER] is
    /// in reference to.
    pub fn to_lle<CoordinateSystem, LLEAngularMeasure>(
        &self,
        obs: &LLE<CoordinateSystem, LLEAngularMeasure>,
    ) -> LLE<CoordinateSystem, LLEAngularMeasure>
    where
        Radians: From<LLEAngularMeasure> + Copy,
        LLEAngularMeasure: From<Radians> + Copy,
        CoordinateSystem: crate::CoordinateSystem<LLEAngularMeasure>,
    {
        let enu: ENU = (*self).into();
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

impl From<AER<Degrees>> for AER<Radians> {
    /// Convert an [AER] in [Degrees] into [Radians]
    fn from(aer: AER<Degrees>) -> Self {
        Self {
            azimuth: aer.azimuth.into(),
            elevation: aer.elevation.into(),
            range: aer.range,
        }
    }
}

impl From<AER<Radians>> for AER<Degrees> {
    /// Convert an [AER] in [Radians] into [Degrees].
    fn from(aer: AER<Radians>) -> Self {
        Self {
            azimuth: aer.azimuth.into(),
            elevation: aer.elevation.into(),
            range: aer.range,
        }
    }
}

// Implement conversions from AER -> ENU; for Radians and Degrees. Again,
// this can't be done easily using generics, so I'm just going to make it
// explicit.

impl From<AER<Radians>> for ENU {
    /// Convert to [ENU] cartesian local tangent plane coordinates from an
    /// [AER] local tangent plane angular measurement in [Radians].
    fn from(aer: AER<Radians>) -> Self {
        let az_rad: Radians = aer.azimuth;
        let el_rad: Radians = aer.elevation;
        let r = Meters::new(aer.range.as_float() * el_rad.as_float().cos());
        ENU {
            east: Meters::new(r.as_float() * az_rad.as_float().sin()),
            north: Meters::new(r.as_float() * az_rad.as_float().cos()),
            up: Meters::new(aer.range.as_float() * el_rad.as_float().sin()),
        }
    }
}

impl From<AER<Degrees>> for ENU {
    /// Convert to [ENU] cartesian local tangent plane coordinates from an
    /// [AER] local tangent plane angular measurement in [Degrees].
    fn from(aer: AER<Degrees>) -> Self {
        let aer: AER<Radians> = aer.into();
        aer.into()
    }
}

/// [ENU] is East, North, Up in Meters.
///
/// East-North-Up are cartesian coordinates rooted at some reference
/// point's local tangent plane. That reference location is not included in
/// the [ENU] struct, so re-contextualizing an [ENU] to a [LLE] or [XYZ]
/// requires knowing the point at which the observation was taken from or
/// was in reference to.
///
/// These measures are cartesian in the local tangent plane, which is to say,
/// increasing "North" will continue straight ahead, but not following the
/// curve of the Earth -- so you'll slowly get further and further away from
/// Earth's surface.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ENU {
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

impl ENU {
    /// Take some observation (`self`), and contextualize it into some
    /// absolute [LLE] given some reference point `obs` that this [AER]
    /// is in reference to.
    pub fn to_lle<CoordinateSystem, LLEAngularMeasure>(
        &self,
        obs: &LLE<CoordinateSystem, LLEAngularMeasure>,
    ) -> LLE<CoordinateSystem, LLEAngularMeasure>
    where
        Radians: From<LLEAngularMeasure> + Copy,
        LLEAngularMeasure: From<Radians> + Copy,
        CoordinateSystem: crate::CoordinateSystem<LLEAngularMeasure>,
    {
        CoordinateSystem::enu_to_lle(obs, self)
    }
}

// Implement conversions from AER -> ENU; for Radians and Degrees. Again,
// this can't be done easily using generics, so I'm just going to make it
// explicit.

impl From<ENU> for AER<Radians> {
    /// Convert to [AER] cartesian local tangent plane angular measurement in
    /// [Radians] from [ENU] local tangent plane coordinates.
    fn from(enu: ENU) -> Self {
        let r = (enu.east.as_float() * enu.east.as_float()
            + enu.north.as_float() * enu.north.as_float())
        .sqrt();

        let tau = std::f64::consts::PI * 2.0;
        Self {
            azimuth: Radians::new(enu.east.as_float().atan2(enu.north.as_float()) % tau).into(),
            elevation: Radians::new(enu.up.as_float().atan2(r)).into(),
            range: Meters::new((r * r + enu.up.as_float() * enu.up.as_float()).sqrt()),
        }
    }
}

impl From<ENU> for AER<Degrees> {
    /// Convert to [AER] cartesian local tangent plane angular measurement in
    /// [Degrees] from [ENU] local tangent plane coordinates.
    fn from(enu: ENU) -> Self {
        let aer: AER<Radians> = enu.into();
        aer.into()
    }
}

/// [CoordinateSystem] is a trait to enable converstion between locations
/// (usually Latitude and Longitude, in the form of [LLE] objects) to absolute
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
    /// Convert [XYZ] to [LLE] in our [CoordinateSystem].
    fn xyz_to_lle(c: &XYZ) -> LLE<Self, AngularMeasure>;

    /// Convert [XYZ] to [ENU], referenced to `refr` ([LLE]) in our
    /// [CoordinateSystem].
    fn xyz_to_enu(refr: &LLE<Self, AngularMeasure>, c: &XYZ) -> ENU;

    /// Convert [ENU] to [XYZ], referenced to `refr` ([LLE]) in our
    /// [CoordinateSystem].
    fn enu_to_xyz(refr: &LLE<Self, AngularMeasure>, c: &ENU) -> XYZ;

    /// Convert [LLE] to [XYZ] in our [CoordinateSystem].
    fn lle_to_xyz(geo: &LLE<Self, AngularMeasure>) -> XYZ;

    /// Convert [ENU] to [LLE], referenced to `refr` ([LLE]) in our
    /// [CoordinateSystem].
    fn enu_to_lle(refr: &LLE<Self, AngularMeasure>, c: &ENU) -> LLE<Self, AngularMeasure> {
        let xyz = Self::enu_to_xyz(refr, c);
        Self::xyz_to_lle(&xyz)
    }

    /// Convert [LLE] to [ENU], referenced to `refr` ([LLE]) in our
    /// [CoordinateSystem].
    fn lle_to_enu(r: &LLE<Self, AngularMeasure>, geo: &LLE<Self, AngularMeasure>) -> ENU {
        let xyz = Self::lle_to_xyz(geo);
        Self::xyz_to_enu(r, &xyz)
    }
}

impl<FromCoordinateSystem, FromAngularMeasure> LLE<FromCoordinateSystem, FromAngularMeasure>
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
    /// that <code>XYZ</code>'s <code>(0.0, 0.0, 0.0)</code> is the same
    /// exact point in space for all the <code>CoordinateSystem</code>s,
    /// although this library <u>will assume they are</u>. As such, I'm not
    /// confident these conversions are implemented properly, and
    /// <b>must not be relied upon for anything approaching important</b>.
    /// </div>
    pub fn translate<ToCoordinateSystem, ToAngularMeasure>(
        &self,
    ) -> LLE<ToCoordinateSystem, ToAngularMeasure>
    where
        Radians: From<ToAngularMeasure> + Copy,
        ToAngularMeasure: From<Radians> + Copy,
        ToCoordinateSystem: CoordinateSystem<ToAngularMeasure>,
    {
        let xyz = FromCoordinateSystem::lle_to_xyz(self);
        ToCoordinateSystem::xyz_to_lle(&xyz)
    }
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
/// very likely case these are [WGS84] coordinates) -- and absolutely won't
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

    let (self_lat, self_lon) = self_lat_lon;
    let (self_lat, self_lon): (Radians, Radians) = (self_lat.into(), self_lon.into());
    let (self_lat, self_lon) = (self_lat.as_float(), self_lon.as_float());

    let (other_lat, other_lon) = other_lat_lon;
    let (other_lat, other_lon): (Radians, Radians) = (other_lat.into(), other_lon.into());
    let (other_lat, other_lon) = (other_lat.as_float(), other_lon.as_float());

    let delta_lat = self_lat - other_lat;
    let delta_lon = self_lon - other_lon;

    let a = (delta_lat / 2.0).sin().powf(2.0)
        + self_lat.cos() * other_lat.cos() * (delta_lon / 2.0).sin().powf(2.0);
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());

    Meters::new(EARTH_RADIUS.as_float() * c)
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
                panic!();
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
}

// vim: foldmethod=marker
