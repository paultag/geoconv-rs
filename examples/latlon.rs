
use geoconv::CoordinateSystem;
use geoconv::{Degrees, Meters, WGS84, LLE};

fn main() {
    let tea_party = LLE{
        latitude: Degrees(42.352211),
        longitude: Degrees(-71.051315),
        elevation: Meters(0.0),
    };
    let georges_island = LLE{
        latitude: Degrees(42.320239),
        longitude: Degrees(-70.929482),
        elevation: Meters(100.0),
    };

    let look = WGS84::lle_to_enu(tea_party, georges_island).to_aer();
    println!("Azimuth: {} deg, Elevation: {} deg, Range: {} meters",
             look.azimuth.0,
             look.elevation.0,
             look.range.0);
}
