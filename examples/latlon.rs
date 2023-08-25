use geoconv::CoordinateSystem;
use geoconv::{Degrees, Meters, LLE, WGS84};

fn main() {
    let tea_party = LLE {
        latitude: Degrees::new(42.352211),
        longitude: Degrees::new(-71.051315),
        elevation: Meters::new(0.0),
    };
    let georges_island = LLE {
        latitude: Degrees::new(42.320239),
        longitude: Degrees::new(-70.929482),
        elevation: Meters::new(100.0),
    };

    let look = WGS84::lle_to_enu(tea_party, georges_island).to_aer();
    println!(
        "Azimuth: {} deg, Elevation: {} deg, Range: {} meters",
        look.azimuth.to_float(),
        look.elevation.to_float(),
        look.range.to_float()
    );
}
