
use geoconv::{Degrees, Meters, LLE};

fn main() {
    let tea_party = LLE{
        latitude: Degrees(42.352211),
        longitude: Degrees(-71.051315),
        elevation: Meters(0.0),
    };
    let georges_island = LLE{
        latitude: Degrees(42.320239),
        longitude: Degrees(-70.929482),
        elevation: Meters(0.0),
    };

    let distance = tea_party.haversine_distance(georges_island)
        .expect("haversine");
    println!("Distance: {} meters", distance.0);
}
