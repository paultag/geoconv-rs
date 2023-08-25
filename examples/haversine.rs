use geoconv::{Degrees, Meters, LLE};

fn main() {
    let tea_party = LLE {
        latitude: Degrees::new(42.352211),
        longitude: Degrees::new(-71.051315),
        elevation: Meters::new(0.0),
    };
    let georges_island = LLE {
        latitude: Degrees::new(42.320239),
        longitude: Degrees::new(-70.929482),
        elevation: Meters::new(0.0),
    };

    let distance = tea_party
        .haversine_distance(georges_island)
        .expect("haversine");
    println!("Distance: {} meters", distance.to_float());
}
