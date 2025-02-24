use geoconv::{Degrees, haversine_distance};

fn main() {
    let tea_party = (Degrees::new(42.352211), Degrees::new(-71.051315));
    let georges_island = (Degrees::new(42.320239), Degrees::new(-70.929482));

    let distance = haversine_distance(tea_party, georges_island);
    println!("Distance: {} meters", distance.as_float());
}
