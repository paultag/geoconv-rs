use geoconv::{Aer, Degrees, Lle, Meters, Wgs84};

fn main() {
    let tea_party = Lle::<Wgs84>::new(
        Degrees::new(42.352211),
        Degrees::new(-71.051315),
        Meters::new(0.0),
    );
    let georges_island = Lle::<Wgs84>::new(
        Degrees::new(42.320239),
        Degrees::new(-70.929482),
        Meters::new(100.0),
    );

    let look: Aer<Degrees> = tea_party.aer_to(&georges_island);
    println!(
        "Azimuth: {} deg, Elevation: {} deg, Range: {} meters",
        look.azimuth.as_float(),
        look.elevation.as_float(),
        look.range.as_float()
    );
}
