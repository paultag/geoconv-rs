use geoconv::{Aer, Degrees, Lle, Meters, Wgs84};

const GOES_16: Lle<Wgs84> = Lle::<Wgs84>::new(
    Degrees::new(0.0),
    Degrees::new(-75.2),
    Meters::new(35_780_200.0),
);

const WHITE_HOUSE: Lle<Wgs84> = Lle::<Wgs84>::new(
    Degrees::new(42.352211),
    Degrees::new(-71.051315),
    Meters::new(0.0),
);

fn main() {
    let look: Aer<Degrees> = WHITE_HOUSE.aer_to(&GOES_16);
    println!(
        "Someone at the White House needs to look at Azimuth: {} deg, Elevation: {} deg to see GOES 16",
        look.azimuth.as_float(),
        look.elevation.as_float(),
    );
}
