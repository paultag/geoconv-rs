use geoconv::{Degrees, Meters, Radians, Wgs72, Wgs84, AER, LLE};

fn main() {
    let tea_party = LLE::<Wgs84>::new(
        Degrees::new(42.352211),
        Degrees::new(-71.051315),
        Meters::new(0.0),
    );
    let somewhere =
        LLE::<Wgs72, Radians>::new(Radians::new(0.3), Radians::new(2.2), Meters::new(10.0));

    // Our goal:
    //
    // Take a Wgs84 Lat/Lon in Degrees,and a Wgs72 Lat/Lon in Radians, and
    // compute the look angle in the local tangent plane in Degrees.

    let look: AER<Degrees> = tea_party.aer_to(&somewhere.translate());

    println!(
        "A/E/R: {:?} {:?} {:?}",
        look.azimuth, look.elevation, look.range
    );
}
