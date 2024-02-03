use geoconv::{Degrees, Meters, Radians, AER, LLE, WGS72, WGS84};

fn main() {
    let tea_party = LLE::<WGS84>::new(
        Degrees::new(42.352211),
        Degrees::new(-71.051315),
        Meters::new(0.0),
    );
    let somewhere =
        LLE::<WGS72, Radians>::new(Radians::new(0.3), Radians::new(2.2), Meters::new(10.0));

    // Our goal:
    //
    // Take a WGS84 Lat/Lon in Degrees,and a WGS72 Lat/Lon in Radians, and
    // compute the look angle in the local tangent plane in Degrees.

    let look: AER<Degrees> = tea_party.aer_to(&somewhere.translate());

    println!(
        "A/E/R: {:?} {:?} {:?}",
        look.azimuth, look.elevation, look.range
    );
}
