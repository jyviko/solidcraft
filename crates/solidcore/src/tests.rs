use super::*;

/// The bound calls on SM_API's primitives, against their closed forms.
#[test]
fn box_end_to_end() {
    let b = Brep::make_box([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]).unwrap();
    assert!((b.volume(1e-7).unwrap() - 6000.0).abs() < 1e-6);
    assert_eq!(b.counts().unwrap(), [6, 12, 8]);
    assert_eq!(b.bounds().unwrap(), ([1.0, 2.0, 3.0], [11.0, 22.0, 33.0]));
    assert!(b.is_manifold_solid().unwrap());
    assert_eq!(b.solid_count().unwrap(), 1);
    let c = b.try_clone().unwrap();
    drop(b);
    assert!((c.volume(1e-7).unwrap() - 6000.0).abs() < 1e-6);
    assert_eq!(current_call(), None);
}

#[test]
fn every_primitive_measures_its_closed_form() {
    use std::f64::consts::PI;
    let o = [1.0, 2.0, 3.0];
    let cases = [
        ("cylinder", Brep::make_cylinder(o, 15.0, 50.0).unwrap(), PI * 15.0 * 15.0 * 50.0),
        ("sphere", Brep::make_sphere(o, 12.0).unwrap(), 4.0 / 3.0 * PI * 12.0f64.powi(3)),
        ("frustum", Brep::make_cone(o, 10.0, 4.0, 20.0).unwrap(), PI * 20.0 / 3.0 * (100.0 + 40.0 + 16.0)),
        ("cone", Brep::make_cone(o, 10.0, 0.0, 20.0).unwrap(), PI * 100.0 * 20.0 / 3.0),
        ("torus", Brep::make_torus(o, 22.0, 4.0).unwrap(), 2.0 * PI * PI * 22.0 * 16.0),
    ];
    for (name, b, exact) in &cases {
        let v = b.volume(MEASURE_ACCURACY).unwrap();
        assert!((v / exact - 1.0).abs() < 1e-6, "{name}: {v} against {exact}");
        assert_eq!(b.solid_count().unwrap(), 1, "{name}");
    }
}

#[test]
fn bad_arguments_are_errors() {
    assert_eq!(Brep::make_box([0.0; 3], [1.0, -1.0, 1.0]).unwrap_err(), Error::Argument);
    assert_eq!(Brep::make_box([f64::NAN, 0.0, 0.0], [1.0; 3]).unwrap_err(), Error::Argument);
    assert_eq!(Brep::make_cone([0.0; 3], 2.0, 2.0, 1.0).unwrap_err(), Error::Argument);
    assert_eq!(Brep::make_torus([0.0; 3], 2.0, 3.0).unwrap_err(), Error::Argument);
    let b = Brep::make_box([0.0; 3], [1.0; 3]).unwrap();
    assert_eq!(b.volume(f64::INFINITY).unwrap_err(), Error::Argument);
}

/// What the kernel cannot do yet is refused by name, never answered.
#[test]
fn unported_calls_are_refused() {
    let b = Brep::make_box([0.0; 3], [1.0; 3]).unwrap();
    let record = solidcraft_geom::brep::Brep {
        version: solidcraft_geom::brep::VERSION,
        tolerance: 1e-6,
        regions: Vec::new(),
        shells: Vec::new(),
        faces: Vec::new(),
        edges: Vec::new(),
        vertices: Vec::new(),
    };
    assert_eq!(Brep::from_record(&record, true).unwrap_err(), Error::NotBuilt);
    assert_eq!(b.to_native().unwrap_err(), export::ExportError::Smlib(Error::NotBuilt));
    assert_eq!(b.boolean(&b, BooleanOp::Union).unwrap_err(), Error::NotBuilt);
    assert_eq!(b.fillet(&[[0.0; 3]], 1e-3, 0.1).unwrap_err(), Error::NotBuilt);
    assert_eq!(b.check(CheckLevel::All).unwrap_err(), Error::NotBuilt);
    assert!(ssi::Surface::plane([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0], [0.0, 1.0]).is_err());
    assert!(Error::NotBuilt.to_string().starts_with("not supported yet"));
    assert_eq!(set_dump_dir(None), Ok(()));
    assert_eq!(set_dump_dir(Some(Path::new("/nowhere"))), Ok(()));
}
