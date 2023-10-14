use std::{sync::Arc, thread};

const INI_INTER: f64 = -5.0;
const FIN_INTER: f64 = -3.5;
const TOLERANCIA: f64 = 10e-6;
const MAX_ITER: u8 = 100;

fn f(x: f64) -> f64 {
    f64::powf(2.0, x) + f64::powf(8.0, x) - 0.05263157895
}

fn biseccion() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut a = INI_INTER;
    let mut b = FIN_INTER;
    for _ in 0..MAX_ITER {
        let p = (a + b) / 2.0;
        let f_p = f(p);

        if f_p.abs() < TOLERANCIA
            || (iteraciones.len() > 0 && (iteraciones.last().unwrap() - f_p).abs() < TOLERANCIA)
        {
            break;
        }

        if f(a) * f_p < 0.0 {
            // la raiz esta entre a y p
            b = p;
        } else {
            // la raiz esta entre p y b
            a = p;
        }
        iteraciones.push(p);
        println!("{p}");
    }

    iteraciones
}

fn punto_fijo() {
    todo!()
}

fn secante() {
    todo!()
}

fn newton_raphson() {
    todo!()
}

fn newton_raphson_modificado() {
    todo!()
}

fn main() {
    biseccion();
    // let mut handler = vec![];

    // let biseccion = thread::spawn(move || biseccion());
    // handler.push(biseccion);

    // let punto_fijo = thread::spawn(move || punto_fijo());
    // handler.push(punto_fijo);

    // let secante = thread::spawn(move || secante());
    // handler.push(secante);

    // let newton_raphson = thread::spawn(move || newton_raphson());
    // handler.push(newton_raphson);

    // let newton_raphson_modificado = thread::spawn(move || newton_raphson_modificado());
    // handler.push(newton_raphson_modificado);

    // let _: Vec<()> = handler.into_iter().flat_map(|x| x.join()).collect();
}
