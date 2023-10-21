use roots::SimpleConvergency;
use std::thread;

const INI_INTER: f64 = 0.0;
const FIN_INTER: f64 = 3.0;
const TOLERANCIA: f64 = 1e-5;
// const TOLERANCIA: f64 = 1e-13;
// const MAX_ITER: u8 = 100;
const MAX_ITER: u8 = 200;

fn f(x: f64) -> f64 {
    f64::powf(2.0, x) + f64::powf(8.0, x) - 19.0
}

fn f_prima(x: f64) -> f64 {
    f64::powf(2.0, x) * f64::ln(2.0) + 3.0 * f64::ln(2.0) * f64::powf(8.0, x)
}

fn f_prima_prima(x: f64) -> f64 {
    f64::powf(
        f64::powf(2.0, x) * f64::log10(2.0) + f64::powf(8.0, x) * f64::log10(8.0),
        2.0,
    )
}

fn biseccion() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut a = INI_INTER;
    let mut b = FIN_INTER;
    for i in 0..MAX_ITER {
        let p = (a + b) / 2.0;
        let f_p = f(p);

        if iteraciones.len() > 0 && (iteraciones.last().unwrap() - p).abs() < TOLERANCIA {
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
        println!("BISECCION --> [p_{i}] = {p}");
    }

    println!("BISECCION = {}", iteraciones.last().unwrap());
    iteraciones
}

fn g(x: f64) -> f64 {
    (19.0 - f64::powf(2.0, x)).log(8.0)
}

fn punto_fijo() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();
    let mut p_n = (INI_INTER + FIN_INTER) / 2.0; // Punto intermedio para la semilla

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = g(p_n);

        iteraciones.push(p_n1);

        println!("PUNTO FIJO --> [p_{i}] = {p_n1}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    println!("PUNTO FIJO = {}", iteraciones.last().unwrap());

    iteraciones
}

fn secante() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n2 = INI_INTER;
    let mut p_n1 = FIN_INTER;

    if f(p_n2) * f(p_n1) > 0.0 {
        println!("No hay cambio de signo. No hay raíz.");
        return iteraciones;
    }

    for i in 0..MAX_ITER {
        let p_n = p_n1 - (f(p_n1) * (p_n1 - p_n2)) / (f(p_n1) - f(p_n2));

        iteraciones.push(p_n);

        println!("SECANTE --> [p_{i}] = {p_n}");

        if (p_n1 - p_n).abs() < TOLERANCIA {
            break;
        }

        p_n2 = p_n1;
        p_n1 = p_n;
    }

    println!("SECANTE = {}", iteraciones.last().unwrap());

    iteraciones
}

fn newton_raphson() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n = 0.5;

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = p_n - f(p_n) / f_prima(p_n);

        iteraciones.push(p_n1);

        println!("NEWTON RAPHSON --> [p_{i}] = {p_n1}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    println!("NEWTON RAPHSON = {}", iteraciones.last().unwrap());

    iteraciones
}

fn newton_raphson_modificado() -> Vec<f64> {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n = 1.0;

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = p_n
            - (f(p_n) * f_prima(p_n))
                / (f64::powf(f_prima(p_n), 2.0) - f(p_n) * f_prima_prima(p_n));

        iteraciones.push(p_n1);

        println!("NEWTON RAPHSON MODIFICADO --> [p_{i}] = {p_n1}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    println!(
        "NEWTON RAPHSON MODIFICADO = {}",
        iteraciones.last().unwrap()
    );
    iteraciones
}

fn main() {
    biseccion();
    punto_fijo();
    secante();
    newton_raphson();
    newton_raphson_modificado();
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

    // for h in handler {
    //     let _ = h.join();
    // }

    let mut convergency = SimpleConvergency {
        eps: TOLERANCIA,
        max_iter: 200,
    };

    let r = roots::find_root_secant(INI_INTER, FIN_INTER, &f, &mut convergency);

    println!("Root found by crate: {}", r.unwrap());
}
