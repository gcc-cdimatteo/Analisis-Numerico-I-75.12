use core::f64;
use plotters::prelude::*;
use roots::SimpleConvergency;
use std::thread;

const INI_INTER: f64 = 0.0;
const FIN_INTER: f64 = 3.0;
// const TOLERANCIA: f64 = 1e-5;
const TOLERANCIA: f64 = 1e-13;
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

fn orden_de_convergencia(iteraciones: Vec<f64>) -> Vec<f64> {
    let mut convergencia: Vec<f64> = Vec::new();

    for i in 0..(iteraciones.len() - 1) {
        if i < 3 {
            convergencia.push(0.0);
        } else {
            let x_n_mas_uno: f64 = iteraciones[i + 1];
            let x_n: f64 = iteraciones[i];
            let x_n_menos_uno: f64 = iteraciones[i - 1];
            let x_n_menos_dos: f64 = iteraciones[i - 2];

            let log_numerador: f64 =
                f64::log10(((x_n_mas_uno - x_n) / (x_n - x_n_menos_uno)).abs());
            let log_denominador: f64 =
                f64::log10(((x_n - x_n_menos_uno) / (x_n_menos_uno - x_n_menos_dos)).abs());

            let alfa: f64 = log_numerador / log_denominador;

            if alfa > TOLERANCIA {
                convergencia.push(alfa);
            } else {
                convergencia.push(*convergencia.last().unwrap());
            }
        }
    }
    convergencia
}

fn constante_asintotica(iteraciones: Vec<f64>, alfa: f64, raiz: f64) -> Vec<f64> {
    let mut constantes: Vec<f64> = Vec::new();

    for i in 0..(iteraciones.len() - 1) {
        if i < 2 {
            constantes.push(0.0);
        } else {
            let x_n_mas_uno: f64 = iteraciones[i + 1];
            let x_n: f64 = iteraciones[i];
            let numerador = (x_n_mas_uno - raiz).abs();
            let denominador = f64::powf((x_n - raiz).abs(), alfa);
            let constante: f64 = numerador / denominador;
            constantes.push(constante);
        }
    }
    constantes
}

fn grafico_convergencia(
    convergencia: Vec<f64>,
    iteraciones: Vec<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    let x_values = iteraciones;
    let x_values_len = x_values.len() as f64;
    let drawing_area = SVGBackend::new("./plot_test.svg", (600, 500)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder
        .margin(10)
        .set_left_and_bottom_label_area_size(20);
    let mut chart_context = chart_builder
        .build_cartesian_2d(0.0..4.0, 0.0..3.0)
        .unwrap();
    chart_context.configure_mesh().draw().unwrap();
    // chart_context
    //     .draw_series(
    //         LineSeries::new(
    //             (0.0..(x_values.len() as f64)).map(|i| (i, convergencia.get(i))),
    //             RED,
    //         )
    //         .point_size(5),
    //     )
    //     // (0..x_values.len()).iter().map(|i| (i, convergencia.get(i)))
    //     .unwrap();
    Ok(())
}

fn main() {
    // let iteraciones: Vec<f64> = biseccion();
    // let iteraciones: Vec<f64> = punto_fijo();
    // let iteraciones: Vec<f64> = secante();
    let iteraciones: Vec<f64> = newton_raphson();
    // let iteraciones: Vec<f64> = newton_raphson_modificado();
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

    // GRAFICOS
    let raiz: f64 = *iteraciones.last().unwrap();
    let convergencia: Vec<f64> = orden_de_convergencia(iteraciones.clone());
    // let _ = grafico_convergencia(convergencia, iteraciones.clone());
    let orden_conv: f64 = *convergencia.last().unwrap();
    let constante_asin: Vec<f64> = constante_asintotica(iteraciones, orden_conv, raiz);

    println!("CONSTANTES: {:?}", constante_asin);

    // let mut convergency = SimpleConvergency {
    //     eps: TOLERANCIA,
    //     max_iter: 200,
    // };

    // let r = roots::find_root_secant(INI_INTER, FIN_INTER, &f, &mut convergency);

    // println!("Root found by crate: {}", r.unwrap());
}
