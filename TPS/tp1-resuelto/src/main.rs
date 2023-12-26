use core::f64;
use plotters::prelude::*;
use roots::SimpleConvergency;
use std::{fmt, thread};

const INI_INTER: f64 = 1.0;
const FIN_INTER: f64 = 1.5;
// const TOLERANCIA: f64 = 1e-5;
const TOLERANCIA: f64 = 1e-13;
// const MAX_ITER: u8 = 100;
const MAX_ITER: u8 = 200;
const GRAPH_RES: (u32, u32) = (1600, 800);

#[derive(Debug)]
enum Metodos {
    Biseccion,
    PuntoFijo,
    Secante,
    NewtonRaphson,
    NewtonRaphsonModificado,
}

impl TryFrom<i32> for Metodos {
    type Error = ();

    fn try_from(v: i32) -> Result<Self, Self::Error> {
        match v {
            x if x == 0 => Ok(Metodos::Biseccion),
            x if x == 1 => Ok(Metodos::PuntoFijo),
            x if x == 2 => Ok(Metodos::Secante),
            x if x == 3 => Ok(Metodos::NewtonRaphson),
            x if x == 4 => Ok(Metodos::NewtonRaphsonModificado),
            _ => Err(()),
        }
    }
}

impl fmt::Display for Metodos {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Metodos::Biseccion => write!(f, "Biseccion"),
            Metodos::PuntoFijo => write!(f, "Punto Fijo"),
            Metodos::Secante => write!(f, "Secante"),
            Metodos::NewtonRaphson => write!(f, "Newton Raphson"),
            Metodos::NewtonRaphsonModificado => write!(f, "Newton Raphson Modificado"),
        }
    }
}

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

fn biseccion() -> (Metodos, Vec<f64>) {
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
        // println!("iteración {i} = {p}");
    }

    // println!("BISECCION = {}", iteraciones.last().unwrap());
    (Metodos::Biseccion, iteraciones)
}

fn g(x: f64) -> f64 {
    (19.0 - f64::powf(2.0, x)).log(8.0)
}

fn punto_fijo() -> (Metodos, Vec<f64>) {
    let mut iteraciones: Vec<f64> = Vec::new();
    let mut p_n = (INI_INTER + FIN_INTER) / 2.0; // Punto intermedio para la semilla

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = g(p_n);

        iteraciones.push(p_n1);

        // println!("iteración {i} = {p_n1}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    // println!("PUNTO FIJO = {}", iteraciones.last().unwrap());
    (Metodos::PuntoFijo, iteraciones)
}

fn secante() -> (Metodos, Vec<f64>) {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n2 = INI_INTER;
    let mut p_n1 = FIN_INTER;

    if f(p_n2) * f(p_n1) > 0.0 {
        println!("No hay cambio de signo. No hay raíz.");
        return (Metodos::Secante, iteraciones);
    }

    for i in 0..MAX_ITER {
        let p_n = p_n1 - (f(p_n1) * (p_n1 - p_n2)) / (f(p_n1) - f(p_n2));

        iteraciones.push(p_n);

        // println!("SECANTE --> [p_{i}] = {p_n}");
        // println!("iteración {i} = {p_n}");

        if (p_n1 - p_n).abs() < TOLERANCIA {
            break;
        }

        p_n2 = p_n1;
        p_n1 = p_n;
    }

    // println!("SECANTE = {}", iteraciones.last().unwrap());

    (Metodos::Secante, iteraciones)
}

fn newton_raphson() -> (Metodos, Vec<f64>) {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n = 0.5;

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = p_n - f(p_n) / f_prima(p_n);

        iteraciones.push(p_n1);

        // println!("NEWTON RAPHSON --> [p_{i}] = {p_n1}");
        // println!("iteración {i} = {p_n}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    // println!("NEWTON RAPHSON = {}", iteraciones.last().unwrap());

    (Metodos::NewtonRaphson, iteraciones)
}

fn newton_raphson_modificado() -> (Metodos, Vec<f64>) {
    let mut iteraciones: Vec<f64> = Vec::new();

    let mut p_n = 1.0;

    iteraciones.push(p_n);

    for i in 0..MAX_ITER {
        let p_n1 = p_n
            - (f(p_n) * f_prima(p_n))
                / (f64::powf(f_prima(p_n), 2.0) - f(p_n) * f_prima_prima(p_n));

        iteraciones.push(p_n1);

        // println!("NEWTON RAPHSON MODIFICADO --> [p_{i}] = {p_n1}");
        // println!("iteración {i} = {p_n}");

        if (p_n - p_n1).abs() < TOLERANCIA {
            break;
        }

        p_n = p_n1;
    }

    // println!(
    //     "NEWTON RAPHSON MODIFICADO = {}",
    //     iteraciones.last().unwrap()
    // );
    (Metodos::NewtonRaphsonModificado, iteraciones)
}

fn orden_de_convergencia(iteraciones: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let mut convergencias: Vec<Vec<f64>> = Vec::new();

    for metodo in iteraciones {
        let mut convergencia: Vec<f64> = Vec::new();

        for i in 0..(metodo.len() - 1) {
            if i < 3 {
                convergencia.push(0.0);
            } else {
                let x_n_mas_uno: f64 = metodo[i + 1];
                let x_n: f64 = metodo[i];
                let x_n_menos_uno: f64 = metodo[i - 1];
                let x_n_menos_dos: f64 = metodo[i - 2];

                let log_numerador: f64 =
                    f64::log10(((x_n_mas_uno - x_n) / (x_n - x_n_menos_uno)).abs());
                let log_denominador: f64 =
                    f64::log10(((x_n - x_n_menos_uno) / (x_n_menos_uno - x_n_menos_dos)).abs());

                if log_denominador.abs() < TOLERANCIA
                    || log_numerador.abs() < TOLERANCIA
                    || log_denominador.is_infinite()
                    || log_numerador.is_infinite()
                {
                    convergencia.push(*convergencia.last().unwrap());
                } else {
                    let alfa: f64 = log_numerador / log_denominador;

                    if alfa.abs() > TOLERANCIA {
                        convergencia.push(alfa);
                    } else {
                        convergencia.push(*convergencia.last().unwrap());
                    }
                }
            }
        }

        convergencias.push(convergencia);
    }

    convergencias
}

fn constante_asintotica(
    iteraciones: &Vec<Vec<f64>>,
    alfa: Vec<f64>,
    raiz: Vec<f64>,
) -> Vec<Vec<f64>> {
    let mut constantes: Vec<Vec<f64>> = Vec::new();

    for (metodo, iteracion) in iteraciones.iter().enumerate() {
        let mut constante: Vec<f64> = Vec::new();

        for i in 0..(iteracion.len() - 1) {
            if i < 2 {
                constante.push(0.0);
            } else {
                let x_n_mas_uno: f64 = iteracion[i + 1];
                let x_n: f64 = iteracion[i];
                let numerador = (x_n_mas_uno - raiz[metodo]).abs();
                let denominador = f64::powf((x_n - raiz[metodo]).abs(), alfa[metodo]);
                if denominador < TOLERANCIA {
                    constante.push(*constante.last().unwrap());
                } else {
                    let cte: f64 = numerador / denominador;
                    if cte.abs() < TOLERANCIA {
                        constante.push(*constante.last().unwrap());
                    } else {
                        constante.push(cte);
                    }
                }
            }
        }

        constantes.push(constante);
    }

    constantes
}

fn grafico_convergencia(
    convergencias: &Vec<Vec<f64>>,
    cant_iteraciones: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let drawing_area = SVGBackend::new("./convergencia.svg", GRAPH_RES).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder
        .margin(10)
        .set_left_and_bottom_label_area_size(20);
    let mut chart_context = chart_builder
        .caption("Gráfico de Convergencia", ("times-new-roman", 40))
        .y_label_area_size(60)
        .x_label_area_size(50)
        .build_cartesian_2d(0f32..((cant_iteraciones) as f32), 0f32..3f32)
        .unwrap();

    chart_context
        .configure_mesh()
        .x_labels(15)
        .y_labels(5)
        .x_desc("número de iteraciones")
        .y_desc("convergencia")
        .axis_desc_style(("times-new-roman", 20))
        .draw()
        .unwrap();

    let colors: Vec<RGBAColor> = vec![
        plotters::style::colors::full_palette::ORANGE_A400.into(), // BISECCION
        plotters::style::colors::full_palette::PINK_A200.into(),   // PUNTO FIJO
        plotters::style::colors::full_palette::DEEPPURPLE_400.into(), // SECANTE
        plotters::style::colors::full_palette::INDIGO.into(),      // NEWTON RAPHSON
        plotters::style::colors::full_palette::TEAL.into(),        // NEWTON RAPHSON MODIFICADO
    ];

    for (metodo, convergencia) in convergencias.iter().enumerate() {
        let eje_x: Vec<u8> = (0..((convergencia.len() - 1) as u8)).collect();

        let color = *colors.get(metodo).unwrap();

        chart_context
            .draw_series(
                LineSeries::new(
                    eje_x.iter().map(|x_i| {
                        (
                            *x_i as f32,
                            *convergencia.get(*x_i as usize).unwrap() as f32,
                        )
                    }),
                    ShapeStyle {
                        color: color.clone(),
                        filled: true,
                        stroke_width: 2,
                    },
                )
                .point_size(5),
            )
            .unwrap()
            .label(format!(
                "{}",
                Metodos::try_from(metodo as i32).unwrap().to_string()
            ))
            .legend(move |(x, y)| {
                Rectangle::new(
                    [(x + 5, y - 5), (x + 15, y + 5)],
                    ShapeStyle {
                        color: color,
                        filled: true,
                        stroke_width: 1,
                    },
                )
            });
    }

    chart_context
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn grafico_constante_asintotica(
    constantes_asintoticas: &Vec<Vec<f64>>,
    cant_iteraciones: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let drawing_area = SVGBackend::new("./constante_asintotica.svg", GRAPH_RES).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder
        .margin(10)
        .set_left_and_bottom_label_area_size(20);

    let mut chart_context = chart_builder
        .caption("Gráfico de Constante Asintótica", ("times-new-roman", 40))
        .y_label_area_size(60)
        .x_label_area_size(50)
        .build_cartesian_2d(0f32..((cant_iteraciones) as f32), 0f32..9f32)
        .unwrap();

    chart_context
        .configure_mesh()
        .x_labels(15)
        .y_labels(5)
        .x_desc("número de iteraciones")
        .y_desc("constante asintótica")
        .axis_desc_style(("times-new-roman", 20))
        .draw()
        .unwrap();

    let colors: Vec<RGBAColor> = vec![
        plotters::style::colors::full_palette::ORANGE_A400.into(), // BISECCION
        plotters::style::colors::full_palette::PINK_A200.into(),   // PUNTO FIJO
        plotters::style::colors::full_palette::DEEPPURPLE_400.into(), // SECANTE
        plotters::style::colors::full_palette::INDIGO.into(),      // NEWTON RAPHSON
        plotters::style::colors::full_palette::TEAL.into(),        // NEWTON RAPHSON MODIFICADO
    ];

    for (metodo, convergencia) in constantes_asintoticas.iter().enumerate() {
        let eje_x: Vec<u8> = (0..((convergencia.len() - 1) as u8)).collect();

        let color = *colors.get(metodo).unwrap();

        chart_context
            .draw_series(
                LineSeries::new(
                    eje_x.iter().map(|x_i| {
                        (
                            *x_i as f32,
                            *convergencia.get(*x_i as usize).unwrap() as f32,
                        )
                    }),
                    ShapeStyle {
                        color: color.clone(),
                        filled: true,
                        stroke_width: 2,
                    },
                )
                .point_size(5),
            )
            .unwrap()
            .label(format!(
                "{}",
                Metodos::try_from(metodo as i32).unwrap().to_string()
            ))
            .legend(move |(x, y)| {
                Rectangle::new(
                    [(x + 5, y - 5), (x + 15, y + 5)],
                    ShapeStyle {
                        color: color,
                        filled: true,
                        stroke_width: 1,
                    },
                )
            });
    }

    chart_context
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn grafico_diferencias_sucesivas(
    iteraciones: &Vec<Vec<f64>>,
    cant_iteraciones: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let drawing_area =
        SVGBackend::new("./diferencias_sucesivas.svg", GRAPH_RES).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder
        .margin(10)
        .set_left_and_bottom_label_area_size(20);

    let mut chart_context = chart_builder
        .caption(
            "Gráfico de Diferencias Sucesivas en Escala Logarítmica",
            ("times-new-roman", 40),
        )
        .y_label_area_size(60)
        .x_label_area_size(50)
        .build_cartesian_2d(0f32..((cant_iteraciones) as f32), -14f32..1f32)
        .unwrap();

    chart_context
        .configure_mesh()
        .x_labels(15)
        .y_labels(5)
        .x_desc("número de iteraciones")
        .y_desc("log10(|Δx|)")
        .axis_desc_style(("times-new-roman", 20))
        .draw()
        .unwrap();

    let colors: Vec<RGBAColor> = vec![
        plotters::style::colors::full_palette::ORANGE_A400.into(), // BISECCION
        plotters::style::colors::full_palette::PINK_A200.into(),   // PUNTO FIJO
        plotters::style::colors::full_palette::DEEPPURPLE_400.into(), // SECANTE
        plotters::style::colors::full_palette::INDIGO.into(),      // NEWTON RAPHSON
        plotters::style::colors::full_palette::TEAL.into(),        // NEWTON RAPHSON MODIFICADO
    ];

    for (metodo, iteracion) in iteraciones.iter().enumerate() {
        let eje_x: Vec<u8> = (1..((iteracion.len() - 1) as u8)).collect();

        let color = *colors.get(metodo).unwrap();

        chart_context
            .draw_series(
                LineSeries::new(
                    eje_x.iter().map(|x_i| {
                        (
                            *x_i as f32,
                            (f64::log10(
                                (iteracion.get((*x_i) as usize).unwrap()
                                    - iteracion.get((*x_i - 1) as usize).unwrap())
                                .abs(),
                            )) as f32,
                        )
                    }),
                    ShapeStyle {
                        color: color.clone(),
                        filled: true,
                        stroke_width: 2,
                    },
                )
                .point_size(5),
            )
            .unwrap()
            .label(format!(
                "{}",
                Metodos::try_from(metodo as i32).unwrap().to_string()
            ))
            .legend(move |(x, y)| {
                Rectangle::new(
                    [(x + 5, y - 5), (x + 15, y + 5)],
                    ShapeStyle {
                        color: color,
                        filled: true,
                        stroke_width: 1,
                    },
                )
            });
    }

    chart_context
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn grafico_error_absoluto(
    iteraciones: &Vec<Vec<f64>>,
    cant_iteraciones: u8,
    raiz_real: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let drawing_area = SVGBackend::new("./error_absoluto.svg", GRAPH_RES).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder
        .margin(10)
        .set_left_and_bottom_label_area_size(20);

    let mut chart_context = chart_builder
        .caption(
            "Gráfico del Error Absoluto en Escala Logarítmica",
            ("times-new-roman", 40),
        )
        .y_label_area_size(60)
        .x_label_area_size(50)
        .build_cartesian_2d(0f32..((cant_iteraciones) as f32), -16f32..1f32)
        .unwrap();

    chart_context
        .configure_mesh()
        .x_labels(15)
        .y_labels(5)
        .x_desc("número de iteraciones")
        .y_desc("log10(|x_candidata - x_real|)")
        .axis_desc_style(("times-new-roman", 20))
        .draw()
        .unwrap();

    let colors: Vec<RGBAColor> = vec![
        plotters::style::colors::full_palette::ORANGE_A400.into(), // BISECCION
        plotters::style::colors::full_palette::PINK_A200.into(),   // PUNTO FIJO
        plotters::style::colors::full_palette::DEEPPURPLE_400.into(), // SECANTE
        plotters::style::colors::full_palette::INDIGO.into(),      // NEWTON RAPHSON
        plotters::style::colors::full_palette::TEAL.into(),        // NEWTON RAPHSON MODIFICADO
    ];

    for (metodo, iteracion) in iteraciones.iter().enumerate() {
        let eje_x: Vec<u8> = (0..((iteracion.len() - 1) as u8)).collect();

        let color = *colors.get(metodo).unwrap();

        chart_context
            .draw_series(
                LineSeries::new(
                    eje_x.iter().map(|x_i| {
                        (
                            *x_i as f32,
                            (f64::log10(
                                (iteracion.get((*x_i) as usize).unwrap() - raiz_real).abs(),
                            )) as f32,
                        )
                    }),
                    ShapeStyle {
                        color: color.clone(),
                        filled: true,
                        stroke_width: 2,
                    },
                )
                .point_size(5),
            )
            .unwrap()
            .label(format!(
                "{}",
                Metodos::try_from(metodo as i32).unwrap().to_string()
            ))
            .legend(move |(x, y)| {
                Rectangle::new(
                    [(x + 5, y - 5), (x + 15, y + 5)],
                    ShapeStyle {
                        color: color,
                        filled: true,
                        stroke_width: 1,
                    },
                )
            });
    }

    chart_context
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn main() {
    let mut iteraciones: Vec<Vec<f64>> = Vec::new();

    let mut handler = vec![];

    let biseccion = thread::spawn(move || biseccion());
    handler.push(biseccion);

    let punto_fijo = thread::spawn(move || punto_fijo());
    handler.push(punto_fijo);

    let secante = thread::spawn(move || secante());
    handler.push(secante);

    let newton_raphson = thread::spawn(move || newton_raphson());
    handler.push(newton_raphson);

    let newton_raphson_modificado = thread::spawn(move || newton_raphson_modificado());
    handler.push(newton_raphson_modificado);

    for h in handler {
        let (met, vec) = h.join().unwrap();
        match met {
            Metodos::Biseccion => iteraciones.insert(0, vec),
            Metodos::PuntoFijo => iteraciones.insert(1, vec),
            Metodos::Secante => iteraciones.insert(2, vec),
            Metodos::NewtonRaphson => iteraciones.insert(3, vec),
            Metodos::NewtonRaphsonModificado => iteraciones.insert(4, vec),
        }
    }

    //// GRAFICOS
    // CONVERGENCIA
    let convergencias: Vec<Vec<f64>> = orden_de_convergencia(&iteraciones);

    let mut max_iter = 0;
    for metodo in &iteraciones {
        if metodo.len() > max_iter {
            max_iter = metodo.len();
        }
    }

    let _ = grafico_convergencia(&convergencias, (max_iter - 1) as u8);

    // CTE ASINTOTICA
    let ultima_convergencia: Vec<f64> = convergencias
        .iter()
        .map(|convergencia| *convergencia.last().unwrap())
        .collect();
    let ultima_raiz: Vec<f64> = iteraciones
        .iter()
        .map(|iteracion| *iteracion.last().unwrap())
        .collect();

    let constantes_asintoticas: Vec<Vec<f64>> =
        constante_asintotica(&iteraciones, ultima_convergencia, ultima_raiz);

    let _ = grafico_constante_asintotica(&constantes_asintoticas, (max_iter - 1) as u8);

    // DIFERENCIAS SUCESIVAS - ESCALA LOGARITMICA
    let _ = grafico_diferencias_sucesivas(&iteraciones, (max_iter - 1) as u8);

    // ERROR ABSOLUTO - ESCALA LOGARITMICA
    let mut convergency = SimpleConvergency {
        eps: TOLERANCIA,
        max_iter: (MAX_ITER as usize),
    };

    let raiz_real =
        roots::find_root_regula_falsi(INI_INTER, FIN_INTER, &f, &mut convergency).unwrap();

    // println!("raiz real: {raiz_real}");

    let _ = grafico_error_absoluto(&iteraciones, (max_iter - 1) as u8, raiz_real);
}
