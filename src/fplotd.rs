use gnuplot::{Figure, Caption, Color, Graph, AxesCommon, PointSymbol, PointSize, LineWidth, Font, Fix, TickOption::{MajorScale, MinorScale}, AutoOption::Auto};

fn fname(mut file: String) -> String
{
    let mut i: u64 = 0;
    while std::path::Path::new(&file).exists()
    {
        i += 1;
        file = format!("eafit_{}.pdf",i);    
    }
    file
}

pub fn plotd(time: &Vec<f64>, data: &Vec<f64>, timeprof: &Vec<f64>, profile: &Vec<f64>, model: &Vec<f64>, savectrl: bool, plotctrl: bool, syms: f64)
{
    let mut fig = Figure::new();
    if std::env::consts::OS.eq("linux") { _ = fig.set_pre_commands(&format!("set term x11")); }
    let ax = fig.axes2d();

    ax.set_x_label("TBJD [days]", &[Font("Arial", 13.0)]);
    ax.set_y_label("Normalized flux", &[Font("Arial", 13.0)]);
    ax.set_legend(Graph(1.0), Graph(1.0), &[], &[]);

    ax.points(time,data,&[Color("#00000"),PointSymbol('O'),PointSize(syms),Caption("Data points")]);
    ax.lines(timeprof,profile,&[Caption("Flare model"), Color("#2ca02c"),LineWidth(3.0)]);
    ax.lines(timeprof,model,&[Caption("Modulated flare model"), Color("#d62728"),LineWidth(3.0)]);
    ax.set_x_range(Fix(time[0]),Fix(time[time.len()-1]));
    ax.set_x_ticks(Some((Auto, 4)), &[MajorScale(1.5), MinorScale(0.75)],&[]);
    ax.set_y_ticks(Some((Auto, 4)), &[MajorScale(1.5), MinorScale(0.75)],&[]);
    if savectrl { _ = fig.save_to_pdf(fname(String::from("eafit_0.pdf")), 10, 6) }
    if plotctrl { _ = fig.show(); }
}