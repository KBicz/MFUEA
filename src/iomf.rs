use crate::functions;

pub fn read_lcprof(path: &str, mut nparams: usize, fitprofile: bool) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, f64, usize, bool, Vec<f64>, usize) 
{
    let profile: Vec<f64>;
    let mut line: Vec<&str>;
    let mut params: Vec<f64> = vec![];
    let datavec: String = std::fs::read_to_string(path).expect("# Unable to read file!");
    let mut datavec: Vec<&str> = datavec.lines().collect();
    let fline: Vec<&str> = datavec[0].split(" ").collect();
    let sline: Vec<&str> = datavec[1].split(" ").collect();
    let (istart, istop, t0): (usize, usize, f64);
    let (mut timeprof, mut newtime): (Vec<f64>, Vec<f64>)  = (vec![], vec![]);

    datavec = datavec.into_iter().filter(|&i| i.trim() != "" && i.chars().next().unwrap() != '#').collect::<Vec<_>>();
    
    let (mut flux, mut model): (f64, f64);
    let (mut time, mut err, mut trend): (Vec<f64>,Vec<f64>,Vec<f64>) = (vec![],vec![],vec![]);

    for i in 0..datavec.len()
    {
        line = datavec[i].split(" ").collect();
        line = line.into_iter().filter(|&i| i.trim() != "").collect::<Vec<_>>();

        time.push((&line[0]).parse().unwrap());
        flux = (&line[1]).parse().unwrap(); flux = flux/1000.0+1.0;
        err.push((&line[2]).parse().unwrap()); err[i] /= 1000.0;
        model = (&line[3]).parse().unwrap(); model = model/1000.0+1.0;
        trend.push(flux/model-1.0);
    }

    (istart, istop) = ((&fline[1]).parse().unwrap(),(&fline[2]).parse().unwrap());
    t0 = {
        if fitprofile { time[istart] }
        else { (&fline[3]).parse().unwrap() }
    };
    for i in istart..istop { newtime.push((time[i]-t0)*24.0*60.0); timeprof.push(time[i]);}

    if fitprofile
    {
        profile = vec![0.0; newtime.len()];
    }
    else 
    {  
        for param in sline.clone() { if !param.eq("#") { params.push(param.parse().unwrap()); } }
        nparams = params.len();
        profile = {
            match nparams
            {
                4 => {functions::profil_wroclawski(&newtime, params[0], params[1], params[2], params[3])},
                7 => {functions::profil_wroclawski_1b(&newtime, params[0], params[1], params[2], params[3], params[4], params[5], params[6])},
                8 => {functions::profil_wroclawski_2b(&newtime, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7])},
                _ => { std::process::exit(1); },
            }
        };
    }   
    
    (time, trend, err, profile, timeprof, t0, istart, fitprofile, newtime, nparams)
}

pub fn read_response(path: &str) -> (Vec<f64>, Vec<f64>)
{
    let mut line: Vec<&str>;
    let data: String = std::fs::read_to_string(path).expect("# Unable to read file!");
    let mut datavec: Vec<&str> = data.lines().collect();
    datavec = datavec.into_iter().filter(|&i| i.trim() != "" && i.chars().next().unwrap() != '#').collect::<Vec<_>>();
    let (mut wave, mut resp): (Vec<f64>, Vec<f64>) = (vec![], vec![]);

    for i in 0..datavec.len()
    {
        line = datavec[i].split(" ").collect();
        line = line.into_iter().filter(|&i| i.trim() != "").collect::<Vec<_>>();

        wave.push((&line[0]).parse().unwrap()); wave[i] *= 10.0;
        resp.push((&line[1]).parse().unwrap());
    }

    (wave, resp)
}