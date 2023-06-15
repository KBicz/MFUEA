mod iomf;
mod errors;
mod fplotd;
mod eacrmf;
mod eamutmf;
mod energymf;
mod eacrossmf;
mod functions;
mod scaledatamf;

use std::env;
use rand::Rng;
use ansi_term::Style;
use statistical::mean;
use std::time::Instant;
use std::process::exit;
use indicatif::{ProgressStyle,ProgressBar};

fn helpf() 
{
    println!("\n Program modulateflare_ea written by K. Bicz for Linux, Mac, and Windows, version of 13 Mar 2023.");
    println!(" Usage: modulateflare <-lc=file> [-cross=f64] [-mut=f64] [-npop=int] [-niter=f64] [-prot=f64]");
    println!("                      [-imin=f64] [-imax=f64] [-lamin=f64] [-lamax=f64] [-lomin=f64] [-lomax=f64]");
    println!("                      [-nparams=4/7/8] [-npts=i64] [-teff=f64] [-rad=f64] [-terr=f64] [-raderr=f64]");
    println!("                      [-logg=f64] [-respf=str] [-frad=f64] [-famp] [--limbd] [--noprof] [--energy]");
    println!("                      [--plot] [--save]\n");
    println!("         option  -lc      : file with the light curve of the event.");
    println!("                 -cross   : crossing value (default cross = 0.1).");
    println!("                 -mut     : mutation value (default mut = 0.0).");
    println!("                 -npop    : number of population (default npop = (nparams+3)*10).");
    println!("                 -niter   : number of iterations with mutation and crossing (default niter = 10000).");
    println!("                 -imin    : minimal inclination angle of the star (default imin = 10).");
    println!("                 -imax    : maximal inclination angle of the star (default imax = 90).");
    println!("                 -lamin   : minimal latitude of the flaring region (default lamin = -90).");
    println!("                 -lamax   : maximal latitude of the flaring region (default lamax = 90).");
    println!("                 -lamin   : minimal longitude of the flaring region (default lamin = 0, east limb).");
    println!("                 -lamax   : maximal longitude of the flaring region (default lamax = 360).");
    println!("                 -nparams : number of parameters to fit profile (default nparams = 4).");
    println!("                 -npts    : make region not the point but area made of npts points (default npts = 100)");
    println!("                 -prot    : rotational period of the star (default prot = 0.2353 days).");
    println!("                 -teff    : effective temperature of the star (default teff = 4885.0 K).");
    println!("                 -terr    : error of the effective temperature of the star (default terr = 128.079 K");
    println!("                 -rad     : radius of the star (default rad = 0.799503 Rsun");
    println!("                 -raderr  : error of the radius of the star (default rad = 0.0540905 Rsun");
    println!("                 -logg    : log(g) of the star (def logg = 4.53552).");
    println!("                 -frad    : radius of the flaring region in degrees (by default it is estimated).");
    println!("                 -famp    : force the new amplitude of the flare (by default it is taken from data).");
    println!("                 -respf   : file with the response function of the observational instrument");
    println!("                 --limbd  : use limb darkening in the fit.");
    println!("                 --noprof : assume that fitted profile from lc file is the final profile of the event.");
    println!("                 --energy : estimate the energy of the flare before and after the correction.");
    println!("                 --plot   : display the results.");
    println!("                 --save   : save the results to file.\n");
    exit(0);
}

fn ea_main(lc: &str, mutat: f64, cross: f64, mut n_pop: usize, niter: usize, plotctrl: bool, savectrl: bool, imin: f64, imax: f64, lamin: f64, lamax: f64, lomin: f64, lomax: f64, prot: f64, teff: f64, tefferr: f64, rad: f64, raderr: f64, energy: bool, respfile: &str, respctrl: bool, mut nparams: usize, fitprofile: bool, npopctrl: bool, npts: usize, mut area: bool, limbdctrl: bool, logg: f64, frad: f64, famp: f64) -> Result<(), &'static str> 
{
    let start: Instant = Instant::now();
    let (time, trend, err): (Vec<f64>,Vec<f64>,Vec<f64>);
    let (mut profile, timeprof, t0, istart, fitprof, newtime): (Vec<f64>, Vec<f64>, f64, usize, bool, Vec<f64>);
    (time, trend, err, profile, timeprof, t0, istart, fitprof, newtime, nparams) = iomf::read_lcprof(lc,nparams,fitprofile);

    let (u1, u2): (f64, f64) = {
        if limbdctrl { functions::limbd_coeffs(logg, teff) }
        else {(0.0, 0.0)}
    };

    let mut t: usize = 10;
    let mut minchi: f64 = 1f64;
    let mut tbar: Vec<f64> = vec![];
    let (mut lat, mut lon): (f64, f64);
    let (mut hlat, mut hlon): (f64, f64);
    let (vmin, vmax): (f64, f64) = (1.0,12.0);
    let mut v: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>);
    let mut u: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>);
    let (angrad, mut theta, mut radscale): (f64, f64, f64);
    let influence: Vec<f64> = vec![1.0/(npts as f64); npts];
    let mut s: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>) = (vec![vec![],vec![]],vec![]);

    for i in 0..timeprof.len() { tbar.push(2.0*functions::PI*(timeprof[i]-t0)/prot) }

    let nmainparams: usize = {
        if fitprof {
            match nparams
            {
                4 => { 7 },
                7 => { 10 },
                8 => { 11 },
                _ => { exit(1) },
            }
        }
        else { 3 }
    };

    if !npopctrl { n_pop = nmainparams*10; }
    let mut indeksy: Vec<f64> = vec![0.0; n_pop];

    if area && frad == 0.0
    {
        if famp == 0.0 { (angrad, area)  = energymf::footrad((trend[istart..istart+newtime.len()]).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)), teff, respfile, respctrl); }
        else { (angrad, area)  = energymf::footrad(famp, teff, respfile, respctrl); }
    }
    else if !area && frad != 0.0 { angrad = frad; area = true; }
    else { angrad = 0.0; }

    for i in 0..n_pop
    {
        lat = rand::thread_rng().gen_range(0f64..1f64)*(lamax-lamin)+lamin;
        lon = rand::thread_rng().gen_range(0f64..1f64)*(lomax-lomin)+lomin;
        s.0[0].push(vec![0.0; npts]);
        s.0[1].push(vec![0.0; npts]);
        for j in 0..npts
        {            
            hlat = 0.0; hlon =  0.0;
            while s.0[0][i].contains(&hlat) && s.0[1][i].contains(&hlon) && area
            { 
                theta = rand::thread_rng().gen_range(0f64..1f64);
                radscale = rand::thread_rng().gen_range(0f64..1f64);
                hlat = scaledatamf::sf64(lat+angrad*radscale.powf(0.5)*(2.0*functions::PI*theta).cos(), lamin, lamax, vmin, vmax);
                hlon = scaledatamf::sf64(lon+angrad*radscale.powf(0.5)*(2.0*functions::PI*theta).sin(), lomin, lomax, vmin, vmax);
            }
            if !area { hlat = lat; hlon = lon; }

            s.0[0][i][j] = scaledatamf::sf64(hlat,lamin,lamax,vmin,vmax); 
            s.0[1][i][j] = scaledatamf::sf64(hlon,lomin,lomax,vmin,vmax); 
        }
    }

    for _ in 2..nmainparams
    {
        for j in 0..n_pop {indeksy[j] = rand::thread_rng().gen_range(0f64..1f64)*(vmax-vmin)+vmin;}
        (s.1).push(indeksy.clone()); 
    }
    drop(indeksy);

    let pb: ProgressBar = ProgressBar::new(niter as u64);
    pb.set_style(ProgressStyle::with_template(&format!("{} {} {{spinner:.green}} [{{elapsed_precise}}] ╢{{bar:45.white/gray}}╟ {{percent}}% [{{eta_precise}}, {{per_sec}}]",Style::new().bold().paint("Fitting").to_string(),Style::new().bold().paint("profile").to_string())).unwrap());

    for _ in 0..niter
    {
        v = eamutmf::ea_mut(mutat,&s.clone(),vmin,vmax,nmainparams,npts,area,lamin,lamax,lomin,lomax,angrad);
        u = eacrossmf::ea_cross(cross,&s.clone(),&v.clone(),vmin,vmax,nmainparams,npts,area,lamin,lamax,lomin,lomax,angrad);
        (s, t, minchi) = eacrmf::ea_cr(&s.clone(),&u.clone(),&tbar,&trend,&err,&profile,istart,imin,imax,lamin,lamax,lomin,lomax,vmin,vmax,fitprof,nmainparams,&newtime,npts,area,&influence,u1,u2);
        
        pb.inc(1);
    }
    pb.finish();

    for i in 0..n_pop
    {
        s.0[0][i] = scaledatamf::scaledata(&s.0[0][i],vmin,vmax,lamin, lamax);
        s.0[1][i] = scaledatamf::scaledata(&s.0[1][i],vmin,vmax,lomin, lomax);
    }
    s.1[0] = scaledatamf::scaledata(&s.1[0],vmin,vmax,imin, imax);

    let mut paramp: Vec<f64> = vec![0.0; nparams];
    if fitprof
    {
        match nparams
        {
            4 => {
                paramp = scaledatamf::re_param(&vec![s.1[1][t],s.1[2][t],s.1[3][t],s.1[4][t]],vmin,vmax);
                profile = functions::profil_wroclawski(&newtime, paramp[0], paramp[1], paramp[2], paramp[3]);
            },
            7 => {
                paramp = scaledatamf::re_param(&vec![s.1[1][t],s.1[2][t],s.1[3][t],s.1[4][t],s.1[5][t],s.1[6][t],s.1[7][t]],vmin,vmax);
                profile = functions::profil_wroclawski_1b(&newtime, paramp[0], paramp[1], paramp[2], paramp[3], paramp[4], paramp[5], paramp[6]); 
            },
            8 => {
                paramp = scaledatamf::re_param(&vec![s.1[1][t],s.1[2][t],s.1[3][t],s.1[4][t],s.1[5][t],s.1[6][t],s.1[7][t],s.1[8][t]],vmin,vmax);
                profile = functions::profil_wroclawski_2b(&newtime, paramp[0], paramp[1], paramp[2], paramp[3], paramp[4], paramp[5], paramp[6], paramp[7]); 
            },
            _ => { exit(1); },
        }
    }

    let duration = start.elapsed();
    println!("Execution time: {:?}", duration);
    println!("chi^2/(N-{}) = {}",3,minchi/( (tbar.len()-3) as f64 ) );

    let para: Vec<f64> = {
        if area { vec![mean(&s.0[0][t]),mean(&s.0[1][t]),s.1[0][t]] }
        else { vec![s.0[0][t][0],s.0[1][t][0],s.1[0][t]] }
    };
    let model: Vec<f64> = {
        if !area { functions::alpha(para[0], para[1], para[2], &tbar, &profile, fitprof, &vec![u1,u2]) }
        else { functions::alpha_area(&s.0[0][t], &s.0[1][t], para[2], &tbar, fitprof, &influence, &profile,&vec![u1,u2]) }
    };

    let laterr: Vec<f64> = errors::latlonerr(&s.0[0], area, para[0]); 
    let lonerr: Vec<f64> = errors::latlonerr(&s.0[1], area, para[1]); 

    println!("> Radius = {} deg",angrad);
    println!("> Latitude = {} + {} - {} deg",para[0],laterr[0],laterr[1]);
    println!("> Longitude = {} + {} - {} deg",para[1],lonerr[0],lonerr[1]);
    println!("> Inclination = {} + {} - {} deg",para[2],(para[2]-s.1[0].iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))).abs(),(para[2]-s.1[0].iter().fold(f64::INFINITY, |a, &b| a.min(b))).abs());
    if fitprof { functions::printparams(&s, nmainparams, &paramp, vmin, vmax); }
    if energy 
    { 
        let eold: (f64, f64) = energymf::shibayama(&timeprof, &model, teff, rad, 0.0, respfile, respctrl);
        let eolderr: f64 = energymf::shibayama(&timeprof, &model, teff+tefferr, rad+raderr, 2000.0, respfile, respctrl).0-eold.0;
        let enew: (f64, f64) = energymf::shibayama(&timeprof, &profile, teff, rad, 0.0, respfile, respctrl);
        let enewerr: f64 = energymf::shibayama(&timeprof, &profile, teff+tefferr, rad+raderr, 2000.0, respfile, respctrl).0-enew.0;
        println!("> Energy : {:e} +/- {:e} -> {:e} +/- {:e}",eold.0,eolderr,enew.0,enewerr); 
    }
    if famp == 0.0 { println!("> Amplitude {} -> {}",(trend[istart..istart+newtime.len()]).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),profile.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))); }
    else {println!("> Amplitude ({}) {} -> {}",(trend[istart..istart+newtime.len()]).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),famp,profile.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))); }

    if savectrl { fplotd::plotd(&time, &trend, &timeprof, &profile, &model, savectrl, false, 0.25); }
    if plotctrl { fplotd::plotd(&time, &trend, &timeprof, &profile, &model, false, plotctrl, 1.0); }
    
    Ok(())
}

fn main() -> Result<(), &'static str> 
{
    let mut v: Vec<&str>;
    let mut n_pop: usize = 50;
    let mut prot: f64 = 0.2354;
    let mut nparams: usize = 4;
    let mut energy: bool = false;
    let argc: usize = env::args().len();
    let mut file: &str = "rozblysk.dat";
    let argv: Vec<String> = env::args().collect();
    let (mut imin, mut imax, mut frad): (f64, f64, f64) = (10.0, 90.0, 0.0);
    let (mut mutat, mut cross, mut niter): (f64, f64, usize) = (0.9, 0.1, 10000);
    let (mut plotctrl, mut savectrl, mut fitprofile): (bool, bool, bool) = (false, false, true);
    let (mut lomin, mut lomax, mut logg, mut famp): (f64, f64, f64, f64) = (0.0, 360.0, 4.53552, 0.0);
    let (mut lamin, mut lamax, mut npts, mut area): (f64, f64, usize, bool) = (-90.0, 90.0, 100, false);
    let (mut teff, mut terr, mut rad, mut raderr): (f64, f64, f64, f64) = (4885.0,128.079,0.799503*6.9634e8,0.0540905*6.9634e8);
    let (mut respfile, mut respctrl, mut npopctrl, mut limbdctrl): (&str, bool, bool, bool) = ("TESS_response.txt", false, false, false);

    if argc < 2 {helpf();}

    for i in 0..argc
    {
        if argv[i].contains("-lc=") {v = argv[i].split("=").collect(); file = v[1];}
        else if argv[i].contains("-cross=") {v = argv[i].split("=").collect(); cross = v[1].parse().unwrap();}
        else if argv[i].contains("-mut=") {v = argv[i].split("=").collect(); mutat = v[1].parse().unwrap();}
        else if argv[i].contains("-npop=") {v = argv[i].split("=").collect(); n_pop = v[1].parse().unwrap(); npopctrl = true;}
        else if argv[i].contains("-niter=") {v = argv[i].split("=").collect(); niter = v[1].parse().unwrap();}
        else if argv[i].contains("-imin=") {v = argv[i].split("=").collect(); imin = v[1].parse().unwrap();}
        else if argv[i].contains("-imax=") {v = argv[i].split("=").collect(); imax = v[1].parse().unwrap();}
        else if argv[i].contains("-lamin=") {v = argv[i].split("=").collect(); lamin = v[1].parse().unwrap();}
        else if argv[i].contains("-lamax=") {v = argv[i].split("=").collect(); lamax = v[1].parse().unwrap();}
        else if argv[i].contains("-lomin=") {v = argv[i].split("=").collect(); lomin = v[1].parse().unwrap();}
        else if argv[i].contains("-lomax=") {v = argv[i].split("=").collect(); lomax = v[1].parse().unwrap();}
        else if argv[i].contains("-nparams=") {v = argv[i].split("=").collect(); nparams = v[1].parse().unwrap();}
        else if argv[i].contains("-npts=") {v = argv[i].split("=").collect(); npts = v[1].parse().unwrap(); area = true;}
        else if argv[i].contains("-prot=") {v = argv[i].split("=").collect(); prot = v[1].parse().unwrap();}
        else if argv[i].contains("-teff=") {v = argv[i].split("=").collect(); teff = v[1].parse().unwrap();}
        else if argv[i].contains("-terr=") {v = argv[i].split("=").collect(); terr = v[1].parse().unwrap();}
        else if argv[i].contains("-rad=") {v = argv[i].split("=").collect(); rad = v[1].parse().unwrap();}
        else if argv[i].contains("-raderr=") {v = argv[i].split("=").collect(); raderr = v[1].parse().unwrap();}
        else if argv[i].contains("-logg=") {v = argv[i].split("=").collect(); logg = v[1].parse().unwrap();}
        else if argv[i].contains("-frad=") {v = argv[i].split("=").collect(); frad = v[1].parse().unwrap();}
        else if argv[i].contains("-famp=") {v = argv[i].split("=").collect(); famp = v[1].parse().unwrap();}
        else if argv[i].contains("-respf=") {v = argv[i].split("=").collect(); respfile = v[1]; respctrl = true;}
        else if argv[i].contains("--limbd") {limbdctrl = true;}
        else if argv[i].contains("--noprof") {fitprofile = false;}
        else if argv[i].eq("--energy") { energy = true; }
        else if argv[i].eq("--plot") {plotctrl = true;}
        else if argv[i].eq("--save") {savectrl = true;}
        else if argv[i].eq("--help") || argv[i].eq("-h")  {helpf();}
    }

    ea_main(file,mutat,cross,n_pop,niter,plotctrl,savectrl,imin,imax,lamin,lamax,lomin,lomax,prot,teff,terr,rad,raderr,energy,respfile,respctrl,nparams,fitprofile,npopctrl,npts,area,limbdctrl,logg,frad,famp)
}