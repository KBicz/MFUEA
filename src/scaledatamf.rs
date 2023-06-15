use rand::Rng;
use statistical::mean;
use crate::functions::PI;

pub fn scaledata(data: &Vec<f64>, arrmin: f64, arrmax: f64, tmin: f64, tmax: f64) -> Vec<f64>
{
    let mut result: Vec<f64> = vec![0f64; data.len()];
    for i in 0..data.len() { result[i] = (data[i]-arrmin)/(arrmax-arrmin)*(tmax-tmin)+tmin; }

    result
}

pub fn sf64(data: f64, arrmin: f64, arrmax: f64, tmin: f64, tmax: f64) -> f64
{
    (data-arrmin)/(arrmax-arrmin)*(tmax-tmin)+tmin
}

pub fn re_param(args: &Vec<f64>, vmin: f64, vmax: f64) -> Vec<f64>
{
    match args.len()
    {
        4 => {
            let mut an = 10f64.powf(sf64(args[0],vmin,vmax,1.0,12.0)-6f64);
            let mut bn = 10f64.powf(sf64(args[1],vmin,vmax,1.0,6.0)-3f64); 
            let mut cn = 10f64.powf(sf64(args[2],vmin,vmax,1.0,6.0)-3f64); 
            let mut dn = 10f64.powf(sf64(args[3],vmin,vmax,1.0,8.0)-5f64); 
            if an < 0f64 {an = 0f64;}
            if bn < 0f64 {bn = 0f64;}
            if cn < 0f64 {cn = 0f64;}
            if dn < 0f64 {dn = 0f64;}
            return vec![an,bn,cn,dn];
        },
        7 => {
            let mut a1n = 10f64.powf(sf64(args[0],vmin,vmax,1.0,12.0)-6f64);
            let mut bn = 10f64.powf(sf64(args[1],vmin,vmax,1.0,6.0)-3f64); 
            let mut c1n = 10f64.powf(sf64(args[2],vmin,vmax,1.0,6.0)-3f64); 
            let mut d1n = 10f64.powf(sf64(args[3],vmin,vmax,1.0,8.0)-5f64); 
            let mut a2n = 10f64.powf(sf64(args[4],vmin,vmax,1.0,12.0)-6f64);
            let mut c2n = 10f64.powf(sf64(args[5],vmin,vmax,1.0,6.0)-6f64);
            let mut d2n = 10f64.powf(sf64(args[6],vmin,vmax,1.0,8.0)-6f64);
            if a1n < 0f64 {a1n = 0f64;}
            if bn < 0f64 {bn = 0f64;}
            if c1n < 0f64 {c1n = 0f64;}
            if d1n < 0f64 {d1n = 0f64;}
            if a2n < 0f64 {a2n = 0f64;}
            if c2n < 0f64 {c2n = 0f64;}
            if d2n < 0f64 {d2n = 0f64;}
            return vec![a1n,bn,c1n,d1n,a2n,c2n,d2n];
        },
        8 => {
            let mut a1n = 10f64.powf(sf64(args[0],vmin,vmax,1.0,12.0)-6f64);
            let mut b1n = 10f64.powf(sf64(args[1],vmin,vmax,1.0,6.0)-3f64); 
            let mut c1n = 10f64.powf(sf64(args[2],vmin,vmax,1.0,6.0)-3f64); 
            let mut d1n = 10f64.powf(sf64(args[3],vmin,vmax,1.0,8.0)-5f64); 
            let mut a2n = 10f64.powf(sf64(args[4],vmin,vmax,1.0,12.0)-6f64);
            let mut b2n = 10f64.powf(sf64(args[5],vmin,vmax,1.0,6.0)-3f64); 
            let mut c2n = 10f64.powf(sf64(args[6],vmin,vmax,1.0,6.0)-3f64); 
            let mut d2n = 10f64.powf(sf64(args[7],vmin,vmax,1.0,8.0)-5f64); 
            if a1n < 0f64 {a1n = 0f64;}
            if b1n < 0f64 {b1n = 0f64;}
            if c1n < 0f64 {c1n = 0f64;}
            if d1n < 0f64 {d1n = 0f64;}
            if a2n < 0f64 {a2n = 0f64;}
            if b2n < 0f64 {b2n = 0f64;}
            if c2n < 0f64 {c2n = 0f64;}
            if d2n < 0f64 {d2n = 0f64;}
            return vec![a1n,b1n,c1n,d1n,a2n,b2n,c2n,d2n];
        },
        _ => std::process::exit(0),
    }
}

pub fn checklat(arr: &Vec<f64>, vmin: f64, vmax: f64, angrad: f64, lamin: f64, lamax: f64) -> Vec<f64>
{
    let mn: f64 = mean(arr);
    let mut result: Vec<f64> = arr.clone();
    let lat: f64 = rand::thread_rng().gen_range(0f64..1f64)*(lamax-lamin)+lamin; 

    if mn < vmin || mn > vmax
    {
        for i in 0..arr.len()
        {
            result[i] = sf64(lat+angrad*rand::thread_rng().gen_range(0f64..1f64)*(2.0*PI*rand::thread_rng().gen_range(0f64..1f64)).cos(), lamin, lamax, vmin, vmax);
        }
    }

    result 
}

pub fn checklon(arr: &Vec<f64>, vmin: f64, vmax: f64, angrad: f64, lomin: f64, lomax: f64) -> Vec<f64>
{
    let mn: f64 = mean(arr);
    let mut result: Vec<f64> = arr.clone();
    let  lon: f64 = rand::thread_rng().gen_range(0f64..1f64)*(lomax-lomin)+lomin; 

    if mn < vmin || mn > vmax
    {
        for i in 0..arr.len()
        {
            result[i] = sf64(lon+angrad*rand::thread_rng().gen_range(0f64..1f64)*(2.0*PI*rand::thread_rng().gen_range(0f64..1f64)).sin(), lomin, lomax, vmin, vmax);
        }
    }

    result 
}