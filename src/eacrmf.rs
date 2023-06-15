use crate::functions;
use crate::scaledatamf::{scaledata,sf64,re_param};

pub fn ea_cr(s: &(Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>), u: &(Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>), time: &Vec<f64>, flux: &Vec<f64>, err: &Vec<f64>, profile: &Vec<f64>, istart: usize, imin: f64, imax: f64, lamin: f64, lamax: f64, lomin: f64, lomax: f64, vmin: f64, vmax: f64, fitprof: bool, nparams: usize, newtime: &Vec<f64>, npts: usize, area: bool, influence: &Vec<f64>,u1: f64, u2: f64) -> ((Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>),usize,f64)
{
    let mut minchi: f64;
    let mut min_t: usize = 100000;
    let mut min_chi: f64 = 100000f64;
    let (mut chi_s, mut chi_u): (f64, f64);
    let (mut model_s, mut model_u): (Vec<f64>,Vec<f64>);
    let mut new_s: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>) = s.clone();
    let (mut profile_s, mut profile_u): (Vec<f64>, Vec<f64>) = (profile.to_vec().clone(), profile.to_vec().clone());
    let (mut pars, mut paru): (Vec<f64>,Vec<f64>);

    for i in 0..s.1[0].len()
    {
        if fitprof
        {
            match nparams
            {
                7 => {
                    pars = re_param(&vec![s.1[1][i],s.1[2][i],s.1[3][i],s.1[4][i]],vmin,vmax);
                    paru = re_param(&vec![u.1[1][i],u.1[2][i],u.1[3][i],u.1[4][i]],vmin,vmax);
                    profile_s = functions::profil_wroclawski(&newtime, pars[0], pars[1], pars[2], pars[3]);
                    profile_u = functions::profil_wroclawski(&newtime, paru[0], paru[1], paru[2], paru[3]);
                },
                10 => {
                    pars = re_param(&vec![s.1[1][i],s.1[2][i],s.1[3][i],s.1[4][i],s.1[5][i],s.1[6][i],s.1[7][i]],vmin,vmax);
                    paru = re_param(&vec![u.1[1][i],u.1[2][i],u.1[3][i],u.1[4][i],u.1[5][i],u.1[6][i],u.1[7][i]],vmin,vmax);
                    profile_s = functions::profil_wroclawski_1b(&newtime, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]);
                    profile_u = functions::profil_wroclawski_1b(&newtime, paru[0], paru[1], paru[2], paru[3], paru[4], paru[5], paru[6]);
                },
                11 => {
                    pars = re_param(&vec![s.1[1][i],s.1[2][i],s.1[3][i],s.1[4][i],s.1[5][i],s.1[6][i],s.1[7][i],s.1[8][i]],vmin,vmax);
                    paru = re_param(&vec![u.1[1][i],u.1[2][i],u.1[3][i],u.1[4][i],u.1[5][i],u.1[6][i],u.1[7][i],u.1[8][i]],vmin,vmax);
                    profile_s = functions::profil_wroclawski_2b(&newtime, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7]);
                    profile_u = functions::profil_wroclawski_2b(&newtime, paru[0], paru[1], paru[2], paru[3], paru[4], paru[5], paru[6], paru[7]);
                }, 
                _ => { std::process::exit(1); },
            }
        }

        if area
        {
            model_s = functions::alpha_area(&scaledata(&s.0[0][i],vmin,vmax, lamin,lamax,),&scaledata(&s.0[1][i],vmin,vmax,lomin,lomax),sf64(s.1[0][i],vmin,vmax,imin,imax),&time,fitprof,&influence,&profile_s,&vec![u1,u2]);
            model_u = functions::alpha_area(&scaledata(&u.0[0][i],vmin,vmax, lamin,lamax,),&scaledata(&u.0[1][i],vmin,vmax,lomin,lomax),sf64(u.1[0][i],vmin,vmax,imin,imax),&time,fitprof,&influence,&profile_u,&vec![u1,u2]);
        }
        else
        {
            model_s = functions::alpha(sf64(s.0[0][i][0],vmin,vmax,lamin,lamax),sf64(s.0[1][i][0],vmin,vmax,lomin,lomax),sf64(s.1[0][i],vmin,vmax,imin,imax),&time,&profile_s,fitprof,&vec![u1,u2]);
            model_u = functions::alpha(sf64(u.0[0][i][0],vmin,vmax,lamin,lamax),sf64(u.0[1][i][0],vmin,vmax,lomin,lomax),sf64(u.1[0][i],vmin,vmax,imin,imax),&time,&profile_u,fitprof,&vec![u1,u2]);    
        }
        
        chi_s = 0.0; chi_u = 0.0;
        for k in 0..time.len()
        {
            chi_s += ((model_s[k]-flux[istart+k])/err[istart+k]).powf(2.0);
            chi_u += ((model_u[k]-flux[istart+k])/err[istart+k]).powf(2.0);
        }

        minchi = chi_s;
        if chi_s > chi_u 
        {
            minchi = chi_u; 
            for j in 0..nparams
            {
                if j < 2
                {
                    if area { for k in 0..npts { new_s.0[j][i][k] = u.0[j][i][k]; } }
                    else { new_s.0[j][i][0] = u.0[j][i][0] }
                }
                else
                {
                    new_s.1[j-2][i] = u.1[j-2][i];
                }
            }
        }
        
        
        if i == 0 || min_chi > minchi
        {
            min_chi = minchi;  
            min_t = i;
        }
    }
    
    (new_s, min_t, min_chi)
}