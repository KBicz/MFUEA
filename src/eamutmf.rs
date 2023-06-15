use rand::Rng;
use crate::scaledatamf;
use statistical::mean;

pub fn ea_mut(f: f64, s: &(Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>), vmin: f64, vmax: f64, nmainparams: usize, npts: usize, area: bool, lamin: f64, lamax: f64, lomin: f64, lomax: f64, angrad: f64) -> (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>)
{
    let ni: usize = s.1[0].clone().len();
    let (mut mean1, mut mean2): (f64, f64);
    let mut v: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>) = s.clone();
    let mut indexes: Vec<Vec<usize>> = vec![vec![0usize; 3]; ni];

    for i in 0..ni { for j in 0..3 {indexes[i][j] = rand::thread_rng().gen_range(0..ni); } }

    for i in 0..ni
    {
        for j in 0..nmainparams
        {
            if j < 2
            {
                if area
                {
                    mean1 = mean(&s.0[j][indexes[i][1]]);
                    mean2 = mean(&s.0[j][indexes[i][2]]);
                    for k in 0..npts { v.0[j][i][k] = s.0[j][i][k] + f * (mean1-mean2); }
                    v.0[j][i] = {
                        if j == 0 { scaledatamf::checklat(&v.0[j][i], vmin, vmax, angrad, lamin, lamax) }
                        else { scaledatamf::checklon(&v.0[j][i], vmin, vmax, angrad, lomin, lomax) }
                    };
                }
                else
                {
                    v.0[j][i][0] = s.0[j][i][0] + f * (s.0[j][indexes[i][1]][0]-s.0[j][indexes[i][2]][0]);
                    if v.0[j][i][0] < vmin || v.0[j][i][0] > vmax {v.0[j][i][0] = rand::thread_rng().gen_range(0f64..1f64)*(vmax-vmin)+vmin;} 
                }
            }
            else
            {
                v.1[j-2][i] = s.1[j-2][i] + f * (s.1[j-2][indexes[i][1]]-s.1[j-2][indexes[i][2]]);
                if v.1[j-2][i] < vmin || v.1[j-2][i] > vmax {v.1[j-2][i] = rand::thread_rng().gen_range(0f64..1f64)*(vmax-vmin)+vmin;}
            }
        }
    }

    v
}