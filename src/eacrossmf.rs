use rand::Rng;
use crate::scaledatamf;

pub fn ea_cross(cr: f64, s: &(Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>), v: &(Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>), vmin: f64, vmax: f64, nmainparams: usize, npts: usize, area: bool, lamin: f64, lamax: f64, lomin: f64, lomax: f64, angrad: f64) -> (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>)
{
    let npop: usize = s.1[0].len();
    let mut u: (Vec<Vec<Vec<f64>>>,Vec<Vec<f64>>) = s.clone();
    let mut indexes: Vec<Vec<f64>> = vec![vec![0.0f64; npop]; nmainparams];

    for j in 0..nmainparams
    {
        for i in 0..npop
        {
            let _noth = rand::thread_rng().try_fill(&mut indexes[j][..]);
            if indexes[j][i] >= cr
            {
                if j < 2
                {
                    if area
                    {
                        for k in 0..npts
                        {
                            u.0[j][i][k] = v.0[j][i][k];
                        }
                        u.0[j][i] = {
                            if j == 0 { scaledatamf::checklat(&u.0[j][i], vmin, vmax, angrad, lamin, lamax) }
                            else { scaledatamf::checklon(&u.0[j][i], vmin, vmax, angrad, lomin, lomax) }
                        };
                    }
                    else
                    {
                        u.0[j][i][0] = v.0[j][i][0];
                        if u.0[j][i][0] < vmin || u.0[j][i][0] > vmax {u.0[j][i][0] = rand::thread_rng().gen_range(0f64..1f64)*(vmax-vmin)+vmin;}
                    }
                }
                else
                {
                    u.1[j-2][i] = v.1[j-2][i];
                    if u.1[j-2][i] < vmin || u.1[j-2][i] > vmax  {u.1[j-2][i] = rand::thread_rng().gen_range(0f64..1f64)*(vmax-vmin)+vmin;}
                }
            }
        }
    }

    u
}