use statistical::mean;

pub fn latlonerr(arr: &Vec<Vec<f64>>, area: bool, val: f64) -> Vec<f64>
{
    let mut values: Vec<f64> = vec![];

    for i in 0..arr.len()
    {
        if area { values.push(mean(&arr[i]));  }
        else { values.push(arr[i][0]); }
    }

    let (min, max): (f64, f64) = ((val - values.iter().fold(f64::INFINITY, |a, &b| a.min(b))).abs(),(val - values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))).abs());
    if min == 0.0 { vec![max,max] }
    else if max == 0.0 { vec![min,min] }
    else { vec![min,max] }

}