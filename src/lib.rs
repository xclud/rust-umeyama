use std::ops::Mul;

use nalgebra::{ArrayStorage, Matrix1x2, Matrix2, Matrix2x1, Matrix3x2};

/// The Kabsch-Umeyama algorithm is a method for finding the optimal translation, rotation,
/// and scaling that aligns two sets of points with minimum root-mean-square deviation (RMSD).
pub fn umeyama<const R: usize>(src: &[(f32, f32); R], dst: &[(f32, f32); R]) -> Matrix3x2<f32> {
    let src_x_sum: f32 = src.iter().map(|v| v.0).sum();
    let src_x_mean = src_x_sum / (R as f32);

    let src_y_sum: f32 = src.iter().map(|v| v.1).sum();
    let src_y_mean = src_y_sum / (R as f32);

    let dst_x_sum: f32 = dst.iter().map(|v| v.0).sum();
    let dst_x_mean = dst_x_sum / (R as f32);

    let dst_y_sum: f32 = dst.iter().map(|v| v.1).sum();
    let dst_y_mean = dst_y_sum / (R as f32);

    let src_demean_s = ArrayStorage(src.map(|v| [v.0 - src_x_mean, v.1 - src_y_mean]));
    let dst_demean_s = ArrayStorage(dst.map(|v| [v.0 - dst_x_mean, v.1 - dst_y_mean]));

    let src_demean = nalgebra::Matrix::from_array_storage(src_demean_s);
    let dst_demean = nalgebra::Matrix::from_array_storage(dst_demean_s);

    let a = std::ops::Mul::mul(dst_demean, &src_demean.transpose()) / (R as f32);
    let svd = nalgebra::Matrix::svd(a, true, true);

    let determinant = a.determinant();

    let mut d = [1f32; 2];

    if determinant < 0.0f32 {
        d[2 - 1] = -1.0f32;
    }

    let mut t = Matrix2::<f32>::identity();
    let s = svd.singular_values;
    let u = svd.u.unwrap();
    let v = svd.v_t.unwrap();

    let rank = a.rank(0.00001f32);

    if rank == 0 {
        panic!("Matrix rank is 0.");
    } else if rank == 2 - 1 {
        //
        if u.determinant() * v.determinant() > 0.0 {
            u.mul_to(&v, &mut t);
        } else {
            let s = d[2 - 1];
            d[2 - 1] = -1f32;
            let dg = Matrix2::<f32>::new(d[0], 0f32, 0f32, d[1]);

            let udg = u.mul(&dg);
            udg.mul_to(&v, &mut t);
            d[2 - 1] = s;
        }
    } else {
        let dg = Matrix2::<f32>::new(d[0], 0f32, 0f32, d[1]);
        let udg = u.mul(&dg);
        udg.mul_to(&v, &mut t);
    }

    let ddd = Matrix1x2::new(d[0], d[1]);
    let d_x_s = ddd.mul(s);

    let var0 = src_demean.remove_row(0).variance();
    let var1 = src_demean.remove_row(1).variance();

    let varsum = var0 + var1;

    let scale = d_x_s.get((0, 0)).unwrap() / varsum;

    let dst_mean = Matrix2x1::<f32>::new(dst_x_mean, dst_y_mean);
    let src_mean = Matrix2x1::<f32>::new(src_x_mean, src_y_mean);
    let t_x_srcmean = t.mul(&src_mean);

    let xxx = scale * t_x_srcmean;
    let yyy = dst_mean - xxx;

    let m31 = yyy.get(0).unwrap();
    let m32 = yyy.get(1).unwrap();

    let m00x22 = t * scale;

    let m11 = m00x22.m11;
    let m21 = m00x22.m21;
    let m12 = m00x22.m12;
    let m22 = m00x22.m22;

    let tt = Matrix3x2::<f32>::new(m11, m12, m21, m22, *m31, *m32);

    return tt;
}

#[cfg(test)]
mod tests {
    use super::*;

    pub const POINTS: [(f32, f32); 5] = [
        (38.2946f32, 51.6963f32),
        (73.5318f32, 51.5014f32),
        (56.0252f32, 71.7366f32),
        (41.5493f32, 92.3655f32),
        (70.7299f32, 92.2041f32),
    ];

    #[test]
    fn estimate() {
        let result = umeyama(&POINTS, &POINTS);
        assert_eq!(result, Matrix3x2::identity());
    }
}
