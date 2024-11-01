use ark_relations::r1cs::Matrix as R1CSMatrix;
use ark_std::One;

use crate::{arith::vec::{hadamard, matrix_vec_mul, vec_scalar_mul, vec_sub}, nova::{Instance, Witness}, utils::BlsScalarField};

pub struct R1CS {
    pub l: usize, // io len
    pub a: R1CSMatrix<BlsScalarField>,
    pub b: R1CSMatrix<BlsScalarField>,
    pub c: R1CSMatrix<BlsScalarField>,
}

impl R1CS {
    pub fn eval_at_z(&self, z: &[BlsScalarField]) -> Vec<BlsScalarField> {
        let az = matrix_vec_mul(&self.a, z);
        let bz = matrix_vec_mul(&self.b, z);
        let cz = matrix_vec_mul(&self.c, z);
        // Multiply Cz by z[-1] (u) here, allowing this method to be reused for
        // both relaxed and plain R1CS.
        let ucz = vec_scalar_mul(&cz, &z[z.len() - 1]);
        let az_bz = hadamard(&az, &bz);
        vec_sub(&az_bz, &ucz)
    }

    pub fn eval_relation(&self, w: &Vec<BlsScalarField>, x: &Vec<BlsScalarField>) -> Vec<BlsScalarField> {
        self.eval_at_z(&[w.to_vec(), x.to_vec(), [BlsScalarField::one()].to_vec()].concat())
    }

    pub fn check_relation(&self, w: &Witness, u: &Instance) {
        let z = [w.w.to_vec(), u.x.to_vec(), [u.u].to_vec()].concat();
        let evals = self.eval_at_z(&z);
        assert_eq!(evals, w.e);
    }
}

#[cfg(test)]
pub mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ff::PrimeField;
    use ark_relations::r1cs::Matrix as R1CSMatrix;
    use ark_std::{cfg_iter, rand::Rng, Zero};

    use crate::arith::vec::is_zero_vec;

    use super::R1CS;

    type BlsScalarField = <Bls12_381 as Pairing>::ScalarField;

    pub fn to_f_vec<F: PrimeField>(z: Vec<usize>) -> Vec<F> {
        cfg_iter!(z)
            .map(| e | F::from(*e as u64))
            .collect()
    }

    pub fn get_test_z<F: PrimeField>(input: usize) -> Vec<F> {
        // z = (W, x, u)
        to_f_vec(vec![
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
            input,                             // io
            1,
        ])
    }

    pub fn get_test_z_split<F: PrimeField>(input: usize) -> (Vec<F>, Vec<F>, F) {
        // z = (W, x, u)
        (
            to_f_vec(vec![
                input * input * input + input + 5, // x^3 + x + 5
                input * input,                     // x^2
                input * input * input,             // x^2 * x
                input * input * input + input,     // x^3 + x
            ]),
            to_f_vec(vec![
                input, // io
            ]),
            F::one(),
        )
    }

    pub fn to_f_matrix<F: PrimeField>(m: Vec<Vec<usize>>) -> R1CSMatrix<F> {
        let mut matrix = Vec::new();
        for m_row in m {
            let mut row = Vec::<(F, usize)>::new();
            for (i, value) in m_row.iter().enumerate() {
                if *value != 0 {
                    row.push((F::from(*value as u64), i));
                }
            }
            matrix.push(row);
        }
        matrix
    }

    pub fn get_test_r1cs() -> R1CS {
        let a = to_f_matrix(vec![
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 0, 1, 0, 1, 0],
            vec![0, 0, 0, 1, 0, 5],
        ]);
        let b = to_f_matrix(vec![
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 1],
            vec![0, 0, 0, 0, 0, 1],
        ]);
        let c = to_f_matrix(vec![
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 0, 1, 0, 0, 0],
            vec![0, 0, 0, 1, 0, 0],
            vec![1, 0, 0, 0, 0, 0],
        ]);
        R1CS { l: 1, a, b, c }
    }

    #[test]
    fn test_eval_r1cs_relation() {
        let mut rng = ark_std::test_rng();
        let r1cs = get_test_r1cs();
        let z = get_test_z::<BlsScalarField>(rng.gen::<u16>() as usize);

        let evals = r1cs.eval_at_z(&z);
        is_zero_vec(&evals);
    }
}
