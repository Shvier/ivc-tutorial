use ark_ff::PrimeField;
use ark_relations::r1cs::Matrix as R1CSMatrix;
use ark_std::cfg_iter;

pub fn matrix_vec_mul<F: PrimeField>(
    m: &R1CSMatrix<F>,
    z: &[F],
) -> Vec<F> {
    // assert_eq!(m[0].len(), z.len());
    cfg_iter!(m)
        .map(| row | row.iter()
            .map(| (value, col_i) | *value * z[*col_i]).sum())
        .collect()
}

pub fn hadamard<F: PrimeField>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    // assert_eq!(a.len(), b.len());
    cfg_iter!(a).zip(b)
        .map(| (a, b) | *a * b)
        .collect()
}

pub fn vec_scalar_mul<F: PrimeField>(
    vec: &[F],
    e: &F,
) -> Vec<F> {
    cfg_iter!(vec)
        .map(| a | *a * e)
        .collect()
}

pub fn vec_add<F: PrimeField>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    cfg_iter!(a).zip(b)
        .map(| (a, b) | *a + b)
        .collect()
}

pub fn vec_sub<F: PrimeField>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    cfg_iter!(a).zip(b)
        .map(| (a, b) | *a - b)
        .collect()
}

pub fn is_zero_vec<F: PrimeField>(vec: &[F]) -> bool {
    cfg_iter!(vec).all(|a| a.is_zero())
}
