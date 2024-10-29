use ark_bls12_381::{Bls12_381, G1Projective};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::PrimeField;
use ark_std::cfg_iter;

use crate::{arith::{r1cs::R1CS, vec::{hadamard, matrix_vec_mul, vec_add, vec_scalar_mul, vec_sub}}, utils::RandomOracle};

type BlsScalarField = <Bls12_381 as Pairing>::ScalarField;

pub struct Witness {
    pub E: Vec<BlsScalarField>,
    pub rE: BlsScalarField,
    pub W: Vec<BlsScalarField>,
    pub rW: BlsScalarField,
}

pub struct Instance {
    E: Vec<BlsScalarField>,
    u: BlsScalarField,
    x: Vec<BlsScalarField>,
}

impl Absorb for Instance {
    fn to_sponge_bytes(&self, dest: &mut Vec<u8>) {
        BlsScalarField::batch_to_sponge_bytes(&self.to_sponge_field_elements_as_vec(), dest);
    }

    fn to_sponge_field_elements<F: PrimeField>(&self, dest: &mut Vec<F>) {
        self.u.to_sponge_field_elements(dest);
        self.x.to_sponge_field_elements(dest);
        // We cannot call `to_native_sponge_field_elements(dest)` directly, as
        // `to_native_sponge_field_elements` needs `F` to be `C::ScalarField`,
        // but here `F` is a generic `PrimeField`.
        // self.cmE
        //     .to_native_sponge_field_elements_as_vec()
        //     .to_sponge_field_elements(dest);
        // self.cmW
        //     .to_native_sponge_field_elements_as_vec()
        //     .to_sponge_field_elements(dest);
    }
}

pub struct Nova {

}

impl Nova {
    fn fold_witness(
        r: BlsScalarField,
        W_i: &Witness,
        w_i: &Witness,
        T: Vec<BlsScalarField>,
    ) {
        let r2 = r * r;
        let E = vec_add(
            &vec_add(&W_i.E, &vec_scalar_mul(&T, &r)), 
            &vec_scalar_mul(&w_i.E, &r2),
        );
        let W: Vec<BlsScalarField> = cfg_iter!(W_i.W).zip(&w_i.W)
            .map(| (a, b) | *a + (r * b))
            .collect();
    }

    fn prove(
        r1cs: &R1CS<BlsScalarField>,
        W_i: &Witness,
        U_i: &Instance,
        w_i: &Witness,
        u_i: &Instance,
    ) {
        let z1: Vec<BlsScalarField> = [W_i.W.to_vec(), U_i.x.to_vec(), vec![U_i.u]].concat();
        let z2: Vec<BlsScalarField> = [w_i.W.to_vec(), u_i.x.to_vec(), vec![u_i.u]].concat();
        let T = Self::compute_T(r1cs, U_i.u, u_i.u, &z1, &z2);

        let r = RandomOracle::<G1Projective, Instance>::get_challenge(U_i, u_i);
        let w = Self::fold_witness(r, W_i, w_i, T);
    }

    fn compute_T(
        r1cs: &R1CS<BlsScalarField>,
        u1: BlsScalarField,
        u2: BlsScalarField,
        z1: &[BlsScalarField],
        z2: &[BlsScalarField],
    ) -> Vec<BlsScalarField> { 
        let (a, b, c) = (r1cs.a.clone(), r1cs.b.clone(), r1cs.c.clone());

        let az1 = matrix_vec_mul(&a, z1);
        let bz1 = matrix_vec_mul(&b, z1);
        let cz1 = matrix_vec_mul(&c, z1);
        let az2 = matrix_vec_mul(&a, z2);
        let bz2 = matrix_vec_mul(&b, z2);
        let cz2 = matrix_vec_mul(&c, z2);

        let az1_bz2 = hadamard(&az1, &bz2);
        let az2_bz1 = hadamard(&az2, &bz1);
        let u1_cz2 = vec_scalar_mul(&cz2, &u1);
        let u2_cz1 = vec_scalar_mul(&cz1, &u2);

        vec_sub(&vec_sub(&vec_add(&az1_bz2, &az2_bz1), &u1_cz2), &u2_cz1)
    }
}
