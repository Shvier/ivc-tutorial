use ark_bls12_381::G1Projective;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_std::{cfg_iter, rand::RngCore, UniformRand, Zero, One};

use crate::{arith::{r1cs::R1CS, vec::{hadamard, matrix_vec_mul, vec_add, vec_scalar_mul, vec_sub}}, utils::{BlsScalarField, RandomOracle}};

pub struct Witness {
    pub e: Vec<BlsScalarField>,
    pub rE: BlsScalarField,
    pub w: Vec<BlsScalarField>,
    pub rW: BlsScalarField,
}

impl Witness {
    pub fn new<const H: bool>(w: Vec<BlsScalarField>, e_len: usize, mut rng: impl RngCore) -> Self {
        let (rW, rE) = if H {
            (
                BlsScalarField::rand(&mut rng),
                BlsScalarField::rand(&mut rng),
            )
        } else {
            (BlsScalarField::zero(), BlsScalarField::zero())
        };

        Self {
            e: vec![BlsScalarField::zero(); e_len],
            rE,
            w,
            rW,
        }
    }
}

pub struct Instance {
    pub e: Vec<BlsScalarField>,
    pub u: BlsScalarField,
    pub w: Vec<BlsScalarField>,
    pub x: Vec<BlsScalarField>,
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

pub struct Nova<
    const H: bool = false,
> {

}

impl<const H: bool> Nova<H> {
    fn new_witness(w: &Vec<BlsScalarField>, e_len: usize, rng: impl RngCore) -> Witness {
        Witness::new::<H>(w.to_vec(), e_len, rng)
    }

    fn new_instance(w: &Witness, x: &Vec<BlsScalarField>) -> Instance {
        Instance { e: w.e.clone(), w: w.w.clone(), u: BlsScalarField::one(), x: x.to_vec() }
    }

    fn fold_witness(
        r: BlsScalarField,
        w_i: &Witness,
        incoming_w_i: &Witness,
        t: &Vec<BlsScalarField>,
    ) -> Witness {
        let r2 = r * r;
        let e = vec_add(
            &vec_add(&w_i.e, &vec_scalar_mul(t, &r)), 
            &vec_scalar_mul(&incoming_w_i.e, &r2),
        );
        let w: Vec<BlsScalarField> = cfg_iter!(w_i.w).zip(&incoming_w_i.w)
            .map(| (a, b) | *a + (r * b))
            .collect();
        Witness { e, rE: BlsScalarField::zero(), w, rW: BlsScalarField::zero() }
    }

    fn fold_instance(
        r: BlsScalarField,
        u_i: &Instance,
        incoming_u_i: &Instance,
        t: &Vec<BlsScalarField>,
    ) -> Instance {
        let r2 = r * r;
        let e = vec_add(
            &u_i.e, 
            &vec_add(
                &vec_scalar_mul(&t, &r), 
                &vec_scalar_mul(&incoming_u_i.e, &r2)
            )
        );
        let u = u_i.u + r * incoming_u_i.u;
        let w = vec_add(&u_i.w, &vec_scalar_mul(&incoming_u_i.w, &r));
        let x = cfg_iter!(u_i.x).zip(&incoming_u_i.x)
            .map(| (a, b) | *a + (r * b))
            .collect();
        Instance { e, u, w, x }
    }

    fn prove(
        r1cs: &R1CS,
        w_i: &Witness,
        u_i: &Instance,
        incoming_w_i: &Witness,
        incoming_u_i: &Instance,
    ) -> (Witness, Instance, Vec<BlsScalarField>) {
        let z1: Vec<BlsScalarField> = [w_i.w.to_vec(), u_i.x.to_vec(), vec![u_i.u]].concat();
        let z2: Vec<BlsScalarField> = [incoming_w_i.w.to_vec(), incoming_u_i.x.to_vec(), vec![incoming_u_i.u]].concat();
        let t = Self::compute_T(r1cs, u_i.u, incoming_u_i.u, &z1, &z2);

        let r = RandomOracle::<G1Projective, Instance>::get_challenge(u_i, incoming_u_i);
        let w = Self::fold_witness(r, w_i, incoming_w_i, &t);
        let u = Self::fold_instance(r, u_i, incoming_u_i, &t);

        (w, u, t)
    }

    fn verify(
        u_i: &Instance,
        incoming_u_i: &Instance,
        t: &Vec<BlsScalarField>,
    ) -> Instance {
        let r = RandomOracle::<G1Projective, Instance>::get_challenge(u_i, incoming_u_i);
        Self::fold_instance(r, u_i, incoming_u_i, t)
    }

    fn compute_T(
        r1cs: &R1CS,
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

#[cfg(test)]
pub mod tests {
    use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, CryptographicSponge};
    use ark_std::test_rng;

    use crate::{arith::r1cs::tests::{get_test_r1cs, get_test_z, get_test_z_split}, utils::poseidon_canonical_config};

    use super::{BlsScalarField, Nova};

    #[test]
    fn test_nova() {
        let r1cs = get_test_r1cs();
        let (w, x, u) = get_test_z_split(2);
        
        let mut rng = ark_std::test_rng();

        let poseidon_config = poseidon_canonical_config::<BlsScalarField>();
        let mut transcript_p = PoseidonSponge::<BlsScalarField>::new(&poseidon_config);
        
        const H: bool = false;

        let mut w_i = Nova::<H>::new_witness(&w, r1cs.a.len(), test_rng());
        let mut u_i = Nova::<H>::new_instance(&w_i, &x);

        for i in  0..4 {
            let (incoming_w, incoming_x, incoming_u) = get_test_z_split(i + 4);
            let incoming_w_i = Nova::<H>::new_witness(&incoming_w, r1cs.a.len(), test_rng());
            let incoming_u_i = Nova::<H>::new_instance(&incoming_w_i, &incoming_x);

            let (folded_witness, u, t) = Nova::<H>::prove(&r1cs, &w_i, &u_i, &incoming_w_i, &incoming_u_i);
            let folded_instance = Nova::<H>::verify(&u_i, &incoming_u_i, &t);

            w_i = folded_witness;
            u_i = folded_instance;
        }

        r1cs.check_relation(&w_i, &u_i);
    }
}
