use std::marker::PhantomData;

use ark_bls12_381::Bls12_381;
use ark_crypto_primitives::sponge::{poseidon::{PoseidonConfig, PoseidonSponge}, Absorb, CryptographicSponge};
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::{BigInteger, PrimeField};

pub type BlsScalarField = <Bls12_381 as Pairing>::ScalarField;

pub struct RandomOracle<
    G: CurveGroup,
    I: Absorb,
> {
   _g: PhantomData<G>,
   _i: PhantomData<I>,
}

impl<G: CurveGroup, I: Absorb> RandomOracle<G, I> {
    pub fn get_challenge(
        U_i: &I,
        u_i: &I,
    ) -> G::ScalarField {
        let poseidon_config = poseidon_canonical_config::<G::ScalarField>();
        let mut sponge = PoseidonSponge::<G::ScalarField>::new(&poseidon_config);

        sponge.absorb(&U_i);
        sponge.absorb(&u_i);
        let bits = sponge.squeeze_bits(128);
        G::ScalarField::from_bigint(BigInteger::from_bits_le(&bits)).unwrap()
    }
}

/// This Poseidon configuration generator agrees with Circom's Poseidon(4) in the case of BN254's scalar field
pub fn poseidon_canonical_config<F: PrimeField>() -> PoseidonConfig<F> {
    // 120 bit security target as in
    // https://eprint.iacr.org/2019/458.pdf
    // t = rate + 1

    let full_rounds = 8;
    let partial_rounds = 60;
    let alpha = 5;
    let rate = 4;

    let (ark, mds) = ark_crypto_primitives::sponge::poseidon::find_poseidon_ark_and_mds::<F>(
        F::MODULUS_BIT_SIZE as u64,
        rate,
        full_rounds,
        partial_rounds,
        0,
    );

    PoseidonConfig::new(
        full_rounds as usize,
        partial_rounds as usize,
        alpha,
        mds,
        ark,
        rate,
        1,
    )
}
