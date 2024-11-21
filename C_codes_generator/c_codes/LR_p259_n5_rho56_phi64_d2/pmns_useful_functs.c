#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

#include "pmns_data.h"
#include "pmns_arith_ops.h"
#include "pmns_useful_functs.h"

//~ --------------------------------------------------------------------

mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];

//Representations of polynomials Pi, used for conversion into the PMNS
//Each Pi is a representation of '(2^CONVBASE_LOG2)^i * alpha^{-1} * phi^2'
int64_t polys_P[NB_COEFF][NB_COEFF] =
{
	{-0x30fd13e2872fL, -0x23c8131daef7b7L, 0x3a2c359b7419b5L, 0x29ee01ed359bfL, 0x117ee91e88f299L},
	{0x147ed787b9950aL, -0x2bc0ab91b51d8dL, 0x1d7b71e045d9c2L, 0xe7004fd984422L, 0x1cc864b75ab960L},
	{0x34f5aabd8100bL, -0x2156caed0eb36eL, -0x3c426ec427b56bL, 0x3f37c9daff2b00L, -0x2a98a19423187L},
	{0x12f34555e67584L, 0x4573c6e35e5d34L, 0x144b3c9d444ca5L, -0xa19d917489bdbL, -0x79ce4851fa7c6L},
	{-0x499c604b3b0c0L, -0x216ac48e1a361fL, 0x13ca4b4f735119L, 0x2d43f9e5e93440L, 0x5d4603d96ee8aL}
};

//~ --------------------------------------------------------------------

//~ Assumes allocation already done for 'rop'.
//~ IMPORTANT : convertion to montgomery domain will be done here, thanks to conv polys 'polys_P'
//~ IMPORTANT : also assumes that the case 'alpha > 1' is already taken into account in conv polys 'polys_P'
void from_int_to_pmns(int64_t *rop, mpz_t op){

	int i;
	mpz_t tmp;
	int128_t tmp_poly[NB_COEFF];
	int128_t tmp_sum[NB_COEFF];

	mpz_init_set(tmp, op);

	for(i=0; i<NB_COEFF; i++){
		rop[i] = 0;
		tmp_sum[i] = 0;
	}

	if(tmp->_mp_size == 0)
		return;

	i = 0;
	while(tmp->_mp_size && (i < NB_COEFF)){
		scalar_mult_lpoly(tmp_poly, polys_P[i++], (tmp->_mp_d[0]) & CONV_MASK);
		add_lpoly(tmp_sum, tmp_sum, tmp_poly);

		mpz_tdiv_q_2exp (tmp, tmp, CONVBASE_LOG2);
	}

	internal_reduction(rop, tmp_sum);

	mpz_clear(tmp);
}

//~ Assumes "rop" already initialized.
//~ IMPORTANT : convertion from montgomery domain will be done here.
//~ IMPORTANT : the case 'alpha > 1' is also taken into account here
void from_pmns_to_int(mpz_t rop, int64_t *op){
	int i;
	mpz_t tmp_sum;
	int64_t tmp_conv[NB_COEFF];

	mpz_init(tmp_sum);

	//~ convertion out of mont domain
	from_mont_domain(tmp_conv, op);

	mpz_set_si(rop, tmp_conv[0]);
	for(i=0; i<POLY_DEG; i++){
		mpz_mul_si(tmp_sum, gama_pow[i], tmp_conv[i+1]);
		mpz_add(rop, rop, tmp_sum);
	}
	mpz_mul_si(rop, rop, E_ALPHA);	// we deal with 'alpha' here
	mpz_mod (rop, rop, modul_p);

	mpz_clear(tmp_sum);
}

//~ ----------------------------------------------------------------------------------------

void init_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_init (gama_pow[i]);
	mpz_init (modul_p);
	mpz_set_str (modul_p, "516748990391900428273385526572212742018522329481256120716365298411794909129283", 10);
	mpz_set_str (gama_pow[0], "516748990390673768726014642635731803952333562134963832412458542806198656004283", 10);
	for(i=1; i<POLY_DEG; i++){
		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);
		mpz_mod (gama_pow[i], gama_pow[i], modul_p);
	}
}

void free_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_clear (gama_pow[i]);
	mpz_clear (modul_p);
}

//~ ----------------------------------------------------------------------------------------

void exact_coeffs_reduction(int64_t *rop, int64_t *op){

	int i;
	int128_t tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++){
		tmp[i] = (int128_t) op[i];
	}

	internal_reduction(rop, tmp);

	mult_mod_poly(rop, rop, polys_P[0]);
}

//~ computes : op/phi
void from_mont_domain(int64_t *rop, int64_t *op){

	int i;
	int128_t tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++){
		tmp[i] = (int128_t) op[i];
	}

	internal_reduction(rop, tmp);
}

//~ ----------------------------------------------------------------------------------------

void copy_poly(int64_t *rop, int64_t *op){
	int i;
	for(i=0; i<NB_COEFF; i++)
		rop[i] = op[i];
}

void print_element(int64_t *poly){
	int i;
	printf("[");
	for (i=0; i<POLY_DEG; i++)
		printf("%2ld, ", poly[i]);
	printf("%2ld]", poly[POLY_DEG]);
}

//~ ----------------------------------------------------------------------------------------

//~ return a positive value if pa > pb, zero if pa = pb, or a negative value if pa < pb.
//~ Important : evaluation is done using the corresponding integers modulo 'p'.
int cmp_poly_evals(int64_t *pa, int64_t *pb){
	int rep;
	mpz_t a, b;
	mpz_inits (a, b, NULL);
	from_pmns_to_int(a, pa);
	from_pmns_to_int(b, pb);
	rep = mpz_cmp (a, b);
	mpz_clears (a, b, NULL);
	return rep;
}

