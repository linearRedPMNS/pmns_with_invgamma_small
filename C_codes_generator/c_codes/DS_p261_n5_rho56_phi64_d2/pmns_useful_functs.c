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
	{0x196d15f6cbed98L, 0x9361d21aa4355L, 0x1e5c7799113222L, -0x1916c1df17ba94L, -0x56c7db9745bc8L},
	{-0x2e466b1c5638acL, 0x1c5c821d8c3b18L, -0x1d92c9511ea4eeL, -0x1661b77800616cL, -0x1b741c351e2ca0L},
	{0x3aef5c71d4e3aaL, 0xe2e410ec41f52L, 0x1824bf27f61027L, 0x825f07367a0dL, -0x1f3c28ae44ace3L},
	{0x67f85467b8bdaL, 0x95ee0c3dfac9bL, 0x1ec2be3d83529bL, -0xb26d24b4b0f05L, 0x76551e2b392b4L},
	{-0x294ef130db6a50L, -0x6738761418e9aL, -0xfb0b42104d3ccL, 0xc4afb945858b0L, 0x1520d58de2580bL}
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
	mpz_set_str (modul_p, "3605246239461134813160563798075998591637758718398896416624541907425448462090271", 10);

	mpz_set_str (gama_pow[0], "3605246239460923457475022654799143484297760323690712621800523331487181246005279", 10);
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

