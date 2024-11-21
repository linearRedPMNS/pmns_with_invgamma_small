#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

#include "gmp_stuff.c"

#include "pmns_data.h"
#include "pmns_arith_ops.h"
#include "pmns_useful_functs.h"

extern mpz_t modul_p; // from "pmns_useful_functs.c"

//~ Compilation command:	make
//~ Execution command:		./main
//~ Cleaning command:		make clean

//~ WARNING: if needed, change appropriately 'pa' and 'pb' type (lines 42 and 43) to avoid a lot of warning and potential errors

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){
	
	srand(time(NULL));	
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);

	int i, nbiter, nb_limbs;
	mpz_t A, B, C, E, F, R2, G, H;
	mpz_inits (A, B, C, E, F, R2, G, H, NULL);
	
	mp_limb_t mip0;
	
	const mp_limb_t *p_limbs;
	mp_limb_t *r2_limbs, *scratch_limbs, *mip_limbs;
	mp_limb_t *a_limbs, *b_limbs, *am1_limbs, *am2_limbs, *bm1_limbs, *bm2_limbs;
	
	int64_t pa[NB_COEFF];
	int64_t pb[NB_COEFF];
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	init_data();
	
	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	mpz_set(E, A);
	
	nb_limbs = mpz_size (modul_p);
	
	p_limbs = mpz_limbs_read (modul_p);   
	
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	mpz_setbit (R2, 2*nb_limbs*8*sizeof(mp_limb_t)); 
	mpz_mod(R2, R2, modul_p);
	r2_limbs = mpz_limbs_modify (R2, nb_limbs); 
	
	a_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	scratch_limbs = (mp_limb_t*) calloc (2*nb_limbs, sizeof(mp_limb_t));

	b_limbs = mpz_limbs_modify (B, nb_limbs);
	copy_limbs(a_limbs, A, nb_limbs);
	copy_limbs(am1_limbs, A, nb_limbs);
	copy_limbs(am2_limbs, A, nb_limbs);
	copy_limbs(bm1_limbs, B, nb_limbs);
	copy_limbs(bm2_limbs, B, nb_limbs);
	
	mpn_binvert (mip_limbs, p_limbs, nb_limbs, scratch_limbs); //clean_limbs(scratch_limbs, 2*nb_limbs);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	from_int_to_pmns(pa, A);
	from_int_to_pmns(pb, B);

	//~ conversion to Montgomery domain (mont par bloc)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	mpn_mont_mul_red_1(bm1_limbs, bm1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	
	//~ conversion to Montgomery domain (mont classique)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	mpn_mont_mul_red_n(bm2_limbs, bm2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	nbiter = 1 << 2;
	
	for (i=0; i<nbiter; i++) {
		mpz_mul (A, A, B);
		mpz_mod (A, A, modul_p);
	}
	
	for (i=0; i<nbiter; i++) {
		mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
	}

	for (i=0; i<nbiter; i++) {
		//~ Montgomery modular multiplication (mont par bloc)
		mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
	}
	
	for (i=0; i<nbiter; i++) {
		//~ Montgomery modular multiplication (mont classique)
		mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
	}
	
	for (i=0; i<nbiter; i++) {
		mult_mod_poly(pa, pa, pb);
	}

	//~ should not modify the value which is represented
	exact_coeffs_reduction(pa, pa);

	from_pmns_to_int(C, pa);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	from_limbs_to_mpz_t(F, a_limbs, nb_limbs);
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am1_limbs[i];
	
	mpn_redc_1 (am1_limbs, scratch_limbs, p_limbs, nb_limbs, mip0);
	
	from_limbs_to_mpz_t(G, am1_limbs, nb_limbs);
	
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am2_limbs[i];
	
	mpn_redc_n (am2_limbs, scratch_limbs, p_limbs, nb_limbs, mip_limbs);
	
	from_limbs_to_mpz_t(H, am2_limbs, nb_limbs);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	printf("\nnbiter = %d\n\n", nbiter);
	gmp_printf("p       : %Zd\n\n", modul_p);
	gmp_printf("A       : %Zd\n", E);
	gmp_printf("B       : %Zd\n\n", B);
	gmp_printf("r_pmns  : %Zd\n", C);
	gmp_printf("r_gmp   : %Zd\n", A);
	gmp_printf("r_lgmp  : %Zd\n", F);
	gmp_printf("r_mcgmp : %Zd\n", H);
	gmp_printf("r_mbgmp : %Zd\n", G);
	
	mpz_clears (A, B, C, E, F, R2, G, H, NULL);
	gmp_randclear(r);
	
	free(a_limbs);
	free(am1_limbs);
	free(am2_limbs);
	free(bm1_limbs);
	free(bm2_limbs);
	free(mip_limbs);
	free(scratch_limbs);
	
	free_data();
	return 0;
}

















