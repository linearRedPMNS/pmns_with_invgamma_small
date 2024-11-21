#include <stdlib.h>
#include <stdint.h>

#include "pmns_data.h"
#include "pmns_arith_ops.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

//~ --------------------------------------------------------------------

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void leftshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ --------------------------------------------------------------------

//~ computes : pa + 2.pb
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + 2*pb[j];
}

//~ computes : pa - 2.pb
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - 2*pb[j];
}

//~ --------------------------------------------------------------------

void add_lpoly(int128_t *rop, int128_t *pa, int128_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_lpoly(int128_t *rop, int64_t *op, uint64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = (int128_t)op[j] * scalar;
}

//~ --------------------------------------------------------------------

//~ Computes: pa*pb mod(E)
//~ Note: E(X) = 2*X^5 + 1
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = 2*((int128_t)pa[0]*pb[0]) - ((int128_t)pa[1]*pb[4] + (int128_t)pa[2]*pb[3] + (int128_t)pa[3]*pb[2] + (int128_t)pa[4]*pb[1]);
	tmp_prod_result[1] = 2*((int128_t)pa[0]*pb[1] + (int128_t)pa[1]*pb[0]) - ((int128_t)pa[2]*pb[4] + (int128_t)pa[3]*pb[3] + (int128_t)pa[4]*pb[2]);
	tmp_prod_result[2] = 2*((int128_t)pa[0]*pb[2] + (int128_t)pa[1]*pb[1] + (int128_t)pa[2]*pb[0]) - ((int128_t)pa[3]*pb[4] + (int128_t)pa[4]*pb[3]);
	tmp_prod_result[3] = 2*((int128_t)pa[0]*pb[3] + (int128_t)pa[1]*pb[2] + (int128_t)pa[2]*pb[1] + (int128_t)pa[3]*pb[0]) - ((int128_t)pa[4]*pb[4]);
	tmp_prod_result[4] = 2*((int128_t)pa[0]*pb[4] + (int128_t)pa[1]*pb[3] + (int128_t)pa[2]*pb[2] + (int128_t)pa[3]*pb[1] + (int128_t)pa[4]*pb[0]);

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and puts the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[5];
	int128_t tmpZero[5];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi)
	tmpQ[0] = (int64_t)op[0]*(-0x4aa2a56cb44550afL) + (int64_t)op[1]*(0x1a49066c20af13f8L) + (int64_t)op[2]*(0x3b65b545c4ef8174L) + (int64_t)op[3]*(0x7e3f5963822b0f6eL) + (int64_t)op[4]*(-0x5b41b64c1e7480fbL);
	tmpQ[1] = tmpQ[0]*(0x8c9ab5d9dca116L) + (int64_t)op[1];
	tmpQ[2] = tmpQ[1]*(0x8c9ab5d9dca116L) + (int64_t)op[2];
	tmpQ[3] = tmpQ[2]*(0x8c9ab5d9dca116L) + (int64_t)op[3];
	tmpQ[4] = tmpQ[3]*(0x8c9ab5d9dca116L) + (int64_t)op[4];

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = -(int128_t)tmpQ[0] + (int128_t)tmpQ[4]*(-0x464d5aecee508bL);
	tmpZero[1] = (int128_t)tmpQ[0]*(0x8c9ab5d9dca116L) - (int128_t)tmpQ[1];
	tmpZero[2] = (int128_t)tmpQ[1]*(0x8c9ab5d9dca116L) - (int128_t)tmpQ[2];
	tmpZero[3] = (int128_t)tmpQ[2]*(0x8c9ab5d9dca116L) - (int128_t)tmpQ[3];
	tmpZero[4] = (int128_t)tmpQ[3]*(0x8c9ab5d9dca116L) - (int128_t)tmpQ[4];

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
}

