#ifndef USEFUL_FUNCTS
#define USEFUL_FUNCTS

#include <stdint.h>
#include <gmp.h>

#include "pmns_data.h"
#include "pmns_arith_ops.h"

void init_data();

void free_data();

void from_int_to_pmns(int64_t *rop, mpz_t op);

void from_pmns_to_int(mpz_t rop, int64_t *op);

void exact_coeffs_reduction(int64_t *rop, int64_t *op);

void from_mont_domain(int64_t *rop, int64_t *op);

void copy_poly(int64_t *rop, int64_t *op);

void print_element(int64_t *poly);

//int64_t equality_check_polys(int64_t *pA, int64_t *pB);

int cmp_poly_evals(int64_t *pa, int64_t *pb);

#endif

