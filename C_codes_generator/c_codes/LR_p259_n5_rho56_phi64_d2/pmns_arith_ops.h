#ifndef POLY_MULT_ADD
#define POLY_MULT_ADD

#include <stdint.h>

#include "pmns_data.h"

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void add_poly(int64_t *rop, int64_t *pa, int64_t *pb);

void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void neg_poly(int64_t *rop, int64_t *op);
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar);
void double_poly_coeffs(int64_t *rop, int64_t *op);
void leftshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos);

void add_lpoly(int128_t *rop, int128_t *pa, int128_t *pb);

void scalar_mult_lpoly(int128_t *rop, int64_t *op, uint64_t scalar);

void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb);

//void square_mod_poly(int64_t *rop, int64_t *pa);

void internal_reduction(int64_t *rop, int128_t *op);

#endif

