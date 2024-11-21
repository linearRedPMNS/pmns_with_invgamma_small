#ifndef STRUCTS_DATA
#define STRUCTS_DATA

#include <stdint.h>

typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 4
#define NB_COEFF 5
#define NB_ADD_MAX 2

#define E_ALPHA 1
#define E_LAMBDA -2

#define RHO_NBITS 56	// note: so 57 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 53

#define CONV_MASK 9007199254740991UL  // = (1 << CONVBASE_LOG2) - 1, for conversion


//Representations of polynomials Pi, used for conversion into the PMNS
//Each Pi is a representation of '(2^CONVBASE_LOG2)^i * alpha^{-1} * phi^2'
extern int64_t polys_P[NB_COEFF][NB_COEFF];

#endif

