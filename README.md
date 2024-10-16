# Friendly PMNS, with $M(X) =$  $\frac{1}{\gamma}$  $X - 1$
This repository contains code to generate Polynomial Modular Number Systems (PMNS) for a friendly class of primes. It also contains an efficient C codes generator. 
The generated PMNS are LinearRed or DoubleSparse, with $M(X) =$  $\frac{1}{\gamma}$  $X - 1$.

The PMNS generator is in the directory 'pmns_generator'.

A standard C codes generator is provied in the directory 'C_codes_generator'. Given a PMNS (obtained from the preceding generator), it generates C codes for all the main operations (including forward and backward conversion to PMNS, addition, multiplication).

In each directory, a file 'EXAMPLE.py' explains how to use the corresponding generator.

Note : To use the PMNS generator, you will need SageMath library which can be found here: http://www.sagemath.org/

