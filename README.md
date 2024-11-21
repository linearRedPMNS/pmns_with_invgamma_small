# Friendly PMNS, with $M(X) =$  $\frac{1}{\gamma}$  $X - 1$
This repository contains code to generate Polynomial Modular Number Systems (PMNS) for a friendly class of primes. It also contains an efficient C codes generator. 

The generated PMNS are LinearRed or DoubleSparse, with $M(X) =$  $\frac{1}{\gamma}$  $X - 1$ and $E(X) = \alpha X^n - \lambda$.
<br>
<br>

The PMNS generator is in the directory 'pmns_generator'.
<br>
<br>

A standard C codes generator is provided in the directory 'C_codes_generator'. <br>
Given a PMNS (obtained from the previous generator), it generates C codes for all the main operations (including forward and backward conversion to PMNS, addition, multiplication). The generated C codes are located in the 'c_codes' subdirectory. <br>
Depending on whether the PMNS is DoubleSparse or LinearRed, the corresponding C code will be in a directory whose name starts with 'DS' (for DoubleSparse PMNS) or 'LR' (for LinearRed PMNS). 
In both cases, the code is compiled with the command: 
> make

The above command also compiles a simple example of main (the file 'main.c') that can be run with the command:
> ./main
<br>

In each directory ('pmns_generator' and 'C_codes_generator'), an 'EXAMPLE.py' file shows how to use the corresponding generator.
<br>
<br>

Note : To use these generators, you will need SageMath library which can be found here: http://www.sagemath.org/

