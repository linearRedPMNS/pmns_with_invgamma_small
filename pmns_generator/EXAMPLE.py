load("primes_and_pmns.sage")

phi_log2 = 64	

p_size = 256	#the modulus (minimum) bit-size we want

#the generation process might find suitable PMNS with modulus bit-size slightly bigger than 'p size', without any penalty on both efficiency and memory requirement. 
#'lv_tolerance' (below) allows to decide whether or not we want such PMNS, with the maximum size we accept.
#i.e.: we accept primes of bit-size <= (p_size + lv_tolerence); if 'lv_tolerence < 0', then no restriction is applied on modulus size.
lv_tolerance = -1	

double_spare = True	#'True' for DoubleSparse PMNS or 'False' for LinearRed PMNS

alpha_max = 4	#maximum absolute value allowed for 'alpha'

lambda_max = 4	#maximum absolute value allowed for 'lambda'

delta = 0		#the number of desired "free" add before a modular multplication  

n = (p_size//phi_log2) + 1	#should be incremented (if the generator takes too long to find a PMNS) as much as necessary, especially if "delta > 0"

nb_pmns = 6	#the desired number of PMNS (if any)


p_list = look_for_good_primes(p_size, lv_tolerance, n, phi_log2, delta, alpha_max, lambda_max, double_spare, nb_pmns)	#THE GENERATOR


#Note (result struct): [p.nbits(), p, n, str(E), t, phi_log2, delta]
#Extended parameters can be chosen/uncommented in the function 'ckeck_poly' 

#-----------------------------------------------------------------------





