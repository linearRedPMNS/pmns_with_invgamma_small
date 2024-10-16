load("primes_and_pmns.sage")

phi_log2 = 64	

p_size = 256	#the modulus (minimum) bit-size we want

lv_tolerence = -1	#we accept primes of bit-size <= (p_size + lv_tolerence); this parameter is ignored if < 0

double_spare = True	#put 'True' for DoubleSparse PMNS

alpha_max = 4	#maximum absolute value allowed for 'alpha'

lambda_max = 4	#maximum absolute value allowed for 'lambda'

delta = 0		#the number of desired "free" add before a modular multplication  

n = (p_size//phi_log2) + 1	#should be incremented (if the generator takes too long to find a PMNS) as much as necessary, especially if "delta > 0"

nb_pmns = 10	#the desired number of PMNS (if any)


p_list = look_for_good_primes(p_size, lv_tolerence, n, phi_log2, delta, alpha_max, lambda_max, double_spare, nb_pmns)	#THE GENERATOR


#Note (result struct): [p.nbits(), p, n, str(E), t, phi_log2, delta]
#Extended parameters can be chosen/uncommented in the function 'ckeck_poly' 

#-----------------------------------------------------------------------





