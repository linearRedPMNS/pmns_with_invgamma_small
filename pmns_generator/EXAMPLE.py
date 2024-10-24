load("primes_and_pmns.sage")

phi_log2 = 64	

p_size = 256	#the (minimum) bit-size we want for the modulus

#the generation process might find suitable PMNS with modulus bit-size slightly bigger than 'p size', without any penalty on both efficiency and memory requirement. 
#'lv_tolerance' (below) allows to decide whether or not we want such PMNS, with the maximum size we accept.
#i.e.: we accept primes of bit-size <= (p_size + lv_tolerence); if 'lv_tolerence < 0', then no restriction is applied on modulus size.
lv_tolerance = -1	

double_spare = True	#'True' for DoubleSparse PMNS or 'False' for LinearRed PMNS

alpha_max = 4	#maximum absolute value allowed for 'alpha'

lambda_max = 4	#maximum absolute value allowed for 'lambda'

delta = 0		#the (minimum) number of desired "free" additions before a modular multplication  

n = (p_size//phi_log2) + 1	#should be incremented (if the generator takes too long to find a PMNS) as much as necessary, especially if "delta > 0" 

nb_pmns = 6	#the desired number of PMNS (if any)


p_list = look_for_good_primes(p_size, lv_tolerance, n, phi_log2, delta, alpha_max, lambda_max, double_spare, nb_pmns)	#THE PMNS GENERATOR


#Note (result struct): [p.nbits(), p, n, str(E), t, phi_log2, delta]

#Extended parameters can be chosen/uncommented in the function 'ckeck_poly' 


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


# ~ Some examples:  (DoubleSparse PMNS)

[261, 3605246239461134813160563798075998591637758718398896416624541907425448462090271, 5, 'X^5 + 2', -18030010715602944, 64, 2]

[385, 58992570678302246419582798577264301612336542685365472869468543189994166969222906573321650316811984155055226880000003, 7, '3*X^7 + 1', 34568881800478720, 64, 1]

[516, 170868453011949448085316123845124502170952331273575992081283315980217240856554505137532688039136697475600968124107411573766474618325591563528039840957342029, 9, '2*X^9 - 1', 326387110422511616, 64, 0]


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


# ~ Some examples:  (LinearRed PMNS)

[259, 516748990391900428273385526572212742018522329481256120716365298411794909129283, 5, '2*X^5 + 1', 39576602574364950, 64, 2]

[384, 29444190713823033278873681323415084508904983637646517819667723087302133133307051023849126398573545314912153522243303, 7, '4*X^7 - 1', 41333260579584285, 64, 1]

[517, 223284846011310584854764074522789473035608291352134700025126499663575997441956415852141217125330921839333699936334766034285376672385829925537109374999999999, 9, '3*X^9 - 1', 206061963501113250, 64, 0]


#-----------------------------------------------------------------------
