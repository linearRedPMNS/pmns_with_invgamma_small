########################################################################

#note: it is not necessary to consider the case 'alpha < 0'
def gen_binomial_polyE(n, alpha_max, lambda_max):
	mid_part = [0]*(n-1)
	lambda_set = list(range(-lambda_max, lambda_max+1))
	lambda_set.remove(0)
	extPol_set = []
	for a in range(1, alpha_max+1):
		for b in lambda_set:
			if gcd(a,b) != 1:	#equivalent (and better) has been already checked
				continue
			bino_E = [b] + mid_part + [a]
			extPol_set.append(bino_E)
			#extPol_set.append(tuple(bino_E))
	extPol_set.sort(key=binoE_sorting_criteria)
	return extPol_set


def binoE_sorting_criteria(ext_pol_coeffs):
	n = len(ext_pol_coeffs) - 1
	return compute_binoE_w(n, ext_pol_coeffs[-1], ext_pol_coeffs[0])


def compute_binoE_w(n, E_alpha, E_lambda):
	a = n * abs(E_alpha)
	b = abs(E_alpha) + (n-1)*abs(E_lambda)
	return max(a,b)

########################################################################

def compute__t_min__t_max(p_size, n, alpha_max, lambda_max, phi_log2, delta):
	a = (1 << (p_size-1)) - alpha_max
	b = (a / lambda_max)**(1/n)
	t_min = ceil(N(b))
	#---------------------------------
	d = 2 * n * ((delta + 1)**2)
	c = (1 << phi_log2) / d
	t_max = floor(N(c)) + 1
	#---------------------------------
	return (t_min, t_max)


def compute__l_max(n, alpha_max, lambda_max, t_max):
	y_max = lambda_max*(t_max**n) + alpha_max
	l_max = ceil(N(log(y_max, 2)))
	return l_max

########################################################################

def check_bound_on_phi(phi_log2, n, E_w, E_alpha, E_lambda, delta, t):
	a = 2 * E_w * ((delta + 1)**2)
	if t%E_alpha == 0:  #a case mentioned in a remark of the paper
		s = t//E_alpha
		b = abs(s*E_lambda)
		rho = max(b, abs(t)) 
	else:
		b = abs(t*E_lambda) + 1
		c = abs(t) + abs(E_alpha)
		rho = max(b, c) - 1
	min_phi = a * (rho - 1)
	phi = 1 << phi_log2
	if phi > min_phi:
		return True
	return False 


def clean_int(y, min_size, max_size, trial_div_limit):
	if N(log(y,2)) <= (min_size-1):
		return 0
	else:
		L = list(y.factor(limit=trial_div_limit))
		p = L[-1][0]
		p_size = N(log(p, 2))
		if (p_size <= (min_size-1)) or (not p.is_prime()):
			return 0
		if (max_size < 0) or (p_size <= max_size):
			return p	
	return 0


def ckeck_poly(abs_t, ext_pol, n, phi_log2, delta, min_size, max_size, tdiv_limit):	
	lmbd = -ext_pol[0]
	alph = ext_pol[-1]
	w = compute_binoE_w(n, alph, lmbd)
	if not check_bound_on_phi(phi_log2, n, w, alph, lmbd, delta, abs_t):
		return []
	R.<X> = ZZ[]
	p_list = []
	for t in {-abs_t, abs_t}:
		y = Integer(lmbd*(t**n) - alph)
		p = clean_int(y, min_size, max_size, tdiv_limit)
		if p == 0:
			continue
		E = R(ext_pol)
		if alph < 0: #if E polys generation process included this case 
			E = -E
		main_params = [p.nbits(), p, n, str(E), t, phi_log2, delta]
		# ~ #--------------------------------------------
		# ~ gmm = (1/t)%p
		# ~ if t%alph == 0:  #a case mentioned in a remark of the paper
			# ~ s = t//alph
			# ~ b = abs(s*lmbd)
			# ~ rho = max(b, abs(t)) 
		# ~ else:
			# ~ b = abs(t*lmbd) + 1
			# ~ c = abs(t) + abs(alph)
			# ~ rho = max(b, c) - 1
		# ~ M = R([-1,t])
		# ~ extended_params = [p.nbits(), p, n, gmm, rho, str(E), str(M), phi_log2, delta]
		# ~ #--------------------------------------------
		p_list.append(main_params)
		print(main_params)
		print()
	return p_list

########################################################################

def look_for_doublesparse(n, phi_log2, delta, sorted_bino_polys, t_min, t_max, min_size, max_size, trial_div_limit, nb_pmns):
	shift_pos = ceil(phi_log2 / 2)
	min_up = t_min >> shift_pos
	max_up = t_max >> shift_pos
	print(f'DoubleSparse PMNS, with (min_up: {min_up}, max_up: {max_up})\n')
	p_list = []
	nb_found = 0
	# ~ for v in range(min_up, max_up):
	for j in range(min_up, max_up):
		v = randint(min_up, max_up)
		abs_t = v << shift_pos
		for ext_pol in sorted_bino_polys:
			ps = ckeck_poly(abs_t, ext_pol, n, phi_log2, delta, min_size, max_size, trial_div_limit)
			if ps != []:
				nb_found += len(ps)
				p_list.append(ps)
			if nb_found >= nb_pmns:
				pmns_list = flatten(p_list, max_level=1)
				return pmns_list[:nb_pmns]
	return flatten(p_list, max_level=1)


def look_for_linearred(n, phi_log2, delta, sorted_bino_polys, t_min, t_max, min_size, max_size, trial_div_limit, nb_pmns):
	print(f'LinearRed PMNS, with (t_min: {t_min}, t_max: {t_max})\n')
	p_list = []
	nb_found = 0
	# ~ for abs_t in range(t_min, t_max):
	for j in range(t_min, t_max):
		abs_t = randint(t_min, t_max)
		for ext_pol in sorted_bino_polys:
			ps = ckeck_poly(abs_t, ext_pol, n, phi_log2, delta, min_size, max_size, trial_div_limit)
			if ps != []:
				nb_found += len(ps)
				p_list.append(ps)
			if nb_found >= nb_pmns:
				pmns_list = flatten(p_list, max_level=1)
				return pmns_list[:nb_pmns]
	return flatten(p_list, max_level=1)
	
########################################################################

#note: assumes 'alpha_max, lambda_max > 0'
def look_for_good_primes(p_size, lv_tolerance, n, phi_log2, delta, alpha_max, lambda_max, double_spare, nb_pmns):
	
	sorted_bino_polys = gen_binomial_polyE(n, alpha_max, lambda_max)
	
	print("\nExternal polynomials generation done: " + str(len(sorted_bino_polys)) + " candidates.\n")
	
	(t_min, t_max) = compute__t_min__t_max(p_size, n, alpha_max, lambda_max, phi_log2, delta)
	
	if t_min > t_max:
		print("No PMNS: t_min > t_max. Please increment the value of n.\n")
		return []
	
	min_size = p_size
	
	if lv_tolerance >= 0:
		max_size = p_size + lv_tolerance
	else:
		max_size = -1
	
	l_max = compute__l_max(n, alpha_max, lambda_max, t_max)
	max_tdiv_2pow = 20	#to increased (or decreased) depending on computer power
	trial_div_limit = 1 << min(l_max, max_tdiv_2pow)
	
	print("					########### STARTING ###########\n")
	
	pmns_found = []
	if double_spare:
		pmns_found = look_for_doublesparse(n, phi_log2, delta, sorted_bino_polys, t_min, t_max, min_size, max_size, trial_div_limit, nb_pmns)
	else:
		pmns_found = look_for_linearred(n, phi_log2, delta, sorted_bino_polys, t_min, t_max, min_size, max_size, trial_div_limit, nb_pmns)
		
	print("\n					########### ENDED ###########")
	
	return pmns_found

########################################################################
