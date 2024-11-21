########################################################################
########################################################################

def get_sign(k):
	if k < 0 :
		return '-'
	else :
		return '+'


def is_power_of_two(k) :
	if k <= 0 :
		return False
	return ((k & (k - 1)) == 0)


def same_sign(id_list):
	s = id_list[0].sign()
	for el in id_list[1:]:
		if el.sign() != s :
			return [False, 1]
	return [True, s]

########################################################################
########################################################################
#~ ---------- for E(X) = alpha*X^n - lambda ----------------------------

#NOTE: assumes 'alpha > 0'
def build_prod_code(n, alpha, lambd, big_int):
	
	c_sign = get_sign(lambd)
	abs_c = abs(lambd)
	
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	tmp_part3 = ['']*n
	
	for i in range(n):
		tmp_part3[i] = 'tmp_prod_result['+str(i)+'] = '
	
	for i in range(n):
		for j in range(n):
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')pa['+ str(i) + ']*pb['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')pa['+ str(i) + ']*pb['+ str(j) + '] +'
	
	if alpha != 1:
		for i in range(n):
			tmp_part1[i] = str(alpha) + '*(' + tmp_part1[i][1:-2] +')'
	else:
		for i in range(n):
			tmp_part1[i] = tmp_part1[i][1:-2]
			
	for i in range(n-1):
		# ~ if abs_c == 1 :
			# ~ tmp_part2[i] = ' ' + c_sign + ' (' + tmp_part2[i][1:-2] +')'
		if lambd == 1:
			tmp_part2[i] = ' + ' + tmp_part2[i][1:-2]
		elif lambd == -1:
			tmp_part2[i] = ' - (' + tmp_part2[i][1:-2] +')'
		else :
			tmp_part2[i] = ' ' + c_sign + ' ' + str(abs_c) + '*(' + tmp_part2[i][1:-2] +')'
	
	result = '\n'
	for i in range(n-1):
		result += "	" + tmp_part3[i] + tmp_part1[i] + tmp_part2[i] + ';\n'
	result +=  "	" + tmp_part3[n-1] + tmp_part1[n-1] + ';\n'
	
	return result

########################################################################
########################################################################
#~ ------------- internal reduction methods ----------------------------

def build_red_int_code(small_int, big_int, n, phi_log2, red_int_mat, neginv_red_int_mat, double_sparse, E_alpha):
	
	#print(f'\n{red_int_mat}\n\n{red_int_mat.det()}\n')  #<-- HERE
	
	result = "	" + small_int
	result += " tmpQ[" + str(n) + "];\n"
	
	result += "	" + big_int
	result += " tmpZero[" + str(n) + "];\n\n"
	
	result += "	//~ computation of : op*neginv_red_int_mat mod(mont_phi)\n"
	
	if double_sparse:
		result += build_redInt_tmpQ_prod_code(n, phi_log2, neginv_red_int_mat, small_int)
	else:
		polM_t = red_int_mat[0][1]
		alpha_div_t = (polM_t%E_alpha == 0)
		result += build_redInt_tmpQ_prod_code_LR(n, phi_log2, neginv_red_int_mat, polM_t, alpha_div_t, small_int)
	
	result += "\n	//~ computation of : tmp_q*red_int_mat\n"
	result += build_redInt_tmpZero_prod_code(n, red_int_mat, big_int)
	
	tmp = ''
	result += "\n	//~ computation of : (op + tmp_zero)/mont_phi\n"
	
	for i in range(n):
		tmp = "	rop[" + str(i) + "] = (op[" + str(i) + "] + tmpZero[" + str(i) + "]) >> PHI_LOG2;\n"
		result += tmp
	
	return result

########################################################################

# IMPORTANT : 'big_int' is assumed big enough so that overflow will not happen on target architecture.
def build_redInt_tmpZero_prod_code(n, red_int_mat, big_int):
	
	tmp_partStart = ['']*n
	for i in range(n):
		tmp_partStart[i] = 'tmpZero['+str(i)+'] = '
	
	prod_list = build_scal_prod(red_int_mat, 1, "tmpZero", big_int=big_int)
	
	result = ""
	for i in range(n):
		if prod_list[i] == '' :
			result += "	" + tmp_partStart[i] + '0;\n'
		else :
			if prod_list[i][1] == '-':
				result += "	" + tmp_partStart[i] + '-' + prod_list[i][3:] + ';\n'
			else:
				result += "	" + tmp_partStart[i] + prod_list[i][3:] + ';\n'
	
	return result


#~ note : jump > 0
def build_scal_prod(mat, jump, call_for, big_int=None, small_int=None):
	
	(n, m) = mat.dimensions()
	tmp = []
	for i in range(m):
		
		tmp_dict = {}
		
		for j in range(n):
			val = Integer(mat[j][i])
			k = abs(val)
			if k == 0 :
				continue
			s = val.sign() 
			try:
				tmp_dict[k].append(s*(j+jump))
			except KeyError :
				tmp_dict[k] = [s*(j+jump)]
		
		if tmp_dict == {} :
			tmp.append("")
		else :
			if call_for == "tmpZero" :
				tmp.append(build_cell_tmpZero(tmp_dict, big_int))
	return tmp


def build_cell_tmpZero(cell_dict, big_int):
	
	ch = ""
	for k in cell_dict.keys() : # keys are > 0
		
		id_list = cell_dict[k]
		
		if len(id_list) == 1 :
			
			idd = id_list[0]
			
			if k == 1 :
				ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd)-1) + "]"
			else :
				ch += " + (" + big_int + ")tmpQ[" + str(abs(idd)-1) + "]*(" + str(hex(k*idd.sign())) + "L)"
			
		elif k == 1 :
			for idd in id_list :
				ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd)-1)+ "]"
		else :
		
			[vr, s] = same_sign(id_list)
			
			if vr :
				ch += " + ((" + big_int + ")tmpQ[" + str(s*id_list[0] - 1)+ "]"
				for idd in id_list[1:] :
					ch += " + (" + big_int + ")tmpQ[" + str(s*idd - 1)+ "]"
				k *= s 
			else :
				ch += " + ("
				
				idd = id_list[0] 
				if idd > 0 :
					ch += "(" + big_int + ")tmpQ[" + str(idd - 1)+ "]"
				else :
					ch += "-(" + big_int + ")tmpQ[" + str(-idd - 1)+ "]"
				for idd in id_list[1:] :
					ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd) - 1)+ "]"
		
			ch += ")*(" + str(hex(k)) + "L)"
	
	return ch

########################################################################

#~ note : elements of 'neg_inv_ri_rep_coeff' are >= 0
def build_redInt_tmpQ_prod_code(n, phi_log2, neginv_red_int_mat, small_int):
	
	tmp_partStart = ['']*n
	for i in range(n):
		tmp_partStart[i] = 'tmpQ['+str(i)+'] = '
	
	prod_list = build_tmpQ_scal_prod(small_int, n, phi_log2, neginv_red_int_mat)
	
	result = ""
	for i in range(n):
		if prod_list[i] == '' :
			result += "	" + tmp_partStart[i] + '0;\n'
		else :
			if prod_list[i][1] == '-':
				result += "	" + tmp_partStart[i] + '-' + prod_list[i][3:] + ';\n'
			else:
				result += "	" + tmp_partStart[i] + prod_list[i][3:] + ';\n'
	
	return result


def build_tmpQ_scal_prod(small_int, n, phi_log2, neginv_red_int_mat):
	phi = 1 << phi_log2
	upB = Integer(phi/2)
	tmp = []
	for i in range(n):
		tmp_dict = {}
		for j in range(n):
			k = Integer(neginv_red_int_mat[j][i]) #note: 0 <= k < phi
			if k == 0 :
				continue
			elif k < upB:
				s = 1
			else:
				s = -1
				k = phi - k
			try:
				tmp_dict[k].append(s*(j+1))
			except KeyError :
				tmp_dict[k] = [s*(j+1)]
		if tmp_dict == {} :
			tmp.append("")
		else :
			tmp.append(build_cell_tmpQ(tmp_dict, small_int))
	return tmp


def build_cell_tmpQ(cell_dict, small_int):
	
	ch = ""
	for k in cell_dict.keys() : # keys are > 0
		
		id_list = cell_dict[k]
		
		if len(id_list) == 1 :
			
			idd = id_list[0]
			
			if k == 1 :
				ch += " " + get_sign(idd) + " (" + small_int + ")op[" + str(abs(idd)-1) + "]"
			else :
				ch += " + (" + small_int + ")op[" + str(abs(idd)-1) + "]*(" + str(hex(k*idd.sign())) + "L)"
			
		elif k == 1 :
			for idd in id_list :
				ch += " " + get_sign(idd) + " (" + small_int + ")op[" + str(abs(idd)-1)+ "]"
		else :
		
			[vr, s] = same_sign(id_list)
			
			ch += " + ("
			
			if vr :
				ch += "(" + small_int + ")op[" + str(s*id_list[0] - 1)+ "]"
				for idd in id_list[1:] :
					ch += " + (" + small_int + ")op[" + str(s*idd - 1)+ "]"
				k *= s 
			else :
				idd = id_list[0] 
				if idd > 0 :
					ch += "(" + small_int + ")op[" + str(idd - 1)+ "]"
				else :
					ch += "-(" + small_int + ")op[" + str(-idd - 1)+ "]"
				for idd in id_list[1:] :
					ch += " " + get_sign(idd) + " (" + small_int + ")op[" + str(abs(idd) - 1)+ "]"
		
			ch += ")*(" + str(hex(k)) + "L)"
	
	return ch

########################################################################

#~ note : elements of 'neg_inv_ri_rep_coeff' are >= 0
def build_redInt_tmpQ_prod_code_LR(n, phi_log2, neginv_red_int_mat, t, alpha_div_t, small_int):
	
	tmp_partStart = ['']*n
	for i in range(n):
		tmp_partStart[i] = 'tmpQ['+str(i)+'] = '
	
	phi = 1 << phi_log2
	
	first_col = list(neginv_red_int_mat.column(0))
	first_res = build_scal_prod_tmpQ_LR(first_col, n, phi, small_int)
	
	result = "	" + tmp_partStart[0] + first_res + ';\n'
	
	cphi_t = centered_t(t, phi)
	
	for i in range(1, n-1):
		result += "	" + tmp_partStart[i] + 'tmpQ['+str(i-1)+']*' + cphi_t
		result += ' + (' + small_int + ')op['+str(i)+'];\n'
	
	if alpha_div_t:
		result += "	" + tmp_partStart[n-1] + 'tmpQ['+str(n-2)+']*' + cphi_t
		result += ' + (' + small_int + ')op['+str(n-1)+'];\n'
	else:
		last_col = list(neginv_red_int_mat.column(-1))
		last_res = build_scal_prod_tmpQ_LR(last_col, n, phi, small_int)
		result += "	" + tmp_partStart[n-1] + last_res + ';\n'

	return result


def centered_t(t, phi):
	a = t % phi
	if a < Integer(phi/2):
		ch = '(' + str(hex(a)) + 'L)'
	else:
		ch =  '(-' + str(hex(phi - a)) + 'L)'
	return ch


def build_scal_prod_tmpQ_LR(col, n, phi, small_int):
	
	phi_c = Integer(phi/2)
	
	tmp_dict = {}
	for j in range(n):
		k = Integer(col[j]) #note: 0 <= k < phi
		if k == 0 :
			continue
		elif k < phi_c:
			s = 1
		else:
			s = -1
			k = phi - k
		try:
			tmp_dict[k].append(s*(j+1))
		except KeyError:
			tmp_dict[k] = [s*(j+1)]
	
	tmp = build_cell_tmpQ(tmp_dict, small_int)
	
	if tmp[1] == '+':
		tmp = tmp[3:]
	else:
		tmp = '-' + tmp[3:]
	
	return tmp

########################################################################
########################################################################
#~ --------------- conversion polynomials P  ---------------------------

def generate_conv_polys_P(small_int_name, p, n, convbase_log2, phi, ri_mat, neg_iri_mat, iG_dom, Zp, ext_pol, R):
	
	conv_polys = compute_conv_polys(p, n, convbase_log2, phi, ri_mat, neg_iri_mat, iG_dom, Zp, ext_pol, R)
	
	result = small_int_name + " polys_P[NB_COEFF][NB_COEFF] =\n"
	result += "{\n"
	
	for i in range(n-1):
		result += "	{" + list__from_int_to_str_hex(list(conv_polys[i])) + "},\n"
	result += "	{" + list__from_int_to_str_hex(list(conv_polys[n-1])) + "}\n"
	
	result += "};\n\n"
	
	return result 


def list__from_int_to_str_hex(llist, el_type='L'):
	res = ''
	for el in llist[:-1]:
		res += hex(el) + el_type + ', '
	res += hex(llist[-1]) + el_type
	return res

########################################################################

#~ computes the polynomials Pi, for conversion: Pi ~ (convbase^i * phi^2 * alpha^{-1})
#~ Note: also works with non monic polynomial E
def compute_conv_polys(p, n, convbase_log2, phi, ri_mat, neg_iri_mat, iG_dom, Zp, ext_pol, R):
	
	phi_pown = phi.powermod(n, p)
	
	alph = list(ext_pol)[-1]
	
	ialpha = (1/alph)%p
	
	ialpha_phi2 = (ialpha * (phi**2)) % p
	
	convbase = 1 << convbase_log2
	
	init_val = (convbase * phi * ialpha) % p
	
	tp = exact_conv_to_pmns(ialpha_phi2, n, phi, ri_mat, neg_iri_mat, iG_dom, Zp, phi_pown)

	rp = exact_conv_to_pmns(init_val, n, phi, ri_mat, neg_iri_mat, iG_dom, Zp, phi_pown)
	
	conv_polys = [tp]
	
	for i in range(1,n):
		
		tp = pmns_bino_mont_mult(tp, rp, n, ri_mat, neg_iri_mat, iG_dom, phi, ext_pol, R)
		
		conv_polys.append(tp)
	
	return conv_polys


def exact_conv_to_pmns(val, n, phi, ri_mat, neg_iri_mat, iG_dom, Zp, phi_pown):
	v = [Integer(Zp(val * phi_pown))] + [0]*(n-1)
	rep = vector(ZZ, v)
	for i in range(n):
		rep = pmns_red_int(rep, n, ri_mat, neg_iri_mat, iG_dom, phi)
	return rep


#~ returns a representation of 'op/phi', WITH centered reduction modulo phi
# ~ Note: assumes that 'op' is of dim 'n'
def pmns_red_int(op, n, ri_mat, neg_iri_mat, iG_dom, phi):
	q = vector(iG_dom, op) * neg_iri_mat
	t = center_vect_coeffs(q, phi) * ri_mat  # to look here
	r = op + t
	return (r/phi)


# ~ IMPORTANT: assumes that 'phi' is even and 'vect' elements are in Zphi
def center_vect_coeffs(vect, phi):
	cmax = Integer(phi/2)
	res = []
	for c in vect:
		c = Integer(c) # note: 0 <= c < phi
		if c < cmax:
			res.append(c)
		else:
			res.append(-phi + c)
	return vector(ZZ, res)

########################################################################
########################################################################

#~ returns a representation of '(op1 * op2 * alpha) / phi'
#~ IMPORTANT: assumes 'E(X) = alph*X^n - lmbd'
def pmns_bino_mont_mult(op1, op2, n, ri_mat, neg_iri_mat, iG_dom, phi, ext_pol, R):
	
	ext_polC = list(ext_pol)
	alph = ext_polC[-1]
	lmbd = -ext_polC[0]
	
	op11 = R(list(op1))
	op22 = R(list(op2))
	
	pprod = list(op11 * op22)
	c = pprod + [0]*(2*n - 1 - len(pprod))
	
	cl = vector(ZZ, c[:n])
	ch = vector(ZZ, c[n:]+[0])
	
	r = (alph * cl) + (lmbd * ch)
	
	return pmns_red_int(r, n, ri_mat, neg_iri_mat, iG_dom, phi)

########################################################################
########################################################################



