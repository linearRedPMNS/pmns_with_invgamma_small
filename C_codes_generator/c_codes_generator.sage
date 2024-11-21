from os import makedirs as create_dir
from shutil import rmtree as delete_code
from shutil import copy as cpfile
load("c_codes_gen_utils.sage")
load("gen_data_file.sage")
load("gen_arith_ops_file.sage")
load("gen_useful_functs_file.sage")

########################################################################

#~ NOTE: polynomial form is P(X) = a_0 + ... + a_n.X^n = [a_0, ..., a_n]

########################################################################

#~ outputs: [small_int_name, unsigned_small_int_name, big_int_name, big_int_new_name]
def choose_target_archi_info(phi_log2):
	if phi_log2 == 8:
		return ["int8_t", "uint8_t", "int16_t", "int16_t"]
	elif phi_log2 == 16:
		return ["int16_t", "uint16_t", "int32_t", "int32_t"]
	elif phi_log2 == 32:
		return ["int32_t", "uint32_t", "int64_t", "int64_t"]
	elif phi_log2 == 64:
		return ["int64_t", "uint64_t", "__int128", "int128_t"]
	else:
		print("ERROR: unspecified info for an architecture with 'phi_log2="+str(phi_log2)+".") 
		print("Please specify it in the function 'choose_target_archi_info'.") 
		print("No code has been generated.")
	return -1

########################################################################

def compute_mat_G(n, t, E_alpha, E_lambda):
	mat = []
	tmp = [-1, t] + [0]*(n-2)
	mat.append(tmp)
	for i in range(n-2):
		tmp = [0] + tmp[:-1]
		mat.append(tmp)
	if t%E_alpha == 0:
		s = t//E_alpha
		tmp = [s*E_lambda] + [0]*(n-2) + [-1]
	else:
		tmp = [t*E_lambda] + [0]*(n-2) + [-E_alpha]
	mat.append(tmp)
	return matrix(ZZ, mat)

def compute_mat_iG(G, iG_dom):
	iG = matrix(iG_dom, G)
	iG = -(iG.inverse())
	return iG

def compute_rho_rhoLog2(t, E_alpha, E_lambda):
	if t%E_alpha == 0:  #a case mentioned in a remark of the paper
		s = t//E_alpha
		b = abs(s*E_lambda)
		rho = max(b, abs(t)) 
	else:
		b = abs(t*E_lambda) + 1
		c = abs(t) + abs(E_alpha)
		rho = max(b, c) - 1
	return (rho, rho.nbits())

def compute_alpha_lambda(ext_pol):
	ext_polC = list(ext_pol)
	E_alpha = ext_polC[-1]
	E_lambda = -ext_polC[0]
	return (E_alpha, E_lambda)

def is_DoubleSparse(t, phi_log2):
	a = ceil(phi_log2/2)
	b = 1 << a
	return (t%b == 0)

########################################################################

#~ IMPORTANT: 
# -'pmns_data' must contain in order: [p.nbits(), p, n, str(E), t, phi_log2, delta]
#Data struct: [p.nbits(), p, n, str(E), t, phi_log2, delta]
# - assumes that E(X) = aX^n - b
def build_pmns_c_codes(pmns_data, p_num=0, pmns_num=0):
		
	p_size = pmns_data[0]
	p = pmns_data[1]
	n = pmns_data[2]
	str_E = pmns_data[3]
	pmns_t = pmns_data[4]
	phi_log2 = pmns_data[5]
	delta = pmns_data[6]
	
	# ~ ----------------------------------------------------------------
	
	target_archi_info = choose_target_archi_info(phi_log2)
	if target_archi_info == -1:
		return -1
	
	small_int = target_archi_info[0]
	unsigned_small_int = target_archi_info[1]
	big_int_name = target_archi_info[2]
	big_int = target_archi_info[3]
	
	# ~ ----------------------------------------------------------------
	
	F = GF(p); R.<X> = ZZ[]
	
	mont_phi = 1 << phi_log2; iG_dom = ZZ.quo(mont_phi)
	
	gmm = (1/pmns_t)%p
	
	ext_pol = R(str_E)
	
	(E_alpha, E_lambda) = compute_alpha_lambda(ext_pol)
	
	(rho, rho_nbits) = compute_rho_rhoLog2(pmns_t, E_alpha, E_lambda)
	
	red_int_mat = compute_mat_G(n, pmns_t, E_alpha, E_lambda)
	
	neginv_red_int_mat = compute_mat_iG(red_int_mat, iG_dom)
	
	convbase_log2 = ceil(log(p,2)/n)
	
	double_sparse = is_DoubleSparse(pmns_t, phi_log2)
	
	# ~ ----------------------------------------------------------------
	
	if double_sparse:
		dir_name = 'DS_'
	else:
		dir_name = 'LR_'
	
	dir_name += 'p' + str(p.nbits()) + '_n'  + str(n) + '_rho'  + str(rho_nbits) + '_phi'  + str(phi_log2) + '_d'  + str(delta)
	if p_num != 0:
		dir_name += '_' + str(p_num)
	if pmns_num != 0:
		dir_name += '_' + str(pmns_num)
	
	dir_path = "c_codes/" + dir_name
	
	try:
		create_dir(dir_path)
	except OSError:  # if this directory already exist
		delete_code(dir_path)
		create_dir(dir_path)
	
	cpfile("constant_files/simple_main.c", dir_path+"/main.c")
	cpfile("constant_files/Makefile", dir_path)
	cpfile("constant_files/gmp_stuff.c", dir_path)
	
	# ~ cpfile("constant_files/main__check_ops.c", dir_path+"/main__check_ops.c")
	# ~ cpfile("constant_files/main__for_nb_cycles.c", dir_path+"/main__for_nb_cycles.c")
	# ~ cpfile("constant_files/gmp_stuff.c", dir_path)
	# ~ cpfile("constant_files/intel_measurement_stuff.c", dir_path)
	
	# ~ ----------------------------------------------------------------
	
	
	build_data_file(dir_path, n, phi_log2, rho_nbits, delta, big_int_name, big_int, small_int, convbase_log2, ext_pol)


	build_arith_ops_h_file(dir_path, small_int, big_int, unsigned_small_int)
	
	build_arith_ops_c_file(dir_path, n, phi_log2, ext_pol, small_int, unsigned_small_int, big_int, red_int_mat, neginv_red_int_mat, double_sparse)
	

	build_useful_functs_h_file(dir_path, small_int)
	
	build_useful_functs_c_file(dir_path, small_int, big_int, p, n, gmm, convbase_log2, phi_log2, red_int_mat, neginv_red_int_mat, iG_dom, F, ext_pol, R)
	
	return;

########################################################################
########################################################################


