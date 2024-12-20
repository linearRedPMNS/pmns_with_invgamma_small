def build_useful_functs_h_file(dir_path, small_int):
	with open(dir_path+"/pmns_useful_functs.h", "w") as f:
		
		f.write("#ifndef USEFUL_FUNCTS\n")
		f.write("#define USEFUL_FUNCTS\n\n")
		
		f.write("#include <stdint.h>\n")
		f.write("#include <gmp.h>\n\n")
		f.write("#include \"pmns_data.h\"\n")
		f.write("#include \"pmns_arith_ops.h\"\n\n")
		
		f.write("void init_data();\n\n")
		f.write("void free_data();\n\n")
		
		f.write("void from_int_to_pmns(" + small_int + " *rop, mpz_t op);\n\n")
		
		f.write("void from_pmns_to_int(mpz_t rop, " + small_int + " *op);\n\n")
		
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("void from_mont_domain(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("void copy_poly(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("void print_element(" + small_int + " *poly);\n\n")
		
		f.write("//" + small_int + " equality_check_polys(" + small_int + " *pA, " + small_int + " *pB);\n\n")
		
		f.write("int cmp_poly_evals(" + small_int + " *pa, " + small_int + " *pb);\n\n")
		
		f.write("#endif\n\n")


def build_useful_functs_c_file(dir_path, small_int, big_int, p, n, gmm, convbase_log2, phi_log2, red_int_mat, neginv_red_int_mat, iG_dom, Zp, ext_pol, R):
	with open(dir_path+"/pmns_useful_functs.c", "w") as f:
		
		f.write("#include <stdio.h>\n")
		f.write("#include <stdlib.h>\n")
		f.write("#include <stdint.h>\n")
		f.write("#include <gmp.h>\n\n")
		f.write("#include \"pmns_data.h\"\n")
		f.write("#include \"pmns_arith_ops.h\"\n")
		f.write("#include \"pmns_useful_functs.h\"\n\n")
		
		f.write("//~ --------------------------------------------------------------------\n\n")
		
		f.write("mpz_t modul_p;\n")
		f.write("mpz_t gama_pow[POLY_DEG];\n\n")
		
		f.write("//Representations of polynomials Pi, used for conversion into the PMNS\n")
		f.write("//Each Pi is a representation of '(2^CONVBASE_LOG2)^i * alpha^{-1} * phi^2'\n")
		f.write(generate_conv_polys_P(small_int, p, n, convbase_log2, (1<<phi_log2), red_int_mat, neginv_red_int_mat, iG_dom, Zp, ext_pol, R))
		
		f.write("//~ --------------------------------------------------------------------\n\n")
		
		f.write("//~ Assumes allocation already done for 'rop'.\n")
		f.write("//~ IMPORTANT : convertion to montgomery domain will be done here, thanks to conv polys 'polys_P'\n")
		f.write("//~ IMPORTANT : also assumes that the case 'alpha > 1' is already taken into account in conv polys 'polys_P'\n")
		f.write("void from_int_to_pmns(" + small_int + " *rop, mpz_t op){\n\n")
		f.write("	int i;\n")
		f.write("	mpz_t tmp;\n")
		f.write("	" + big_int + " tmp_poly[NB_COEFF];\n")
		f.write("	" + big_int + " tmp_sum[NB_COEFF];\n\n")
		f.write("	mpz_init_set(tmp, op);\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++){\n")
		f.write("		rop[i] = 0;\n")
		f.write("		tmp_sum[i] = 0;\n")
		f.write("	}\n\n")
		f.write("	if(tmp->_mp_size == 0)\n")
		f.write("		return;\n\n")
		f.write("	i = 0;\n")
		f.write("	while(tmp->_mp_size && (i < NB_COEFF)){\n")
		f.write("		scalar_mult_lpoly(tmp_poly, polys_P[i++], (tmp->_mp_d[0]) & CONV_MASK);\n")
		f.write("		add_lpoly(tmp_sum, tmp_sum, tmp_poly);\n\n")
		f.write("		mpz_tdiv_q_2exp (tmp, tmp, CONVBASE_LOG2);\n")
		f.write("	}\n\n")
		f.write("	internal_reduction(rop, tmp_sum);\n\n")
		f.write("	mpz_clear(tmp);\n")
		f.write("}\n\n")
		
		f.write("//~ Assumes \"rop\" already initialized.\n")
		f.write("//~ IMPORTANT : convertion from montgomery domain will be done here.\n")
		f.write("//~ IMPORTANT : the case 'alpha > 1' is also taken into account here\n")
		f.write("void from_pmns_to_int(mpz_t rop, " + small_int + " *op){\n")
		f.write("	int i;\n")
		f.write("	mpz_t tmp_sum;\n")
		f.write("	" + small_int + " tmp_conv[NB_COEFF];\n\n")
		f.write("	mpz_init(tmp_sum);\n\n")
		f.write("	//~ convertion out of mont domain\n")
		f.write("	from_mont_domain(tmp_conv, op);\n\n")
		f.write("	mpz_set_si(rop, tmp_conv[0]);\n")
		f.write("	for(i=0; i<POLY_DEG; i++){\n")
		f.write("		mpz_mul_si(tmp_sum, gama_pow[i], tmp_conv[i+1]);\n")
		f.write("		mpz_add(rop, rop, tmp_sum);\n")
		f.write("	}\n")
		f.write("	mpz_mul_si(rop, rop, E_ALPHA);	// we deal with 'alpha' here\n")  
		f.write("	mpz_mod (rop, rop, modul_p);\n\n")
		f.write("	mpz_clear(tmp_sum);\n")
		f.write("}\n\n")
	
		f.write("//~ ----------------------------------------------------------------------------------------\n\n")
		
		f.write("void init_data(){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_init (gama_pow[i]);\n")
		f.write("	mpz_init (modul_p);\n")
		f.write("	mpz_set_str (modul_p, \"" + str(p) + "\", 10);\n")
		f.write("	mpz_set_str (gama_pow[0], \"" + str(gmm) + "\", 10);\n")
		f.write("	for(i=1; i<POLY_DEG; i++){\n")
		f.write("		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);\n")
		f.write("		mpz_mod (gama_pow[i], gama_pow[i], modul_p);\n")
		f.write("	}\n")
		f.write("}\n\n")
		
		f.write("void free_data(){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_clear (gama_pow[i]);\n")
		f.write("	mpz_clear (modul_p);\n")
		f.write("}\n\n")
		
		f.write("//~ ----------------------------------------------------------------------------------------\n\n")
		
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op){\n\n")
		f.write("	int i;\n")
		f.write("	" + big_int + " tmp[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++){\n")
		f.write("		tmp[i] = (" + big_int + ") op[i];\n")
		f.write("	}\n")
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n\n")
		f.write("	mult_mod_poly(rop, rop, polys_P[0]);\n")
		f.write("}\n\n")
		
		f.write("//~ computes : op/phi\n")
		f.write("void from_mont_domain(" + small_int + " *rop, " + small_int + " *op){\n\n")
		f.write("	int i;\n")
		f.write("	" + big_int + " tmp[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++){\n")
		f.write("		tmp[i] = (" + big_int + ") op[i];\n")
		f.write("	}\n")
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write("}\n\n")
		
		f.write("//~ ----------------------------------------------------------------------------------------\n\n")
		
		f.write("void copy_poly(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		rop[i] = op[i];\n")
		f.write("}\n\n")

		f.write("void print_element(" + small_int + " *poly){\n")
		f.write("	int i;\n")
		f.write("	printf(\"[\");\n")
		f.write("	for (i=0; i<POLY_DEG; i++)\n")
		f.write("		printf(\"%2ld, \", poly[i]);\n")
		f.write("	printf(\"%2ld]\", poly[POLY_DEG]);\n")
		f.write("}\n\n")
		
		f.write("//~ ----------------------------------------------------------------------------------------\n\n")
		
		f.write("//~ return a positive value if pa > pb, zero if pa = pb, or a negative value if pa < pb.\n")
		f.write("//~ Important : evaluation is done using the corresponding integers modulo 'p'.\n")
		f.write("int cmp_poly_evals(" + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int rep;\n")
		f.write("	mpz_t a, b;\n")
		f.write("	mpz_inits (a, b, NULL);\n")
		f.write("	from_pmns_to_int(a, pa);\n")
		f.write("	from_pmns_to_int(b, pb);\n")
		f.write("	rep = mpz_cmp (a, b);\n")
		f.write("	mpz_clears (a, b, NULL);\n")
		f.write("	return rep;\n")
		f.write("}\n\n")

















