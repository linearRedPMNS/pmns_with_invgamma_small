
def build_arith_ops_h_file(dir_path, small_int, big_int, unsigned_small_int):
	with open(dir_path+"/pmns_arith_ops.h", "w") as f:
		
		f.write("#ifndef POLY_MULT_ADD\n")
		f.write("#define POLY_MULT_ADD\n\n\n")
		
		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n\n")
		
		f.write("void double_add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void double_sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar);\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos);\n\n")
		
		f.write("void add_lpoly(" + big_int + " *rop, " + big_int + " *pa, " + big_int + " *pb);\n\n")
		f.write("void scalar_mult_lpoly(" + big_int + " *rop, " + small_int + " *op, " + unsigned_small_int + " scalar);\n\n")
		
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa);\n\n")	
		
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("void from_mont_domain(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write(small_int + " equality_check_polys(" + small_int + " *pA, " + small_int + " *pB);\n\n")
		
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op);\n\n")
		
		f.write("#endif\n\n")


def build_arith_ops_c_file(dir_path, n, phi_log2, ext_pol, small_int, unsigned_small_int, big_int, red_int_mat, neginv_red_int_mat, double_sparse):
	with open(dir_path+"/pmns_arith_ops.c", "w") as f:
		
		f.write("#include \"pmns_arith_ops.h\"\n\n\n")
		
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + pb[j];\n")
		f.write("}\n\n")

		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] - pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ --------------------------------------------------------------------\n\n")
	
		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = -op[j];\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << nb_pos;\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << 1;\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = scalar * op[j];\n")
		f.write("}\n\n")
		
		f.write("//~ --------------------------------------------------------------------\n\n")
		
		f.write("//~ computes : pa + 2.pb\n")
		f.write("void double_add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + 2*pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ computes : pa - 2.pb\n")
		f.write("void double_sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] - 2*pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ --------------------------------------------------------------------\n\n")
		
		f.write("void add_lpoly(" + big_int + " *rop, " + big_int + " *pa, " + big_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void scalar_mult_lpoly(" + big_int + " *rop, " + small_int + " *op, " + unsigned_small_int + " scalar){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = (" + big_int + ")op[j] * scalar;\n")
		f.write("}\n\n")
			
		f.write("//~ --------------------------------------------------------------------\n\n")
			
		f.write("//~ Computes: pa*pb mod(E)\n")
		f.write("//~ Note: E(X) = " + str(ext_pol) + "\n")
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n\n")
		f.write("	" + big_int + " tmp_prod_result[NB_COEFF];\n")
		f.write(build_prod_code(n, list(ext_pol)[-1], -list(ext_pol)[0], big_int) + "\n")
		f.write("	internal_reduction(rop, tmp_prod_result);\n")
		f.write("}\n\n")
		
		f.write("//~ performs the internal reduction on 'op' and PUTS the result in 'rop'\n")
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op){\n\n")
		f.write(build_red_int_code(small_int, big_int, n, phi_log2, red_int_mat, neginv_red_int_mat, double_sparse))
		f.write("}\n\n")
		
		f.write("//~ --------------------------------------------------------------------\n\n")
		
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
		
	return;

# ~ --------------------------------------------------------------------






































