
def build_data_file(dir_path, Zp, p, n, phi_log2, rho_nbits, delta, big_int_name, big_int, small_int, convbase_log2, ext_pol, ri_mat, neg_iri_mat, R, iM_dom):
	
	with open(dir_path+"/pmns_data.h", "w") as f:
		f.write("#ifndef STRUCTS_DATA\n")
		f.write("#define STRUCTS_DATA\n\n\n")
		
		if big_int_name != big_int:
			f.write("typedef " + big_int_name + " " + big_int + ";\n\n")
		
		f.write("#define PHI_LOG2 " + str(phi_log2) + "\n")
		f.write("#define POLY_DEG " + str(n-1) + "\n")
		f.write("#define NB_COEFF " + str(n) + "\n")
		f.write("#define NB_ADD_MAX " + str(delta) + "\n\n")
		
		f.write("#define RHO_NBITS " + str(rho_nbits) + "	// note: so " + str(rho_nbits +  1) + " bits are used to store each coeficient of an element\n\n")
		
		f.write("#define CONVBASE_LOG2 " + str(convbase_log2) + "\n\n")
		
		f.write("#define CONV_MASK " + str((1<<convbase_log2) - 1) + "UL  // = (1 << CONVBASE_LOG2) - 1, for conversion\n\n")
		
		f.write("\n\n")
		f.write("//Representations of polynomials Pi, used for conversion into the PMNS\n")
		f.write("//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2\n")
		f.write(generate_conv_polys_P(small_int, p, n, convbase_log2, (1<<phi_log2), ri_mat, neg_iri_mat, iM_dom, Zp, ext_pol, R))
		
		f.write("\n")
		f.write("mpz_t modul_p;\n\n")
		f.write("mpz_t gama_pow[POLY_DEG];\n\n")
			
		f.write("#endif\n\n")
		
		
	


