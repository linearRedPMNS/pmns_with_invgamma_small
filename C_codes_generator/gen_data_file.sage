
def build_data_file(dir_path, n, phi_log2, rho_nbits, delta, big_int_name, big_int, small_int, convbase_log2, ext_pol):
	
	with open(dir_path+"/pmns_data.h", "w") as f:
		f.write("#ifndef STRUCTS_DATA\n")
		f.write("#define STRUCTS_DATA\n\n")
		
		f.write("#include <stdint.h>\n\n")
		
		if big_int_name != big_int:
			f.write("typedef " + big_int_name + " " + big_int + ";\n\n")
		
		f.write("#define PHI_LOG2 " + str(phi_log2) + "\n")
		f.write("#define POLY_DEG " + str(n-1) + "\n")
		f.write("#define NB_COEFF " + str(n) + "\n")
		f.write("#define NB_ADD_MAX " + str(delta) + "\n\n")
		
		ext_polC = list(ext_pol)
		f.write("#define E_ALPHA " + str(ext_polC[-1]) + "\n")
		f.write("#define E_LAMBDA " + str(-ext_polC[0]) + "\n\n")
		
		f.write("#define RHO_NBITS " + str(rho_nbits) + "	// note: so " + str(rho_nbits +  1) + " bits are used to store each coeficient of an element\n\n")
		
		f.write("#define CONVBASE_LOG2 " + str(convbase_log2) + "\n\n")
		
		f.write("#define CONV_MASK " + str((1<<convbase_log2) - 1) + "UL  // = (1 << CONVBASE_LOG2) - 1, for conversion\n\n")
		
		f.write("\n")
		f.write("//Representations of polynomials Pi, used for conversion into the PMNS\n")
		f.write("//Each Pi is a representation of '(2^CONVBASE_LOG2)^i * alpha^{-1} * phi^2'\n")
		f.write("extern " + str(small_int) + " polys_P[NB_COEFF][NB_COEFF];\n\n")
			
		f.write("#endif\n\n")
		
		

