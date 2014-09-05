feelpp_add_application( seme-acoustic 
		SRCS seme-acoustic.cpp 
		CFG seme-acoustic.cfg 
		GEO
		acoustic3d.geo
		acoustic.geo
		acoustic-vecpropagation.geo
		INCLUDE_IN_ALL )

feelpp_add_application( seme-acoustic-fe 
		SRCS seme-acoustic-fe.cpp 
		CFG seme-acoustic-fe.cfg  
		GEO
		acoustic3d.geo
		acoustic.geo
		acoustic-vecpropagation.geo
		INCLUDE_IN_ALL )
