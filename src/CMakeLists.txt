
include_directories(./)

#############################################################
# SEME acoustic applications
#############################################################
feelpp_add_application( seme-acoustic-fd 
		SRCS seme-acoustic-fd.cpp 
		CFG seme-acoustic-fd.cfg 
		GEO
		acoustic3d.geo
		acoustic-vecpropagation.geo
		INCLUDE_IN_ALL )
feelpp_add_application( seme-acoustic-fe 
		SRCS seme-acoustic-fe.cpp 
		CFG seme-acoustic-fe.cfg  
		GEO
		acoustic3d.geo
		acoustic-vecpropagation.geo
		INCLUDE_IN_ALL )



#############################################################
# acoustic model library
#############################################################
add_library( cerema_acoustic_model
  SHARED
  acousticmodel.cpp
  )
target_link_libraries( cerema_acoustic_model
  feelpp ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES}
  )
#############################################################
# acoustic applications
#############################################################
feelpp_add_application( cerema_acoustic_simulation
  SRCS main_acoustic_simulation.cpp
  LINK_LIBRARIES cerema_acoustic_model
  INCLUDE_IN_ALL
  )



#############################################################
# transport applications
#############################################################
feelpp_add_application( transport-equation
  SRCS transport-equation.cpp
  )
