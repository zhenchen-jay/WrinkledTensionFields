# SPQR lib usually requires linking to a blas and lapack library.
# It is up to the user of this module to find a BLAS and link to it.

# SPQR lib requires Cholmod, colamd and amd as well. 
# FindCholmod.cmake can be used to find those packages before finding spqr

if (SPQR_INCLUDES AND SPQR_LIBRARIES)
  set(SPQR_FIND_QUIETLY TRUE)
endif (SPQR_INCLUDES AND SPQR_LIBRARIES)

find_path(SPQR_INCLUDES
  NAMES
  SuiteSparseQR.hpp
  PATHS
  $ENV{SPQRDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  include
  suitesparse
  include/suitesparse
  ufsparse
)

find_library(SPQR_LIBRARIES 
  NAMES
  spqr 
  PATHS
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )

if(SPQR_LIBRARIES)

  find_library(SUITESPARSE_LIBRARY 
  NAMES
  SuiteSparse suitesparseconfig
  PATHS $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if (SUITESPARSE_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${SUITESPARSE_LIBRARY})
  endif()

  find_library(CHOLMOD_LIBRARY 
  NAMES
  cholmod 
  PATHS 
  $ENV{SPQRDIR} $ENV{CHOLMOD_LIBDIR} $ENV{CHOLMODDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${CHOLMOD_LIBRARY})
  endif()
  
  find_library(METIS_LIBRARY 
  NAMES
  metis
  metis.so.5 
  PATHS 
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${METIS_LIBRARY})
  endif()
  
  find_library(COLAMD_LIBRARY 
  NAMES
  colamd 
  PATHS 
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${COLAMD_LIBRARY})
  endif()
  
  find_library(CCOLAMD_LIBRARY 
  NAMES
  ccolamd 
  PATHS 
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${CCOLAMD_LIBRARY})
  endif()
  
  find_library(AMD_LIBRARY 
  NAMES
  amd 
  PATHS 
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${AMD_LIBRARY})
  endif()
  
  find_library(CAMD_LIBRARY 
  NAMES
  camd 
  PATHS 
  $ENV{SPQRDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
  )
  if(CHOLMOD_LIBRARY)
    set(SPQR_LIBRARIES ${SPQR_LIBRARIES} ${CAMD_LIBRARY})
  endif()
  
endif(SPQR_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPQR DEFAULT_MSG SPQR_INCLUDES SPQR_LIBRARIES)

mark_as_advanced(SPQR_INCLUDES SPQR_LIBRARIES)