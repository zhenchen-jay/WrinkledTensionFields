# - Try to find the NASOQ library
# Once done this will define
#
#  NASOQ_FOUND - system has LIBIGL
#  NASOQ_INCLUDE_DIRS - the NASOQ include directories
find_path(NASOQ_ROOT_DIR QP/nasoq.h
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/../tools/nasoq
        ${CMAKE_SOURCE_DIR}/nasoq
        ${CMAKE_SOURCE_DIR}/../nasoq
        ${CMAKE_SOURCE_DIR}/../../nasoq
        ${CMAKE_SOURCE_DIR}/../external/nasoq
        ${CMAKE_SOURCE_DIR}/external/nasoq
        /usr
        /usr/local
        /usr/local/igl/nasoq
)

if(NASOQ_ROOT_DIR)
	set(NASOQ_FOUND TRUE)
	list(APPEND NASOQ_INCLUDE_DIRS
			${NASOQ_ROOT_DIR}/symbolic/
			${NASOQ_ROOT_DIR}/common/
			${NASOQ_ROOT_DIR}/ldl/
			${NASOQ_ROOT_DIR}/matrixMatrix/
			${NASOQ_ROOT_DIR}/matrixVector/
			${NASOQ_ROOT_DIR}/linear_solver/
			${NASOQ_ROOT_DIR}/gmres/
			${NASOQ_ROOT_DIR}/QP/
			${NASOQ_ROOT_DIR}/triangularSolve/
            ${NASOQ_ROOT_DIR}/eigen_interface/
	)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NASOQ
    "\nNASOQ not found"
    NASOQ_INCLUDE_DIRS)
mark_as_advanced(NASOQ_INCLUDE_DIR)
