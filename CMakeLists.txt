cmake_minimum_required(VERSION 3.0.0)

# add_compile_options(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# set(CMAKE_CXX_COMPILER "/usr/local/bin/g++11")

project("Thin Shell Benchmarks")
# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# build type 
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


### Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    message("Setting gcc-specific options")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()


  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  #SET(CMAKE_CXX_FLAGS_RELEASE        "-mcpu=apple-m1 -O3 -DNDEBUG") # uncomment if you use M1 cpu
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)

else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# polyscope
message("build with polyscope")
include(polyscope)

#libigl
message("build with libigl")
option(LIBIGL_GLFW                "Build target igl::glfw"                OFF)
option(LIBIGL_PNG                 "Build target igl::png"                 OFF)
option(LIBIGL_IMGUI               "Build target igl::imgui"               OFF)
option(LIBIGL_OPENGL              "Build target igl::opengl"              OFF)
option(LIBIGL_PREDICATES          "Build target igl::predicates"          ON)
option(LIBIGL_COPYLEFT_COMISO     "Build target igl_copyleft::comiso"     ON)
option(LIBIGL_COPYLEFT_CORE       "Build target igl_copyleft::core"       ON)
option(LIBIGL_COPYLEFT_TETGEN     "Build target igl_copyleft::tetgen"     ON)
option(LIBIGL_COPYLEFT_CGAL       "Build target igl_copyleft::cgal"       OFF)
option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" ON)
include(libigl)

#tbb
set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC ON CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC_PROXY ON CACHE BOOL " " FORCE)
set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
option(BUILD_WITH_TBB        "Enable TBB"         ON)
include(tbb)

# NASOQ, refer https://github.com/cheshmi/nasoq-example.git
# set NASOQ properties such BLAS Type, etc
set(NASOQ_BLAS_BACKEND "OpenBLAS" CACHE STRING "BLAS implementation for NASOQ to use")

if (MSVC)
 set(NASOQ_USE_CLAPACK  ON CACHE STRING "Use CLAPACK as the LAPACK implementaiton")
endif()

add_subdirectory(external/nasoq)

if(NOT DEFINED NASOQ_ROOT)
    set(NASOQ_ROOT_DIR "external/nasoq") #assumed the example folder is within nasoq
else()
    set(NASOQ_ROOT_DIR "${NASOQ_ROOT}")
endif()

# Set include and library directories
set(NASOQ_INC ${NASOQ_ROOT_DIR}/include)
set(NASOQ_EIGEN_INC ${NASOQ_ROOT_DIR}/eigen_interface/include)
set(NASOQ_LIB ${NASOQ_ROOT_DIR}/lib/)
include_directories(${NASOQ_INC} ${NASOQ_EIGEN_INC})
link_directories(${NASOQ_LIB})

#suite-sparse
message("build with suite-sparse")
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

#CLI11
include(CLI11)

# nlohmann json
set(JSON_INC external/json)
include_directories(${JSON_INC})

# spectra
include(spectra)

# Binary executables
add_subdirectory(apps)


