cmake_minimum_required(VERSION 3.1)
project(shellbenchmark)

set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_LINKER_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
  set(MSVC_RUNTIME "dynamic")
  set(BLA_STATIC ON)
endif()

set(LIBIGL_EIGEN_VERSION 3.3.7 CACHE STRING "Eigen version")

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)


# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL REQUIRED)
message("found libigl at ${LIBIGL_INCLUDE_DIR}")

option(NASOQ "Use NASOQ" ON)
if(NASOQ)
  find_package(MKL REQUIRED)
  find_package(NASOQ REQUIRED)
  include_directories(${MKL_INCLUDE_DIR} ${NASOQ_INCLUDE_DIRS})  
  message("found MKL at ${MKL_INCLUDE_DIR}")
  message("found NASOQ at ${NASOQ_INCLUDE_DIRS}")
endif()


# Use TBB
# MUST COME BEFFORE SUITESPARSE!

option(WTF_WITH_TBB        "Enable TBB"         ON)

if(WTF_WITH_TBB)
  message("-- Building with TBB.")
  set(TBB_FIND_RELEASE_ONLY TRUE)
  if(MSVC)
	  find_package(TBB REQUIRED)
	  include_directories(${TBB_INCLUDE_DIRS})    
  else()
    find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
  endif()
endif()
if(MSVC)
  add_compile_definitions(__TBB_NO_IMPLICIT_LINKAGE=1)
endif()
  

find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS} ${METIS_INCLUDE_DIR})

# use gurobi
find_package(GUROBI)
if(NOT GUROBI_FOUND)
  message("Unable to find Gurobi. Make sure that the env GUROBI_HOME is set")
  set (GUROBI_INCLUDE_DIRS " C:/gurobi901/win64/include")
  set (GUROBI_CXX_LIBRARY "C:/gurobi901/win64/lib/gurobi_c++md2019.lib")
  set (GUROBI_LIBRARY "C:/gurobi901/win64/lib/gurobi90.lib")
else()
  message("Found Gurobi in : " ${GUROBI_INCLUDE_DIRS})
  message("Found Gurobi cxx lib in : " ${GUROBI_CXX_LIBRARY})
  message("Found Gurobi lib in : " ${GUROBI_LIBRARY})
  include_directories(${GUROBI_INCLUDE_DIRS})
endif()

# use json, TODO: switch to use json.hpp provided by nlohmann
add_subdirectory(external/jsoncppWrapper)
include_directories(external/jsoncppWrapper/include)

# Binary executables
add_subdirectory(SolverApp)