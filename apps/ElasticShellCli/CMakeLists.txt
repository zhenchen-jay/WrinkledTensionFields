################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.1)
project(ElasticShellCli)

# Add your project files
# Create an executable
file(GLOB SRCFILES 
*.cpp
../../*.cpp ../../*.h ../../*.hpp
../../SecondFundamentalForm/*.cpp 
../../SecondFundamentalForm/*.h 
../../ElasticShell/*.h
../../ElasticShell/*.cpp
../../TFWShell/*.h
../../TFWShell/*.cpp
../../MeshLib/*.h
../../MeshLib/*.cpp
../../Optimization/*.h
../../Optimization/*.cpp
../../OtherEnergies/*.h
../../OtherEnergies/*.cpp
../../MIPWrapper/*.h
../../MIPWrapper/*.cpp
../../Collision/*.h
../../Collision/*.cpp
../../ComplexLoop/*.h
../../ComplexLoop/*.cpp
../../external/SecStencils/*.h
../../external/SecStencils/*.cpp
)
    
add_executable(
        ${PROJECT_NAME}_bin
        ${SRCFILES}
        )

set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD_REQUIRED ON)

# Link settings
target_include_directories(${PROJECT_NAME}_bin PUBLIC ${spectra_SOURCE_DIR}/include/Spectra)
target_link_libraries( ${PROJECT_NAME}_bin
        nasoq
        TBB::tbb 
        tbb_static 
        tbbmalloc_static
        igl::core 
        igl_restricted::triangle 
        igl_copyleft::comiso 
        ${SUITESPARSE_LIBRARIES}
        CLI11::CLI11)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME}_bin OpenMP::OpenMP_CXX)
endif()