################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.1)
project(WTF_viewer)

# Add your project files
file(GLOB SRCFILES 
*.cpp
../../*.cpp ../../*.h ../../*.hpp
../../SecondFundamentalForm/*.cpp 
../../SecondFundamentalForm/*.h 
../../external/cppoptlib/*.h 
../../external/cppoptlib/linesearch/*.h 
../../external/cppoptlib/solver/*.h 
../../external/eigengurobi/*.h 
../../external/eigengurobi/*.cpp
../../collisionDetectionWrapper/*.h
../../collisionDetectionWrapper/*.cpp
../../WTFShell/*.h
../../WTFShell/*.cpp
../../ElasticShell/*.h
../../ElasticShell/*.cpp)
	  

add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::triangle igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::comiso ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${METIS_LIBRARY} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY} ${TBB_LIBRARIES} json)

# Use C++14
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME}_bin OpenMP::OpenMP_CXX)
endif()