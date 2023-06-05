
if(TARGET igl::core)
return()
endif()

include(FetchContent)
FetchContent_Declare(
libigl
GIT_REPOSITORY https://github.com/libigl/libigl.git
GIT_TAG 7b6cc27
)
FetchContent_MakeAvailable(libigl)