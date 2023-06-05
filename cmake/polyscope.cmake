include(FetchContent)
FetchContent_Declare(
        polyscope
        GIT_REPOSITORY "https://github.com/nmwsharp/polyscope.git"
        GIT_TAG "4c173af5685aa2a446efe5fb3d87fb7f6d26a274"
)
FetchContent_MakeAvailable(polyscope)