# Wrinkled Tension Field

This code implements the paper "Fine Wrinkling on Coarsely-Meshed Thin Shells". (Only tested with Linux)

## Elastic Material Model (benchmarks)

The elastic material model is wrappered from [libshell](https://github.com/evouga/libshell), where we implement a parallel version using [TBB](https://github.com/wjakob/tbb). See the files in "ElasticShell" for details

## Wrinkled Tension Field (WTF) Model

As mentioned in the paper, we seperate the solver into two parts: base mesh extraction using tension field thoery (TFT), and wrinkle variable computation. The former one is implemented inside "ElasticShell". The later is inside "WTFShell".

## Dependencies
- [libigl](https://libigl.github.io/) (used as submodule)
- [Gurobi](https://www.gurobi.com/downloads/?campaignid=193283256&adgroupid=51266130904&creative=419644944624&keyword=gurobi&matchtype=p&gclid=Cj0KCQjwit_8BRCoARIsAIx3Rj6JdxrDRsUeWNRjj8ABmFg40kehVwvHoTsi28UxLeOqa8GhYTQU9usaAucxEALw_wcB) (Free for academic use).
- [NASOQ](https://nasoq.github.io/) (used as submodule, which requires the [MKL](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html))
- [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)
- [TBB](https://github.com/wjakob/tbb)
- [Json](https://github.com/open-source-parsers/jsoncpp) (wrappered in "external/jsoncppWrapper")
- [cppoptlib](https://github.com/PatWie/CppNumericalSolvers/tree/v2/include/cppoptlib) (a wrapper version in inside "external/cppoptlib", where we build our SQP solver with NASOQ)

## To download and build

Compile this project using the standard cmake routine:

    git clone --recurse-submodules https://github.com/csyzzkdcz/WrinkledTensionFields.git
    mkdir build
    cd build
    cmake ..
    make -j4

This procedure will build:
 - Elastic simulator program (benchmarks);
 - WTF similator program, which compute the wrinkle variables (amplitude and frequency) given the TFT base mesh.

## Example Program

After build, there will be two programs (Elatic_viewer_bin, WTF_viewer_bin) inside "ElasticGui" and "WTFGui" respectively.
- To run elastic benchmarks, try with
```bash
   ./SolverApp/ElasticGui/Elastic_viewer_bin ../Sims/disc_elastic/disc_elastic.json
```
  This will generate libigl viewer. Then click "Optimize ElasticShell" to run the simulation process starting from the current configuration. Click "Reset ElasticShell" to reset the current state to initial state (where in the disc case, it is the simulated TFT mesh).
- To run wrinkled tension field model, try with
```bash
   ./SolverApp/WTFGui/WTF_viewer_bin ../Sims/disc/disc_WTF.json
```
To see the upsample wrinkled mesh, check on the "Visualize Wrinkles" button. This will render the upsampled wrinkled mesh in red. 

## Compiling on Windows

Due to poor interoperation of the Eigen library with the MSVC compiler, Release mode compilation of the derivative code on Windows can take forever (over 8 hours). To solve this issue add EIGEN_STRONG_INLINE=inline to your preprocessor macros when building libshell.

## More data

You can find more complicate data https://www.dropbox.com/sh/b33crodfwjbf982/AABhaIjJHhiIK4HeFMedpcgoa?dl=0