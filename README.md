# Wrinkled Tension Field

This code implements the paper "Fine Wrinkling on Coarsely-Meshed Thin Shells".

## Elastic Material Model (benchmarks)

The elastic material model is wrappered from [libshell](https://github.com/evouga/libshell), where we implement a parallel version using [TBB](https://github.com/wjakob/tbb). See the files in "ElasticShell" for details

## Wrinkled Tension Field (WTF) Model

As mentioned in the paper, we seperate the solver into two parts: base mesh extraction using tension field thoery (TFT), and wrinkle variable computation. The former one is implemented inside "ElasticShell". The later is inside "WTFShell".


## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This procedure will build:
 - Elastic simulator program (benchmarks);
 - WTF similator program, which compute the wrinkle variables (amplitude and frequency) given the TFT base mesh.

## Dependencies
- [libigl](https://libigl.github.io/)
- [Gurobi](https://www.gurobi.com/downloads/?campaignid=193283256&adgroupid=51266130904&creative=419644944624&keyword=gurobi&matchtype=p&gclid=Cj0KCQjwit_8BRCoARIsAIx3Rj6JdxrDRsUeWNRjj8ABmFg40kehVwvHoTsi28UxLeOqa8GhYTQU9usaAucxEALw_wcB) (Free for academic use).
- [NASOQ](https://nasoq.github.io/)
- [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (with and without [OpenMP](https://computing.llnl.gov/tutorials/openMP/))
- [TBB](https://github.com/wjakob/tbb)
- [Json](https://github.com/nlohmann/json)
- [cppoptlib](https://github.com/PatWie/CppNumericalSolvers/tree/v2/include/cppoptlib) (a wrapper version in inside "external/cppoptlib", where we build our SQP solver with NASOQ)

## Compiling on Windows

Due to poor interoperation of the Eigen library with the MSVC compiler, Release mode compilation of the derivative code on Windows can take forever (over 8 hours). To solve this issue add EIGEN_STRONG_INLINE=inline to your preprocessor macros when building libshell.

## Example Program

After build, there will be two programs (Elatic_viewer_bin, WTF_viewer_bin) inside "ElasticGui" and "WTFGui" respectively.
- To run elastic benchmarks, try with
```bash
   ./Elatic_viewer_bin ../../Sims/disc_elastic/disc_elastic.json
```
  This will generate libigl viewer. Then click "Optimize ElasticShell" to run the simulation process starting from the current configuration. Click "Reset ElasticShell" to reset the current state to initial state (where in the disc case, it is the simulated TFT mesh).
- To run wrinkled tension field model, try with
```bash
   ./WTF_viewer_bin ../../Sims/disc/disc.json
```
To see the upsample wrinkled mesh, check on the "Visualize Wrinkles" button. This will render the upsampled wrinkled mesh in red. 