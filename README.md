# Fine Wrinkling on Coarsely-Meshed Thin Shells

This code implements the paper [Fine Wrinkling on Coarsely-Meshed Thin Shells](https://dl.acm.org/doi/10.1145/3462758). 

## Elastic Material Model (benchmarks)

The elastic material model is wrappered from [libshell](https://github.com/evouga/libshell), where we implement a parallel version using [TBB](https://github.com/wjakob/tbb). See the files in "ElasticShell" for details

## Wrinkled Tension Field (TFW) Model

As mentioned in the paper, we seperate the solver into two parts: base mesh extraction using tension field thoery (TFT), and wrinkle variable computation. The former one is implemented inside "ElasticShell". The later is inside "TFWShell".

## Dependencies
- [Libigl](https://github.com/libigl/libigl.git)
- [Polyscope](https://github.com/nmwsharp/polyscope.git)
- [TBB](https://github.com/wjakob/tbb.git)
- [NASOQ](https://nasoq.github.io/) 
- [Json](https://github.com/nlohmann/json.git) (We include the json.hpp file in the external/json)
- [Suite Sparse](https://people.engr.tamu.edu/davis/suitesparse.html)

All the dependencies are solved by Fetcontent, except Suite Sparse (See below for instruction) and NASOQ (included as a submodule). 


## Build with Suite-Sparse
For linux, you should use 
```
sudo apt-get update -y
sudo apt-get install -y libsuitesparse-dev
```

For macOS, this can be done with [Homebrew](https://brew.sh/):
```
brew install suite-sparse
```

For windows, please follow the guidence provided in [suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows). If you have anaconda installed, this part can be solved by:
```
conda install -c conda-forge suitesparse
```

## To download and build
To build this project, you can use the following commands
```
    git clone --recursive https://github.com/zhenchen-jay/WrinkledTensionFields.git
    cd WrinkledTensionField
    mkdir build
    cd build
    cmake ..
    make -j4
```
**important:** Make sure you have `--recursive` when you clone the repository.

## To run
You can use the following commands to run the program:
```
./bin/TFWShellCli_bin -i ../data/TFW/dress/dress2_TFW.json -o ../data/TFW/dress/results -r -u 3
```
The detailed parameters are listed below:

* `-i/--input`: path for input json file (`.json`), which contains the problem setup (see below for details).
* `-o/--output`: the folder for output results (ampliude, frequency and wrinkled mesh).
* `-r/--reinitialization`: reinitialize the amplitude and frequency based on the current strain.
* `-u/--upsampleTimes`:  the upsampling (subdivision) times to recover the final wrinkled mesh.
* `-n/--numIter`: Number of iterations, default is 1000.
* `-g/--gradTol`: The tolerance for gradient norm termination, default is 1e-6.
* `-x/--xTol`: The tolerance of variable update termination, default is 0.
* `-f/--fTol`: The tolerance of function update termination, default is 0.
* `-q/--quiet`: The flag to turn off printing the optimization log, default is false.

## JSON file
The setup json file has the following format:

* `rest_mesh`: The .obj file for the rest (coarse) mesh.
* `base_mesh`: The .obj file for the base (coarse) mesh.
* `poisson_ratio`: The Poisson's ratio
* `thickness`: Thickness (m)
* `youngs_modulus`: Young's modulus (Pa)
* `amp_path`: The .txt file to store the current/final amplitude file (per vertex value)
* `phi_path`: The .txt file to store the current/final phase file (per vertex value)
* `dphi_path`: The .txt file to store the current/final frequency file (per edge value)
* `clamp`:  The flag to indicate whether we keep the clamped vertices clamped (force the ampltidue to be zero), default is true
* `clamped_DOFs`: The .dat file stores all the clamped vertices. (vid, x, y, z).
* `nasoq_eps`: The tolerance for NASOQ solver, default is 1e-6
* `obs_mesh`: The .obj file for the obstacle file.
* `quadpoints`: The number of quadrature point used to evaluate the integration, default is 3
* `sff_type`: The discrete second fundamental form type, default is midedgeTan (see [libsell](https://github.com/evouga/libshell) for details)
* `rest_flat`: The flag for flat (zero second fundamental form) rest mesh.

## Compiling Issues

* It may take a while to compile on Windows, see [libsell](https://github.com/evouga/libshell) for possible solutions. Beside, you may encounter the miss of `libopenblas.dll` issue. In this case, please copy the `libopenblas.dll` inside the `build/_deps/comiso-src/ext/OpenBLAS-v0.2.14-Win64-int64/bin` and paste it into `build/bin/Release`. We also include pre-compiled executable files under the folder `WindowsExe/Release`.

* When compiling in MacOS with C++17, you may encounter this issue: 
```
build/_deps/comiso-src/ext/gmm-4.2/include/gmm/gmm_domain_decomp.h:84:2: error: ISO C++17 does not allow 'register' storage class specifier [-Wregister]
```
To solve this, please remove the replace `register double` by `double`.  

## Input coase mesh
Although our TFWShellCli_bin can take any coarse mesh as input, as stated in our paper, the tension field results give the best final wrinkled mesh. To turn the any coarse shape into tension field result, you can try with [libsell](https://github.com/evouga/libshell). We also provide an executable (ElasticShellCli_bin). Please refer the README.md file under apps/ElasticShellCli for details.

## GUI
We also offer GUI versions of these two executable programs: `ElaticShellGui_bin` and `TFWShellGui_bin`. The inputs for these GUI versions are identical to those required for the command line counterparts.

## More data

You can find more complicate data https://www.dropbox.com/sh/b33crodfwjbf982/AABhaIjJHhiIK4HeFMedpcgoa?dl=0
