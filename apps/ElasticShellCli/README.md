# Elastic Shell Simulator
This folder constains an executable program for (quasi-)static elatic simulation. The implementation is based on [libshell](https://github.com/evouga/libshell), and adapted for the paper "Fine Wrinkling on Coarsely-Meshed Thin Shells".

# Run
Assume that you are under the folder build, then you can try with:
```
./bin/ElasticShellCli_bin -i ../data/elastic/disk/disk_elastic.json -o ../data/elastic/disk/results
```
The detailed parameters are listed below:

* `-i/--input`: path for input json file (`.json`), which contains the problem setup (see below for details).
* `-o/--output`: the folder for output results (ampliude, frequency and wrinkled mesh).
* `-n/--numIter`: Number of iterations, default is 1000.
* `-g/--gradTol`: The tolerance for gradient norm termination, default is 1e-6.
* `-x/--xTol`: The tolerance of variable update termination, default is 0.
* `-f/--fTol`: The tolerance of function update termination, default is 0.
* `-q/--quiet`: The flag to turn off printing the optimization log, default is false.
* `-p/--randomPerturb`: Add some random perturbation, default is 0 (not add)

## JSON file
* `rest_mesh`: The .obj file for the rest mesh.
* `init_mesh`: The .obj file for the initial simulated mesh.
* `cur_mesh`: The .obj file for the current simulated mesh.
* `obs_mesh`: The .obj file for the obstacle mesh.
* `curedge_DOFs`: The .txt file storing current edge DOFs (see [libshell](https://github.com/evouga/libshell) for details)
* `poisson_ratio`: The Poisson's ratio
* `thickness`: Thickness (m)
* `youngs_modulus`: Young's modulus (Pa)
* `density`: The material density (kg/m^3)
* `clamped_DOFs`: The .dat file stores all the clamped vertices. (vid, x, y, z).
* `collision_eta`: The effective region to consider the collision: d(p1, obs) < eta,
* `collision_penalty`: The collision stiffness
* `pressure`: The pressure stiffness
* `gravity`: Gravity (or constant external force) 
* `stretching_type`: The elastic strain type: "NeoHookean", "tensionField" or "StVK"
* `bending_type`: The bending type: "elastic" or "quadratic"
* `sff_type`: The way to discretize the Second Fundamental Form, see [libshell](https://github.com/evouga/libshell) for details
* `rest_flat`: Whether the rest shape is flat (if this flag turns on, we will set the rest second fundamental form to 0)
* `perturb`: The magnitude to perturb the intial guess
* `frame_frequency`: Save the intermidate results every "frame_frequency" iteration. 
* `max_stepsize`: The maximum step size for line search
* `num_interpolation`: The number of quasi-static step (gradually move the clamped vertices to the target position) 

## Collision
We use a really native collsion model, and may slow down the simulation. For more advanced collsion model, please refer [C-IPC](https://github.com/ipc-sim/Codim-IPC) repo for details.