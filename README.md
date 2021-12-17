

# Efficient Differentiable Simulation of Articulated Bodies

*Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin*

<!-- [[Project]](https://gamma.umd.edu/researchdirections/mlphysics/diffsim/) -->
[[Paper]](http://proceedings.mlr.press/v139/qiao21a/qiao21a.pdf)
[[Video]](https://icml.cc/virtual/2021/poster/9049)
[[Slides]](https://icml.cc/media/icml-2021/Slides/9049.pdf)
[[Code]](https://github.com/YilingQiao/diffarticulated)

## Setup
0. This project is still work-in-progress. It can be built with gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04), cmake 3.17.3, and CUDA 10.2 (only required for MBPO).
1. Clone this repo and setup the python environment.
```bash
git clone git@github.com:YilingQiao/diffarticulated.git
cd diffarticulated
conda env create -f gpu-env.yml
conda activate diffarti_36
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
```

2. Compile this library.
```bash
mkdir build
cd build
cmake ..
make -j
cd ..
python setup.py install
```
If you do not have the `uuid` lib, you may install it by `apt-get install uuid-dev`. We will add `uuid` to the `third_party` directory shortly.

3. Run the examples
## Examples
### Start visualization
We use `meshcat` to visualize the simulation. Please run the following command to start a session ***before*** runing the simulation.
```bash
meshcat-server --open
```
### Run the simulation code
An exmaple of solving inverse physics using the differentiable physics is `python/exp_throw.py`. Please run
```bash
python python/exp_throw.py
```
and it will optimize the control input in 60 episodes.
### Simulation results
1. When you want to export the simulation results in one frame, you can call `world.sync_visual_meshcat(sim_step)` in the python code. An example can be found in `python/exp_throw.py`
2. By default, simulation records are stored in `meshcat/his/`. This path can be customized by modifying  the `sync_visual_transforms()` function in `examples/meshcat_urdf_visualizer.h`.
3.  To convert the simulation records into obj meshes, please run 
```bash
cd python/results
cp -r ../../meshcat/* ./
python convert_history.py
cd ../..
```
Before generating new records for conversion, please clean the cache in `meshcat/his` and `meshcat/obj`

### Learn to throw a marble
To run the demo of throwing a marble,
```bash
python python/exp_throw.py
```
<div align="center">
<img width="300px" src="https://github.com/YilingQiao/linkfiles/raw/master/icml21/throw.gif"> 
</div>


### Learn to hit a golf ball
```bash
python python/exp_punch.py
```
<div align="center">
<img width="300px" src="https://github.com/YilingQiao/linkfiles/raw/master/icml21/punch.gif"> 
</div>


### Learn the frictional coefficient
```bash
python python/exp_car.py
```
<div align="center">
<img width="300px" src="https://github.com/YilingQiao/linkfiles/raw/master/icml21/car.gif"> 
</div>

## Enhance RL with differentiable physics
We place the RL code in another seperate [repo](https://github.com/YilingQiao/diffarti_mbpo).  
1. Clone the submodule for RL.
```bash
git submodule init
git submodule update
```
2. Install the packages for MBPO.
```bash
cd diffarti_mbpo/
git checkout master
pip install -e viskit
pip install -e .
```
3. Run the experiments.
### Policy Enhancement
In `diffarti_mbpo/`, run the command
```bash
chmod +x ./run_mbpo.sh
./run_mbpo.sh 8_7 pendulumours
```
<div align="center">
<img width="300px" src="https://github.com/YilingQiao/linkfiles/raw/master/icml21/pen.gif"> 
</div>

### Sample Enhancement
In `diffarti_mbpo/`, run the command
```bash
./run_mbpo.sh final_3 antours
```
<div align="center">
<img width="300px" src="https://github.com/YilingQiao/linkfiles/raw/master/icml21/ant.gif"> 
</div>

## Our Related Repos
Differentiable Soft Body Dynamics [Code](https://github.com/YilingQiao/diff_fem) [Paper](http://vladlen.info/publications/differentiable-simulation-soft-multi-body-systems/)
*Differentiable Simulation of Soft Multi-body Systems. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (Neurips 2021)*

Differentiable Articulated Body Dynamics [Code](https://github.com/YilingQiao/diffarticulated) [Paper](https://arxiv.org/abs/2109.07719)
*Efficient Differentiable Simulation of Articulated Bodies. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (ICML 2021)*

Differentiable Dynamics for Rigid Body and Cloth Coupling [Code](https://github.com/YilingQiao/diffsim) [Paper](https://arxiv.org/abs/2007.02168)
*Scalable Differentiable Physics for Learning and Control. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (ICML 2020)*

Differentiable Cloth Dynamics [Code](https://github.com/williamljb/DifferentiableCloth) [Paper](https://www.cs.umd.edu/~liangjb/docs/NIPS2019.pdf)
*Differentiable Cloth Simulation for Inverse Problems. Junbang Liang, Ming C. Lin, Vladlen Koltun. (NeurIPS 2019)*

## Acknowledgments
Thanks for the great open-source project [tiny-differentiable-simulatior](https://github.com/google-research/tiny-differentiable-simulator). This repo is derived from tinydiffsim. You might want to try the original tinydiffsim for its templatized auto-differentiation.

## Bibtex
```
@inproceedings{Qiao2021Efficient,
author  = {Qiao, Yi-Ling and Liang, Junbang and Koltun, Vladlen and Lin, Ming C.},
title  = {Efficient Differentiable Simulation of Articulated Bodies},
booktitle = {ICML},
year  = {2021},
}
```
