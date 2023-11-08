# just_math

Just Math - A collection of pure math demos.

The goal of Just Math is to provide both a visual example and code demonstrations of specific concepts in computer graphics, mathematics, simulation and AI. This repository consists of 'libmin', a **minimal utility library** for math and graphics, along with **several math examples** for various concepts. 

Copyright 2007-2022 (c) Quanta Sciences, Rama Hoetzlein, <a href="http://ramakarl.com">htpp://ramakarl.com</a>. MIT License.<br>
Contact: ramahoetzlein@gmail.com

## Sample Gallery

<div style="display:flex">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_3ddda.JPG" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_basis.JPG" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_bp.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_cells.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_deform.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_invk.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_quatsquad.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_raycast.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_trajectories.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_wangtiles.jpg" width="200">
<img src="https://github.com/ramakarl/just_math/blob/main/gallery/img_wangtiles3d.jpg" width="200">
</div>

## Just Math Samples

Each sample in Just Math demonstrates a specific concept in code and visually.
The samples provided are briefly described:
- 3DDDA - 3D Differential Analyzer. March through a volume to identify voxels on a line.
- Basis - Orthonormal bases. Transformation from one space to another using bases.
- BeliefProp - Belief propagation. Demonstrated on a 3D grid domain with speedups.
- Cells - Cellular membrane simulation. Simulated with physics using circles for cells.
- Deform - 3D spatial deformations, including bending, twisting and folding.
- InvK - Inverse Kinematics using quaternions. Demo of robot and human arm IK.
- QuatSquad - Quaternion Squad. A C1 continuous method for interpolating orientations.
- QuatTrajectory - Trajectory interpolation of both position and orientation,
using B-Splines, Bezier Curves, and Catmull-Rom splines for position. Slerp or Squad for orientation.
- Raycast - Rendering of a 3D volume with an opacity-based volume integral, on CPU.
- Voxelizer - Voxelization of triangle into a volume, using several methods.
- WangTiles - Sampling of spatial distribution functions with scale invariance.
- WangTiles3D - Alternative demo of Wang Tiles for 3D geometry instancing over a density map landscape.
- Flightsim - *NEW* 2023. Simple force-driven flight simulator.

## How to Build
Platforms:<br>
Win10, Visual Studio 2019 - definitely<br>
Win10/11, VS{other} - probably<br>
Linux - yes, libmin and math_belief_prop only<br>
Mac - unknown<br>
Dependencies: OpenGL only, CUDA is optional (flag at cmake time)<br>

Cmake build options should default to BUILD_OPENGL=ON, BUILD_CUDA=off, BUILD_CONSOLE=off.<br>
Keep these settings. CUDA and/or Console mode are not yet well supported.

**Step 1)** Cmake and build Libmin. <br><br>
Windows: `cmake -S \just_math\libmin -B \build\libmin`<br><br>
Linux: `cmake -DBUILD_OPENGL=OFF`
 (you may need to create a libmin\bin folder and copy liblibmin.so into libmin\bin\libmind.so and libmin\bin\libmin.so)<br><br>

The binary (build) path should be outside of the source \just_math folder.<br>
You must successfully build libmin before proceeding to step 2.<br>

**Step 2)** Cmake and build sample. <br><br>
Windows: `cmake -S \just_math\raycast -B \build\raycast`<br><br>
Linux: `cmake -DBUILD_OPENGL=OFF -DBUILD_CONSOLE=ON -DLIBMIN_ROOT_DIR={abspath}`<br><br>

The binary (build) path should be outside of the source \just_math folder.<br>
Build and run the sample.<br>


## Contributions
I am interested in building a community around simple, well documented, math codes, in pure C/C++ for CPU (no shaders), with interactive graphical demos (not just youtube videos) that are MIT/BSD Licensed. If you have similar interests contact me at: Rama Hoetzlein, ramahoetzlein@gmail.com

## License

MIT License <br>
Copyright 2007-2022 (c) Quanta Sciences, Rama Hoetzlein, ramakarl.com<br>

The Just Math samples are MIT Licensed.<br>

Libmin is MIT Licensed with contributions from other BSD and MIT licensed sources.<br>

Contact: Rama Hoetzlein at ramahoetzlein@gmail.com



