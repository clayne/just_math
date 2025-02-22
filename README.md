# just_math

Just Math - A collection of pure math demos.

The goal of Just Math is to provide both a visual example and code demonstrations of specific concepts in computer graphics, mathematics, simulation and AI. 

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
- ANN - Artificial Neural Network learning the sine function.
- Basis - Orthonormal bases. Transformation from one space to another using bases.
- Bilinear Patch - Fast raytracing of a bilinear patch (curved quad defined by 4x points).
- Cells - Cellular membrane simulation. Simulated with physics using circles for cells.
- Deform - 3D spatial deformations, including bending, twisting and folding.
- Gcoder - Generation of CNC toolpath g-code from depth images.
- InvK - Inverse Kinematics using quaternions. Demo of robot and human arm IK.
- Obj Materials - Reads and renders 3D meshes (.obj) with their material (.mtl) definitions.
- QuatSquad - Quaternion Squad. A C1 continuous method for interpolating orientations.
- QuatTrajectory - Trajectory interpolation of both position and orientation,
using B-Splines, Bezier Curves, and Catmull-Rom splines for position. Slerp or Squad for orientation.
- Raycast - Rendering of a 3D volume with an opacity-based volume integral, on CPU.
- Voxelizer - Voxelization of triangle into a volume, using several methods.
- WangTiles - Sampling of spatial distribution functions with scale invariance.
- WangTiles3D - Alternative demo of Wang Tiles for 3D geometry instancing over a density map landscape.

## Breakout Model Synthesis (BMS) - math_belief_prop
Abram Connelly and Rama Hoetzlein, 2023<br>
<sub>Background: Breakout Model Synthesis (BMS) was a collaboration by Abram Connelly and Rama Hoetzlein in 2023 to explore the constraint based tiling generation (CBTG) problem. We found a shared interested when Abram started exploring the Belief Propagation (BP) algorithm in javascript and introduced me to it. I (Rama) developed a C++ port of BP and WFC (Wave Function Collapse) so that we could work together in a common framework - with the first commit to math_belief_prop in Sept 2022. After we became comfortable with data structures, optimizations and linux/win builds, we began to explore the problem further. Abram developed test 2D tilesets for pacman (pm), supermario (smb), roads and pillmortal (pm) in Feb & March 2023, while Rama tried a wavefront approach to BP which turned out to be no better than WFC. At the end of March 2023, we realized that BP and BP-W were very slow and memory intensive compared to WFC, so much so that it didn't justify the marginal improvement in success rate. I tried some random experiments with LSA, localized simulated annealing, but this was also slower and had more difficulty converging than WFC.	Since it was fast, one could just run WFC multiple times to improve the success rate. Only MMS (Merrill's Model Synthesis) could handle larger maps but it requires a solved state as an initial condition. We continued to meet over coffee in summer of 2023 to share ideas, mutually throwing around ideas about hierarchical, block-based or local optimization. I was skeptical about the increased code complexity of a hierarchical approach, versus the speed of WFC, and stopped my own experiments. Abram continued to believe a solution was possible. There was a commit gap of four months a from April to July. Then in late July, Abram committed the first Block-WFC and then on Aug 27, 2023 the first Breakout Model Synthesis (BMS) code. This breakthru gave new	life to the project and we commenced testing. I added visualization fixes in Sept 2023, while Abram created a 3D stair tileset. We found that BMS was a significant improvement	over both WFC and MMS, especially on larger maps. We could quickly solve Pacman and others maps up to 256^2. BMS still had issues with certainhighly constrainted tilesets. Abram was doing correlation radius experiments by Sept 14, 2023, and I developed a harder 3D Railway tileset. The hierarchical, block approach to localized constraint tiling generation demonstrated that it could be valuable technical contribution, and I am glad that Abram introduced the problem, revitalized the BP algorithm	which was (and still is) largely unknown and, despite my doubts, pushed to make BMS into a working technique for large maps. This repo was never intended to be a final repository and is maintained here for historic reasons.</sub>
<br><br>
Commit history for BMS can be <a href="https://github.com/ramakarl/just_math/commits/main/math_belief_prop">found here.</a><br>
BMS code can be found in <a href="https://github.com/ramakarl/just_math/tree/main/math_belief_prop">math_belief_prop</a>, which also includes BP and WFC.<br>
Abram's more recent & updated work, Punch Out Model Synthesis, can be <a href="https://zzyzek.github.io/PunchOutModelSynthesisPaper">found here</a> and on <a href="https://arxiv.org/abs/2501.14786">arXiv here.</a><br><br>

## News & Updates

**UPDATES**: <br>
Dec 5, 2023 - Build steps confirmed working with latest libmin repo and Visual Studio 2019 in Dec 2023 for all samples. Removed <a href="https://github.com/ramakarl/libmin">libmin</a> from this repo and put into its own separate repository.<br>
Sep 16,2023 - Obj materials sample added.<br>
Aug 12,2023 - G-coder sample added.<br>
Mar 13,2023 - Bilinear Patch intersection sample added.<br>
Feb 7, 2023 - Artificial Neural Network sample added. PS: For fun I wrote the original wikipedia article for <a href="https://en.wikipedia.org/wiki/Tensor_(machine_learning)">Tensor (machine learning)</a>.<br>

## How to Build
<br>

**Platforms:**
- Win10, Visual Studio 2019/2022 - definitely<br>
- Win10/11, VS{other} - probably<br>
- Linux - should work, post any issues<br>
<br>

**Dependencies:**
- <a href="https://github.com/ramakarl/libmin">Libmin</a> - minimal utilitiy libary for graphics.<br>
- OpenGL <br>
- CUDA is optional (flag at cmake time)<br><br>

**Step 1) Build Libmin** <br>
Cmake and build Libmin from <a href="https://github.com/ramakarl/libmin">here</a> <br>
Libmin repo: <a href="https://github.com/ramakarl/libmin">https://github.com/ramakarl/libmin</a><br>
The binary (build) path should be outside of the \libmin repo folder as follows.<br>
Windows: <br>
1.1) `cmake CMakeLists.txt -B..\build\libmin {options}`<br>
1.2) Open and compile the generated \build\libmin\libmin.sln in Visual Studio 2019+<br>
Linux: <br>
1.1) `cmake CMakeLists.txt -B../build/libmin {options}`<br>
1.2) `make ../build/libmin`<br>
By default, you should not need to specify any additional options for just_math samples.<br>
Options (multiple may be specified):<br>
-DBUILD_OPENGL=true/false - for interactive apps, required for just_math (default=true)<br>
-DBUILD_GLEW=true/false - for interactive apps, required for just_math (default=true)<br>
-DBUILD_CUDA=true/false - for GPU-based apps (default=false)<Br>
-DBUILD_OPENSSL=true/false - for secure network apps (default=false)<br>
-DBUILD_BCRYPT=true/false - for secure network apps (default=false)<br>
-DUSE_NVTX - for Nvidia Nsight timeline-based profiling, cpu & gpu (default=false)<br>
-DUSE_PROFILE_NET - profile network code (default=false)<br><br>

**You must successfully build libmin before proceeding to step 2**.<br>

**Step 2) Build sample**<br>
Cmake and build a sample.<br>
The binary (build) path should be outside of the source \just_math folder as follows (-B option).<br>
Windows: <br>
2.1) `cmake -S \just_math\{sample_name} -B \build\{sample_name} -DLIBMIN_PATH=\build\libmin`<br>
2.2) Open and compile the generated \build\{sample_folder}\{sample_name}.sln in Visual Studio 2019+<br>
Linux: <br>
2.1) `cmake -S \just_math\{sample_name} -B \build\{sample_name} -DLIBMIN_PATH=/build/libmin`<br>
2.2) `make \build\{sample_name}`<br>
All samples are apps that make use of libmin.<br>
Provide the Libmin location:<br>
- Specify the installed path of libmin as LIBMIN_PATH during cmake. eg. `-DLIBMIN_PATH=/usr/local/libmin/`<br>
- If that doesn't work, also specify the path to the source repository with LIBMIN_REPO. eg. `-DLIBMIN_REPO=/libmin`<br>
Build and run the sample.<br>

## Contributions
I am interested in building a community around simple, well documented, math codes, in pure C/C++ for CPU (no shaders), with interactive graphical demos (not just youtube videos) that are MIT/BSD Licensed. If you have similar interests contact me at: Rama Hoetzlein, ramahoetzlein@gmail.com

## License & Copyright info
MIT License.<br>
Copyright 2007-2024 (c) Quanta Sciences & Rama Hoetzlein: 3DDDA, ANN, Basis, Bilinear_Patch, Cells, Deform, Displace, Gcoder, Invk, Obj_Materials, QuadSquad, QuatTraj, Raycast, Voxelizer<br>
Copyright 2022-2023 (c) Abram Connelly and Rama Hoetzlein: math_belief_prop (incl. BMS, WFC, BP, LSA)<br>
The Just Math samples are MIT Licensed.<br>
Libmin is MIT Licensed with contributions from other BSD and MIT licensed sources.<br>
Contact: Rama Hoetzlein at ramahoetzlein@gmail.com



