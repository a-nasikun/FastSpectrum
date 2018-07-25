## Fast Approximation of Eigenproblem of Laplace-Beltrami Operator from Triangle Mesh

We _approximate_ the lowest part of the _eigenvalues_ and _eigenfunctions_ of **Laplace-Beltrami operator**, to have faster speed, less storage, and quicker mapping to world coordinates. Our general idea is to construct a **subspace** (from the original space) that can approximate the _low frequency (smooth)_ functions on the manifold. Then we will solve the eigenproblem _restricted to such subspace_.

Here is main construct of our approximation scheme:
1. Construct a subspace  
   1. sampling 
   2. radial basis functions 
   3. partition of unity  
2. Solve the eigenproblem in subpace  
   1. low-dim eigenproblem
   2. lift the eigenfunction

Please note that this program requires **Eigen**, **LibIGL**, and an **eigensolver** (either **CUDA's cuSOLVER** or **Matlab's eigs** function).

## How to Compile and Run Our Implementation

We use **CMake** for easy installation/deployment on your machine, by creating a CMakeFiles that hopefully do (almost) everything to make our code runs in your machine. 

1. **Download** required libraries.
   1. Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
   2. LibIGL: http://libigl.github.io/libigl/ (I'm using the version with commit ID=274dc2de4535e7a91aac782d088c0ee14c6d3f1e) 
   5. GLFW: http://www.glfw.org/download.html   
   6. CUDA (with cuSolver): https://developer.nvidia.com/cuda-zone
2. Set **environment variables** related to those libraries
   1. EIGEN_DIR = /path/to/eigen/main/directory
   2. IGL_DIR = /path/to/libigl/main/directory
   5. GLFW_DIR = /path/to/GLFW/main/directory
3. Make sure you have [CMake](https://cmake.org/ "Download CMake") in your machine
4. **Download/Clone** our repository
5. Use CMake to **set-up** our implementation in your favorite IDE/system. 

Having followed those steps, our program **should compile and run** well. If you have any problem compiling it, please drop me an **email** at a.nasikun[at]tudelft.nl.

## SGP 2018 (Symposium on Geometry Processing 2018)

This work was presented at [SGP 2018](https://sgp2018.sciencesconf.org/resource/page/id/10#s4 "SGP 2018") at Paris TelecomTech, organized by EuroGraphics and ACM SIGGRAPH, and will be published at [Computer Graphics Forum (CGF) Vol.37(5)](https://onlinelibrary.wiley.com/loi/14678659 "CGF Volume 37 Issue 5"). You can read the paper from [our website](https://graphics.tudelft.nl/ahmad-nasikun/ "Fast Approximation of Laplace-Beltrami Eigenproblems") for a more detailed explation of our approximation algorithm. 
