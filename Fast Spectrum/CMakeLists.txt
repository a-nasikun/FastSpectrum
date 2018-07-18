cmake_minimum_required (VERSION 3.8 FATAL_ERROR)
set(PROJECT_NAME FastSpectrum)

project (${PROJECT_NAME} LANGUAGES CXX CUDA)
find_package(CUDA QUIET REQUIRED)

set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)

#Defining the system varibles where the include files reside
include_directories($ENV{IGL_DIR} 
					$ENV{EIGEN_DIR}
					$ENV{OPENMESH_DIR}/include/
					$ENV{MATLAB_DIR}/include/
)

#Defining the system varibles where the link libraries reside					
link_directories($ENV{MATLAB_DIR}/lib/win64/microsoft/
				 $ENV{OPENMESH_DIR}/lib/
				 $ENV{GLFW_DIR}/lib-vc2015/			#Change according to your Visual Studio
)	

#Stating what files to include in the project				 
add_executable(${PROJECT_NAME}
	Main.cpp 
	Basis.h 
	Basis.cpp 
	EigenSolver.h 
	EigenSolver.cpp 
	Sampling.h 
	Sampling.cpp 
	Utility.h 
	Utility.cpp 
	Visualization.h 
	Visualization.cpp
)

target_link_libraries(
	${PROJECT_NAME}
	OpenGL32.lib; 			
	glew32s.lib;
	glfw3.lib; 				glfw3dll.lib;
	cublas.lib; 			cudart.lib; 
	cusolver.lib; 			cusparse.lib;
	libMatlabEngine.lib; 	libMatlabDataArray.lib; 
	libeng.lib; 			libut.lib; 
	libmex.lib; 			libmx.lib;
	OpenMeshCore.lib; 		OpenMeshCored.lib; 
	OpenMeshTools.lib; 		OpenMeshToolsd.lib;
)	

if(MSVC)	
	add_definitions(-D_USE_MATH_DEFINES)
endif()

#To activate OpenMP
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

set(CMAKE_SUPPRESS_REGENERATION true)