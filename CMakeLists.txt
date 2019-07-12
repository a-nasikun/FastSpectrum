cmake_minimum_required (VERSION 3.10 FATAL_ERROR)
project(FastSpectrum LANGUAGES CXX CUDA)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
cmake_policy(SET CMP0076 NEW)

# OpenMP, OpenGL
# #################################################################################################
find_package(OpenMP)
find_package(OpenGL)

# EIGEN
# #################################################################################################
add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE libigl/external/eigen)
if(MSVC)
	target_sources(eigen INTERFACE libigl/external/eigen/debug/msvc/eigen.natvis)
endif()

# GLAD
# #################################################################################################
add_subdirectory(libigl/external/glad)

# GLFW
# #################################################################################################
option(GLFW_BUILD_EXAMPLES "" OFF)
option(GLFW_BUILD_TESTS    "" OFF)
option(GLFW_BUILD_DOCS     "" OFF)
option(GLFW_INSTALL        "" OFF)
option(GLFW_VULKAN_STATIC  "" OFF)
option(GLFW_USE_OSMESA     "" OFF)
option(GLFW_USE_HYBRID_HPG "" OFF)
option(GLFW_USE_WAYLAND    "" OFF)
option(GLFW_USE_MIR        "" OFF)
add_subdirectory(libigl/external/glfw)

# IMGUI
# #################################################################################################
add_subdirectory(libigl/external/imgui)

# IGL
# #################################################################################################
add_library(libigl INTERFACE)
target_include_directories(libigl INTERFACE libigl/include)

# Fast spectrum
# #################################################################################################
add_executable(FastSpectrum
	"Fast Spectrum/Main.cpp"
	"Fast Spectrum/FastSpectrum.h"
	"Fast Spectrum/FastSpectrum.cpp"
	"Fast Spectrum/FastSpectrumGUI.h"
	"Fast Spectrum/Basis.h"
	"Fast Spectrum/Basis.cpp"
	"Fast Spectrum/EigenSolver.h"
	"Fast Spectrum/EigenSolver.cpp"
	"Fast Spectrum/Sampling.h"
	"Fast Spectrum/Sampling.cpp"
	"Fast Spectrum/Utility.h"
	"Fast Spectrum/Utility.cpp"
	"Fast Spectrum/Visualization.h"
	"Fast Spectrum/Visualization.cpp")
target_link_libraries(FastSpectrum
	PUBLIC libigl
	PUBLIC eigen
	PUBLIC glad
	PUBLIC glfw
	PUBLIC imgui
	PUBLIC OpenGL::GL
	PUBLIC OpenMP::OpenMP_CXX
	PUBLIC cudart
	PUBLIC cublas
	PUBLIC cusolver
	PUBLIC cusparse)
set_target_properties(FastSpectrum PROPERTIES
	VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
	RUNTIME_OUTPUT_DIRECTORY $<1:${CMAKE_CURRENT_LIST_DIR}/bin>)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT FastSpectrum)
if(MSVC)
	target_compile_definitions(FastSpectrum PUBLIC _USE_MATH_DEFINES)
endif()
