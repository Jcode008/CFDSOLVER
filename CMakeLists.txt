cmake_minimum_required(VERSION 3.10)
project(CFDSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

include_directories(${CMAKE_SOURCE_DIR}/include)

# Main CFD solver executable (exclude test files)
add_executable(CFDSolver 
	src/main.cpp
	src/grid.cpp
	src/field.cpp
	src/solver.cpp
	src/utils.cpp
	src/io.cpp
)

# Grid generator test executable
add_executable(TestMesh
	src/test_mesh.cpp
	src/mesh_generator.cpp
	src/ogrid_generator.cpp
)

# Curvilinear solver test executable
add_executable(TestCurvilinear
	src/test_curvilinear.cpp
	src/curvilinear_grid.cpp
	src/curvilinear_field.cpp
	src/curvilinear_solver.cpp
	src/io.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	target_link_libraries(CFDSolver PRIVATE OpenMP::OpenMP_CXX)
endif()

if (MSVC)
	target_compile_options(CFDSolver PRIVATE
		/MP
		$<$<BOOL:${OpenMP_CXX_FOUND}>:/openmp>
		$<$<CONFIG:Release>:/O2>
		$<$<CONFIG:RelWithDebInfo>:/O2>
		$<$<CONFIG:Release>:/fp:fast>
		$<$<CONFIG:RelWithDebInfo>:/fp:fast>
	)
	target_link_options(CFDSolver PRIVATE
		$<$<CONFIG:Release>:/LTCG>
		$<$<CONFIG:RelWithDebInfo>:/LTCG>
	)
else()
	target_compile_options(CFDSolver PRIVATE
		$<$<CONFIG:Release>:-O3>
		$<$<CONFIG:RelWithDebInfo>:-O3>
		$<$<CONFIG:Release>:-ffast-math>
		$<$<CONFIG:RelWithDebInfo>:-ffast-math>
		$<$<CONFIG:Release>:-march=native>
		$<$<CONFIG:RelWithDebInfo>:-march=native>
	)
endif()
