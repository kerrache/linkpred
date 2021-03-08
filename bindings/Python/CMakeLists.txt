cmake_minimum_required (VERSION 3.0)
project ( _LinkPredPython )
# The version number.
# The project version number.
set(_LinkPredPython_VERSION_MAJOR   1   CACHE STRING "Project major version number.")
set(_LinkPredPython_VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(_LinkPredPython_VERSION_PATCH   0   CACHE STRING "Project patch version number.")
mark_as_advanced(_LinkPredPython_VERSION_MAJOR _LinkPredPython_VERSION_MINOR _LinkPredPython_VERSION_PATCH)


find_package (Python COMPONENTS Interpreter Development REQUIRED)

message("${Python_INCLUDE_DIRS}")
message("${Python_LIBRARIES}")

# Check build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -ansi -fno-strict-aliasing")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

# We want to see all details
set( CMAKE_VERBOSE_MAKEFILE on )

# Compiler options
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/LinkPredPythonConfig.hpp.in"
  "${PROJECT_BINARY_DIR}/include/LinkPredPythonConfig.hpp"
  )

# Add the binary tree to the search path for include files so that we will find _LinkPredPythonConfig.hpp
include_directories("${PROJECT_BINARY_DIR}")

# Should we use MKL?
OPTION (LINKPRED_WITH_MKL "Use Intel Math Kernel Library ?" OFF) 
if (LINKPRED_WITH_MKL)
	message("MKL enabled")
	set (EXTRA_LIBS ${EXTRA_LIBS} mkl_rt)
endif ()

OPTION (LINKPRED_WITH_OPENMP "Use OpenMP ?" ON)
IF(LINKPRED_WITH_OPENMP)
	FIND_PACKAGE(OpenMP)
	message("OpenMP enabled")
	IF(OPENMP_FOUND)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (EXTRA_LIBS ${EXTRA_LIBS} gomp)
	ENDIF()
ENDIF()

OPTION (LINKPRED_WITH_MPI "Use MPI ?" ON)
IF(LINKPRED_WITH_MPI)
	message("MPI enabled")
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)
		include_directories("${MPI_C_INCLUDE_PATH}")
		add_definitions(-DOMPI_SKIP_MPICXX)
		link_libraries("${MPI_C_LIBRARIES}")
		SET (EXTRA_LIBS ${EXTRA_LIBS} ${MPI_C_LIBRARIES})
	ENDIF()
ENDIF()

# Should we use mlpack?
OPTION (LINKPRED_WITH_MLPACK "Use mlpack Library ?" ON) 
if (LINKPRED_WITH_MLPACK)
	message("mlpack enabled")
	set (EXTRA_LIBS ${EXTRA_LIBS} mlpack armadillo)
endif ()

# Python and JNI includes
include_directories("${Python_INCLUDE_DIRS}")
		
# The include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

# Adding Libraries
set (EXTRA_LIBS ${EXTRA_LIBS} gsl gslcblas LinkPred ${Python_LIBRARIES})

# Add target library 
add_library(_LinkPredPython SHARED
	LinkPredPython.cpp
	../../src/simp/predictor.cpp
	../../src/simp/evaluator.cpp
)
target_link_libraries (_LinkPredPython ${EXTRA_LIBS})
set_target_properties(_LinkPredPython PROPERTIES VERSION 1)
# Remove the prefix lib from the targte. This is necessary fro SWIG genarted python code
SET_TARGET_PROPERTIES(_LinkPredPython PROPERTIES PREFIX "")


##############################################################################

# Add the install targets
if(WIN32)
	install(TARGETS _LinkPredPython
	RUNTIME DESTINATION lib)
else()
	install(TARGETS _LinkPredPython
	LIBRARY DESTINATION lib)
endif() 

# Install headers
INSTALL (
	DIRECTORY ${CMAKE_SOURCE_DIR}/include/
	DESTINATION include
	FILES_MATCHING 
		PATTERN "*.h" 
		PATTERN "*.hpp"
)
INSTALL (
	DIRECTORY ${CMAKE_BINARY_DIR}/include/
	DESTINATION include
	FILES_MATCHING 
		PATTERN "*.h" 
		PATTERN "*.hpp"
)

# Build documentation
#find_package(Doxygen)
#option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)" ${DOXYGEN_FOUND})
#if(BUILD_DOCUMENTATION)
#	ADD_SUBDIRECTORY( doc )
#endif()

