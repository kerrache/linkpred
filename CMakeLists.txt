cmake_minimum_required (VERSION 3.0)
project ( LinkPred)
# The version number.
# The project version number.
set(LinkPred_VERSION_MAJOR   1   CACHE STRING "Project major version number.")
set(LinkPred_VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(LinkPred_VERSION_PATCH   0   CACHE STRING "Project patch version number.")
mark_as_advanced(LinkPred_VERSION_MAJOR LinkPred_VERSION_MINOR LinkPred_VERSION_PATCH)

# Check build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -ansi")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

# We want to see all details
set( CMAKE_VERBOSE_MAKEFILE on )

# Compiler options
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/LinkPredConfig.hpp.in"
  "${PROJECT_BINARY_DIR}/include/LinkPredConfig.hpp"
  )

# Add the binary tree to the search path for include files so that we will find LinkPredConfig.hpp
include_directories("${PROJECT_BINARY_DIR}")

# Should we use MKL?
OPTION (WITH_MKL "Use Intel Math Kernel Library ?" OFF) 
if (WITH_MKL)
	message("MKL enabled")
	set (EXTRA_LIBS ${EXTRA_LIBS} mkl_rt)
endif ()

OPTION (WITH_OPENMP "Use OpenMP ?" ON)
IF(WITH_OPENMP)
	FIND_PACKAGE(OpenMP)
	message("OpenMP enabled")
	IF(OPENMP_FOUND)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (EXTRA_LIBS ${EXTRA_LIBS} gomp)
	ENDIF()
ENDIF()

OPTION (WITH_MPI "Use MPI ?" ON)
IF(WITH_MPI)
	message("MPI enabled")
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)
		include_directories("${MPI_C_INCLUDE_PATH}")
		add_definitions(-DOMPI_SKIP_MPICXX)
		link_libraries("${MPI_C_LIBRARIES}")
		SET (EXTRA_LIBS ${EXTRA_LIBS} ${MPI_C_LIBRARIES})
	ENDIF()
ENDIF()


# Enable building examples
OPTION( LINKPRED_WITH_EXAMPLES "Building LinkPred examples" ON )

# The include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

# Adding Libraries
set (EXTRA_LIBS ${EXTRA_LIBS} gsl gslcblas)

# Add target library 
add_library(LinkPred SHARED
	src/utils/utils.cpp
	src/utils/log.cpp
	src/numerical/plfit/mt.cpp
	src/numerical/plfit/platform.cpp
	src/numerical/plfit/zeta.cpp
	src/numerical/plfit/plfit.cpp
	src/numerical/plfit/stats.cpp
	src/numerical/plfit/gss.cpp
	src/numerical/plfit/error.cpp
	src/numerical/plfit/rbinom.cpp
	src/numerical/plfit/options.cpp
	src/numerical/plfit/sampling.cpp
	src/numerical/plfit/kolmogorov.cpp
	src/numerical/plfit/lbfgs.cpp
	src/numerical/linear/vec.cpp
	src/numerical/linear/gfmatrix.cpp
	src/numerical/cg/cgdescent.cpp
	src/numerical/logistic/logisticregresser.cpp
	src/numerical/logistic/logregcg.cpp
	src/core/unetwork.cpp
	src/core/dnetwork.cpp
	src/core/graphtraversal.cpp
	src/core/netdistcalculator.cpp
	src/core/dijkstra.cpp
	src/predictors/upatpredictor.cpp
	src/predictors/ujidpredictor.cpp
	src/predictors/uhpipredictor.cpp
	src/predictors/uenspredictor.cpp
	src/predictors/ucnepredictor.cpp
	src/predictors/ucrapredictor.cpp
	src/predictors/ukabpredictor.cpp
	src/predictors/ushppredictor.cpp
	src/predictors/uralpredictor.cpp
	src/predictors/uhrgpredictor.cpp
	src/predictors/ucstpredictor.cpp
	src/predictors/uhdipredictor.cpp
	src/predictors/uhrgpredictor/hrgds.cpp
	src/predictors/uhrgpredictor/rbtree.cpp
	src/predictors/uhrgpredictor/dendro_pr.cpp
	src/predictors/ulhnpredictor.cpp
	src/predictors/ufbmpredictor.cpp
	src/predictors/uhyppredictor.cpp
	src/predictors/usaipredictor.cpp
	src/predictors/usnspredictor.cpp
	src/predictors/ulcppredictor.cpp
	src/predictors/uadapredictor.cpp
	src/predictors/urndpredictor.cpp
	src/predictors/ukabpredictor/kablambdacg.cpp
	src/predictors/ukabpredictor/kabwcg.cpp
	src/predictors/usoipredictor.cpp
	src/predictors/usbmpredictor/tools.cpp
	src/predictors/usbmpredictor/graph.cpp
	src/predictors/usbmpredictor/conflict.cpp
	src/predictors/usbmpredictor/datastruct.cpp
	src/predictors/usbmpredictor/graph_simp.cpp
	src/predictors/usbmpredictor/layout.cpp
	src/predictors/usbmpredictor/recommend.cpp
	src/predictors/usbmpredictor/graph_pr.cpp
	src/predictors/usbmpredictor/bipartite.cpp
	src/predictors/usbmpredictor/models.cpp
	src/predictors/usbmpredictor/matrix.cpp
	src/predictors/usbmpredictor/missing.cpp
	src/predictors/usbmpredictor/modules.cpp
	src/predictors/usbmpredictor.cpp
	src/predictors/usumpredictor.cpp
	src/predictors/upoppredictor.cpp
	src/predictors/unedpredictor.cpp
	src/predictors/upndpredictor.cpp
	src/predictors/dcnepredictor.cpp
	src/predictors/dadapredictor.cpp
	src/predictors/dsaipredictor.cpp
	src/predictors/djidpredictor.cpp
	src/predictors/dsoipredictor.cpp
	src/predictors/dhdipredictor.cpp
	src/predictors/dhpipredictor.cpp
	src/predictors/dlhnpredictor.cpp
	src/predictors/dpatpredictor.cpp
	src/predictors/dlcppredictor.cpp
	src/perf/networkmanipulator.cpp
	src/perf/perfevaluator.cpp
)
target_link_libraries (LinkPred ${EXTRA_LIBS})

# Add the install targets
if(WIN32)
	install(TARGETS LinkPred
	RUNTIME DESTINATION lib)
else()
	install(TARGETS LinkPred
	LIBRARY DESTINATION lib)
endif() 

set_target_properties(LinkPred PROPERTIES VERSION 1)

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

# Build examples
IF( LINKPRED_WITH_EXAMPLES )
	ADD_SUBDIRECTORY( examples )
ENDIF( LINKPRED_WITH_EXAMPLES )

# Build apps
ADD_SUBDIRECTORY( apps )

# Build documentation
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if(BUILD_DOCUMENTATION)
	ADD_SUBDIRECTORY( doc )
endif()

