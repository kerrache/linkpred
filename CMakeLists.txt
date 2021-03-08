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

# Should we use Armadillo?
OPTION (LINKPRED_WITH_ARMADILLO "Use Armadillo Library ?" ON) 
if (LINKPRED_WITH_ARMADILLO)
	message("Armadillo enabled")
	set (EXTRA_LIBS ${EXTRA_LIBS} armadillo)
endif ()

# The include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

# Adding Libraries
set (EXTRA_LIBS ${EXTRA_LIBS} gsl gslcblas)

# Add target library 
add_library(LinkPred SHARED
	src/core/dnetwork/dnetwork.cpp
	src/core/unetwork/unetwork.cpp
	src/graphalg/encoders/deepwalk/deepwalk.cpp
	src/graphalg/encoders/hmsm/hmsm.cpp
	src/graphalg/encoders/largevis/largevis.cpp
	src/graphalg/encoders/lem/lem.cpp
	src/graphalg/encoders/line/line.cpp
	src/graphalg/encoders/lle/lle.cpp
	src/graphalg/encoders/matfact/matfact.cpp
	src/graphalg/encoders/matfact/matfactcg.cpp
	src/graphalg/encoders/node2vec/node2vec.cpp
	src/graphalg/shortestpaths/dijkstra.cpp
	src/graphalg/shortestpaths/netdistcalculator.cpp
	src/graphalg/traversal/graphtraversal.cpp
	src/ml/classifiers/ffn/ffn.cpp
	src/ml/classifiers/linearsvm/linearsvm.cpp
	src/ml/classifiers/logistic/logisticregresser.cpp
	src/ml/classifiers/logistic/logregcg.cpp
	src/ml/classifiers/naivebayes/naivebayes.cpp
	src/ml/classifiers/rndclassifier/rndclassifier.cpp
	src/ml/simmeasures/cosinesim.cpp
	src/ml/simmeasures/dotprod.cpp
	src/ml/simmeasures/l1sim.cpp
	src/ml/simmeasures/l2sim.cpp
	src/ml/simmeasures/lpsim.cpp
	src/ml/simmeasures/pearson.cpp
	src/numerical/cg/cgdescent.cpp
	src/numerical/linear/gfmatrix.cpp
	src/numerical/linear/vec.cpp
	src/numerical/mds/logmdscg.cpp
	src/numerical/mds/mds.cpp
	src/numerical/plfit/error.cpp
	src/numerical/plfit/gss.cpp
	src/numerical/plfit/kolmogorov.cpp
	src/numerical/plfit/lbfgs.cpp
	src/numerical/plfit/mt.cpp
	src/numerical/plfit/options.cpp
	src/numerical/plfit/platform.cpp
	src/numerical/plfit/plfit.cpp
	src/numerical/plfit/rbinom.cpp
	src/numerical/plfit/sampling.cpp
	src/numerical/plfit/stats.cpp
	src/numerical/plfit/zeta.cpp
	src/perf/networkmanipulator.cpp
	src/perf/perfevaluator.cpp
	src/predictors/directed/dadapredictor.cpp
	src/predictors/directed/dcnepredictor.cpp
	src/predictors/directed/dhdipredictor.cpp
	src/predictors/directed/dhpipredictor.cpp
	src/predictors/directed/djidpredictor.cpp
	src/predictors/directed/dlcppredictor.cpp
	src/predictors/directed/dlhnpredictor.cpp
	src/predictors/directed/dpatpredictor.cpp
	src/predictors/directed/dpstpredictor.cpp
	src/predictors/directed/dsaipredictor.cpp
	src/predictors/directed/dsoipredictor.cpp
	src/predictors/undirected/uadapredictor.cpp
	src/predictors/undirected/ucnepredictor.cpp
	src/predictors/undirected/ucrapredictor.cpp
	src/predictors/undirected/ucstpredictor.cpp
	src/predictors/undirected/ueclpredictor.cpp
	src/predictors/undirected/uesmpredictor.cpp
	src/predictors/undirected/ufbmpredictor.cpp
	src/predictors/undirected/uhdipredictor.cpp
	src/predictors/undirected/uhpipredictor.cpp
	src/predictors/undirected/uhrgpredictor.cpp
	src/predictors/undirected/uhrgpredictor/dendro_pr.cpp
	src/predictors/undirected/uhrgpredictor/hrgds.cpp
	src/predictors/undirected/uhrgpredictor/rbtree.cpp
	src/predictors/undirected/uhyppredictor.cpp
	src/predictors/undirected/ujidpredictor.cpp
	src/predictors/undirected/ukabpredictor.cpp
	src/predictors/undirected/ulcppredictor.cpp
	src/predictors/undirected/ulhnpredictor.cpp
	src/predictors/undirected/unedpredictor.cpp
	src/predictors/undirected/upatpredictor.cpp
	src/predictors/undirected/upstpredictor.cpp
	src/predictors/undirected/uralpredictor.cpp
	src/predictors/undirected/urndpredictor.cpp
	src/predictors/undirected/usaipredictor.cpp
	src/predictors/undirected/usbmpredictor.cpp
	src/predictors/undirected/usbmpredictor/bipartite.cpp
	src/predictors/undirected/usbmpredictor/conflict.cpp
	src/predictors/undirected/usbmpredictor/datastruct.cpp
	src/predictors/undirected/usbmpredictor/graph_pr.cpp
	src/predictors/undirected/usbmpredictor/graph_simp.cpp
	src/predictors/undirected/usbmpredictor/graph.cpp
	src/predictors/undirected/usbmpredictor/layout.cpp
	src/predictors/undirected/usbmpredictor/matrix.cpp
	src/predictors/undirected/usbmpredictor/missing.cpp
	src/predictors/undirected/usbmpredictor/models.cpp
	src/predictors/undirected/usbmpredictor/modules.cpp
	src/predictors/undirected/usbmpredictor/recommend.cpp
	src/predictors/undirected/usbmpredictor/tools.cpp
	src/predictors/undirected/ushppredictor.cpp
	src/predictors/undirected/usoipredictor.cpp
	src/predictors/undirected/usumpredictor.cpp
	src/simp/evaluator.cpp
	src/simp/predictor.cpp
	src/utils/log.cpp
	src/utils/utils.cpp
)
target_link_libraries (LinkPred ${EXTRA_LIBS})
set_target_properties(LinkPred PROPERTIES VERSION 1)
##############################################################################

# Add the install targets
if(WIN32)
	install(TARGETS LinkPred
	RUNTIME DESTINATION lib)
else()
	install(TARGETS LinkPred
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

# Build the bindings
OPTION( LINKPRED_WITH_BINDINGS "Building LinkPred bindings" ON )
IF( LINKPRED_WITH_BINDINGS )
	ADD_SUBDIRECTORY( bindings )
ENDIF( LINKPRED_WITH_BINDINGS )

# Enable building examples
OPTION( LINKPRED_WITH_EXAMPLES "Building LinkPred examples" ON )
IF( LINKPRED_WITH_EXAMPLES )
	ADD_SUBDIRECTORY( examples )
ENDIF( LINKPRED_WITH_EXAMPLES )

# Build documentation
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if(BUILD_DOCUMENTATION)
	ADD_SUBDIRECTORY( doc )
endif()

