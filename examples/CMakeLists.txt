# Example folders
SET(EXAMPLE_DIRS
	./unetwork/core/datastructures
	./unetwork/core/maps
	./unetwork/core/traversal
	./unetwork/core/network
	./unetwork/core/shortest-paths
	./unetwork/distributed
	./unetwork/parallel
	./unetwork/predictors
	./unetwork/perf/perfevaluator
	./unetwork/perf/networkmanip
	./unetwork/perf/perfmeasures
	./dnetwork/core/network
	./dnetwork/core/traversal
	./dnetwork/core/maps
	./dnetwork/core/shortest-paths
	./dnetwork/perf/networkmanip
	./dnetwork/perf/perfevaluator
	./dnetwork/predictors
)

# Compiling examples
FOREACH(DIR ${EXAMPLE_DIRS} )
	FILE( GLOB SOURCES ./${DIR}/*.cpp )
	SET( SOURCE_PATH ${DIR} )
	STRING( REPLACE "/" "_" DIR "${DIR}" )
	STRING( REPLACE "._" "" DIR "${DIR}" )
	FOREACH( SRC ${SOURCES} )
		GET_FILENAME_COMPONENT( EXEC_NAME ${SRC} NAME_WE )
		SET( CURR_EXE examples_${DIR}_${EXEC_NAME} )
		ADD_EXECUTABLE( ${CURR_EXE} ${SRC} )
		TARGET_LINK_LIBRARIES(
			${CURR_EXE}
			LinkPred
		)
	ENDFOREACH( SRC ${SOURCES} )
ENDFOREACH( DIR ${EXAMPLE_DIRS} )

