# Example folders
SET(EXAMPLE_DIRS
	./unetwork/core/ds
	./unetwork/core/network
	./unetwork/core/maps
	./unetwork/graphalg/traversal
	./unetwork/graphalg/shortest-paths
	./unetwork/graphalg/encoders
	./unetwork/distributed
	./unetwork/parallel
	./unetwork/predictors
	./unetwork/perf/perfevaluator
	./unetwork/perf/networkmanip
	./unetwork/perf/perfmeasures
	./dnetwork/core/network
	./dnetwork/core/maps
	./dnetwork/graphalg/traversal
	./dnetwork/graphalg/shortest-paths
	./dnetwork/perf/networkmanip
	./dnetwork/perf/perfevaluator
	./dnetwork/predictors
	./simp/CPP/evaluator
	./simp/CPP/predictor
	./ml/classifiers
	./ml/simmeasures
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

