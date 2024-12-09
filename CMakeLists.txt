cmake_minimum_required(VERSION 2.8)

project(flow_and_transport_2d)

option(USE_MPI "Compile with MPI support" OFF)

set(SOURCES main.cpp
	        Process/process.cpp
			Process/process_flow2d.cpp
			Process/process_advection2d.cpp
            ProcessData/processdata.cpp
			ProcessData/processdata_flow2d.cpp
			ProcessData/processdata_advection2d.cpp
            Discretization/discretization.cpp
			Discretization/fv_diffusion2d.cpp
			Discretization/fv_diffusion2d_tpfa.cpp
#		Utils/utils.cpp
)

add_executable(main ${SOURCES})

find_package(INMOST REQUIRED)
if(NOT INMOST_FOUND)
    message("INMOST not found!")
endif()

#set(INMOST_LIBRARIES ../../INMOST/build-debug/inmost)
include_directories(${INMOST_INCLUDE_DIRS})
message(${INMOST_INCLUDE_DIRS})
message("\n===\n")
message(${INMOST_LIBRARIES})
message("\n===\n")
target_link_libraries(main ${INMOST_LIBRARIES})

if(USE_MPI)
    message("Dealing with MPI")
    find_package(MPI REQUIRED)
    if(NOT MPI_FOUND)
        message("MPI not found!")
    endif()

    message((${MPI_INCLUDE_PATH}))
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(main ${MPI_CXX_LIBRARIES})

    if(MPI_LINK_FLAGS)
        set_target_properties(main PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif()
