cmake_minimum_required(VERSION 3.11)

project(molDyn)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(MOLDYN_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/types.h
        ${CMAKE_CURRENT_SOURCE_DIR}/atoms.h
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/lj_direct_summation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/bwgl.h
        ${CMAKE_CURRENT_SOURCE_DIR}/xyz.h
        ${CMAKE_CURRENT_SOURCE_DIR}/kin_energy.h
        ${CMAKE_CURRENT_SOURCE_DIR}/berendsen_thermostat.h
        ${CMAKE_CURRENT_SOURCE_DIR}/set_cubic_lattice.h
        ${CMAKE_CURRENT_SOURCE_DIR}/neighbors.h
        ${CMAKE_CURRENT_SOURCE_DIR}/lj.h
        ${CMAKE_CURRENT_SOURCE_DIR}/gupta.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_support.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domain.h
        ${CMAKE_CURRENT_SOURCE_DIR}/force_on_whisker.h
        )

set(MOLDYN_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lj_direct_summation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/bwgl.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kin_energy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/berendsen_thermostat.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/set_cubic_lattice.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lj.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/gupta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/domain.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/force_on_whisker.cpp
        )

#set(MOLDYN_executable_SRCS
#        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
#        )

add_subdirectory(tests)
add_subdirectory(milestones)

# Simulation code
add_executable(molDyn ${MOLDYN_SRCS} ${MOLDYN_HDRS} ${MOLDYN_executable_SRCS})
target_link_libraries(molDyn PUBLIC Eigen3::Eigen)
