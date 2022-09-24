add_executable(08 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(08 ${CMAKE_SOURCE_DIR})
target_link_libraries(08 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})