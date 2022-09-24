add_executable(05 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(05 ${CMAKE_SOURCE_DIR})
target_link_libraries(05 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})