add_executable(06 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(06 ${CMAKE_SOURCE_DIR})
target_link_libraries(06 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})