add_executable(07 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(07 ${CMAKE_SOURCE_DIR})
target_link_libraries(07 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})