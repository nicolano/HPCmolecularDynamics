add_executable(09 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(09 ${CMAKE_SOURCE_DIR})
target_link_libraries(09 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})