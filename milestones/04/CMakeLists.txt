add_executable(04 ${MOLDYN_SRCS} ${MOLDYN_HDRS} main.cpp)
include_directories(04 ${CMAKE_SOURCE_DIR})
target_link_libraries(04 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})