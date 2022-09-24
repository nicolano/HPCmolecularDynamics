include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(MOLDYN_tests_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lj_direct_summation_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/berendsen_thermostat_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/neighbors_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/gupta_test.cpp
        )

# Tests
include_directories(molDyn_tests ${CMAKE_SOURCE_DIR})
add_executable(molDyn_tests ${MOLDYN_SRCS} ${MOLDYN_HDRS} ${MOLDYN_tests_SRCS})
target_link_libraries(molDyn_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(molDyn_tests)

# For tests that do not use GTest
# add_test(NAME <test name> COMMAND <test executable>)
