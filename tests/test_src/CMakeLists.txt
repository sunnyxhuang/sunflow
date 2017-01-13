set(gtest_SOURCE_DIR
        ../test_lib/googletest)
set(ximulator_SOURCE_DIR
        ../..)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(ximulator_test
        ximulator_test.cc)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(${ximulator_SOURCE_DIR})

target_link_libraries(ximulator_test gtest gtest_main)
target_link_libraries(ximulator_test ximulator)