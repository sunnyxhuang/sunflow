cmake_minimum_required(VERSION 3.6)
project(Ximulator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.cc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_executable(ximulator_main ${SOURCE_FILES})
target_link_libraries(ximulator_main ximulator)

set(LIB_SOURCE_FILES
        src/coflow.cc
        src/coflow.h
        src/edmond.cc
        src/edmond.h
        src/events.cc
        src/events.h
        src/global.cc
        src/global.h
        src/hopcroft.cc
        src/hopcroft.h
        src/scheduler.cc
        src/scheduler.h
        src/schedulerAalo.cc
        src/schedulerBvN.cc
        src/schedulerEdmond.cc
        src/schedulerOptc.cc
        src/schedulerSolstice.cc
        src/schedulerSunflow.cc
        src/schedulerTMS.cc
        src/schedulerVarys.cc
        src/sunflowSolver.cc
        src/sunflowSolver.h
        src/trafficGen.cc
        src/trafficGen.h
        src/util.cc
        src/util.h)

add_library(ximulator STATIC ${LIB_SOURCE_FILES})

add_subdirectory(tests)