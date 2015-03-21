cmake_minimum_required(VERSION 3.1)
project(deepflow_release1_0_1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    deepflow.c
    image.c
    image.h
    io.c
    io.h
    opticalflow.c
    opticalflow.h
    opticalflow_aux.c
    opticalflow_aux.h
    README
    solver.c
    solver.h)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
add_executable(deepflow_release1_0_1 ${SOURCE_FILES})

target_link_libraries(deepflow_release1_0_1 jpeg m png)

