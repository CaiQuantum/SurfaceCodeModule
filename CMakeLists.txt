cmake_minimum_required(VERSION 3.8)
project(ToricCode)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp ErrorFunc.h HelperFunc.h)

#find_package(blossom5-v2.05 REQUIRED)
add_library(blossom5-v2.05
        blossom5-v2.05/MinCost/MinCost.cpp
        blossom5-v2.05/MinCost/MinCost.h
        blossom5-v2.05/block.h
        blossom5-v2.05/LCA.h
        blossom5-v2.05/misc.cpp
        blossom5-v2.05/PerfectMatching.h
        blossom5-v2.05/PMduals.cpp
        blossom5-v2.05/PMexpand.cpp
        blossom5-v2.05/PMimplementation.h
        blossom5-v2.05/PMinit.cpp
        blossom5-v2.05/PMinterface.cpp
        blossom5-v2.05/PMmain.cpp
        blossom5-v2.05/PMrepair.cpp
        blossom5-v2.05/PMshrink.cpp
        blossom5-v2.05/PQ.h
        blossom5-v2.05/timer.h)

add_executable(ToricCode ${SOURCE_FILES})
target_include_directories(ToricCode PRIVATE blossom5-v2.05)

target_link_libraries(ToricCode blossom5-v2.05)