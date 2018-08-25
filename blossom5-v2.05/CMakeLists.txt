cmake_minimum_required(VERSION 3.8)

#PROJECT(blossom5 CXX)
set(CMAKE_CXX_STANDARD 11)
IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
ENDIF()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -stdlib=libc++")

add_library(blossom
MinCost/MinCost.cpp
MinCost/MinCost.h
block.h
LCA.h
misc.cpp
PerfectMatching.h
PMduals.cpp
PMexpand.cpp
PMimplementation.h
PMinit.cpp
PMinterface.cpp
PMmain.cpp
PMrepair.cpp
PMshrink.cpp
PQ.h
timer.h)
