cmake_minimum_required(VERSION 3.9)
project(SurfaceCodeModule)

set(CMAKE_CXX_STANDARD 11)
IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
ENDIF()

if (APPLE)
    set(PYTHON_LIBRARY /Users/fish/anaconda3/lib/libpython3.6m.dylib)
    #    set(PYTHON_LIBRARY /Users/fish/anaconda3/lib/python3.6)
    set(PYTHON_INCLUDE_DIR /Users/fish/anaconda3/include/python3.6m)
    #    set(PYTHON_INCLUDE_DIR /Users/fish/anaconda3/bin/python3.6m)
    #    include_directories(/Users/fish/anaconda3/bin/python3.6m)
endif()


set(SOURCE_FILES main.cpp ErrorFunc.h HelperFunc.h ErrorFunc.cpp)

get_filename_component(PARENT_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)

set(BLOSSOM_LIBRARY blossom5-v2.05/lib/libblossom.a)


add_executable(SurfaceCode ${SOURCE_FILES})

get_filename_component(PARENT_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
target_include_directories(SurfaceCode PRIVATE blossom5-v2.05)
include_directories(${PYTHON_INCLUDE_DIR})

set(BLOSSOM_LIBRARY blossom5-v2.05/lib/libblossom.a)
target_link_libraries(SurfaceCode ${PYTHON_LIBRARY})
target_link_libraries(SurfaceCode ${BLOSSOM_LIBRARY})





#set(SOURCE_FILES main.cpp ErrorFunc.h HelperFunc.h ErrorFunc.cpp)
#
#get_filename_component(PARENT_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
#
#set(BLOSSOM_LIBRARY blossom5-v2.05/lib/libblossom.a)
#
#
##find_package(PythonLibs REQUIRED)
#
#add_library(SurfaceCode SHARED ${SOURCE_FILES})
#if (APPLE)
#    set(PYTHON_LIBRARY /Users/fish/anaconda3/lib/libpython3.6m.dylib)
#    #    set(PYTHON_LIBRARY /Users/fish/anaconda3/lib/python3.6)
#    set(PYTHON_INCLUDE_DIR /Users/fish/anaconda3/include/python3.6m)
#    #    set(PYTHON_INCLUDE_DIR /Users/fish/anaconda3/bin/python3.6m)
#    #    include_directories(/Users/fish/anaconda3/bin/python3.6m)
#    set_property(TARGET SurfaceCode PROPERTY PREFIX "")
#    set_property(TARGET SurfaceCode PROPERTY OUTPUT_NAME "SurfaceCode.so")
#    set_property(TARGET SurfaceCode PROPERTY SUFFIX "")
#endif()
#
#target_include_directories(SurfaceCode PRIVATE blossom5-v2.05)
#include_directories(${PYTHON_INCLUDE_DIR})
#
#target_link_libraries(SurfaceCode ${PYTHON_LIBRARY})
#target_link_libraries(SurfaceCode ${BLOSSOM_LIBRARY})
#IF(CMAKE_COMPILER_IS_GNUCXX)
#    ADD_DEFINITIONS("-Wall")
#ENDIF()
