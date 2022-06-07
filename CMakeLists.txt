# - Basic CMake setup
# Check version meets your requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.1)
project(TBeta)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# enable threading
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

#Executables
add_executable(test1.exe examples/test1.cpp)
target_include_directories(test1.exe PRIVATE ${PROJECT_SOURCE_DIR})

#Executables
add_executable(test2.exe examples/test2.cpp)
target_include_directories(test2.exe PRIVATE ${PROJECT_SOURCE_DIR})

#Executables
add_executable(testbg.exe examples/generatorExample.cpp)
target_include_directories(testbg.exe PRIVATE ${PROJECT_SOURCE_DIR})

