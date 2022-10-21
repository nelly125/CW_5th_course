cmake_minimum_required(VERSION 3.21)
project(Task_1_Karman)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_XCODE_SCHEME_WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_executable(Task_1_Karman main.cpp solver.cpp solver.hpp)
