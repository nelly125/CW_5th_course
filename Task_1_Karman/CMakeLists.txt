cmake_minimum_required(VERSION 3.23)
project(Task_1_Karman)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_executable(Task_1_Karman main.cpp solver.cpp solver.hpp)
