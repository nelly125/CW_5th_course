cmake_minimum_required(VERSION 3.21)
project(Task_1_Karman)

set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_XCODE_SCHEME_WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
add_compile_options(-mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format)

add_executable(Task_1_Karman main.cpp solver.cpp solver.hpp system/system_helpers.cpp system/system_helpers.hpp ../Task_3/main.cpp)
