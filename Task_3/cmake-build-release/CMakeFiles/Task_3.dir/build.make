# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2022.3.2/bin/cmake/linux/x64/bin/cmake

# The command to remove a file.
RM = /opt/clion-2022.3.2/bin/cmake/linux/x64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nelly/Documents/MSU/CW_5th_course/Task_3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/Task_3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Task_3.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Task_3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Task_3.dir/flags.make

CMakeFiles/Task_3.dir/main.cpp.o: CMakeFiles/Task_3.dir/flags.make
CMakeFiles/Task_3.dir/main.cpp.o: /home/nelly/Documents/MSU/CW_5th_course/Task_3/main.cpp
CMakeFiles/Task_3.dir/main.cpp.o: CMakeFiles/Task_3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Task_3.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Task_3.dir/main.cpp.o -MF CMakeFiles/Task_3.dir/main.cpp.o.d -o CMakeFiles/Task_3.dir/main.cpp.o -c /home/nelly/Documents/MSU/CW_5th_course/Task_3/main.cpp

CMakeFiles/Task_3.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Task_3.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nelly/Documents/MSU/CW_5th_course/Task_3/main.cpp > CMakeFiles/Task_3.dir/main.cpp.i

CMakeFiles/Task_3.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Task_3.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nelly/Documents/MSU/CW_5th_course/Task_3/main.cpp -o CMakeFiles/Task_3.dir/main.cpp.s

# Object files for target Task_3
Task_3_OBJECTS = \
"CMakeFiles/Task_3.dir/main.cpp.o"

# External object files for target Task_3
Task_3_EXTERNAL_OBJECTS =

Task_3: CMakeFiles/Task_3.dir/main.cpp.o
Task_3: CMakeFiles/Task_3.dir/build.make
Task_3: CMakeFiles/Task_3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Task_3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Task_3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Task_3.dir/build: Task_3
.PHONY : CMakeFiles/Task_3.dir/build

CMakeFiles/Task_3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Task_3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Task_3.dir/clean

CMakeFiles/Task_3.dir/depend:
	cd /home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nelly/Documents/MSU/CW_5th_course/Task_3 /home/nelly/Documents/MSU/CW_5th_course/Task_3 /home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release /home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release /home/nelly/Documents/MSU/CW_5th_course/Task_3/cmake-build-release/CMakeFiles/Task_3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Task_3.dir/depend

