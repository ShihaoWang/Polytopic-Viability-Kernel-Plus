# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/cmake/bin/cmake

# The command to remove a file.
RM = /opt/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build

# Include any dependencies generated for this target.
include CMakeFiles/MyApp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MyApp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MyApp.dir/flags.make

CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o: CMakeFiles/MyApp.dir/flags.make
CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o: ../src/HJB_VI_Plus.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o -c /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/src/HJB_VI_Plus.cpp

CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/src/HJB_VI_Plus.cpp > CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.i

CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/src/HJB_VI_Plus.cpp -o CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.s

# Object files for target MyApp
MyApp_OBJECTS = \
"CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o"

# External object files for target MyApp
MyApp_EXTERNAL_OBJECTS =

MyApp: CMakeFiles/MyApp.dir/src/HJB_VI_Plus.cpp.o
MyApp: CMakeFiles/MyApp.dir/build.make
MyApp: CMakeFiles/MyApp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MyApp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MyApp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MyApp.dir/build: MyApp

.PHONY : CMakeFiles/MyApp.dir/build

CMakeFiles/MyApp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MyApp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MyApp.dir/clean

CMakeFiles/MyApp.dir/depend:
	cd /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build /home/shihao/Desktop/Polytopic-Viability-Kernel-Plus/build/CMakeFiles/MyApp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MyApp.dir/depend

