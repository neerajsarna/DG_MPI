# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/neerajsarna/sciebo/DG_MPI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/neerajsarna/sciebo/DG_MPI

# Utility rule file for debug.

# Include the progress variables for this target.
include CMakeFiles/debug.dir/progress.make

CMakeFiles/debug:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/neerajsarna/sciebo/DG_MPI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Switch CMAKE_BUILD_TYPE to Debug"
	/Applications/CMake.app/Contents/bin/cmake -DCMAKE_BUILD_TYPE=Debug /Users/neerajsarna/sciebo/DG_MPI
	/Applications/CMake.app/Contents/bin/cmake --build /Users/neerajsarna/sciebo/DG_MPI --target all

debug: CMakeFiles/debug
debug: CMakeFiles/debug.dir/build.make

.PHONY : debug

# Rule to build all files generated by this target.
CMakeFiles/debug.dir/build: debug

.PHONY : CMakeFiles/debug.dir/build

CMakeFiles/debug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/debug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/debug.dir/clean

CMakeFiles/debug.dir/depend:
	cd /Users/neerajsarna/sciebo/DG_MPI && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI/CMakeFiles/debug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/debug.dir/depend

