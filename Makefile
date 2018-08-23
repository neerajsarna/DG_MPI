# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/Applications/CMake.app/Contents/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/Applications/CMake.app/Contents/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/neerajsarna/sciebo/DG_MPI/CMakeFiles /Users/neerajsarna/sciebo/DG_MPI/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/neerajsarna/sciebo/DG_MPI/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named debug

# Build rule for target.
debug: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 debug
.PHONY : debug

# fast build rule for target.
debug/fast:
	$(MAKE) -f CMakeFiles/debug.dir/build.make CMakeFiles/debug.dir/build
.PHONY : debug/fast

#=============================================================================
# Target rules for targets named release

# Build rule for target.
release: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 release
.PHONY : release

# fast build rule for target.
release/fast:
	$(MAKE) -f CMakeFiles/release.dir/build.make CMakeFiles/release.dir/build
.PHONY : release/fast

#=============================================================================
# Target rules for targets named 2x3v_heated_cavity_Adp.out

# Build rule for target.
2x3v_heated_cavity_Adp.out: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 2x3v_heated_cavity_Adp.out
.PHONY : 2x3v_heated_cavity_Adp.out

# fast build rule for target.
2x3v_heated_cavity_Adp.out/fast:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build
.PHONY : 2x3v_heated_cavity_Adp.out/fast

2x3v_heated_cavity_Adp.o: 2x3v_heated_cavity_Adp.cc.o

.PHONY : 2x3v_heated_cavity_Adp.o

# target to build an object file
2x3v_heated_cavity_Adp.cc.o:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/2x3v_heated_cavity_Adp.cc.o
.PHONY : 2x3v_heated_cavity_Adp.cc.o

2x3v_heated_cavity_Adp.i: 2x3v_heated_cavity_Adp.cc.i

.PHONY : 2x3v_heated_cavity_Adp.i

# target to preprocess a source file
2x3v_heated_cavity_Adp.cc.i:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/2x3v_heated_cavity_Adp.cc.i
.PHONY : 2x3v_heated_cavity_Adp.cc.i

2x3v_heated_cavity_Adp.s: 2x3v_heated_cavity_Adp.cc.s

.PHONY : 2x3v_heated_cavity_Adp.s

# target to generate assembly for a file
2x3v_heated_cavity_Adp.cc.s:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/2x3v_heated_cavity_Adp.cc.s
.PHONY : 2x3v_heated_cavity_Adp.cc.s

src/ic_bc_base.o: src/ic_bc_base.cc.o

.PHONY : src/ic_bc_base.o

# target to build an object file
src/ic_bc_base.cc.o:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/ic_bc_base.cc.o
.PHONY : src/ic_bc_base.cc.o

src/ic_bc_base.i: src/ic_bc_base.cc.i

.PHONY : src/ic_bc_base.i

# target to preprocess a source file
src/ic_bc_base.cc.i:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/ic_bc_base.cc.i
.PHONY : src/ic_bc_base.cc.i

src/ic_bc_base.s: src/ic_bc_base.cc.s

.PHONY : src/ic_bc_base.s

# target to generate assembly for a file
src/ic_bc_base.cc.s:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/ic_bc_base.cc.s
.PHONY : src/ic_bc_base.cc.s

src/run_problem.o: src/run_problem.cc.o

.PHONY : src/run_problem.o

# target to build an object file
src/run_problem.cc.o:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/run_problem.cc.o
.PHONY : src/run_problem.cc.o

src/run_problem.i: src/run_problem.cc.i

.PHONY : src/run_problem.i

# target to preprocess a source file
src/run_problem.cc.i:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/run_problem.cc.i
.PHONY : src/run_problem.cc.i

src/run_problem.s: src/run_problem.cc.s

.PHONY : src/run_problem.s

# target to generate assembly for a file
src/run_problem.cc.s:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/run_problem.cc.s
.PHONY : src/run_problem.cc.s

src/solve_system_SS_adaptive.o: src/solve_system_SS_adaptive.cc.o

.PHONY : src/solve_system_SS_adaptive.o

# target to build an object file
src/solve_system_SS_adaptive.cc.o:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/solve_system_SS_adaptive.cc.o
.PHONY : src/solve_system_SS_adaptive.cc.o

src/solve_system_SS_adaptive.i: src/solve_system_SS_adaptive.cc.i

.PHONY : src/solve_system_SS_adaptive.i

# target to preprocess a source file
src/solve_system_SS_adaptive.cc.i:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/solve_system_SS_adaptive.cc.i
.PHONY : src/solve_system_SS_adaptive.cc.i

src/solve_system_SS_adaptive.s: src/solve_system_SS_adaptive.cc.s

.PHONY : src/solve_system_SS_adaptive.s

# target to generate assembly for a file
src/solve_system_SS_adaptive.cc.s:
	$(MAKE) -f CMakeFiles/2x3v_heated_cavity_Adp.out.dir/build.make CMakeFiles/2x3v_heated_cavity_Adp.out.dir/src/solve_system_SS_adaptive.cc.s
.PHONY : src/solve_system_SS_adaptive.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... debug"
	@echo "... release"
	@echo "... 2x3v_heated_cavity_Adp.out"
	@echo "... 2x3v_heated_cavity_Adp.o"
	@echo "... 2x3v_heated_cavity_Adp.i"
	@echo "... 2x3v_heated_cavity_Adp.s"
	@echo "... src/ic_bc_base.o"
	@echo "... src/ic_bc_base.i"
	@echo "... src/ic_bc_base.s"
	@echo "... src/run_problem.o"
	@echo "... src/run_problem.i"
	@echo "... src/run_problem.s"
	@echo "... src/solve_system_SS_adaptive.o"
	@echo "... src/solve_system_SS_adaptive.i"
	@echo "... src/solve_system_SS_adaptive.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

