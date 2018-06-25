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

# Include any dependencies generated for this target.
include CMakeFiles/YourLib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/YourLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/YourLib.dir/flags.make

CMakeFiles/YourLib.dir/src/solve_system.cc.o: CMakeFiles/YourLib.dir/flags.make
CMakeFiles/YourLib.dir/src/solve_system.cc.o: src/solve_system.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/neerajsarna/sciebo/DG_MPI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/YourLib.dir/src/solve_system.cc.o"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/YourLib.dir/src/solve_system.cc.o -c /Users/neerajsarna/sciebo/DG_MPI/src/solve_system.cc

CMakeFiles/YourLib.dir/src/solve_system.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/YourLib.dir/src/solve_system.cc.i"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/neerajsarna/sciebo/DG_MPI/src/solve_system.cc > CMakeFiles/YourLib.dir/src/solve_system.cc.i

CMakeFiles/YourLib.dir/src/solve_system.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/YourLib.dir/src/solve_system.cc.s"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/neerajsarna/sciebo/DG_MPI/src/solve_system.cc -o CMakeFiles/YourLib.dir/src/solve_system.cc.s

CMakeFiles/YourLib.dir/src/solve_system.cc.o.requires:

.PHONY : CMakeFiles/YourLib.dir/src/solve_system.cc.o.requires

CMakeFiles/YourLib.dir/src/solve_system.cc.o.provides: CMakeFiles/YourLib.dir/src/solve_system.cc.o.requires
	$(MAKE) -f CMakeFiles/YourLib.dir/build.make CMakeFiles/YourLib.dir/src/solve_system.cc.o.provides.build
.PHONY : CMakeFiles/YourLib.dir/src/solve_system.cc.o.provides

CMakeFiles/YourLib.dir/src/solve_system.cc.o.provides.build: CMakeFiles/YourLib.dir/src/solve_system.cc.o


# Object files for target YourLib
YourLib_OBJECTS = \
"CMakeFiles/YourLib.dir/src/solve_system.cc.o"

# External object files for target YourLib
YourLib_EXTERNAL_OBJECTS =

libYourLib.dylib: CMakeFiles/YourLib.dir/src/solve_system.cc.o
libYourLib.dylib: CMakeFiles/YourLib.dir/build.make
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/lib/libdeal_II.g.8.4.1.dylib
libYourLib.dylib: /usr/lib/libbz2.dylib
libYourLib.dylib: /usr/lib/libz.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libparmetis.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libmetis.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtrilinoscouplings.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpiro.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librol.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_muelu.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_ifpack2.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_amesos2.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_tpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_sacado.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librythmos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-adapters.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-interface.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocathyra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocaepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocalapack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libloca.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxlapack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnox.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libintrepid.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteko.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosbelos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosaztecoo.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosamesos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosml.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosifpack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2-adapters.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazitpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libModeLaplace.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasaziepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazi.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkomplex.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos2.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelostpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelosepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libml.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan2.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen_extras.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-xpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-epetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libaztecoo.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libisorropia.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liboptipack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra-sup.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyratpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetraext.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyracore.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetraext.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraext.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrainout.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkostsqr.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrakernels.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassiclinalg.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassicnodeapi.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassic.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtriutils.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libglobipack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libshards.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetra.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libsacado.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librtop.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscomm.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscompat.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosremainder.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosnumerics.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscomm.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosparameterlist.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscore.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkosalgorithms.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscontainers.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscore.dylib
libYourLib.dylib: /usr/lib/liblapack.dylib
libYourLib.dylib: /usr/lib/libblas.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi_cxx.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libparpack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libarpack.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5_hl.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBO.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBool.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBRep.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKernel.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFeat.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFillet.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG2d.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG3d.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomAlgo.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomBase.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKHLR.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKIGES.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMath.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMesh.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKOffset.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKPrim.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKShHealing.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEP.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPAttr.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPBase.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTL.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKTopAlgo.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKXSBase.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libp4est.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libsc.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/slepc-2c065dd/lib/libslepc.dylib
libYourLib.dylib: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libpetsc.dylib
libYourLib.dylib: CMakeFiles/YourLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/neerajsarna/sciebo/DG_MPI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libYourLib.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/YourLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/YourLib.dir/build: libYourLib.dylib

.PHONY : CMakeFiles/YourLib.dir/build

CMakeFiles/YourLib.dir/requires: CMakeFiles/YourLib.dir/src/solve_system.cc.o.requires

.PHONY : CMakeFiles/YourLib.dir/requires

CMakeFiles/YourLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/YourLib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/YourLib.dir/clean

CMakeFiles/YourLib.dir/depend:
	cd /Users/neerajsarna/sciebo/DG_MPI && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI /Users/neerajsarna/sciebo/DG_MPI/CMakeFiles/YourLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/YourLib.dir/depend

