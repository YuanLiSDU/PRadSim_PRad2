# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/liyuan/PRad/PRadSim_PRad2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/liyuan/PRad/PRadSim_PRad2

# Include any dependencies generated for this target.
include digitization/CMakeFiles/PRadRec.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include digitization/CMakeFiles/PRadRec.dir/compiler_depend.make

# Include the progress variables for this target.
include digitization/CMakeFiles/PRadRec.dir/progress.make

# Include the compile flags for this target's objects.
include digitization/CMakeFiles/PRadRec.dir/flags.make

digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o: digitization/CMakeFiles/PRadRec.dir/flags.make
digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o: digitization/PRadRec.cc
digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o: digitization/CMakeFiles/PRadRec.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/liyuan/PRad/PRadSim_PRad2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o"
	cd /home/liyuan/PRad/PRadSim_PRad2/digitization && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o -MF CMakeFiles/PRadRec.dir/PRadRec.cc.o.d -o CMakeFiles/PRadRec.dir/PRadRec.cc.o -c /home/liyuan/PRad/PRadSim_PRad2/digitization/PRadRec.cc

digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PRadRec.dir/PRadRec.cc.i"
	cd /home/liyuan/PRad/PRadSim_PRad2/digitization && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/liyuan/PRad/PRadSim_PRad2/digitization/PRadRec.cc > CMakeFiles/PRadRec.dir/PRadRec.cc.i

digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PRadRec.dir/PRadRec.cc.s"
	cd /home/liyuan/PRad/PRadSim_PRad2/digitization && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/liyuan/PRad/PRadSim_PRad2/digitization/PRadRec.cc -o CMakeFiles/PRadRec.dir/PRadRec.cc.s

# Object files for target PRadRec
PRadRec_OBJECTS = \
"CMakeFiles/PRadRec.dir/PRadRec.cc.o"

# External object files for target PRadRec
PRadRec_EXTERNAL_OBJECTS =

digitization/PRadRec: digitization/CMakeFiles/PRadRec.dir/PRadRec.cc.o
digitization/PRadRec: digitization/CMakeFiles/PRadRec.dir/build.make
digitization/PRadRec: /home/liyuan/root/install/lib/libCore.so
digitization/PRadRec: /home/liyuan/root/install/lib/libImt.so
digitization/PRadRec: /home/liyuan/root/install/lib/libRIO.so
digitization/PRadRec: /home/liyuan/root/install/lib/libNet.so
digitization/PRadRec: /home/liyuan/root/install/lib/libHist.so
digitization/PRadRec: /home/liyuan/root/install/lib/libGraf.so
digitization/PRadRec: /home/liyuan/root/install/lib/libGraf3d.so
digitization/PRadRec: /home/liyuan/root/install/lib/libGpad.so
digitization/PRadRec: /home/liyuan/root/install/lib/libROOTDataFrame.so
digitization/PRadRec: /home/liyuan/root/install/lib/libTree.so
digitization/PRadRec: /home/liyuan/root/install/lib/libTreePlayer.so
digitization/PRadRec: /home/liyuan/root/install/lib/libRint.so
digitization/PRadRec: /home/liyuan/root/install/lib/libPostscript.so
digitization/PRadRec: /home/liyuan/root/install/lib/libMatrix.so
digitization/PRadRec: /home/liyuan/root/install/lib/libPhysics.so
digitization/PRadRec: /home/liyuan/root/install/lib/libMathCore.so
digitization/PRadRec: /home/liyuan/root/install/lib/libThread.so
digitization/PRadRec: /home/liyuan/root/install/lib/libMultiProc.so
digitization/PRadRec: /home/liyuan/root/install/lib/libROOTVecOps.so
digitization/PRadRec: digitization/CMakeFiles/PRadRec.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/liyuan/PRad/PRadSim_PRad2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable PRadRec"
	cd /home/liyuan/PRad/PRadSim_PRad2/digitization && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PRadRec.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
digitization/CMakeFiles/PRadRec.dir/build: digitization/PRadRec
.PHONY : digitization/CMakeFiles/PRadRec.dir/build

digitization/CMakeFiles/PRadRec.dir/clean:
	cd /home/liyuan/PRad/PRadSim_PRad2/digitization && $(CMAKE_COMMAND) -P CMakeFiles/PRadRec.dir/cmake_clean.cmake
.PHONY : digitization/CMakeFiles/PRadRec.dir/clean

digitization/CMakeFiles/PRadRec.dir/depend:
	cd /home/liyuan/PRad/PRadSim_PRad2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/liyuan/PRad/PRadSim_PRad2 /home/liyuan/PRad/PRadSim_PRad2/digitization /home/liyuan/PRad/PRadSim_PRad2 /home/liyuan/PRad/PRadSim_PRad2/digitization /home/liyuan/PRad/PRadSim_PRad2/digitization/CMakeFiles/PRadRec.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : digitization/CMakeFiles/PRadRec.dir/depend

