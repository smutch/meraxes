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
CMAKE_COMMAND = /fred/oz113/yqin/3rd_party/lib/python3.6/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /fred/oz113/yqin/3rd_party/lib/python3.6/site-packages/cmake/data/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /fred/oz025/yqin/bitbucket/meraxes/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /fred/oz025/yqin/bitbucket/meraxes/src

# Include any dependencies generated for this target.
include CMakeFiles/mlog.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mlog.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mlog.dir/flags.make

CMakeFiles/mlog.dir/mlog/mlog.c.o: CMakeFiles/mlog.dir/flags.make
CMakeFiles/mlog.dir/mlog/mlog.c.o: mlog/mlog.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/mlog.dir/mlog/mlog.c.o"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mlog.dir/mlog/mlog.c.o   -c /fred/oz025/yqin/bitbucket/meraxes/src/mlog/mlog.c

CMakeFiles/mlog.dir/mlog/mlog.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mlog.dir/mlog/mlog.c.i"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /fred/oz025/yqin/bitbucket/meraxes/src/mlog/mlog.c > CMakeFiles/mlog.dir/mlog/mlog.c.i

CMakeFiles/mlog.dir/mlog/mlog.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mlog.dir/mlog/mlog.c.s"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /fred/oz025/yqin/bitbucket/meraxes/src/mlog/mlog.c -o CMakeFiles/mlog.dir/mlog/mlog.c.s

# Object files for target mlog
mlog_OBJECTS = \
"CMakeFiles/mlog.dir/mlog/mlog.c.o"

# External object files for target mlog
mlog_EXTERNAL_OBJECTS =

lib/libmlog.a: CMakeFiles/mlog.dir/mlog/mlog.c.o
lib/libmlog.a: CMakeFiles/mlog.dir/build.make
lib/libmlog.a: CMakeFiles/mlog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library lib/libmlog.a"
	$(CMAKE_COMMAND) -P CMakeFiles/mlog.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mlog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mlog.dir/build: lib/libmlog.a

.PHONY : CMakeFiles/mlog.dir/build

CMakeFiles/mlog.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mlog.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mlog.dir/clean

CMakeFiles/mlog.dir/depend:
	cd /fred/oz025/yqin/bitbucket/meraxes/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles/mlog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mlog.dir/depend

