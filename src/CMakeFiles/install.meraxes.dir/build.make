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

# Utility rule file for install.meraxes.

# Include the progress variables for this target.
include CMakeFiles/install.meraxes.dir/progress.make

CMakeFiles/install.meraxes:
	/fred/oz113/yqin/3rd_party/lib/python3.6/site-packages/cmake/data/bin/cmake -DBUILD_TYPE= -DCOMPONENT=bin -P /fred/oz025/yqin/bitbucket/meraxes/src/cmake_install.cmake

install.meraxes: CMakeFiles/install.meraxes
install.meraxes: CMakeFiles/install.meraxes.dir/build.make

.PHONY : install.meraxes

# Rule to build all files generated by this target.
CMakeFiles/install.meraxes.dir/build: install.meraxes

.PHONY : CMakeFiles/install.meraxes.dir/build

CMakeFiles/install.meraxes.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/install.meraxes.dir/cmake_clean.cmake
.PHONY : CMakeFiles/install.meraxes.dir/clean

CMakeFiles/install.meraxes.dir/depend:
	cd /fred/oz025/yqin/bitbucket/meraxes/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles/install.meraxes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/install.meraxes.dir/depend

