CMAKE_BUILD_TYPE ?= Debug
CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE)
CMAKE_EXTRA_FLAGS ?=
CMAKE_GENERATOR ?= $(shell (type ninja > /dev/null 2>&1 && echo "Ninja") || \
    echo "Unix Makefiles")

ifeq (,$(BUILD_TOOL))
  ifeq (Ninja,$(CMAKE_GENERATOR))
    ifneq ($(shell cmake --help 2>/dev/null | grep Ninja),)
      BUILD_TOOL := ninja
    else
      # User's version of CMake doesn't support Ninja
      BUILD_TOOL = $(MAKE) --no-print-directory
      CMAKE_GENERATOR := Unix Makefiles
    endif
  else
    BUILD_TOOL = $(MAKE) --no-print-directory
  endif
endif

ifneq ($(VERBOSE),)
  # Only need to handle Ninja here.  Make will inherit the VERBOSE variable.
  ifeq ($(CMAKE_GENERATOR),Ninja)
    VERBOSE_FLAG := -v
  endif
endif

BUILD_CMD = $(BUILD_TOOL) $(VERBOSE_FLAG)

# WEIRD MAKEFILE THING HERE...
# This target **MUST** exist and cannot be called 'default'!
init: all

all: build/.ran-cmake
	$(BUILD_CMD) -C build

build/.ran-cmake:
	mkdir -p build && cd build && cmake .. -G '$(CMAKE_GENERATOR)' $(CMAKE_FLAGS) $(CMAKE_EXTRA_FLAGS)
	touch $@

cmake:
	rm -f build/.ran-cmake
	$(MAKE) build/.ran-cmake

cleandist:
	rm -rf build

check: build/compile_commands.json
	mkdir -p cppcheck_report
	cppcheck --enable=all --project=$^ --xml-version=2 2> cppcheck_report/cppcheck_report.xml
	/home/smutch/3rd_party/cppcheck/src/htmlreport/cppcheck-htmlreport --title=Meraxes --source-dir=./ --file=cppcheck_report/cppcheck_report.xml --report-dir=cppcheck_report

% ::
	@$(BUILD_CMD) -C build $@

.PHONY: init all cmake cleandist check
