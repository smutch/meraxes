.PHONY: default
default: build  ## Make the code. Pass `make` flags with MFLAGS.
	cmake --build build -- $(MFLAGS)

build:  ## Configure the build. Pass `cmake` flags with CMAKE_OPTS.
	cmake -B build -S . ${CMAKE_OPTS} ${CMAKE_BUILD_TYPE}
	ln -sf build/compile_commands.json

debug: CMAKE_BUILD_TYPE = -DCMAKE_BUILD_TYPE=Debug
debug: build default  ## Configure and build with -DCMAKE_BUILD_TYPE=Debug

.PHONY: clean
clean:  ## Clean up the build
	rm -rf build
	rm -f compile_commands.json

help:  ## Prints help for targets with comments
	@cat $(MAKEFILE_LIST) | grep -E '^[a-zA-Z_-]+:.*?## .*$$' | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
