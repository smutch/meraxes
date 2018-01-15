.PHONY: default config cleanbuild

ifneq ("$(wildcard local.cmake)", "")
INIT_CACHE = "-Clocal.cmake"
endif

default:
	@if [ ! -d "build" ]; then echo $(INIT_CACHE); cmake . -Bbuild $(INIT_CACHE) $(CMAKE_OPTS); fi
	@make -C build --no-print-directory $(MFLAGS)

config:
	@cmake . -Bbuild $(CMAKE_OPTS)

cleanbuild:
	@rm -r build

% :: build/Makefile
	@make -C build --no-print-directory $@ $(MFLAGS)
