.PHONY: default config cleanbuild

default:
	@if [ ! -d "build" ]; then cmake . -Bbuild $(CMAKE_OPTS); fi
	@+make -C build --no-print-directory $(MFLAGS)

config:
	@cmake . -Bbuild $(CMAKE_OPTS)

cleanbuild:
	@rm -r build

% :: build/Makefile
	@+make -C build --no-print-directory $@ $(MFLAGS)
