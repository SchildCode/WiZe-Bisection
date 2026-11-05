FC = gfortran
FFLAGS = -O3 -Wall -Wextra -ffree-line-length-none

# Disable built-in rules that could treat source files as outputs.
MAKEFLAGS += --no-builtin-rules
.SUFFIXES:

ifeq ($(OS),Windows_NT)
    EXE = bin/run_me.exe
else
    EXE = bin/run_me
endif

.PHONY: all test clean
all: $(EXE)

$(EXE): bin/main.o
	$(FC) $(FFLAGS) $< -o $@

bin/main.o: src/main.f90 | bin
	$(FC) $(FFLAGS) -Jbin -Ibin -c $< -o $@

bin:
	mkdir -p bin

test: $(EXE)
	./$(EXE)

clean:
	rm -f bin/*.o bin/*.mod bin/run_me bin/run_me.exe