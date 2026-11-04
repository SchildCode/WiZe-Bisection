FC = gfortran
FFLAGS = -O3 -Wall -Wextra 
all: run_me # Standardmål når du skriver bare 'make'
run_me: src/main.f90 # Regel for å bygge den kjørbare filen direkte fra hovedprogrammet
	$(FC) $(FFLAGS) src/main.f90 -o run_me
clean: # Slett midlertidige filer vedmake 'make clean'
	rm -f *.o *.mod run_me run_me.exe