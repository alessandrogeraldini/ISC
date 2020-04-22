# Makefile for ISC
# Note: -g adds debug symbols.

# Alessandro's laptop:
EXTRA_COMPILE_FLAGS = -g -O2 -Wall -I /opt/local/include/gsl
EXTRA_LINK_FLAGS = -Wall -ffast-math -lgsl

#EXTRA_COMPILE_FLAGS = -g -O2 -Wall -I /opt/local/include
#EXTRA_LINK_FLAGS = -Wall -ffast-math -L /opt/local/lib -lgsl

MAG: main.o magfield.o linetodata.o RungeKutta.o findcentre.o iota_Poincare.o linalg2x2.o magfield_Dommaschk.o
	cc $(EXTRA_LINK_FLAGS) -o MAG main.o magfield.o linetodata.o RungeKutta.o findcentre.o iota_Poincare.o linalg2x2.o magfield_Dommaschk.o

main.o: main.c 
	cc $(EXTRA_COMPILE_FLAGS) -c -o main.o main.c

magfield.o: magfield.c 
	cc $(EXTRA_COMPILE_FLAGS) -c -o magfield.o magfield.c

magfield_Dommaschk.o: magfield_Dommaschk.c 
	cc $(EXTRA_COMPILE_FLAGS) -c -o magfield_Dommaschk.o magfield_Dommaschk.c

linetodata.o: linetodata.c
	cc $(EXTRA_COMPILE_FLAGS) -c -o linetodata.o linetodata.c

RungeKutta.o: RungeKutta.c
	cc $(EXTRA_COMPILE_FLAGS) -c -o RungeKutta.o RungeKutta.c

findcentre.o: findcentre.c
	cc $(EXTRA_COMPILE_FLAGS) -c -o findcentre.o findcentre.c

iota_Poincare.o: iota_Poincare.c
	cc $(EXTRA_COMPILE_FLAGS) -c -o iota_Poincare.o iota_Poincare.c

linalg2x2.o: linalg2x2.c
	cc $(EXTRA_COMPILE_FLAGS) -c -o linalg2x2.o linalg2x2.c
