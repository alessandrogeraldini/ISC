# Makefile for ISC

MAG: main.o magfield.o linetodata.o RungeKutta.o findcentre.o iota_Poincare.o linalg2x2.o magfield_Dommaschk.o
	cc -O2 -Wall -ffast-math -lgsl -o MAG main.o magfield.o linetodata.o RungeKutta.o findcentre.o iota_Poincare.o linalg2x2.o magfield_Dommaschk.o

main.o: main.c 
	cc -g -Wall -I/opt/local/include/gsl -c -o main.o main.c

magfield.o: magfield.c 
	cc -g -Wall -I/opt/local/include/gsl -c -o magfield.o magfield.c

magfield_Dommaschk.o: magfield_Dommaschk.c 
	cc -g -Wall -I/opt/local/include/gsl -c -o magfield_Dommaschk.o magfield_Dommaschk.c

linetodata.o: linetodata.c
	cc -g -Wall -I/opt/local/include/gsl -c -o linetodata.o linetodata.c

RungeKutta.o: RungeKutta.c
	cc -Wall -I/opt/local/include/gsl -c -o RungeKutta.o RungeKutta.c

findcentre.o: findcentre.c
	cc -Wall -I/opt/local/include/gsl -c -o findcentre.o findcentre.c

iota_Poincare.o: iota_Poincare.c
	cc -Wall -I/opt/local/include/gsl -c -o iota_Poincare.o iota_Poincare.c

linalg2x2.o: linalg2x2.c
	cc -Wall -I/opt/local/include/gsl -c -o linalg2x2.o linalg2x2.c
