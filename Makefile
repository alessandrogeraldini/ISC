MAG: main.o magfield.o extract_coils.o linetodata.o RungeKutta.o
	gcc -O2 -Wall -ffast-math -o MAG main.o magfield.o extract_coils.o linetodata.o RungeKutta.o

main.o: main.c 
	gcc -Wall -I/opt/local/include/gsl -c -o main.o main.c

magfield.o: magfield.c 
	gcc -Wall -I/opt/local/include/gsl -c -o magfield.o magfield.c

linetodata.o: linetodata.c
	gcc -Wall -I/opt/local/include/gsl -c -o linetodata.o linetodata.c

RungeKutta4.o: RungeKutta.c
	gcc -Wall -I/opt/local/include/gsl -c -o RungeKutta.o RungeKutta.c
