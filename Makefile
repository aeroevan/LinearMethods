FC=mpif90 -Wall -Wextra -fbounds-check -O3 -mtune=native -march=native
#FC=gfortran -Wall -Wextra -fbounds-check -O3 -mtune=native -march=native
#FC=mpif90 -fopenmp -O3 -mtune=native -march=native

PROGRAMS=project prob10 prob18 test_lu test_qr


all: $(PROGRAMS)

clean:
	rm -f *.o *.mod
	rm -f $(PROGRAMS)

test_lu: test_lu.o linearmethods.o
	$(FC) -o test_lu test_lu.o linearmethods.o

test_lu.o: test_lu.f90
	$(FC) -c test_lu.f90

test_qr: test_qr.o linearmethods.o
	$(FC) -o test_qr test_qr.o linearmethods.o

test_qr.o: test_qr.f90
	$(FC) -c test_qr.f90

prob18: prob18.o linearmethods.o
	$(FC) -o prob18 prob18.o linearmethods.o

prob18.o: prob18.f90
	$(FC) -c prob18.f90

linearmethods.o: linearmethods.f90
	$(FC) -c linearmethods.f90

prob10.o: prob10.f90 linearmethods.o
	$(FC) -c prob10.f90

prob10: prob10.o
	$(FC) -o prob10 prob10.o linearmethods.o

project.o: project.f90 linearmethods.o
	$(FC) -c project.f90

project:project.o linearmethods.o
	$(FC) -o project project.o linearmethods.o
