implicit:
	gfortran inverse.f90 1d_heat_implicit_euler.f90
	./a.out
explicit:
	gfortran 1d_heat_explicit_euler.f90
	./a.out
clean:
	rm a.out *~
