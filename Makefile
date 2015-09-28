default:
	gfortran parabolic.f90
	./a.out
clean:
	rm a.out *~
