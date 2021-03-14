INCLU_HEAD=${HOME}/local/cal_hit_probability_electron_laser/myHeadfile

libBeam.a:
	g++ -c -o beam.o beam.cpp `root-config --libs --cflags` -I${INCLU_HEAD}
	ar rcs libBeam.a beam.o
