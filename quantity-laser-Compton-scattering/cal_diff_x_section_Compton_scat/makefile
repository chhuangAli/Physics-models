INCLU_HEAD=${HOME}/local/cal_hit_probability_electron_laser/myHeadfile
INCLU_LIB=${HOME}/local/cal_hit_probability_electron_laser/myLibrary

SRC=cal_diff_Compton_scat.cpp

OBJS = $(SRC:.c=.o)

GOrun: $(SRC)
#	g++ -c -o beam.o ${INCLU_LIB}/$(FUNC) `root-config --libs --cflags` -I${INCLU_HEAD} -L${INCLU_LIB}
	g++ -o GOrun $(SRC) `root-config --libs --cflags` -I${INCLU_HEAD} -L${INCLU_LIB} -lBeam
clean:
	rm GOrun
