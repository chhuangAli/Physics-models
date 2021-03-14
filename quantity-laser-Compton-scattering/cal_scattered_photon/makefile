INCLU_HEAD=${HOME}/local/cal_hit_probability_electron_laser/myHeadfile
INCLU_LIB=${HOME}/local/cal_hit_probability_electron_laser/myLibrary

SRC2=cal_scat_photon_headon.cpp
SRC=cal_num_scat_photons.cpp
OBJ=cal_diff_Compton_scat.cpp


OBJS = $(SRC:.c=.o)
OBJS += $(OBJ:.c=.o)	

GOrun: $(SRC)
	g++ -c -o cal_diff_Compton_scat.o $(OBJ) `root-config --libs --cflags` -I${INCLU_HEAD} -L${INCLU_LIB} -lBeam
	g++ -o GOrun $(OBJS) `root-config --libs --cflags` -I${INCLU_HEAD} -L${INCLU_LIB} -lBeam
GOrun2: $(SRC2)
	g++ -o GOrun $(SRC2) `root-config --libs --cflags` -I${INCLU_HEAD} -L${INCLU_LIB} -lBeam
clean:
	rm GOrun
