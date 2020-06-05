BIN=${CZTNOISECLEAN}/bin/
SRC=${CZTNOISECLEAN}/src/
CFITSIO=-lcfitsio
compile_gcc=-I${as1czt}/include/ -L${as1czt}/lib/ ${CFITSIO}  -I${as1czt}/lib/pil  -lpil -lm

all:
	gcc -g ${SRC}cztdataselection.c $(compile_gcc) -o ${BIN}cztdataselection
	gcc -g ${SRC}cztnoisypixclean.c $(compile_gcc) -o ${BIN}cztnoisypixclean 
	gcc -g ${SRC}cztsuperbunchclean.c $(compile_gcc) -o ${BIN}cztsuperbunchclean 
	gcc -g ${SRC}cztheavybunchclean.c $(compile_gcc) -o ${BIN}cztheavybunchclean 
	gcc -g ${SRC}cztflickpixclean.c $(compile_gcc) -o ${BIN}cztflickpixclean
	gcc -g ${SRC}czteventsep.c $(compile_gcc) -o ${BIN}czteventsep
