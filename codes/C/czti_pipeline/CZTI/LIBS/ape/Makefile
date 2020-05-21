
INSTALLPATH=../../
INCLUDE=ape cfortran.h
INCLUDEPATH= -I .
OST :=$(shell uname)
OS=$(strip $(OST))
# $(shell uname)	
all : libape_2.so
	cp -r ${INCLUDE} ${INSTALLPATH}/include
	cp pquery2 ${INSTALLPATH}/czti/bin
libape_2.so : 
ifeq ($(OS),Darwin)
	gcc -fPIC -dynamiclib  -o libape_2.8.dylib *.c -I .
	gcc pquery2.c -I "." -o pquery2 -L . -lape_2.8
	cp libape_2.8.dylib ${INSTALLPATH}/lib
else
	gcc -fPIC -shared  -o libape_2.8.so *.c -I . 
	gcc pquery2.c -I "." -o pquery2 -L . -lape_2.8
	cp libape_2.8.so ${INSTALLPATH}/lib
endif

clean:
	rm ape*.o  libape_2.8.so
#objects :

#OS := $(shell uname)
#all:
#	gcc  -c -fPIC *.c -I .

#ifeq ($(OS),Darwin)
#	gcc  -c -fPIC *.c -I .
#else
#	gcc  -c -fPIC *.c -I . 
#	gcc -shared -Wl,-soname,libape_2.8.so -o libape_2.8.so *.o
#	gcc -fPIC -dynamiclib  -o libape_2.8.dylib *.c
  	# check for Linux and run other commands
#endif




#pil_error.o : pil_error.c ${INCLUDE}
#	gcc -fPIC  -c pil_error.c ${INCLUDE}

#pil_sys.o : pil_sys.c
#	gcc -fPIC  -c pil_sys.c ${INCLUDE}

#clean:
#	rm pil.o pil_error.o pil_sys.o libpil.a
