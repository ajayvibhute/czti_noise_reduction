CPP=g++
FLAGS=-w -c -Wall
INCLUDE=-I../include -I../../include
SOURCES= utils.cpp Exposure.cpp  ExpMap.cpp alglibinternal.cpp  alglibmisc.cpp  ap.cpp linalg.cpp cztstring.cpp cztmatrix.cpp validations.cpp jpeg_handling.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
MODULELIBDIR=../lib
MODULELIB=$(MODULELIBDIR)/libcommon.a
LIB=-L../../lib -lcfitsio -lpil -lm -lglog -lpthread -ljpeg
RM=rm
MK=mkdir -p

all:	$(SOURCES) $(MODULELIB)

$(MODULELIB): $(OBJECTS)
	$(MK) $(MODULELIBDIR)
	ar cr $(MODULELIB) $(OBJECTS) 

.cpp.o: 
	$(CPP) $(FLAGS) $< -o $@ $(INCLUDE) 

clean:
	$(RM) $(OBJECTS)
