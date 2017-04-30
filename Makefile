include ./Defs.inc

SOURCESCPP = htfeti.cpp

include ./Make.inc


all:  
	(cd src; make) 
	$(LD) $(LDOPT) -o main htfeti.o \
	src/*.o $(SHARED_COMPILER)
	-mv main build/

clean:
	-rm *.o main
	(cd src; make clean)
