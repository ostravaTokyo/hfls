include ./Defs.inc

SOURCESCPP = htfeti.cpp

include ./Make.inc


all:  
	(cd src; make) 
	$(LD) $(LDOPT) -o main htfeti.o \
	src/*.o $(SHARED_COMPILER) \
	$(LIB_METIS)
	-mkdir -p build
	-mv main build/

clean:
	-rm *.o main
	(cd src; make clean)

tar:
	tar czvf hfls.tar.gz htfeti.cpp src/*.[ch]pp  Defs.inc Make.inc Makefile src/Makefile
