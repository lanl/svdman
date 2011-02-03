CPP = g++
CFLAGS = -O2 -I./F2CLIBS

all: svdman

clean:
	rm *.o *.a svdman makepng F2CLIBS/*.o CLAPACK/*.o

svdman: svdman.o svdtool.o libclasvd.a libF77.a svd.h
	$(CPP) $(CFLAGS) -o svdman svdman.o svdtool.o libclasvd.a libF77.a -lm

svdtool.o:svdtool.cc svd.h
	$(CPP) $(CFLAGS) -c svdtool.cc

svdman.o:svdman.cc svd.h
	$(CPP) $(CFLAGS) -c svdman.cc 

libclasvd.a:
	cd CLAPACK; make

libF77.a:
	cd F2CLIBS; make

