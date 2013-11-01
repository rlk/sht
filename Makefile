CC=   $(firstword $(wildcard /usr/local/bin/gcc /usr/bin/gcc))
CPPC= $(firstword $(wildcard /usr/local/bin/g++ /usr/bin/g++))

OPTS= -Wall -g -m64 -Wno-unknown-pragmas
#OPTS= -Wall -O3 -fopenmp -m64 -fstrict-aliasing -fno-exceptions -fno-rtti -march=athlon64

UTIL= ../util3d

LIBS= -L/usr/local/lib -lpng -ljpeg -ltiff
INCS= -I$(UTIL)

#-------------------------------------------------------------------------------

all : shtrans shimage sherror shascii

shtrans : shtrans.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

sherror : sherror.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

shimage : shimage.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

shascii : shascii.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

#-------------------------------------------------------------------------------

shtrans.o : Makefile sh.hpp
sherror.o : Makefile sh.hpp
shimage.o : Makefile sh.hpp
shascii.o : Makefile sh.hpp

#-------------------------------------------------------------------------------

%.o : %.cpp
	$(CPPC) $(OPTS) $(INCS) -o $@ -c $<

%.o : %.c
	$(CC) $(OPTS) $(INCS) -o $@ -c $<

clean :
	rm -f shtrans sherror shimage *.o

#-------------------------------------------------------------------------------

# Test to confirm that round trips at several frequencies.

test0 : sherror
	./sherror -n 8
	./sherror -n 16
	./sherror -n 32
	./sherror -n 64
	./sherror -n 128
	./sherror -n 256

# Test to confirm that the a round trip reproduces the input.

test1 : shtrans
	./shtrans -i -osyn.tif test.tif
	./shtrans    -oana.tif syn.tif
