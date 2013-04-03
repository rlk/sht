UTIL= ../util3d

CC=   /usr/local/bin/gcc
CPPC= /usr/local/bin/g++
OPTS= -Wall -O3 -fopenmp -m64 -fstrict-aliasing -fno-exceptions -fno-rtti -march=athlon64
LIBS= -lpng -ljpeg -ltiff
INCS= -I$(UTIL)

#-------------------------------------------------------------------------------

all : shtrans shimage sherror

shtrans : shtrans.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

sherror : sherror.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

shimage : shimage.o $(UTIL)/image.o
	$(CPPC) $(OPTS) $(INCS) -o $@ $^ $(LIBS)

#-------------------------------------------------------------------------------

shtrans.o : Makefile sh.hpp
sherror.o : Makefile sh.hpp
shimage.o : Makefile sh.hpp

#-------------------------------------------------------------------------------

%.o : %.cpp
	$(CPPC) $(OPTS) $(INCS) -o $@ -c $<

%.o : %.c
	$(CC) $(OPTS) $(INCS) -o $@ -c $<

clean :
	rm -f shtrans sherror shimage *.o

