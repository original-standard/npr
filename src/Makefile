INCLIST = -I../include  
CXX=g++
CFLAGS=  -O2 -O2 -msse -msse2 -msse3
CXXFLAGS=   -O2 -msse -msse2 -msse3 -D_X_SITES=64 -D_Y_SITES=64 -D_Z_SITES=64 -D_T_SITES=64 -DWILSON
ASFLAGS= 



LIBLIST =\

#
#  targets
#



.SUFFIXES:
.SUFFIXES:.o .C .S .c

CCSRC :=$(wildcard *.C)
SSRC :=$(wildcard *.S)

CCOBJ=$(CCSRC:.C=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))


all:$(OBJS)

.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS)  $(INCLIST) -c $<

clean:
	rm -f *.dat *.o */*.o  $(BIN)
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog
