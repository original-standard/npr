INCLIST = -I../include  
CXX=g++
CFLAGS= -g -O2 -O2 -msse -msse2 -msse3
CXXFLAGS=  -g -O2 -msse -msse2 -msse3
ASFLAGS= 
DFLAGS +=  -DUSE_SSE -DUSE_QMP -DUSE_QIO



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
