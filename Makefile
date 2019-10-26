
INCLIST = -I./include 

CXX=g++
CFLAGS= -O2 -msse -msse2 -msse3
CXXFLAGS= -O2 -msse -msse2 -msse3 ${INCLIST}
ASFLAGS= 



LIBLIST = $(wildcard src/*.C)

#
#  targets
#

MAIN :=$(wildcard *.C)

all:$(MAIN:.C=.x)

.SUFFIXES:
.SUFFIXES:.cpp .x .o .C .S .c

CCSRC :=$(wildcard src/*.C)
SSRC :=$(wildcard *.S)

CCOBJ=$(CCSRC:.C=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))



.C.x:$(CCSRC:.C=.o)
	make -C src
	$(CXX) -c $^ $(CXXFLAGS) $(LDFLAGS) $(INCLIST)
	@echo OBJS = $(OBJS_SRC)
	$(CXX) $(^:.C=.o) $(OBJS_SRC) $(CXXFLAGS)  $(LDFLAGS) -o $@

clean:
	rm -f *.dat *.o */*.o  $(BIN)
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog
	rm *.x
