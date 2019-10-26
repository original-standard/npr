# $Id: Makefile.users.in,v 1.6 2011/03/24 16:20:52 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------




#
# include list for the Columbia code
#
INCLIST = -I./include 

CXX=g++
CFLAGS= -g -O2 -O2 -msse -msse2 -msse3
CXXFLAGS= -g -O2 -msse -msse2 -msse3 ${INCLIST}
ASFLAGS= 

#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#


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
