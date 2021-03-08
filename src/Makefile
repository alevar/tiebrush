INC := ../include
INCDIRS := -I. -I${INC}
GDIR := ${INC}/gclib

CXX  := $(if $(CXX),$(CXX),g++)

CXXFLAGS  := -std=c++11 -g ${INCDIRS}
LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

HTSLIB := ${INC}/htslib/libhts.a
BWIGLIB := ${INC}/libBigWig/libBigWig.a

LIBS    :=  -lz -llzma -lbz2 -pthread 

DMACH := $(shell ${CXX} -dumpmachine)

#ifneq (, $(findstring linux, $(DMACH)))
#  # these static libraries must be installed for this to work:
#  LDFLAGS += -static-libgcc -static-libstdc++
#endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o \
   ./GSam.o 
TBOBJS := ${OBJS} ./tmerge.o 

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

.PHONY : all

all static : tiebrush tiecov

tiebrush: $(TBOBJS) tiebrush.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${HTSLIB} ${LIBS}

tiecov: $(OBJS) tiecov.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${HTSLIB} ${BWIGLIB} ${LIBS} -lcurl



.PHONY : clean 

clean:
	${RM} tiebrush tiecov tiecov.o tiebrush.o $(TBOBJS)
	${RM} core.*
