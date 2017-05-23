#-*- Makefile -*-
## Choose compiler: gfortran,ifort (g77 not supported, F90 constructs in use!)
COMPILER=gfortran
FC=$(COMPILER)
## Choose analysis
# NNLO, pheno, pheno_2, release
ANALYSIS=release
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77=gfortran -fno-automatic
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -fbounds-check #-Wimplicit-interface 
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow #,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb 
endif

ifeq ("$(COMPILER)","ifort")
F77 = ifort
#CXX = icpc
#LIBS = -limf
#FFLAGS =  -check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
endif

ifdef DEBUG
#FPE=-ffpe-trap=invalid,zero,overflow #,underflow
OPT=-O0
endif

ifeq ("$(COMPILER)","gfortran")
FFLAGS +=-ffixed-line-length-132
else
FFLAGS +=-extend-source
endif


PWD=$(shell pwd)
OBJ=$(PWD)/obj-$(COMPILER)

# The following assumed that NNLOPS is a subdirectory of HJ 
# if it is not then this ath needs to be adjusted 
BOXmaindir=$(PWD)/../../../
HJdir=$(BOXmaindir)/HJ/

#Needed only for V1
VPATH=$(PWD):$(BOXmaindir)/Version-pre2-1/progress/:$(BOXmaindir)/Version-pre2-1/:$(BOXmaindir):$(BOXmaindir)/HJ:$(OBJ)


INCLUDE0=$(BOXmaindir)/Version-pre2-1/include/
INCLUDE1=$(HJdir)/
INCLUDE2=$(BOXmaindir)/include/

FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2)  

INCLUDE=$(INCLUDE0)*.h $(INCLUDE1)*.h $(INCLUDE2)*.h 


# ifeq ("$(ANALYSIS)","pheno")
# ##To include Fastjet configuration uncomment the following lines. 
# FASTJET_CONFIG=$(shell which fastjet-config)
# LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
# FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
# PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-pheno.o jet_finder.o fastjetfortran.o boostrot.o
# endif

ifeq ("$(ANALYSIS)","pheno_2")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-pheno_2.o jet_finder.o fastjetfortran.o boostrot.o
endif

ifeq ("$(ANALYSIS)","release")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-release.o jet_finder-release.o fastjetfortran.o boostrot.o
endif



INCLUDE =$(wildcard $(INCLUDE0)/*.h $(INCLUDE1)/*.h $(INCLUDE2)/*.h)

%.o: %.f $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.F $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)


MINNLO= minnlo.o \
        cernroutines.o auxiliary.o opencount.o powheginput.o \
        lhefread.o lhef_readwrite.o newunit.o \
        get_hdamp.o \
        genclust_kt.o swapjet.o miscclust.o ptyrap.o r.o

# cernroutines: contains pwhg_exit, needed by pwhg_bookhist-multi.f
# auxiliary: needed by minnlo.f
# opencount: needed by minnlo.f
# powheginput: needed by cernroutines.f and pwhg_analysis_driver.o
# lhefread / lhef_readwrite: needed by minnlo.f
# newunit: needed by cernroutines.f and pwhg_bookhist-multi.f
# jet_finder / get_hdamp: needed by minnlo.f and by analysis
# genclust_kt: needed by jet_finder.f
# swapjet / miscclust: needed by genclust_kt.f
# ptyrap: needed in several places
# r: needed by miscclust.f

MINNLO+=$(PWHGANAL)

minnlo: $(MINNLO) 
	$(FF) $(patsubst %.o,$(OBJ)/%.o,$(MINNLO)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@


# target to cleanup
.PHONY: clean veryclean $(OBJ)

$(OBJ):
	if ! [ -d $(OBJ) ] ; then mkdir $(OBJ) ; fi

clean:
	rm -f minnlo ; \
        cd $(OBJ) ; rm -f *.o 


veryclean:
	cd $(OBJ) ; \rm *


