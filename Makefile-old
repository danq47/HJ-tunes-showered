#-*- Makefile -*-
## Choose compiler: gfortran,ifort (g77 not supported, F90 constructs in use!)
COMPILER=gfortran
FC=$(COMPILER)
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
## ANALYSIS: none, BnJ, HnJ, alt, NNLOPS
##           the first 4 analyses require the FASTJET package, that has to be 
##           installed separately (see below)
ANALYSIS=NZ
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77=gfortran -fno-automatic -ffixed-line-length-none
# -finit-real=nan -ffpe-trap=invalid,zero,overflow
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif

ifeq ("$(COMPILER)","g77")
F77= g77 -fno-automatic 
#FFLAGS= -Wall -ffortran-bounds-check
## For floating point exception trapping  uncomment the following 
#FPEOBJ=trapfpe.o
OPT=-O3
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif


ifeq ("$(COMPILER)","ifort")
F77 = ifort
#CXX = icpc
#LIBS = -limf
#FFLAGS =  -check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
#OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
endif

ifdef DEBUG
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
OPT=-O0
endif

PWD=$(shell pwd)
OBJ=$(PWD)/obj-$(COMPILER)
WDNAME=$(shell basename $(PWD))
HJ=$(PWD)
VPATH= ./:../:./Madlib/:./MODEL/:./DHELAS/:./Topmass/:$(OBJ)/


INCLUDE1=$(HJ)
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=$(HJ)/MCFM_Include 
INCLUDE4=$(HJ)/Topmass 
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) -I$(INCLUDE4)

LIBFILES=$(shell  for dir in $(HJ)/Madlib $(HJ)/MODEL $(HJ)/DHELAS ; do cd $$dir ; echo *.[fF] ' ' | sed 's/[fF] /o /g' ; cd .. ; done  )

INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
#LHAPDF_CONFIG=~/Pheno/PDFpacks/lhapdf-5.8.4-$(FC)/bin/lhapdf-config
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
# LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

ifeq ("$(ANALYSIS)","HnJ")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-HnJ.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","NZ")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-NZ.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","tmp")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-tmp.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","fintest")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-fintest.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif



ifeq ("$(ANALYSIS)","BnJ")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-BnJ.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","alt")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-alt-H2J.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","NNLOPS")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-pheno_2.o fastjetfortran.o \
         genclust_kt.o miscclust.o ptyrap.o r.o swapjet.o jet_finder.o  \
         auxiliary.o get_hdamp.o
endif

%.o: %.f $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.F $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o\
     real.o matchcolour.o $(PWHGANAL)

MCFMFILES=i2mcfm_2_POWHEG.o set_interface_MCFM.o i2mcfm_fill.o \
	gg_hg_vi.o gg_hg_eval_v.o hjetfill.o dotem.o set_epinv.o coupling.o 


PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	reshufflemoms.o setlocalscales.o bmass_in_minlo.o validflav.o \
        mint_upb.o minlo_checks.o pwhgreweight.o opencount.o \
        ubprojections.o sigcollsoft.o sigvirtual.o \
        Born_tm.o  bornmatrix_tm.o  init_couplings-topmass.o nloutils.o   \
        ew.o EWgint.o decay-util.o HPL_full.o \
        $(PDFPACK) $(USER)  $(FPEOBJ) $(MCFMFILES) $(H2JVIRTFILES) $(H2JVFILES)

#LIBDIRMG=.
#LINKMGLIBS =  -L$(LIBDIRMG)  -lmadgraph -lmodel -ldhelas3 

#MADLIBS=libdhelas3.a libmadgraph.a libmodel.a

# Get SVN info for SVN version stamping code
$(shell ../svnversion/svnversion.sh>/dev/null)

# target to generate LHEF output
pwhg_main:$(PWHG) libfiles.a
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(OBJ)/libfiles.a $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@


NNLOPSREWEIGHTER=nnlopsreweighter.o cernroutines.o powheginput.o lhefread.o \
                 newunit.o boostrot.o
NNLOPSREWEIGHTER+=pwhg_bookhist-multi.o fastjetfortran.o

nnlopsreweighter: $(NNLOPSREWEIGHTER) 
	$(FF) $(patsubst %.o,$(OBJ)/%.o,$(NNLOPSREWEIGHTER)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

libfiles.a: $(LIBFILES)
	cd $(OBJ) ; \rm libfiles.a ; ar cru libfiles.a $(LIBFILES)

LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o $(PWHGANAL)	\
     lhefread.o newunit.o pwhg_analysis_driver.o bra_ket_subroutines.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o	\
	powheginput.o $(PWHGANAL) lhefread.o	\
	pdfdummies.o opencount.o $(FPEOBJ) 

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o		\
	$(PWHGANAL) lhefread.o newunit.o pdfdummies.o  bra_ket_subroutines.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBSFASTJET)  $(STATIC) -o $@


DELTASUD=integral-deltasud.o pdfcalls.o cernroutines.o newunit.o \
         pwhg_bookhist-multi.o random.o  $(PDFPACK)

deltasud: $(DELTASUD)
	$(FF)  $(patsubst %,$(OBJ)/%,$(DELTASUD)) $(LIBS) -o $@
# target to cleanup
.PHONY: clean veryclean $(OBJ)

$(OBJ):
	if ! [ -d $(OBJ) ] ; then mkdir $(OBJ) ; fi

clean:
	rm -f pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA-lhef ; \
        cd $(OBJ) ; rm -f $(PWHG) $(HERWIG) $(PYTHIA) $(LHEF) $(NNLOPSREWEIGHTER) 


veryclean:
	cd $(OBJ) ; \rm *

# Dependencies of SVN version stamp code
pwhg_main.o: svn.version
lhefwrite.o: svn.version

ifeq ("$(COMPILER)","gfortran")
XFFLAGS +=-ffixed-line-length-132
else
XFFLAGS +=-extend-source
endif
