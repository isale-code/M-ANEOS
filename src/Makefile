FC = gfortran
FCFLAGS = -O3 
AR = ar
ARFLAGS = -cru
LDFLAGS = 
LIBRARY = libaneos.a
PROGRAM = aneos
PROGOBJECTS = ANEOSTEST.o EOS.o
LIBOBJECTS = ANDATA.o \
	ANDEBY.o \
	ANE2PH.o \
	ANEHPP.o \
	ANEI3.o \
	ANEINI.o \
	ANELSM.o \
	ANEOS1.o \
	ANEOS2.o \
	ANEOSD.o \
	ANEOS.o \
	ANEOSI.o \
	ANEOSS.o \
	ANEOSV.o \
	ANEVAL.o \
	ANHUG.o \
	ANION1.o \
	ANION2.o \
	ANMARK.o \
	ANMAXW.o \
	ANN1AS.o \
	ANN1VL.o \
	ANN2AS.o \
	ANNDPR.o \
	ANPHAS.o \
	ANPHTR.o \
	ANPRTR.o \
	ANSIOF.o \
	ANSMFT.o \
	ANTOMF.o \
	ANUEOS.o \
	ANUSET.o \
	ANWARN.o \
	ANZRTR.o

%.o: %.f
	$(FC) -c $(FCFLAGS) $<

default: $(PROGRAM)

$(LIBRARY): $(LIBOBJECTS) 
	rm -f $(LIBRARY)
	$(AR) -cr $(LIBRARY) $(LIBOBJECTS)

$(PROGRAM): $(PROGOBJECTS) $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(PROGRAM) $(PROGOBJECTS) $(LIBRARY)

clean:
	rm -f $(LIBOBJECTS) $(PROGOBJECTS)

allclean:
	rm -f $(PROGRAM) $(LIBRARY)
	make clean
