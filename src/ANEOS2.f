C
C     start of aneos set up routines
C
      SUBROUTINE ANEOS2 (IGK,NUM,ITAPE,IZETL)
C
C***********************************************************************
C
C     SET UP FOR ANEOS PACKAGE
C
C     DIMENSIONS ARE SET FOR MAXMAT     EQUATIONS OF STATE
C                            30*MAXMAT  ELEMENTS
C                            100*MAXMAT TWO-PHASE BOUNDARY POINTS (VAPOUR)
C                            1000*MAXMAT TWO-PHASE BOUNDARY POINTS (MELT)
C     INPUTS
C
C       IGK = 1   DEFINE NEW MATERIALS
C           = 2   WRITE RESTART FILE ITAPE
C           = 3   READ RESTART FILE ITAPE
C           = 4   TURN RADIATION OPACITY CALCULATION OFF
C
C       FOR IGK = 1
C       NUM = NUMBER OF MATERIALS IN PROBLEM
C
C       IZETL = ARRAY OF NUM EOS NUMBERS
C               ONLY DEFINE HERE IF IZETL(I) < 0
C
C   IMPORTANT ARRAYS:
C     UI = CONTAINS DATA READ FROM INPUT FILE ANEOS.INPUT, LENGTH NINPUT
C     DIN = SCRATCH COPY OF UI, MODIFIED DURING SETUP
C     ACK = COEFFICIENT ARRAY TRANSMITTED TO ANEOS1 FOR COMPUTATION,
C                 LENGTH = 99 ENTRIES FOR EACH MATERIAL. CONSTRUCTED
C                 FROM DATA IN UI AND DIN.
C     ZZS = ATOMIC NUMBER OF EACH ELEMENT, 30 FOR EACH MATERIAL
C     COT = ATOMIC FRACTION OF EACH ELEMENT, 30 FOR EACH MATERIAL
C
C************************************** 4/90 version slt ***************
C************************************** 5/06 version hjm ***************
C
C   EQUATION OF STATE TYPE NUMBERS:  array element ACK(IT+30) = DIN(30) = UI(2)
C
C      0.   Solid-gas without electronic terms or detailed liquid-vapor region
C      1.   Solid-gas with electronic terms but without detailed liquid-vapor
C      2.   Gas only with electronic terms
C      3.   Same as 0., but with detailed liquid-vapor region
C      4.   Same as 1., but with detailed liquid-vapor region
C      5.   HPP of immediately preceeding input file
C     -1.   Low temperature solid with linear shock-particle velocity relation
C     -2.   Ideal gas
C     -4.   User-supplied EoS, contained in ANUEOS.f
C
C************************************** 10/10 version hjm ***************
C
C   MELT CURVE TABLE
C
C       MELT_CURVE_TABLE = 1   1000 point table storing melt curve is constructed
C                              here and stored in common /MELTTRANS/ for use in ANEOSV
C                              Constructing the melt curve in this way enables combined
C                              use of solid-solid phase transition and the melt transition
C       MELT_CURVE_TABLE = 0   Melt curve is computed using ANELSM (as before);
C                              this method can be used in conjunction with the solid-solid
C                              phase transition but is not recommended (leads to a pathology
C                              in the mixed phase region)
C
C       Note that provision has also been added for storing one (or more) solid-solid
C       phase curves as tables, but this is not yet functional.
C
C
C************************************** 09/13 version gsc/hjm ***********

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IZETL(21)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10) !maximum number of materials
      PARAMETER (MAXHPP=3)  !maximum number of hpp phases per material
      PARAMETER(NINPUT=48)  !NINPUT must be a multiple of 8!
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /MELTTRANS/ TMC(1000*MAXMAT),RSMC(1000*MAXMAT)
     1 ,RLMC(1000*MAXMAT),NUMTMC(MAXMAT),MELT_CURVE_TABLE
      COMMON /HPPTRANS/ THPP(1000*MAXMAT)
     2 ,RLPP(1000*MAXMAT),RHPP(1000*MAXMAT),NUMHPP(MAXMAT),NUMPHASE
      COMMON /ANZB/ ZB(92), ZRAT
      COMMON /FNES/ CMLT7,ACK46
      COMMON /TPMEOS/ NPORFE(MAXMAT)
      COMMON /MIXCON/  RHOMPH(MAXMAT),RHOMPL(MAXMAT),IPEOS
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     following aneos1 variables are for aneos2 edit only
      COMMON /ANEDIS/ GAMMA,PSI,THETA
      CHARACTER*50 URTIT
      DIMENSION SZ23(23),SC23(23),DIN(NINPUT),UI(NINPUT)
C
C     kexprt is flag for extended print
C     = 0, no print      =1, print
C
      COMMON /KEXARP/ KEXPRT
C
      CHARACTER*8 VERANE
      COMMON /NAMEVERS/ URTIT,VERANE
C
      COMMON /ANEST/ NOEOSU
      SAVE   /FILEOS/,/MIXCON/,/ANES/,/FNES/,/ANEDIS/,/KEXARP/,
     & /ANZB/,/TPMEOS/,/ANEST/
      DIMENSION IZETC(MAXMAT)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (FIVE=5.D0)
      PARAMETER (SIX=6.D0)
      PARAMETER (SEVEN=7.D0)
      PARAMETER (EIGHT=8.D0)
      PARAMETER (CNINE=9.D0)
      PARAMETER (TEN=10.D0)
      PARAMETER (BIGNUM=1.D50)
      PARAMETER (PIE=3.1415926536D0)
      PARAMETER (ATHIRD=1.D0/3.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (FLGJWL=0.0252525D0)
      PARAMETER (QCL1=0.95D0)       !default value for C54, eq. 3.33
      PARAMETER (QCL2=0.8D0)        !default minimum density = QCL2 * ref. density
      PARAMETER (QCL3=5.48D12)      ! 16*sigma/3, constant in eq 7.8
      PARAMETER (QCL4=144.D0)       ! = 144, constant in eq 7.15
      PARAMETER (QCL5=0.1D0)        !constant in equation 7.16
      PARAMETER (QCL6=1.D-4)        !minimum value for C42 in eq 7.15
      PARAMETER (QCL7=1.66026D-24)  ! 1/Avogadro Number
      PARAMETER (QCL8=1.5D0)        ! 3/2
      PARAMETER (QCL9=4.36050D-42)  ! hPlanck^2/(2*Pi*kBoltzman) in eV
      PARAMETER (QCL10=200.D0)      ! = 200, constant in eq 3.32
      PARAMETER (QCL11=27.D0)       ! = 27, constant in eq 3.32
      PARAMETER (QCL12=0.99D0)
      PARAMETER (QCL13=1.1D0)
      PARAMETER (QCL14=0.9D0)
      PARAMETER (QCL15=1.00001D0)
      PARAMETER (QCL16=1.005D0)
      PARAMETER (QCL17=1.D-9)       !bisection convergence for loop at 330
      PARAMETER (QCL18=6.6252D-27)  !Planck constant
      PARAMETER (QCL19=20.D0)
      PARAMETER (QCL20=9.1084D-28)
      PARAMETER (QCL21=0.9D0)
      PARAMETER (QCL22=4.80288D-10)
      PARAMETER (QCL23=18.D0)
      PARAMETER (QCL24=11.D0)
      PARAMETER (QCL25=12.D0)
      PARAMETER (QCL26=15.D0)
      PARAMETER (QCL27=1.5D0)
      PARAMETER (QCL28=0.9999D0)    !decrement factor for iterations
      PARAMETER (QCL29=10000.D0)
      PARAMETER (QCL30=0.002D0)
      PARAMETER (QCL31=1.D-10)
      PARAMETER (QCL32=1.00001D0)
      PARAMETER (QCL33=1.D-4)       !tolerance for GammaHat=Gamma0+Tgamma/3
      PARAMETER (QCL34=1.D-2)
      PARAMETER (QCL35=1.117D12)    !constant in equation 5.59 for Hf
      PARAMETER (QCL36=0.95D0)      !default for RHOliq/RHOsolid at triple pt.
      PARAMETER (QCL37=0.2D0)       !critical density estimate, fraction of RHO ref.
      PARAMETER (QCL38=0.8D0)       !upper limit to C52, upper melt region low limit
      PARAMETER (QCL39=1.02D0)      !log temperature increment for melt search
      PARAMETER (QCL40=5.D-5)       !convergence criterion for 1.-RL/RS
      PARAMETER (QCL41=100.D0)      !upper limit for melt density, QCL41*RHO0
      PARAMETER (QCL42=1.05D0)      !multiplication factor for curve A slope
      PARAMETER (QCL43=0.99D0)      !log decrement for Tmin search
      PARAMETER (QCL44=200.D0)
      PARAMETER (QCL45=1.D12)       !JWL high explosive EOS
      PARAMETER (QCL46=1.2D0)       !JWL high explosive EOS
      PARAMETER (QCL47=1.01D0)      !increment for critical temperature
      PARAMETER (QCL48=1.0001D0)    !min difference of solid/melt density for test
      PARAMETER (QCL49=0.06D0)      !log increment for melt/solid table
      PARAMETER (QCL50=0.05D0)      !log increment for cold density table
      PARAMETER (QCL51=25.D0)       !density step for solid-solid phase transition
      PARAMETER (QCL52=0.9120108D0) !log decrement for cold compression curve
      PARAMETER (QCL53=1.D-5)
      PARAMETER (QCL54=1.D-6)
      PARAMETER (QCL55=ONE-QCL54)
      PARAMETER (QCL56=-1.11D22)
CHJM
      PARAMETER (QCL57=1.195D18)       !constant for diatomic CXX
      PARAMETER (QCL58=1.602D-12)      !ergs/eV conversion
      PARAMETER (QCL59=1.574960995D1)  !constant for diatomic CXX
      PARAMETER (QCL60=1.595D18)       !constant for triatomic CXX
      PARAMETER (QCL61=1.145695D3)     !constant for triatomic CXX
CHJM
C
C     atomic numbers for high z elements - Pu to No
C
      DIMENSION BIGZAS(10)
      DATA BIGZAS / 237.122,239.127,241.132,243.137,245.143,
     &              246.146,253.164,255.170,256.174,253./
C
      DATA IZ,IKPN/1,1/
      DATA IT/0/
      DATA IEXPRT/1/  !initialize extended print flag
      CALL ANZRTR     !zero all counters
C
C     package version number
C
      VERANE='09/30/13'
      MELT_CURVE_TABLE = 0  ! Default is no melt-curve table
      NUMPHASE = 1          ! Default is only one solid phase
C
C    define some conversion factors for ease of readers:
C      tconv = converts temperature in eV to K
C      sconv = converts entropy in erg/gm-eV to J/kg-K
C      econv = converts energy in erg/gm to J/kg
C      dconv = converts density in gm/cc to Kg/m**3
C      pconv = converts pressure in dyne/cm to GPa
C
      TCONV = 1.160564D4
      SCONV = 8.617D-9
      ECONV = 1.0D-4
      DCONV = 1.0D3
      PCONV = 1.0D-10
C
      DO IQ=1,MAXMAT  !initialize porosity flags for all materials
        NPORFE(IQ)=0
      ENDDO
      IPSQTS=1  !pointer to current value of sqrt of temperature array
      IF(IGK.EQ.4) THEN  !turn off radiation opacity calculation
        CALL ANEOSS (IGK,ITAPE)  !sets NOTRAD = 0
        GO TO 2323  !return
      ELSE IF (IGK.EQ.1) THEN  !new material is being defined
        NOTRAD=1    !turn radiation opacity flag on for new materials
      END IF
      ZRAT=0.000045D0
      IPEOS=1
C
C     BRANCH FOR RESTART READ/WRITE
      GO TO (10,1240,1250), IGK
C
C=======================================================================
C     DEFINE NEW MATERIALS
C
C     LOOP FOR EACH MATERIAL *******************************************
C
   10 IAMBAD=0
      TGAM=ZERO
      IF(NUM.GT.MAXMAT) THEN
        CALL ANMARK('ANEOS2 REQUEST MORE MATERIALS THAN DIMENSIONS')
        STOP
      END IF
      DO IQ=1,NUM
        IZETC(IQ)=MIN(0,IZETL(IQ))
      ENDDO
      DO I=1,92
        S=I
        SAVER(I)=LOG(S+HALF)  !log of (atomic number + 1/2) in array SAVER
      ENDDO
C
C     start big loop
C
      DO 1220 IQ=1,NUM
      KUWAY=0
      CALL ANZRTR   !zero all counters
      IF (IZETL(IQ).GE.0) GO TO 1220   !skip to next material
      KEXPRT=IEXPRT
      RHUGMX=100000.  !maximum density for Hugoniot iteration
      NOEOSU=IZETL(IQ)
C
C     READ USER INPUT
C
      ZZS(IZ)=ZERO
      DO  JJ=1,NINPUT   !initialize UI array for input
        UI(JJ)=ZERO
      ENDDO
C
      CALL ANEOSI(KINP,KLST,ISE,URTIT,ISETAB,IZI,RHUG,THUG,UI,SZ23,SC23)
C
C      UI(NINPUT) carries data read from input file.  Write header to
C      output file (this is the first line in the file)
C
      WRITE(KLST,1440) ISE,ISETAB,IZI,VERANE,URTIT,RHUG,THUG
C
C     Check for materials out of range of definition
C
      IF(ISE.GE.0) GO TO 30
      IF(ISE.LT.-MAXMAT) GO TO 30
        DO 20 JJ=1,NUM
        IF(IZETC(JJ).EQ.ISE) THEN
           MAT=IZETC(JJ)
           IZETC(JJ)=ABS(MAT)
           NOEOSU=MAT
           GO TO 40
        END IF
   20   CONTINUE
   30 CONTINUE
      CALL ANWARN(1)
      WRITE(KLST,1450) ISE,(IZETC(JJ),JJ=1,NUM)
      CALL ANMARK('ANEOS DATA MISSING')
      IAMBAD=1
      GO TO 1220   !skip to next material
C
C     Looks OK.  Proceed with definition
C
   40 MAT=-MAT        !this is material # MAT lose the - sign
      LOCSV(MAT)=IT   !pointer to this material in ACK array
      LOCKP(MAT)=IKPN !pointer to start of vapor/liq table
C
C    Get data for UI(NINPUT) from stored table, if requested
C
      IF (ISETAB.NE.0) CALL ANDATA (IZ,ISETAB,UI)
C
   50 DO I=1,NINPUT    !initialize DIN
        DIN(I)=ZERO
      END DO
      IF (UI(4).LE.ZERO) UI(4)=298.D0 !default to 298 K if ref. temp. le 0
      DIN(28)=UI(1)        !number of elements in material
      DIN(30)=UI(2)        !type of equation of state
      DIN(11)=UI(3)        !reference density (cgs)
      DIN(12)=UI(4)/TCONV  !reference temperature, converted to eV
      DIN(20)=UI(5)        !reference pressure (cgs)
      DIN(15)=UI(7)        !reference gruneisen coefficient, eq. 4.11
      DIN(25)=UI(8)/TCONV  !reference debye temperature, converted to eV
      DIN(24)=UI(10)/THREE !constant C24, equation 4.11
      DIN(10)=UI(11)       !evaporation energy, Evap
      DIN(18)=UI(12)       !melting temperature (or energy if .lt. 0)
      IF(UI(12).GT.ZERO) DIN(18)=UI(12)/TCONV  !convert K to eV
      DIN(23)=UI(17)       !lowest allowed solid density
      DIN(1)=UI(18)        !D1               ##
      DIN(2)=UI(19)        !D2               ##
      DIN(7)=UI(20)        !D3     Parameters for solid-solid phase
      DIN(39)=UI(21)       !D4          Transition
      DIN(40)=UI(22)       !D5               ##
      ACK(IT+54)=ZERO
      ACK(IT+46)=ZERO
      ACK46=ACK(IT+46)
        DO 64 I=1,16          !copy UI(25 to 40) to ACK(55 to 70)
   64   ACK(IT+I+54)=UI(I+24) !this is data from 'extended' input format
C
C     save energy & entropy shifts - zero during set up
C
      EMOVER=UI(34)
      SMOVER=UI(35)
C
      ACK(IT+64)=ZERO
      ACK(IT+65)=ZERO
      DO 66 I=71,99     !initialize ACK 71 to 99
   66 ACK(IT+I)=ZERO
C
C    process input 6 --> B0 if positive, S0 if negative
C
      IF (UI(6)) 90,70,80
   70 IF ((NINT(ABS(DIN(30))).EQ.2).OR.
     &    (NINT(DIN(30)).EQ.-4)) GO TO 110 !don't use for types 2, -4
      IAMBAD=36        !B0 is zero--set error flag
      CALL ANWARN(1)   !warn the user
      WRITE(KLST,1260) 6,UI(6),5.E11
      UI(6)=5.E11      !put B0=5.E11 and go on
C
C     setup for use of eq 3.9 for computation of B00 and Gamma0
C
   80 DIN(21)=UI(6)    !DIN(21) is B0
      TGAM=UI(9)       !TGAM parameter
      BOOT=ZERO        !flag for eq 3.9 computation of gamma
      GO TO 100
C
C    setup for computation using slope and intercept of hugoniot
C    curve in linear shock velocity--particle velocity approximation
C
   90 BOOT=-UI(6)      !BOOT is intercept of hugoniot, S0
      GAM=UI(9)        !TGAMMA in this case is slope, S1
C
C     insert data for critical point modification scheme, 3.33, Appendix E
C
  100 IF (UI(13).EQ.ZERO) GO TO 110   !don't modify critical point
      S=UI(14)           !C54 parameter for CP modification, eq 3.33
      IF ((S.LE.ZERO).OR.(S.GT.QCL1)) S=QCL1  !default value of C54
      ACK(IT+54)=S       !C54 in eq 3.33
C
  110 ACK(IT+53)=UI(13)  !C53 in eq 3.33
C
C   insert data on melt transition
C
      ACK(IT+43)=UI(23)  !heat of fusion, sec V-1
      ACK(IT+44)=UI(24)  !ratio of liquid to solid density at melt temp.
      ACK(IT+46)=ZERO    !density of solid at melt temp.
      ACK46=ACK(IT+46)   !flag for melt transition (presently none)
C
C    check equation type designation DIN(30) = UI(2), equation type
C
      IF (NINT(DIN(30)).EQ.5) NUMPHASE = NUMPHASE + 1  ! Count number of solid phases
      IF ((ISETAB.EQ.0).OR.(NINT(DIN(30)).EQ.2).OR.(IZI.LT.0).OR.
     1    (IZI.GT.4).OR.(DIN(30).LT.ZERO)) GO TO 120
      DIN(30)=IZI
  120 DIN(31)=IZ
      IF (NINT(DIN(30)).GE.-4) GO TO 140
  130 CONTINUE           !bad EOS type number
      IAMBAD=38
      CALL ANWARN(1)     !issue a warning
      WRITE(KLST,1260) 2,DIN(30),ZERO
      DIN(30)=ZERO       !and go on
      GO TO 160
  140 IF (NINT(DIN(30)).GT.5) GO TO 130
C
C     check for possible presence of melt transition, and insert
C     constants for eq 5.24 if there is one, for all EOS types except 2
C
      IF ((NINT(DIN(30)).GE.0).AND.(NINT(DIN(30)).NE.2)) THEN
         IF((ACK(IT+57).EQ.ZERO).OR.(ACK(IT+58).EQ.ZERO).OR.
     1                (ACK(IT+59).EQ.ZERO)) THEN
           ACK(IT+57)=0.3D0      !constants alpha, beta, gamma
           ACK(IT+58)=0.1D0      !for melt energy, eq. 5.24
           ACK(IT+59)=0.2D0
         END IF
         ACK(IT+71)=ONE          !add melt terms unless T < Tm
         CMLT7=ACK(IT+71)        !flag for inclusion of melt transition
         ACK(IT+72)=ZERO         !flag set for fast iteration
         ! Construct a tabular melt-curve if melt transition and solid-solid
         ! transition are selected at the same time.
         IF (UI(18).GT.ZERO.AND.ACK(IT+43).NE.ZERO) THEN
            MELT_CURVE_TABLE = 1
            WRITE(KLST,1365)
         END IF
      END IF
C
C    set ionization model flag for EOS types 1, 2, 4 and 5
C
      IF ((NINT(DIN(30)).EQ.1).OR.(NINT(DIN(30)).EQ.2).OR.
     &       (NINT(DIN(30)).EQ.4).OR.(NINT(DIN(30)).EQ.5)) THEN
        IF(ACK(IT+63).EQ.ZERO) THEN
          ACK(IT+63)=ONE   !use Saha computation
        ELSE
          ACK(IT+63)=-ONE  !use Thomas-Fermi table
        END IF
      END IF
C
C    Set default minimum density for EOS types other than 2:
C    if DIN(23) = UI(17).eq.0, minimum density equals
C    a fraction QCL2 of the reference density DIN(11).
C
      IF (DIN(23).LE.ZERO.AND.ABS(NINT(DIN(30))).NE.2)
     1       DIN(23)=QCL2*DIN(11)
C
      IF (DIN(25).EQ.ZERO.AND.ABS(NINT(DIN(30))).NE.2) DIN(25)=0.025D0
      IF(NOTRAD.EQ.0) GO TO 150
C
C     Test for presence of thermal conduction inputs
C
      ACK(IT+41)=UI(16)   !temperature dependence of thermal conduction coeff.
      IF (UI(15).LE.ZERO) GO TO 150
      DIN(22)=QCL3/UI(15) !C22 in equation 7.8
      ACK(IT+42)=MAX(QCL6,MIN(DIN(22)/QCL4,QCL5))  !C42 in eq 7.15
      GO TO 160
  150 DIN(22)=ZERO         !no phonon conduction
      ACK(IT+42)=QCL5      !constant in equation 7.16
c
c    Echo contents of UI to output file
c
  160 I2=0
        DO 162 J1=1,NINPUT  !find location of last nonzero entry in UI array
        IF(UI(J1).NE.ZERO) I2=J1 
  162   CONTINUE
      I2=(MAX(I2,24)+2)/3   !compute indices to format list in 3 columns
      J2=3*I2
        DO 170 I=1,I2
        WRITE(KLST,1470) (J1,UI(J1),J1=I,J2,I2)
  170   CONTINUE
C
      IF(DIN(11).LE.ZERO) THEN  !negative or zero reference density
        IAMBAD=38
        CALL ANWARN(1)          !issue warning to user
        WRITE(KLST,1260) 3,DIN(11),ONE
        DIN(11)=ONE             !default to density 1
      END IF
C
      IF(DIN(28).LE.ZERO) THEN  !zero or negative number of elements
        IAMBAD=39
        CALL ANWARN(1)          !warn user
        WRITE(KLST,1260) 1,DIN(28),ONE
        DIN(28)=ONE             !default to 1 element
      END IF
C
C     Process input data on elemental abundances
C
      J1=DIN(28)                !number of elements to process
      S=ZERO
      IZI=IZ+J1-1
      IF (ZZS(IZ).EQ.ZERO) THEN
         DO 172 I=IZ,IZI
         ZZS(I)=SZ23(I-IZ+1)   !data from input file is in SZ23
  172    COT(I)=SC23(I-IZ+1)   !data from input file is in SC23
      END IF
         DO 180 I=IZ,IZI
         IF (COT(I).GT.ZERO) GO TO 180  !if COT .lt. 0, input is 
         IKK=ZZS(I)                     !atomic weight fraction, so
         IKK=(IKK*(IKK+1))/2            !convert to atomic number fraction
         COT(I)=-COT(I)/EIP(IKK)
  180    S=S+COT(I)                     !total abundance
      DIN(29)=ZERO
      DIN(26)=ZERO
      S1=ZERO
        DO 200 I=IZ,IZI                 !loop over elements present
        COT(I)=COT(I)/S                 !normalized atomic number fraction
        DIN(26)=DIN(26)+ZZS(I)*COT(I)   !construct mean atomic number
        IKK=ZZS(I)
        IKJ=IKK+(IKK*(IKK+1))/2
        IF (IKK.GE.1.AND.IKK.LE.92) GO TO 190
        IF ((ACK(IT+63).LT.ZERO).AND.(IKK.LE.102)) THEN
           TEIP=BIGZAS(IKK-92)  
           GO TO 191
        END IF
        CALL ANWARN(1)        !cannot use saha method for Z .gt. 92
        WRITE(KLST,1490) IKK  !warn user that this will not work
        CALL ANMARK('IONIZATION POTENTIALS')
        IAMBAD=3
        GO TO 1220                  !skip to next material
  190   TEIP=EIP(IKJ-IKK)           !Atomic weight of species I
  191   DIN(29)=DIN(29)+COT(I)*TEIP !mean atomic weight
  200   S1=S1+COT(I)*TEIP*QCL7
c
c    Compute total number of atoms N0 (eq. 6.10)
c      stored in DIN(27), then to ACK(IT+27)
c      FNI(I) is number of atoms/gm of element I
c      QCL7 is 1/Avagadro Number
c   
      DIN(27)=ZERO
        DO 210 I=IZ,IZI
        FNI(I)=COT(I)/S1
  210   DIN(27)=DIN(27)+FNI(I)
c
c    Echo data on elements, abundances to output file
c      ZZS(I) = atomic number of element I
c      COT(I) = number fraction of element I
c      FNI(I) = atoms/gm of element I
c
      WRITE(KLST,1480)  !skip a line
      I1=IZ
      I2=DIN(28)        !number of elements in material
c
c    write DIN(28) lines in output file summarizing element abundances
c
        DO 890 I=1,I2
        WRITE(KLST,1520) I,ZZS(I1),I,COT(I1),I,FNI(I1)
  890   I1=I1+1
c
      IF(IAMBAD.NE.0) GO TO 1220    !skip to next material if bad flag set
      IF(NINT(DIN(30)).EQ.-4) THEN  !type -4 only
      DO 2001 I=1,99         !Initialize ACK array
 2001   ACK(IT+I)=ZERO
        ACK(IT+11)=UI(3)         !reference density 
        ACK(IT+12)=UI(4)/TCONV   !reference temperature
        ACK(IT+26)=DIN(26)       !mean atomic number
        ACK(IT+27)=DIN(27)       !number of atoms/gm
        ACK(IT+28)=UI(1)         !number of elements in material
        ACK(IT+29)=DIN(29)       !mean atomic weight
        ACK(IT+30)=DIN(30)       !EOS type switch
        ACK(IT+31)=DIN(31)       !internal storage location
        WRITE(KLST,2004)
 2004   FORMAT(/,'  Defining user supplied eos')
        CALL ANUSET(ISE,UI,ACK(IT+1))   !call user-defined eos routine
        KUWAY=1
        GO TO 8790           !proceed to final checks for type -4
      END IF
C
      DO I=26,31     !copy DIN to ACK for locations 26 to 31
        ACK(IT+I)=DIN(I)
      END DO
C
      IF(DIN(15).LE.ZERO) THEN   !reference gamma is zero or negative
        IF (NINT(DIN(30)).NE.2) THEN
          IAMBAD=15
          CALL ANWARN(1)         !warn the user
          WRITE(KLST,1260) 7,DIN(15),ONE
          DIN(15)=ONE            !reset gamma to 1.
        END IF
      END IF
C
      IF (NINT(DIN(30)).LE.-1) GO TO 8700   !eos type -1
      IF (NINT(DIN(30)).EQ.2) GO TO 240     !eos type 2
C
      IF(ACK(IT+61).LE.-ONE) THEN   !gamma interpolation paramter is bad
        IAMBAD=31
        CALL ANWARN(1)              !warn the user
        WRITE(KLST,1260) 31,ACK(IT+61),ZERO
        ACK(IT+61)=ZERO             !do not use interpolation method
      END IF
C
      IF(DIN(24).LT.ZERO) THEN      !C24 cannot be less than zero
        IAMBAD=33
        CALL ANWARN(1)              !warn the user
        WRITE(KLST,1260) 10,THREE*DIN(24),ZERO
        DIN(24)=ZERO                !reset C24 to zero
      END IF
C
C    See if user Evap input is in ev/atom. If so, convert to ev/gm
C
      IF (DIN(10).LT.ZERO) DIN(10)=-DIN(10)*DIN(27)*BOLTS
C
      IF (BOOT.LE.ZERO) GO TO 240  !skip hugoniot calc of Tgamma
C
C    Compute B0 and Tgamma in eq 3.10 from slope, intercept of linear hughniot.
C     BOOT is intercept S0 here
C
      S1=THREE*DIN(27)*BOLTS*DIN(12)*DIN(15)**2
      DIN(21)=DIN(11)*(BOOT**2-S1)               !B0 is defined
      S2=DIN(11)*BOOT**2/DIN(21)
      S2=S2*(TWO*GAM-ONE-(DIN(15)-TWO)*(HALF-HALF/S2))
      S3=S1*DIN(11)/DIN(15)
      S4=S3*(ONE+TWO*S2)+DIN(21)
      S5=S4**2-EIGHT*S2*S3**2
      IF (S5.GT.ZERO) THEN
        S6=HALF*(S4+SQRT(S5))
        S6=DIN(21)*S6/(S6-S3)**2
        S5=ONE-S6
        IF (ABS(S5).GT.QCL5) S6=ONE-SIGN(QCL5,S5) !QCL5 = 0.1
      ELSE
        S6=ONE
      END IF
      TGAM=THREE*(S6*S2-DIN(15))               !Tgamma is defined
      WRITE(KLST,1600) DIN(21),TGAM
 1600 FORMAT(/'  HUGONIOT INPUT DATA YIELDS B0 = ',1PD15.6,
     & ' TGAMMA = ',D15.6)
C
C       Compute constant C13 of equation 4.20 for free energy of a
C       monatomic gas
C
  240 S1=ZERO
      DO I=IZ,IZI
        IKK=ZZS(I)
        IKK=(IKK*(IKK+1))/2
        S=EIP(IKK)*QCL7  !mass of atom I
        S1=S1+LOG(FNI(I)/(DIN(27)*(DIN(27)*S)**QCL8))*FNI(I)/DIN(27)
      END DO
      DIN(13)=QCL9*DIN(27)**(FIVE/THREE)*EXP(TWO*S1/THREE)
CHJM
C       Read in data for molecular clusters:
C        UI(36) = Number of atoms in an atomic cluster
C        UI(37) = Binding energy in eV (skip input if zero)
C        UI(38) = Number of rotational degrees of freedom
C        UI(39) = Molecular radius in cm
C        UI(40) = Number of vibrational degrees of freedom
C        UI(41) = Vibrational "Debye temperature"
C 
C      NB:  DIN(27) = number of atoms/gm
C           DIN(29) = mean molecular weight, gm/mole
C
      TPWR=ZERO
      IF(UI(36).EQ.TWO) THEN           !diatomic molecular clusters
        TPWR=HALF*(THREE-UI(38))
        CXX=(QCL57*DIN(29))**TPWR/(QCL59*DIN(27))
        RPWR=ONE
      ELSE IF (UI(36).EQ.THREE) THEN   !triatomic molecular clusters
        TPWR=THREE-UI(38)*HALF
        CXX=(QCL60*DIN(29))**TPWR/(QCL61*DIN(27)**2)
        RPWR=TWO
      ENDIF
	IF(UI(38)*UI(39).NE.ZERO) CXX=CXX/UI(39)**UI(38) !molecular radius term
	ACK(IT+79)=CXX
      ACK(IT+80)=UI(37)
      ACK(IT+81)=TPWR          !exponent of temperature for Y
      ACK(IT+84)=UI(40)
      ACK(IT+85)=UI(41)/TCONV
      ACK(IT+86)=UI(36)        !number of atoms in a molecular cluster
      ACK(IT+87)=RPWR          !exponent of density for Y
C
C    Transfer data for Mie potential
C       UI(42) = flag for Mie treatment
C       UI(43) = power law for potential
C
      ACK(IT+83)=UI(43)
      IF(UI(42).EQ.ONE) ACK(IT+82)=ONE
C
CHJM
      IKK=0                              !loop counter for B00, P0 iteration
      GAM=DIN(15)+TGAM/THREE             !DIN(15) is reference gamma, eq 3.15
      IF (NINT(DIN(30)).EQ.2) GO TO 410  !gas only with electronic terms, skip to later
      IF (NINT(DIN(30)).GE.0) THEN       !all other forms
C
C      compute and check coefficients for computation of gamma in expanded states
C
        S3=ACK(IT+61)+ONE  !C61 = interpolation paramter, usually 0
        DIN(14)=DIN(25)    !constant C14 in equation 4.26  DIN(25) = theta0
     &         *EXP(QCL8*S3-(TWO-HALF*ACK(IT+60))*DIN(15))/DIN(11)**S3
        DIN(16)=(S3-(TWO-ACK(IT+60))*DIN(15))/DIN(11)**2    !C16 in eq. 4.24
        DIN(17)=((THREE-ACK(IT+60))*DIN(15)-TWO*S3)/DIN(11) !C17 in eq. 4.24
        IF(ABS(DIN(16)).GT.QCL17) THEN   !QCL17 = 1.E-9 ...a small number
          S=-DIN(17)/(TWO*DIN(16))
          IF(S.GE.ZERO) THEN
            S1=S3-DIN(17)*DIN(17)/(FOUR*DIN(16))
C
            IF(S1.LE.ZERO) THEN   !negative GAMMA at density S
              IAMBAD=32
              CALL ANWARN(1)      !warn user of possible trouble
              WRITE(KLST,377) S1,S
  377         FORMAT(' INPUTS YIELD GAMMA OF',1PE13.6,' AT RHO',E13.6,/,
     &        ' ABORT POSSIBLE')
            END IF
C
          END IF
        END IF
      END IF
C
C
C  The following large block of spaghetti code is SLT's multi-loop to compute
C  the coefficients for the  cold compressed presssure and energy. 
C  It involves a lot of fiddiling to check on the right range of
C  gamma and to find values of B00 and rho00 that give the correct reference
C  pressure P0 and bulk modulus B0 in the reference state.  This method
C  makes essential use of the Morse potential form for expanded states, because
C  the reference state with rho0 < rho00 is, by definition, an expanded state.
C
C  Rather than try to alter these hundreds of lines of code, HJM chose the messy
C  route of duplicating the code section for the case of the Mie-type
C  potential.  The second derivative of pressure with respect to eta is *not*
C  contrained to be continuous at eta = 1 (even SLT did not *require* this!).
C
C  This modification has little affect on the high pressure and temperature
C  behavior of the code, but it does improve the code's behavior near the
C  reference pressue and temperature 
C
C
      IF(ACK(IT+82).EQ.ONE) GO TO 555    !jump to special evaluation for Mie potential
C
C     The following block of code is SLT's unmodified 1990 version (except for lots of
C     comments added!) for the Morse potential in expanded states
C
C    loop to assure proper upper limit on cold pressure, see eq 3.32
C    this computation depends crucially on the form of the Morse potential
C
      I=0         !initialize iteration counter
      S3=GAM      !save old gamma value
      SPS=1.D6    !SPS LIMITS POTENTIAL RANGE IF POSSIBLE = 1 bar
  260 S=ONE-DIN(21)/(DIN(11)*DIN(10)*GAM**2)  !argument of Sqrt in eq 3.27
      IF (S.LE.ZERO) GO TO 280      !beware of negative argument
      S=SQRT(S)                     !take the Sqrt if S .gt. 0.
      S1=DIN(21)/(QCL10*SPS*GAM*S)  !QCL10 = 200.
      S1=LOG(S1)
      S2=QCL11*GAM*(ONE-S)  !Exponent in eq 3.32.  QCL11 = 27.
      IF (S2.GE.S1) GO TO 280
      I=I+1         !increment counter
      IF (I.GT.400) GO TO 270  !no convergence after 400 steps
      GAM=QCL12*GAM !decrease gamma by 1%.  QCL12 = 0.99
      GO TO 260     !loop again
  270 GAM=S3        !convergence failure, use old gamma
  280 DFB1=GAM      !iteration succeeded, save new gamma in DFB1
C
  290 S=DIN(13)*DIN(12)*(DIN(11)**ATHIRD/DIN(25))**2
      SPS=S
      I=0           !initialize gamma loop counter
      S1=DIN(20)-DIN(11)*((THREE*DIN(15)+S)/
     1    (ONE+S))*DIN(27)*DIN(12)*BOLTS
      IKK=IKK+1     !update iteration loop counter
      IF (IKK.EQ.2) GO TO 310    !call to ANEOS1 on *second* iteration
      IF (DIN(15).EQ.ONE) GO TO 300
      S=ONE+((TWO*DIN(15)-ONE)**2-TWO)*S1/DIN(21)
      DIN(3)=S**2+FOUR*DIN(15)*(DIN(15)-ONE)*(ONE-TWO*S1/DIN(21))**2
      DIN(3)=DIN(21)*(SQRT(DIN(3))-S)*HALF/(DIN(15)-ONE)  !B00
      GO TO 320
  300 DIN(3)=((DIN(21)-TWO*S1)**2)/(DIN(21)-S1)           !B00
      GO TO 320
  310 SQTS(1)=SQRT(DIN(12))
C
C   Call ANEOS1 at reference temp DIN(12), density DIN(11).
C   Output S6 is derivative DPDR. Called only after elastic
C   constants are defined below, on *second* iteration IKK = 2
C
      CALL ANEOS1 (DIN(12),DIN(11),S,S2,S3,S4,S5,S6,IT)
C
      S6=MAX(QCL14,MIN(QCL13,(DIN(21)/(DIN(11)*S6))))

C      QCl14 = 0.9, QCL13 = 1.1; S6 = B0/(rho0 * dP/drho), if it is
C      within these bounds
      DIN(3)=S6*DIN(3)                     !adjust B00 = S6*B00
  320 GAM=DFB1
C
C     start of GAMMA iteration loop
C
  330 S2=DIN(3)/(DIN(11)*DIN(10)*GAM**2) !factor in eqs. 3.27 and 3.28
C                                        !S2 = B00/(rho0*Evap*GAMMA**2)
      IF(S2.LE.ZERO) THEN
        CALL ANWARN(1)
        WRITE(KLST,332) S2,DIN(3),DIN(10),DIN(11),GAM,IKK
  332   FORMAT(' GAM LOOP',1P5E13.6,I5,/,' ANEOS MATERIAL ABORT')
        IAMBAD=7
        GO TO 1220    !skip to next material
      END IF
C
C     next calculation changed 6/87 to give user more control over
C     expanded state cold curve - will give the same result as old
C     model if no user input
C
      IF(ACK(IT+55).LE.ZERO) ACK(IT+55)=ONE
      ACK(IT+55)=MAX(QCL53,MIN(ACK(IT+55),QCL55))
      ACK(IT+56)=MIN((HALF*ACK(IT+55)),MAX(ACK(IT+56),QCL54))
      IF((S2.LE.ACK(IT+55)).AND.(S2.GE.ACK(IT+56))) GO TO 350
      S3=GAM*SQRT(QCL15*S2/ACK(IT+55))
      S4=GAM*SQRT(S2/ACK(IT+56))
      GAM=MAX(S3,MIN(S4,GAM))
      IF (I.GT.15) GAM=GAM*QCL16
      I=I+1           !increment loop counter
C
      IF (I.GT.40) THEN  !fails to converge
         CALL ANWARN(1)  !warn user
         CALL ANMARK('BAD ANEOS GENERATION / GAM')
         IAMBAD=8
         GO TO 1220      !skip to next material
      END IF
C      
      IF (IKK.GE.2) GO TO 330  !back to start of GAMMA loop
C
      S1=DIN(20)-DIN(11)*((THREE*DIN(15)+SPS)/(ONE+SPS))
     1   *DIN(27)*DIN(12)*BOLTS
      IF (GAM.EQ.ONE) GO TO 340
      S=ONE+((TWO*GAM-ONE)**2-TWO)*S1/DIN(21)
      DIN(3)=S**2+FOUR*GAM*(GAM-ONE)*(ONE-TWO*S1/DIN(21))**2
      DIN(3)=DIN(21)*(SQRT(DIN(3))-S)*HALF/(GAM-ONE)
      GO TO 330    !loop back to 330
  340 DIN(3)=((DIN(21)-TWO*S1)**2)/(DIN(21)-S1)  !B00
      GO TO 330    !loop back to 330
  350 CONTINUE     !successful completion of GAMMA loop
C     
C    Iteration loop to adjust coefficients for Morse potential.
C    Loop uses bisection between upper and lower limits S4 and S3.
C    S5 is rho0/rho00 upon exit from loop
C 
      S3=ONE     !upper limit of rho0/rho00 for bisection
      S4=.8D0    !lower limit of rho0/rho00 for bisection
  360 S5=HALF*(S3+S4)           !rho0/rho00
      S6=SQRT(ONE-S2*S5)        !argument of SQRT in eqs. 3.27 and 3.28
      DIN(5)=THREE*GAM*(ONE+S6) !C5 in eq 3.20
      DIN(6)=THREE*GAM*(ONE-S6) !C6 in eq 3.20
      S6=SIX*GAM*S6
      DIN(4)=S5**(-ATHIRD)      !C4 in eq 3.20
      IF (S3-S4.LE.QCL17) GO TO 390  !successful exit
      S6=DIN(5)*(ONE-DIN(4))
      PCP1=DIN(6)*(ONE-DIN(4))
      S6=ONE-THREE*DIN(3)*(EXP(S6)-EXP(PCP1))/(DIN(4)**2*S1*S6)
      IF (S6) 370,390,380       !exit if S6 = 0.
  370 S4=S5
      GO TO 360     !loop back to 360
  380 S3=S5
      GO TO 360     !loop back to 360
  390 CONTINUE      !successful completion
C
C     compute coefficients C32 to C37 for compressed pressure, energy
C     formulae, eqs 3.4 and 3.17
C
      DIN(19)=DIN(11)/S5       !density rho00
      S3=DIN(5)*(ONE-DIN(4))
      S4=DIN(6)*(ONE-DIN(4))
      DIN(4)=S1/(S5**(TWO/THREE)*(EXP(S3)-EXP(S4)))  !C4 in eq 3.20
C
C    C32 in equation 3.5
C
  400 DIN(32)=(THREE*QCL18**2/(QCL19*QCL20*PIE))*
     1  (PIE/THREE)**(ATHIRD)*(DIN(19)*DIN(26)*DIN(27))**(FIVE/THREE)
C
C    C33 in equation 3.6
C
      DIN(33)=(PIE*QCL20/QCL21)*(QCL22/QCL18)**2*
     1  (QCL23*DIN(26)**(ATHIRD)/FIVE+QCL24/(QCL25*PIE**2*DIN(26))
     2  **(ATHIRD))/(TWO*DIN(19)*DIN(27))**(ATHIRD)
C
C    C34 in equation 3.12
C
      S2=DIN(33)                  !C33
      S=DIN(32)*EXP(-S2)          !C32*EXP(-C33)
      GAMHAT=DIN(15)+TGAM/THREE   !GammaHat in eq 3.15
      DIN(34)=S*(SIX+THREE*DIN(33)+HALF*DIN(33)**2)-CNINE*DIN(3)*GAMHAT
C
C    C35 in equation 3.13
C
      DIN(35)=THREE*DIN(3)*(SIX*GAMHAT+ONE)
     1   -S*(QCL26+SEVEN*DIN(33)+DIN(33)**2)
C
C    C36 in equation 3.14
C
      DIN(36)=S*(TEN+FOUR*DIN(33)+HALF*DIN(33)**2)
     1   -THREE*DIN(3)*(THREE*GAMHAT+ONE)
C
C     Compute C37 for compressed energy expression, eq 3.18
C
      S1=EXP(-S2)
      CALL ANEI3 (S2,S1,S)   !third exponential integral
      DIN(37)=THREE*DIN(32)*S+DIN(34)+QCL27*DIN(35)+THREE*DIN(36)
C
C    Transfer data to define first 40 locations of ACK
C
  410   DO 420 I=1,40
  420   ACK(IT+I)=DIN(I)
C
C   more loop control.  If IKK = 1, loop back to 290 and call ANEOS1
C   using the just-computed values of pressure coefficients.
C   If IKK .gt. 1, then go to 440 and compute setup for P0, B00 loop
C   If IKK .lt. 1, skip to 540 and exit B00, rho00 loop
C
      IF (IKK-1) 540,290,440
  440 IKK=0
      S1=DIN(11)/DIN(19)  !current rho0/rho00
      S2=DIN(3)           !current B00
      DFB1=ACK(IT+4)      !save C4 for Morse or L-J potentials
      DFB2=ACK(IT+53)     !save C53 for cold potential modification
      ACK(IT+53)=ZERO
      ACK(IT+4)=ZERO
      SQTS(1)=SQRT(DIN(12))  !Sqrt of reference temperature
C
C    call ANEOS1 for reference temperature, density with
C    C4 and C53 temporarily zeroed--no cold pressure contribution
C    output S8 is DPDR, exclusive of cold pressure.
C
      CALL ANEOS1 (DIN(12),DIN(11),S9,S3,S4,S5,S6,S8,IT)
C
      S8=S8*DIN(11)    !S8 is now DPD(eta) at rho0 for thermal terms only
      ACK(IT+4)=DFB1   !restore C4
      ACK(IT+53)=DFB2  !restore C53
C
C     Start of loop for double newtonian iteration for P0 and B00
C     described by equations 3.36 and 3.37.  B00 and rho00 are adjusted
C     until P0 and B0 are matched at the reference density rho0
C
  450 IKK=IKK-1        !decrement IKK at start of loop
C
      IF(IKK.LT.-9990) THEN          !trouble--more than 9990 iterations
        WRITE(KLST,537) IKK,GAM,S4,S5,S6,S1,S3,S7
  537   FORMAT(' error RHO00/B00',i6,1p8e13.6)
        IF(IKK.LT.-10000) GO TO 538  !give up after 10000 iterations
      END IF
C
      S=-ONE          !flag for loop control
      BOOT=QCL28*S2   !QCL28=0.9999, S2 = current B00
      ETAOT=S1        !current rho0/rho00
      GO TO 500
  460 PCO1=S3         !current energy E
      PCP1=S7         !current CV
      BOOT=S2         !current B00
      ETAOT=QCL28*S1  !decrement rho0/rho00
      S=ZERO          !flag for loop control
      GO TO 500
  470 PCO2=S3         !current energy E
      PCP2=S7         !current CV
      ETAOT=S1        !current rho0/rho00
      S=ONE           !flag for loop control
      GO TO 500
  480 DFB1=QCL29*(S3-PCO1)/S2  !QCL29 = 10,000
      DFB2=QCL29*(S7-PCP1)/S2
      DFN1=QCL29*(S3-PCO2)/S1
      DFN2=QCL29*(S7-PCP2)/S1+(DIN(21)-S8)/S1**2
      S3=S3+S9-DIN(20)        !DIN(20) is target pressure P0
      S7=S7-(DIN(21)-S8)/S1   !DIN(21) it target modulus B0
      S=DFN1*DFB2-DFN2*DFB1   !loop control variable
      IF(S.EQ.ZERO) GO TO 538     !iteration failue
      DFB1=(S7*DFB1-S3*DFB2)/S
      DFB2=(S3*DFN2-S7*DFN1)/S
      IF (ABS(DFB1).GT.QCL30*S1) DFB1=SIGN(QCL30*S1,DFB1) !QCL30 = .002
      IF (ABS(DFB2).GT.QCL30*S2) DFB2=SIGN(QCL30*S2,DFB2)
      IF (S3.LT.ZERO.AND.IKK.GT.-200) GO TO 490
      IF (ABS(DFB1).GT.QCL31*S1) GO TO 490  !QCL31 = 1.e-10
      IF (ABS(DFB2).LE.QCL31*S2) GO TO 530  !successful completion
  490 S1=S1+DFB1     !new rho0/rho00
      S2=S2+DFB2     !new B00
      IF((S1.LE.ZERO).OR.(S2.LE.ZERO)) GO TO 538  !iteration failure
      GO TO 450      !loop back to 450
  500 S4=BOOT*ETAOT/(DIN(11)*DIN(10)*GAM**2)  !factor in eq 3.27
      IF(IKK.GT.-100) THEN  !fewer than 100 iterations
        S5=ACK(IT+55)
        S6=ACK(IT+56)
      ELSE         !more than 100 iterations, try fresh start
        S5=QCL55   !QCL55 = 1 - QCL54, just shy of 1
        S6=QCL54   ! = 1.E-6, very small number
      END IF
C
      IF(S4.GT.S5) THEN
        GAM=SQRT(BOOT*ETAOT/(DIN(11)*DIN(10)*S5))
        S4=MIN(QCL55,BOOT*ETAOT/(DIN(11)*DIN(10)*GAM**2))
      ELSEIF(S4.LT.S6) THEN
        GAM=SQRT(BOOT*ETAOT/(DIN(11)*DIN(10)*S6))
        S4=MIN(QCL55,BOOT*ETAOT/(DIN(11)*DIN(10)*GAM**2))
      END IF
C
      S4=SQRT(ONE-S4)         !Sqrt in eqs 3.27 and 3.28
      S5=THREE*GAM*(ONE+S4)   !C5 in eq 3.27
      S6=THREE*GAM*(ONE-S4)   !C6 in eq 3.28
      S4=BOOT/(TWO*GAM*S4)    !C4 in eq 3.29
      DFN1=ETAOT**ATHIRD
      DFN2=EXP(S5*(ONE-ONE/DFN1))
      DFB1=EXP(S6*(ONE-ONE/DFN1))
      DFB2=DFN1**2
      S3=S4*(DFN2-DFB1)*(ETAOT/DFN1)
      S7=S4*((TWO*DFN1+S5)*DFN2-(TWO*DFN1+S6)*DFB1)/(THREE*DFB2)
      IF (S) 460,470,480   !loop back, depending on S
C
C    loop successfully completed.  Use new B00, rho00 to compute constants
C    for cold pressure constants for both compressed and expanded states. 
C
  530 DIN(3)=S2           !new B00
      DIN(19)=DIN(11)/S1  !new rho00
      DIN(4)=S4           !new C4
      DIN(5)=S5           !new C5
      DIN(6)=S6           !new C6
      GO TO 400           !compute new compressed state constants, C32, etc.
C
C     Iteration failure to converge on B0, P0
C
  538 CONTINUE
      CALL ANWARN(1)
      WRITE(KLST,539) IKK,S1,S2,DFB1,DFB2,GAM,S
  539 FORMAT(/,' RHO00-B00 ITERATION FAILURE',I10,/,6D16.7,/,
     1' CHECK REFERENCE POINT')
      IAMBAD=6
C
  540 CONTINUE   !successful exit from looping.  Constants C32, etc. defined
      IF(NINT(DIN(30)).EQ.2) GO TO 550   !type 2 jump
      ACK(IT+75)=TGAM   !Tgamma saved
      ACK(IT+76)=TGAM
      IF(GAM.NE.ZERO)
     & ACK(IT+77)=ONE-((ACK(IT+5)-ACK(IT+6))/(SIX*GAM))**2
      IF (ABS(DIN(15)+TGAM/THREE-GAM).LT.QCL33) GO TO 550  !QLC33 = 1.E-4
      ACK(IT+75)=THREE*(GAM-DIN(15))  !backup calculation of Tgamma
      GO TO 550
C
CHJM *****************************************************************************
C     The following code section is implemented when the Mie potential
C     option is employed.  The Morse potential code above is short-circuited.
C
  555 CONTINUE
C
C     copy first 40 locations from DIN to ACK array
C
      DO I=1,40
        ACK(IT+I)=DIN(I)
      END DO
C
C    FIRST CUT AT REFERENCE DENSITY--MUST BE DEFINED BEFORE CALL TO ANEOS1
C
      ACK(IT+19) = 1.1D0*ACK(IT+11)  
C
C    Begin search for rho00, B00 by computing the derivative of P with respect
C    to density for the thermal part of the eos at the reference density.  This
C    thermal contribution to B0 will not change with later iterations on rho00
C
      DFB1=ACK(IT+4)      !save C4 for Morse or L-J potentials
      DFB2=ACK(IT+53)     !save C53 for cold potential modification
      ACK(IT+53)=ZERO
      ACK(IT+4)=ZERO
C
C    call ANEOS1 for reference temperature, density with
C    C4 and C53 temporarily zeroed--no cold pressure contribution
C    output DPDRTH, exclusive of cold pressure, and PTH is thermal pressure only.
C
      CALL ANEOS1 (DIN(12),DIN(11),PTH,S3,S4,S5,S6,DPDRTH,IT)
C
      ACK(IT+4)=DFB1        !restore C4
      ACK(IT+53)=DFB2       !restore C53
      BTH=DPDRTH*DIN(11)    !BTH is DPD(eta) at rho0 for thermal terms only
      C1=ACK(IT+20)-PTH     !C1 is Pref - Pthermal
      C2=ACK(IT+21)-BTH     !C2 is Bref-Bthermal
      ALJEXP=ACK(IT+83)     !input power law, a
      EVAP=ACK(IT+10)       !vaporization energy, Evap
C
C    Make sure that the reference pressure is less than the thermal pressure.
C    If not, Pcold is positive, not negative, and no solution can be found!
C
      IF(C1.GT.ZERO) THEN
        WRITE(KLST,1620) ACK(IT+20),PTH
        STOP
 1620 FORMAT(/' REFERENCE PRESSURE = ',1PD12.4,
     & ' IS LARGER THAN THERMAL PRESSURE = ',D12.4/
     & ' NO SOLUTION FOR THE COLD REFERENCE STATE IS POSSIBLE!'/
     & ' EITHER LOWER REFERENCE PRESSURE OR RAISE REFERENCE',
     &' TEMPERATURE')
      ENDIF
C
C    compute the constants needed for the solution for rho00 and B00
C
      D1=-ALJEXP*C1/C2
      D2=(C2-ALJEXP*C1)/(ACK(IT+11)*EVAP*(ALJEXP-ONE))
      GMIN=(ALJEXP-ONE)/D2
      ALPHA=D1*D2/((ONE+D1)*ALJEXP)
C
C    first compute an approximate solution to this problem, then
C    complete the search with a bisection algorithm
C
      GSTAR=ONE+ALPHA
      GL=GMIN     !minimum value to search for solution
      GH=100.D0   !maximum value to search for solution
      ITER=0
C
C    the following two write statements are for debugging a failure
C    to converge--uncomment them to see the course of convergence
C
C      WRITE(KLST,1651)
C 1651 FORMAT(//'      GSTAR           GL              GH',
C     & '          TEST'/)
 556  TEST=(GSTAR-GMIN)*DLOG(GSTAR)-GSTAR*DLOG(ONE+ALPHA*(GSTAR-GMIN))
      ITER=ITER+1
      IF(TEST.LT.ZERO) THEN
        GL=GSTAR
      ELSE
        GH=GSTAR
      END IF
      DG=HALF*(GH-GL)
      IF(DABS(DG).LE.1.D-12) GO TO 557
      IF(ITER.GT.100) THEN
        WRITE(KLST,*) 'L-J ITERATION FAILED TO CONVERGE'
        STOP
      END IF
      GSTAR=GL+DG
C      WRITE(KLST, 1650) GSTAR,GL,GH,TEST
C1650  FORMAT(1P4D15.6)
      GO TO 556
  557 CONTINUE
C
C   Use this solution to compute the parameters describing the cold
C   compression curve.
C
      SFACT=GSTAR*D2
      IF((ONE+SFACT).LT.ALJEXP) THEN
        WRITE(KLST,1630) ONE+SFACT, ALJEXP
1630  FORMAT(/' COLD PRESSURE EXPONENT B = ',1PD12.4,' IS LESS THAN',
     & ' INPUT EXPONENT A = ',D12.4/
     & ' NO SOLUTION IS POSSIBLE.  POSSIBLY EVAP IS TOO SMALL'/)
        STOP
      ENDIF
      DIN(19)=ACK(IT+11)*GSTAR**(ONE/SFACT)         !rho00
      DIN(3)=DIN(19)*EVAP*(ALJEXP-ONE)*SFACT        !B00
      DIN(4)=DIN(3)/(ONE+SFACT-ALJEXP)              !C in expression for Pcold
      DIN(5)=ONE+SFACT                              !exponent b
      DIN(6)=ACK(IT+83)                             !exponent a
C
C     check that the vapor pressure vanishes fast enough at low density--use SLT
C     criterion that Pcold < 1 bar at eta = 1.d-4
C
      ELIM=1.D-4
      PCLIM=DIN(4)*(ELIM**DIN(5)-ELIM**DIN(6))
      IF(DABS(PCLIM).GT.1.D6) THEN
        WRITE(KLST,1640) PCLIM,ELIM
1640  FORMAT(/'  **WARNING**'/
     &' HIGH PRESSURE ENCOUNTERED IN COLD, EXPANDED STATE, PC = ',
     & 1PD12.4,' AT ETA = ',D12.4/' TRY INCREASING MIE',
     &' POTENTIAL EXPONENT'/)
      ENDIF
C
C    Once rho00 = DIN(19), TGAM and B00 = DIN(3) are computed, the constants
C    for the compressed pressure, eq 3.4, and energy, eq 3.17, can be evaluated
C
C    C32 in equation 3.5
C
      DIN(32)=(THREE*QCL18**2/(QCL19*QCL20*PIE))*
     1  (PIE/THREE)**(ATHIRD)*(DIN(19)*DIN(26)*DIN(27))**(FIVE/THREE)
C
C    C33 in equation 3.6
C
      DIN(33)=(PIE*QCL20/QCL21)*(QCL22/QCL18)**2*
     1  (QCL23*DIN(26)**(ATHIRD)/FIVE+QCL24/(QCL25*PIE**2*DIN(26))
     2  **(ATHIRD))/(TWO*DIN(19)*DIN(27))**(ATHIRD)
C
C    C34 in equation 3.12
C
      S2=DIN(33)                  !C33
      S=DIN(32)*EXP(-S2)          !C32*EXP(-C33)
      GAMHAT=DIN(15)+TGAM/THREE   !GammaHat in eq 3.15
      DIN(34)=S*(SIX+THREE*DIN(33)+HALF*DIN(33)**2)-CNINE*DIN(3)*GAMHAT
C
C    C35 in equation 3.13
C
      DIN(35)=THREE*DIN(3)*(SIX*GAMHAT+ONE)
     1   -S*(QCL26+SEVEN*DIN(33)+DIN(33)**2)
C
C    C36 in equation 3.14
C
      DIN(36)=S*(TEN+FOUR*DIN(33)+HALF*DIN(33)**2)
     1   -THREE*DIN(3)*(THREE*GAMHAT+ONE)
C
C     Compute C37 for compressed energy expression, eq 3.18
C
      S1=EXP(-S2)
      CALL ANEI3 (S2,S1,S)   !third exponential integral
      DIN(37)=THREE*DIN(32)*S+DIN(34)+QCL27*DIN(35)+THREE*DIN(36)
C
C             END OF MIE POTENTIAL MODIFICATION
CHJM ************************************************************************
C
C    Call to ANPHTR to add a solid-solid phase transition
C
  550 CALL ANPHTR (DIN,TGAM)
C
C    Transfer data from DIN to define first 40 locations of ACK
C
      DO 560 I=1,40
  560 ACK(IT+I)=DIN(I)
C
C    If DIN(18) = Melt Temperature is GT 0, then go on -> 630
C
      IF (DIN(18).GT.ZERO) GO TO 630

C    Otherwise if LT 0, it is the energy to the melting point
C    at zero pressure from the reference state, so compute
C    the melt temperature from this energy
C
      IF (ABS(NINT(DIN(30))).EQ.2) GO TO 630
      S7=ONE
  570 SPS=DIN(12)       !entry point for re-evalutaion of melt temperature
      S=DIN(11)
      SQTS(1)=SQRT(DIN(12))
      CALL ANEOS1 (SPS,S,S1,S2,S3,S4,S5,S6,IT)
      S9=S7*S2-DIN(18)
      IF (S9) 610,580,580
  580 SPS=SPS+QCL34     !add QCL34 = 0.001 eV to reference temperature
      S8=S2
      IF (SPS.GT.ONE) GO TO 610
      JJ=0
  590 JJ=JJ+1
      IF (JJ.GT.1000) GO TO 610
      SQTS(1)=SQRT(SPS)
      CALL ANEOS1 (SPS,S,S1,S2,S3,S4,S5,S6,IT)
      IF (ABS(S1).LE.TEN) GO TO 600
      S5=S1/S6
      IF (ABS(S5).GT.QCL34*S) S5=SIGN(QCL34*S,S5)
      S=S-S5
      GO TO 590
  600 IF (S2-S9) 580,620,620
  610 CONTINUE
      CALL ANWARN(1)
      WRITE(KLST,1390) SPS,S,S9,S2,S1,S5,S6,JJ
      CALL ANMARK('ANEOS MELT TEMPERATURE ERROR')
      IAMBAD=4
      GO TO 1220    !skip to next material
  620 DIN(18)=((S9-S8)*SPS+(S2-S9)*(SPS-QCL34))/(S2-S8)
      ACK(IT+18)=DIN(18)
C
C    end of melt temperature computation.  Now
C    compute constants for free energy of liquid phase,
C    C43, C44, C45, Equation 5.22
C
  630 IF (NINT(ACK(IT+30)).EQ.2) ACK(IT+43)=ZERO  !material type 2
      IF (ACK(IT+43)) 640,830,650              !ACK(IT+43) is Hf on entry
  640 ACK(IT+43)=ACK(IT+18)*QCL35/ACK(IT+29)   !estimate Hf from eq. 5.59
  650 IF (ACK(IT+44).EQ.ZERO) ACK(IT+44)=QCL36 !ACK(IT+44) is liq/solid
C                                               density ratio at triple pt.

  660 I1=0                 !initialize error flag for ANMAXW
      ACK(IT+49)=ZERO      !lower bound on temperature
      S2=ZERO              !S2 must be defined, but 0. has no effect in ANMAXW
      S1=ZERO              !use RHO00 as initial density guess in ANMAXW
      ACK(IT+47)=BIGNUM    !liquid density at triple point
      IF (NINT(ACK(IT+30)).LT.2) S2=-ONE  !forces P=0 for mat. < 2 in ANMAXW
      IF (NINT(ACK(IT+30)).EQ.5) S2=-ONE  !forces P=0 for HPP in ANMAXW
      S=ACK(IT+18)         !melt temperature Tm
      GAM=ACK(IT+43)       !heat of fusion Hf
      RCT(MAT)=QCL37*ACK(IT+19)  !critical density estimate from ref. density
C
C      Call to maxwell construction of vapor/liquid states
C      at temperature S = Tm
C
      CALL ANMAXW (S,S1,S2,IT,MAT,I1)
C
      IF (ACK(IT+43).EQ.ZERO) GO TO 850  !path for material types < 2
      PCO2=ACK(IT+72)                    !save melt iteration flag
      ACK(IT+72)=ONE                     !force slow melt iteration
      IF (I1.GE.0) GO TO 680             !successful vapor/liquid construction
C
C       The following block of code is a last-ditch attempt to salvage
C       a successful equation of state by raising the melt temperature
C
  670 CONTINUE
      CALL ANWARN(0)      !unable to include melt transiton
      WRITE(KLST,1270)    !will continue without
      ACK(IT+72)=PCO2     !restore melt iteration flag
      IKK=0               !initialize error flag for ANMAXW
      ACK(IT+43)=ZERO     !zero out Hf
      ACK(IT+46)=ZERO     !suppress addition of melt terms
      ACK46=ACK(IT+46)    !ditto
      S2=-ONE             !forces ANMAXW to perform evaluation at P = 0.
C
C      Call to maxwell construction of vapor/liquid states at temperature S
C
      CALL ANMAXW (S,S1,S2,IT,MAT,IKK)
C
      IF (IKK.LT.0) GO TO 830  !still have convergence failure--drop melt
      SQTS(1)=SQRT(S)
      CALL ANEOS1 (S,S1,S4,S5,S6,S7,S8,S9,IT)
      CALL ANWARN(0)        !try inceasing melt temperature S
      WRITE(KLST,1280) S    !warn user of attempt
      DIN(18)=-S5-GAM       !add Hf to energy from ref. to melt temp.
      S7=ZERO
      GO TO 570             !back to compute a new melt temperature
C
C      check liquid and solid densities at triple point, adjust if needed
C
  680 IF (ACK(IT+44).LT.ZERO) ACK(IT+44)=-ACK(IT+44)/S1
      IF (ACK(IT+44).GT.ONE)   !liquid is denser than solid, so adjust it
     &        ACK(IT+44)=ONE/(ONE+S1*(ACK(IT+44)-ONE))
      S2=ACK(IT+44)*S1         !solid density
      IF (S2.LT.S1) GO TO 690  !liquid density < solid density, so go on
      WRITE(KLST,1290) S2,S1   !warn user that liquid > = solid density
      GO TO 670                !try raising Tm
C
C       compute values of constants C43, C44, C45 for melt transition
C
  690 SQTS(1)=SQRT(S)
C
C      Get thermodynamic properties at Tm and density of liquid
C      and solid at triple point to define data for eqs. 5.32-5.35
C
C       note:   S4,   S5,   S6 are Pressure, Energy, Entropy of solid
C       note: DFB1, DFB2, BOOT are Pressure, Energy, Entropy of liquid
C
      CALL ANEOS1 (S,S1,  S4,  S5,  S6,S7,S8,S9,IT) !solid properties
      CALL ANEOS1 (S,S2,DFB1,DFB2,BOOT,S7,S8,S9,IT) !liquid properties
C
      S5=ACK(IT+43)+S5-DFB2                         !Em, equation 5.33
      S6=S6-BOOT+(ACK(IT+43)+S4*(ONE/S2-ONE/S1))/S  !Sm, equation 5.32
      S4=S4-DFB1                                    !Pm, equation 5.35
      ACK(IT+46)=S1           !solid density at triple point
      ACK46=ACK(IT+46)
      ACK(IT+47)=S2           !liquid density at triple point
      S1=S4/S2
      S8=S6*S
      S9=(S1+(TWO*ACK(IT+57)-ACK(IT+58))*S8-ACK(IT+58)*S5)/
     1    (ACK(IT+59)-ACK(IT+58))
      S8=S5+S8-S9
C          use intermediate quantities to compute final results:
      ACK(IT+43)=-TWO*SQRT(S)*S6/S2**ACK(IT+57)  !C43, eq. 5.27
      ACK(IT+44)=S8/S2**ACK(IT+58)               !C44, eq. 5.25 and 5.26
      ACK(IT+45)=S9/S2**ACK(IT+59)               !C45, eq. 5.25 and 5.26
C
C       Use these results to determine boundary lines for the melt region
C
      ACK(IT+48)=BIGNUM   !upper limit to melt temperature
      ACK(IT+50)=BIGNUM   !slope of boundary line A
      ACK(IT+51)=BIGNUM   !slope of boundary line B
      DFN1=ZERO           !initial slope of line A
      ETAOT=ZERO          !will carry value for bottom of upper melt region
      ACK(IT+52)=ZERO     !minimum value for bottom of upper melt region
      PCO1=MIN(TEN*S,QCL38)  !estimate of bottom of upper melt region
C
C       loop to determine boundaries of upper and lower melt regions
C
      S3=S                ! start search at Tm
      PHI=S2              ! First guess at solid melt density
      I2=0                ! Zero ANELSM failure counter
      DO I1=1,501
        S3=QCL39*S3       ! increment S3 by factor QCL39 = 1.02
        IKK=1             ! set operation flag for ANELSM
        SQTS(1)=SQRT(S3)
C
C             get density of liquid and solid at temperature S3
C             S1 = RLIQUID, BOOT = RSOLID
C
        CALL ANELSM (S3,PHI,S1,BOOT,IT,IKK)
C
C       If we are creating a tabular melt curve, failure of ANELSM at this
C       point may not be catastrophic, so continue tentatively, keeping 
C       track of the number of failures. Otherwise, raise the melt temperature
C       and retry. . .
C
        IF (IKK.LT.0) THEN 
           IF (MELT_CURVE_TABLE.EQ.1) THEN
              IF (KEXPRT.NE.0) WRITE(KLST,*)
     $             "Suspected convergence failure in ANELSM."//
     $             " Trying to skip this temperature: ",S3*TCONV
              I2=I2+1           ! count number of ANELSM failures
              CYCLE
           ELSE
              GO TO 670         ! serious error--go back for another try
           END IF
        ENDIF

        PHI=S1        ! Revise guess of solid melt density for next iteration
        DFN2=(S3-S)/(S1-S2)          ! slope of liquid line to liquid triple point
        IF (DFN2.GT.DFN1) DFN1=DFN2  ! keep track of maximum slope
C
C      Usually in lower melt region if temp < lower limit of upper region.
C      However, if the density is less than the reference density, then
C      we are in a hot, low density solid region that is processed here.
C
        IF((S3.LE.PCO1).AND.(BOOT.GT.ACK(IT+11))) THEN !in lower melt region
          S7=(S3-S)/(BOOT-ACK(IT+11)) !slope of solid line to reference point
          IF (S7.LT.ACK(IT+51)) ACK(IT+51)=S7  !keep track of minimum slope
          ETAOT=S3     !save current temperature. Highest value will be C48
        END IF
C
C       We are done if RS/RL-1 is close enough to QCL40 = 5 x 10^-5,
C       OR, if RL is more than QCL41 (= 100) x reference density
C
        IF((BOOT/S1-ONE.LE.QCL40).OR.(S1.GT.QCL41*ACK(IT+11))) THEN
          ACK(IT+48)=S3            !record maximum temperature
          GO TO 750                !done
        END IF
C
      END DO
C        
      CALL ANWARN(0)                 !more than 500 iterations
      WRITE(KLST,1300) S3,S1,BOOT    !print high temperature melt error
      GO TO 670                      !go back for another try
C
C     Skip to here if the above loop to define the melt curve boundaries
C     was successful. . .
C
  750 CONTINUE
      IF (MELT_CURVE_TABLE.EQ.1.AND.I2.GT.0) THEN
         CALL ANWARN(0)
         WRITE(KLST,*)"ANELSM convergence failures in melt"//
     $        "boundary loop: ",I2
         WRITE(KLST,*)"Try reducing step size of it. in ANELSM"
         WRITE(KLST,*)"to improve convergence (QCI4)."
      ENDIF
      ACK(IT+50)=QCL42*DFN1     !increase max slope slightly (QCL42=1.05)
      IF (ETAOT.GT.ZERO) THEN
        ACK(IT+52)=ETAOT             !lower boundary of upper melt region
        ACK(IT+51)=QCL28*ACK(IT+51)  !decrease min slope slightly (QCL28=0.9999)
      END IF
C
C       begin a new iteration downward from Tm to determine Tmin
C
      S1=BIGNUM
      S3=S        !start at Tm
C
      DO I1=1,501
        S3=QCL43*S3  !decrement temperature by QCL43 = 0.99
        IKK=1        !set control flag
        S8=S1        !S8 is the last RS determined
        IF (S8.LT.ACK(IT+47)) IKK=3  !flag if S8 < liq triple pt. density
        SQTS(1)=SQRT(S3)
C          get solid and liquid density for T < Tm, in tensile region
        CALL ANELSM (S3,S2,BOOT,S1,IT,IKK)
C
        IF (IKK.LT.0) GO TO 790          !serious error occurred
        IF (S1.GT.ACK(IT+23)) THEN  !RS still above minimum density
          DFN1=S3                   !save this temperature
          CYCLE
        END IF
  780   ACK(IT+49)=S3             !successful convergence
        GO TO 820                 !jump out to next phase
  790   IF (IKK.NE.-3) GO TO 670  !fatal error, go back for another try
        CALL ANWARN(0)            !warn user of error condition
        WRITE(KLST,1310) ACK(IT+23),S8  !which we will try to fix by
        ACK(IT+23)=S8                   !increasing minimum density
        S3=DFN1    !reset temperature to last successful value, then
        GO TO 780  !go back to set this as minimum temperature, and leave
      END DO
C  
      CALL ANWARN(0) !if we are here, we did not coverge after 500 iterations
      WRITE(KLST,1320) S3,S1,BOOT  !write a warning message
      GO TO 670                    !and try once more
  820 CONTINUE        !done with melt regions
      ACK(IT+72)=PCO2 !restore iteration flag for ANELSM
      GO TO 870       !all done with melt--skip to final summary
C
C      Termination of melt transition section.  If no melt transiton
C      start at 830 and set ACK entries 43 to 52, inclusive, to zero.
C
  830 DO I=43,52
        ACK(IT+I)=ZERO
      END DO
C
C     For material types .le. 2, this is the end.  For other
C     types, loop back to melt transiton section
C
      IF (NINT(ACK(IT+30))-2) 870,870,660
  850 IF (I1.LT.0) GO TO 860
      ACK(IT+47)=S1
      GO TO 870
  860 ACK(IT+47)=ACK(IT+11)
  870 IF (NINT(DIN(30))+1) 8790,8780,8790
C
C     LINEAR SHOCK-PARTICLE VELOCITY EOS (TYPE -1)
C     OR IDEAL GAS EOS (TYPE -2)
C
 8700 ACK(IT+43)=ZERO
      ACK(IT+42)=ZERO
      DIN(10)=ZERO
      DIN(32)=ABS(UI(6))
      DIN(33)=UI(9)
      DIN(4)=UI(11)
      IF(DIN(4).LE.ZERO) DIN(4)=THREE*DIN(27)*BOLTS
      IF(NINT(DIN(30)).LE.-2) GO TO 8740
C     SOLID (TYPE -1)
      DIN(25)=ZERO
      DIN(20)=ZERO
      DIN(6)=-200.D0
      IF(DIN(15).LE.ZERO) DIN(15)=ONE
      DIN(5)=DIN(11)*DIN(15)*DIN(4)
      DIN(21)=DIN(11)*DIN(32)**2
      DIN(36)=ONE/(DIN(11)+DIN(11))
      DIN(35)=1.D10
      DIN(3)=DIN(4)*DIN(12)
      IF (DIN(18).LE.ZERO) DIN(18)=DIN(12)-DIN(18)/DIN(4)
      DIN(1)=ZERO
        DO 8820 I=1,40
 8820   ACK(IT+I)=DIN(I) 
      ACK(IT+48)=UI(13)       !constant C53 to move critical point
      ACK(IT+37)=UI(14)
      ACK(IT+44)=UI(18)
      ACK(IT+45)=UI(19)
      ACK(IT+53)=UI(20)
      ACK(IT+54)=UI(21)
      CALL ANN1AS(IT)
      IF(DIN(33).GT.ZERO) THEN
        SPS=DIN(33)-0.99
        IF(SPS.GT.ZERO) THEN
          RHUGMX=DIN(11)*DIN(33)/SPS
        END IF
      ELSE
          RHUGMX=MIN(ACK(IT+37),FIVE*ACK(IT+11))
      END IF
      GO TO 8790
 8810 CONTINUE
      IF (DIN(18).LE.ZERO) DIN(18)=DIN(12)-DIN(18)/DIN(4)
      PCO1=UI(18)
      PCO2=UI(19)
      BOOT=UI(20)
      ETAOT=UI(21)
      DIN(1)=ZERO
      GO TO 550
C
C     IDEAL GAS (TYPE -2)
C
 8740 IF(UI(11).LE.ZERO) DIN(4)=HALF*DIN(4)
      DIN(33)=ZERO
      DIN(32)=ZERO
      DIN(5)=DIN(15)*DIN(4)
      DIN(3)=ONE+DIN(15)
      DIN(6)=LOG(DIN(12))-DIN(15)*LOG(DIN(11))
      DIN(20)=UI(5)*DIN(25)**3
      IF(NINT(DIN(30)).EQ.-2) GO TO 550
C
C     LLL JWL HIGH EXPLOSIVE EOS
C     SCALED INPUT NUMBERS ASSUMED - SEE UCID-16189
C
      S1=DIN(11)
      IF(UI(11).LE.ZERO) DIN(4)=TWO*DIN(4)
      DIN(5)=S1*UI(20)
      DIN(6)=S1*UI(21)
      DIN(32)=QCL45*UI(18)
      DIN(33)=QCL45*UI(19)
      BOOT=(DIN(32)*EXP(-DIN(5)/S1)+DIN(33)*EXP(-DIN(6)/S1)
     1     +DIN(15)*DIN(4)*DIN(12)*S1)/(S1**(DIN(15)+ONE))
      DFN1=S1
      DFN2=QCL46*S1
 8800 DIN(19)=HALF*(DFN1+DFN2)
      DIN(35)=(DIN(32)*EXP(-DIN(5)/DIN(19))+DIN(33)*
     1        EXP(-DIN(6)/DIN(19)))/(DIN(19)**(DIN(15)+ONE))
      IF(DFN2-DFN1.LE.QCL31*DIN(19)) GO TO 8806
      IF(DIN(35)-BOOT) 8802,8806,8804
 8802 DFN1=DIN(19)
      GO TO 8800
 8804 DFN2=DIN(19)
      GO TO 8800
 8806 DIN(35)=DIN(35)/DIN(15)
      DIN(3)=DIN(35)*DIN(19)**DIN(15)-DIN(32)*EXP(-DIN(5)/DIN(19))
     1/DIN(5)-DIN(33)*EXP(-DIN(6)/DIN(19))/DIN(6)
C
C     MODIFICATION OF JWL FORM FOR VERY HIGH DENSITIES
C     AND MIXED-PHASE CONSTRUCTION
C
      DIN(10)=QCL26*DIN(11)
      DIN(24)=3.2D12
      S2=DIN(32)*EXP(-DIN(5)/DIN(10))
      S3=DIN(33)*EXP(-DIN(6)/DIN(10))
      S4=DIN(10)**DIN(15)
      S5=S2+S3-DIN(10)*S4*DIN(15)*DIN(35)
      S6=S2/DIN(5)+S3/DIN(6)-S4*DIN(35)+DIN(3)
      S2=DIN(10)**ATHIRD
      DIN(13)=(S5-DIN(24)*DIN(10)*S2*S2)/(DIN(10)*S2)
      DIN(14)=S6-THREE*(HALF*DIN(24)*S2+DIN(13))*S2
      ACK(IT+53)=DIN(10)
      ACK(IT+54)=DIN(13)
      DIN(10)=DIN(3)
      DIN(13)=7.77777D-5
      ACK(IT+42)=0.1D0
      DIN(30)=FOUR
      IF(DIN(22).EQ.ZERO) THEN
        DIN(30)=THREE
        ACK(IT+63)=ZERO
      ELSEIF(DIN(22).GT.ZERO) THEN
        ACK(IT+63)=ONE
      ELSE
        ACK(IT+63)=-ONE
      END IF
      DIN(22)=ZERO
      DIN(25)=FLGJWL
      WRITE(KLST,8808) DIN(30)
 8808 FORMAT (/,'  JWL form changed to type',F3.0,
     1          ' for mixed-phase construction',/)
      GO TO 8810
 8780 ACK(IT+35)=DIN(18)
      ACK(IT+44)=PCO1
      ACK(IT+45)=PCO2
      ACK(IT+53)=BOOT
      ACK(IT+54)=ETAOT
C
C----------------------
C     ALL FORM HERE
C----------------------
C
 8790 CONTINUE
C       Watch out for negative reference temperature or density
      IF ((ACK(IT+12).LE.ZERO).OR.(ACK(IT+11).LE.ZERO)) GO TO 900
C
C    call ANEOSV for reference conditions
C
      CALL ANEOSV(1,ACK(IT+12),ACK(IT+11),MAT,S1,S2,S3,S4,S5,S6,
     1    DIN(1),DIN(2),I,EXTRAA,EXTRAB,EXTRAC)
      DIN(3)=ACK(IT+11)*S6   !rho0 * DPdR = DPD(eta) = B0
      S9=S5/(S4*ACK(IT+11))  !Gamma
C
C      write out reference point conditions
C
      WRITE(KLST,1530) ACK(IT+12),ACK(IT+11),S1,S2,S5,S4,S3,S6,DIN(3),
     1    DIN(2),S9
C
  900 SPS=ACK(IT+30)
      IF(KUWAY.EQ.1) GO TO 920
C
C    Compute critical point and liquid/vapor phase curve for later
C    use by ANEOS1.  ANPHAS then echos this curve to the output file
C
      CALL ANPHAS (MAT,IT,IKPN)
C
      IF (SPS-ACK(IT+30)) 910,912,910
  910 IF (ACK(IT+46).LE.ZERO) GO TO 920
      WRITE(KLST,1330) MAT    !recalculation of material MAT
      UI(2)=ACK(IT+30)
      GO TO 50
  912 IF(TCT(MAT).LE.ZERO) GO TO 920   !critical temperature .le. 0
      BOOT=QCL47*TCT(MAT)
      PCO2=ACK(IT+11)
      PCO1=QCL34*PCO2
      SQTS(1)=SQRT(BOOT)
        DO 914 I1=1,99
        CALL ANEOS1 (BOOT,PCO2,S1,S2,S3,S4,S5,S6,IT)
        IF (S6.LT.ZERO) GO TO 916
  914   PCO2=PCO2-PCO1
      GO TO 920
  916 CONTINUE
      CALL ANWARN(0)
      WRITE(KLST,1380)      !eos is not properly defined!
      WRITE(KLST,1382)      !possibly because of multiple critical points
  920 LOCKPL(MAT)=IKPN-1    !pointer to end of vapor/liq table
      IF (TCT(MAT).GE.TWO) THEN  !critical temperature is more than 2 eV
        CALL ANWARN(0)
        WRITE(KLST,1384)    !warn user that critical temperature is too high
      END IF
C     Evap in ev/atom
      IF(ACK(IT+30).GE.ZERO) THEN
        ACK(IT+78)=ACK(IT+10)/(ACK(IT+27)*BOLTS)
        IF(ACK(IT+78).GE.CNINE) THEN
          CALL ANWARN(0)
          WRITE(KLST,924)
  924     FORMAT(' Evap is probably too large')
        END IF
      END IF
C
C     list computational array--start by printing material #
C
      WRITE(KLST,1570) MAT
C     restore energy & entropy shifts for print
      ACK(IT+64)=EMOVER
      ACK(IT+65)=SMOVER
      I2=0
        DO 882 I=1,99
        IF(ACK(IT+I).NE.ZERO) I2=I
  882   CONTINUE
C
C    Echo contents of ACK array to output
C
      I2=(MAX(I2,10)+2)/3
      J1=3*I2
        DO 880 I=1,I2
        WRITE(KLST,1510) (I1,ACK(IT+I1),I1=I,J1,I2)
  880   CONTINUE
      IF(KUWAY.EQ.1) GO TO 1210
C
C       write out pointers to beginning and end of vapor/liquid-solid
C       phase curve storage array.
C
      WRITE(KLST,1540) MAT,LOCKP(MAT),MAT,LOCKPL(MAT)
C     remove energy & entropy shifts again
      ACK(IT+64)=ZERO
      ACK(IT+65)=ZERO
C
C       WRITE TABLE OF SOLID/MELT PHASE CURVE TO OUTPUT FILE
C         and, at the same time, check that convergence is good for the
C         entire phase curve from temperature C49 (min) to C48 (max).
C
      IF (ACK(IT+46).LE.ZERO) GO TO 1190    !Skip this if negative solid melt density
C
C      Write out parameters of region bounding melt/solid phase curve
C
      WRITE(KLST,1345)  UI(23),ACK(IT+49)*TCONV,ACK(IT+18)*TCONV,
     1 ACK(IT+52)*TCONV,ACK(IT+48)*TCONV,ACK(IT+50)*TCONV/DCONV,
     2 ACK(IT+51)*TCONV/DCONV,ACK(IT+23)*DCONV,ACK(IT+46)*CONV,
     3 ACK(IT+47)*DCONV

      ! =========================================================================
      ! =========================================================================
      ! New code to construct table containing solidus and liquidus
      ! 
      ! Note: ACK(IT+49) Is the low-temperature limit of the melt curve
      !       ACK(IT+18) Is the triple point (zero pressure melting)
      !       ACK(IT+48) Is the highest temperature allowed
      !
      ! The tables are:
      !       TMC(IT+1->1000): Temperature along the melting curves
      !       RSMC(IT+1->1000): Density along the solidus
      !       RLMC(IT+1->1000): Density along the liquidus
      !
      ! The tables are constructed so that 100 points lie below the melt
      ! temperature at zero pressure and 900 points lie above; the temperature 
      ! interval is linear. . .
      ! =========================================================================
      IF (MELT_CURVE_TABLE.EQ.1) THEN 
         I1=0
         ACK(IT+71)=ZERO
         CMLT7=ACK(IT+71)
         ACK(IT+72)=ONE         ! Force ANELSM to use slow iteration
         WRITE(KLST,1365)       ! Add note that melt curve table created
         IF(KEXPRT.NE.0) WRITE(KLST,1370) ! write header for melt curve output
         BOOT=ACK(IT+49)        ! define initial temp. for melt curves
         DTBELOW=(ACK(IT+18)-ACK(IT+49))/99.D0 ! 99 linear steps between Tmin and Tm
         DTABOVE=(ACK(IT+48)-ACK(IT+18))/900.D0 ! 900 linear steps between Tm and Tmax
         GEORAT=(ACK(IT+48)/ACK(IT+18))**(1.D0/900.D0)
         I3=0
         I4=0
         S7=ACK(IT+23)          ! Initial guess of density is ro_min
         DO WHILE (BOOT.LT.ACK(IT+48)) ! Loop over temperature until upper bound
            I3=I3+1
            I2=2                ! Raise debugging flag for ANELSM
            SQTS(1)=SQRT(BOOT)
            CALL ANELSM (BOOT,S7,S1,S2,IT,I2) ! get density of liquid (S1), solid (S2)
            S7=S2               ! update density guess with solid density
            
            ! Flag the table entries just below and above the solid-solid phase
            ! transition. . .
            IF (S2.LE.ACK(IT+1)*ACK(IT+19)) ITLS1 = I3 ! Solid density below lpp transition
            IF (S1.LT.ACK(IT+2)*ACK(IT+19)) ITLS2 = I3+1 ! Liquid density above hpp transition

            IF (I2.GE.0) THEN   ! no error from ANELSM
                                ! Store results in table   
               I4=I4+1
               IF (KEXPRT.NE.0) THEN
                  IF (MOD(I4-1,10).EQ.0) 
     $                 WRITE(KLST,1375) BOOT*TCONV,S1*DCONV,S2*DCONV,I2
               END IF
               TMC(IT+I4)=BOOT  ! Temperature on melt curves
               RSMC(IT+I4)=S2   ! Solid density
               RLMC(IT+I4)=S1   ! Liquid (melt) density
            ELSE
               I1=I1+1          ! Error counter
               IF (KEXPRT.NE.0) WRITE(KLST,*)
     $              "Suspected convergence failure in ANELSM."//
     $              " Trying to skip this temperature: ",BOOT*TCONV
            END IF
 
            ! Increment the temperature differently above and below the triple point
            IF (I3 < 100) THEN  ! Temp. below the melt temperature
               BOOT=BOOT+DTBELOW
            ELSEIF (I3 > 100) THEN ! Temp. above the melt temperature
               !BOOT=BOOT+DTABOVE
               BOOT=ACK(IT+18)*GEORAT**(I3-100)
            ELSE
               BOOT=ACK(IT+18)  ! The melt temperature
            END IF

         END DO
         NUMTMC(IT+1)=I4-1
         WRITE(KLST,*)"Number of points in melt table: ",NUMTMC(IT+1)

         ! Report error if EoS is not properly defined
         IF (I1.GE.1) THEN
            WRITE(KLST,*)
            WRITE(KLST,*)"Had to skip points in the melt curve table"
            WRITE(KLST,*)"Try reducing step size of it. in ANELSM"
            WRITE(KLST,*)"to improve convergence (QCI4)."
            WRITE(KLST,*)"Skipped entries in melt table:  ",I1
            I1 = 1
         ENDIF

         ! If a solid-solid phase transition exists, 
         ! fix the doubly-mixed region with a linear interpolation
         IF (ACK(IT+1).GT.ZERO) THEN
            WRITE(KLST,*)
            WRITE(KLST,*) "Overwriting melt curve table across "//
     $           "melt-transition between solid/liquid densities "//
     $           "from: ",RSMC(IT+ITLS1),RLMC(IT+ITLS1),
     $           "to: ",RSMC(IT+ITLS2),RLMC(IT+ITLS2)
            DO I3 = ITLS1,ITLS2
               PHI = (TMC(IT+I3)-TMC(IT+ITLS1))/
     $              (TMC(IT+ITLS2)-TMC(IT+ITLS1))
               RSMC(IT+I3) = RSMC(IT+ITLS1) + 
     $              PHI*(RSMC(IT+ITLS2)-RSMC(IT+ITLS1))
               RLMC(IT+I3) = RLMC(IT+ITLS1) + 
     $              PHI*(RLMC(IT+ITLS2)-RLMC(IT+ITLS1))
            END DO
         END IF

      END IF
      ! End of new code to contruct table containing melt curve information
      ! =========================================================================
      ! =========================================================================
C
      DFN1=ACK(IT+72)        !save current state of ANELSM flag
      I1=0                   !initialize error flag
      I4=1                   !initialize melt curve table index
      DFN2=ZERO              !initialize line number counter
      PCP1=ZERO              !will contain material number if error occurs
      ACK(IT+71)=ZERO        !add melt terms
      CMLT7=ACK(IT+71)       !enable addition of melt terms in ANEOS1
      ACK(IT+72)=ONE         !set slow iteration flag for ANSLEM = 1
      WRITE(KLST,1350)       !write header for melt curve output
      BOOT=ACK(IT+49)        !C49 = min initial temperature for melt curve
C
  930 CONTINUE               !start loop for printing liquid phase table
      I2=2                   !set process flag to return # iterations
      S7=1.D10               !huge initial density ensures solid result
      SQTS(1)=SQRT(BOOT)
C
C     If not using tabular melt curve, then get density of solid S2, liquid S1 using ANELSM
C     Otherwise, loop through the table. . .
C
      IF (MELT_CURVE_TABLE.EQ.0) THEN
         CALL ANELSM (BOOT,S7,S1,S2,IT,I2)
      ELSE
         BOOT=TMC(IT+I4)
         S1=RLMC(IT+I4)
         S2=RSMC(IT+I4)
         I4=I4+10        ! Probe every tenth table entry
      ENDIF
C
      IF (I2.LT.0) THEN       !ANELSM fatal error condition
        I1=1                  !set error flag
        GO TO 1110            !process error
      END IF
      ACK(IT+71)=-ONE         !do not add melt terms
      CMLT7=ACK(IT+71)        !disable addition of melt terms in ANEOS1
      SQTS(1)=SQRT(BOOT)
C
C      use ANEOS1 to compute solid phase properties at
C      density S2 (from previous call to ANELSM), temperature BOOT
C
      CALL ANEOS1 (BOOT,S2,S3,S4,S5,S6,S7,S8,IT)  !solid phase properties
C
      S9=S4-BOOT*S5+S3/S2     !solid gibbs free energy, E-T*S+P/RHO
      S3S=S3                  !solid pressure
      S4S=S4                  !solid energy
      S5S=S5                  !solid entropy
      S9S=S9                  !solid gibbs energy
      ACK(IT+71)=ZERO         !add melt terms
      CMLT7=ACK(IT+71)        !enable addition of melt terms in ANEOS1
C
C      use ANEOS1 to compute melt phase properties at
C      density S1 (from previous call to ANELSM), temperature BOOT
C
      CALL ANEOS1 (BOOT,S1,S3,S4,S5,S6,S7,S8,IT)  !liquid phase properties
C
      S9=S4-BOOT*S5+S3/S1     !liquid gibbs free energy, E-T*S+P/RHO
C
C      WRITE OUT LINE IN TABLE DESCRIBING MELT/SOLID PHASE CURVE
C
      WRITE(KLST,1360) BOOT*TCONV,S1*DCONV,S2*DCONV,S3*PCONV,S3S*PCONV,
     & S4*ECONV,S4S*ECONV,S5*SCONV,S5S*SCONV,S9*ECONV,S9S*ECONV,I2

C
C          CHECK FAST ITERATION
C
C       This is a very comprehensive check that evaluates the material
C       state at nine different densities for each temperature.
C       The densities are chosen in three groups: The first is certain to
C       be in the liquid field, the second between the solid and liquid
C       phase curves (mixed state) and the third in the solid field.  In
C       each group a density is chosen close to the low density limit, 
C       halfway between the low and high density limits and just below the
C       high density limit.  If ANEOSV, using ANELSM in the fast mode, fails
C       to return the expeced phase designation an error message is printed.
C          This check is performed for each temperature that is printed
C       in the table.  If more temperatures are desired, decrease the
C       value of parameter QCL49.
C
C       For melt region phase boundary definitions see Thompson and
C       Lauson 1972, section V-1.
C
      IF (MELT_CURVE_TABLE.EQ.0) THEN
         ACK(IT+72)=ZERO        ! Use fast iteration if no table defined
      END IF
      DO I2=1,3
         IF(I2.EQ.1) THEN       !liquid group density limits
            PCO1=ACK(IT+47)+(BOOT-ACK(IT+18))/ACK(IT+50) !use curve A lhs for T > Tm
            IF(BOOT.EQ.ACK(IT+18)) CYCLE !don't bother for T = Tm
            IF(BOOT.LT.ACK(IT+18)) PCO1=ACK(IT+23) !use lhs = RHOMIN for T < Tm
            PCO2=S1             !rhs limit is liquid
         ELSEIF(I2.EQ.2) THEN   !mixed phase group density limits
            PCO1=S1             !lhs is RL
            PCO2=S2             !rhs is RS
         ELSEIF(I2.EQ.3) THEN   !solid group density limits
            PCO1=S2             !lhs is RS
            PCO2=S2+TEN         !rhs is RS plus 10 gm/cc
         END IF
         IF (QCL48*PCO1.GE.PCO2) CYCLE !upper and lower limits too close
         DO IKJ=1,3             !loop over three points within the density ranges
            IF(IKJ.EQ.1) PCP2=0.99D0 !left fraction of density range
            IF(IKJ.EQ.2) PCP2=HALF !middle fraction of density range
            IF(IKJ.EQ.3) PCP2=0.01D0 !right fraction of density range
            PCP2=PCP2*PCO1+(ONE-PCP2)*PCO2 !density for evaluation
            IF (PCP2.LE.ACK(IT+23)) CYCLE !skip if it is below minimum
C     
C        call ANEOSV in fast execution mode to make sure it works
C        JJ returns the phase code number: JJ = 4(solid), 5(mixed), 6(melt)
C        BOOT is the input temperature, PCP2 is the density.
C
            CALL ANEOSV(1,BOOT,PCP2,MAT,S3,S4,S5,S6,S7,S8,S9,SPS,JJ,
     &           EXTRAA,EXTRAB,EXTRAC)
C
            IF((I2.EQ.1).AND.(JJ.EQ.6)) CYCLE !success
            IF((I2.EQ.2).AND.(JJ.EQ.5)) CYCLE !success
            IF((I2.EQ.3).AND.(JJ.EQ.4)) CYCLE !success
C              wrong result!  write out an error message
            CALL ANWARN(0)
            WRITE(KLST,1340) MAT,PCP2,BOOT,JJ,IKJ,I2,PCO1,PCO2
            PCP1=MAT            !save material number
         END DO
      END DO
C
      ACK(IT+72)=ONE  !back to slow iteration for next property evaluation
 1110 DFN2=DFN2+ONE   !line count. Also, come here on error
C
C       choose new temperature for next line of output
C
      IF(DFN2.LT.THREE) THEN      !three lines for T in tensional region
        BOOT=BOOT+(ACK(IT+18)-ACK(IT+49))/THREE
      ELSEIF(DFN2.EQ.THREE) THEN  !single line for T = Tm (fourth line)
        BOOT=ACK(IT+18)
      ELSEIF(DFN2.GT.THREE) THEN  !logarithmic spacing for T > Tm to limit
        BOOT=BOOT*(ACK(IT+48)/ACK(IT+18))**QCL49
      END IF
C
C     Loop back to next temperature or entry in the table, until table is
C     exhausted or maximum temperature is reached. . .
C
      IF (BOOT.LT.ACK(IT+48).AND.I4.LT.1000) GO TO 930
C
      ACK(IT+72)=MAT                   !error has occurred, save material number
      IF(PCP1.EQ.ZERO) ACK(IT+72)=DFN1 !restore initial state of flag 72
      IF (I1.EQ.1) THEN       !error condition was been found in full table
         CALL ANWARN(0)
         WRITE(KLST,1380)     !eos is not properly defined!
      END IF

 1190 CONTINUE  !on to check solid-solid phase transition

      ! =========================================================================
      ! =========================================================================
      ! New code to construct table containing LPP and HPP
      ! 
      ! The tables are:
      !       THPP(IT+1->1000): Temperature along the phase transition
      !       RLPP(IT+1->1000): Density along the low-pressure phase boundary
      !       RHPP(IT+1->1000): Density along the high-pressure phase boundary
      !
      ! =========================================================================
      IF (DIN(30).EQ.5) THEN  ! Current material is a HPP of previous material 
         I1=0
!         ACK(IT+71)=ZERO
!         CMLT7=ACK(IT+71)
!         ACK(IT+72)=ONE         ! Force ANELSM to use slow iteration
         WRITE(KLST,1371)       ! Add note that melt curve table created
         IF(KEXPRT.NE.0) WRITE(KLST,1372) ! write header for melt curve output
         BOOT=1.D0/TCONV        ! Start phase boundary at 100 K
         GEORAT=(10./BOOT)**(1.D0/1000.D0) ! common ratio between table temps.
         I3=0
         IMAT=(IT-99)/99+1      ! Material number of LPP
         I4=(IMAT-1)*1000       ! Index of LPP material.
         S7=ACK(IT+11)          ! Initial guess of density is ro_min
         S1=ACK(IT-99+11)       ! Zero-kelvin density of LPP (previous material)
         S2=ACK(IT+11)*2.       ! Zero-kelvin density of HPP
c$$$         IF (IT.EQ.0) THEN
c$$$            ACK(IT+8)=0.D0
c$$$         ELSEIF (IT.EQ.99) THEN
c$$$            ACK(IT+8)=12.D9
c$$$         ELSEIF (IT.EQ.198) THEN
c$$$            ACK(IT+8)=.6D10
c$$$         END IF
         DO WHILE (BOOT.LT.10.) ! Loop over temperature until upper bound
            I3=I3+1
            I2=2                ! Raise debugging flag for ANEHPP
            SQTS(1)=SQRT(BOOT)
            CALL ANEHPP (BOOT,S7,S1,S2,IT,I2,ERROR_P,ERROR_G) ! get density of LPP (S1), HPP (S2)
            S7=S1               ! update density guess with LPP density
            !write(*,*)IMAT,BOOT*TCONV,S7,S1,S2,IT,I2
            !read(*,*)
            IF (I2.GE.0 .AND. S1.GT.1.D-6) THEN   ! no error from ANEHPP
                                ! Store results in table   
               I4=I4+1
               IF (KEXPRT.NE.0) THEN
                  IF (MOD(I4-1,10).EQ.0) 
     $                 WRITE(KLST,1375) BOOT*TCONV,S1*DCONV,S2*DCONV,I2,
     $                 ERROR_P*PCONV,ERROR_G*ECONV
               END IF
               THPP(I4)=BOOT  ! Temperature on phase boundary
               RLPP(I4)=S1    ! LPP density
               RHPP(I4)=S2    ! HPP density
            ELSE
               I1=I1+1          ! Error counter
               IF (KEXPRT.NE.0) WRITE(KLST,*)
     $              "Suspected convergence failure in ANEHPP."//
     $              " Trying to skip this temperature: ",BOOT*TCONV
            END IF
 
            ! Increment the temperature geometrically
            BOOT=1.D0*GEORAT**(I3)/TCONV

         END DO
         NUMHPP(IMAT)=I4-1
         WRITE(KLST,*)"Number of points in solid-phase table: ",I4-1

         ! Report error if EoS is not properly defined
         IF (I1.GE.1) THEN
            WRITE(KLST,*)
            WRITE(KLST,*)"Had to skip points in the sol-sol table"
            WRITE(KLST,*)"Try reducing step size of it. in ANEHPP"
            WRITE(KLST,*)"to improve convergence (QCI4)."
            WRITE(KLST,*)"Skipped entries in LPP/HPP table:  ",I1
            I1 = 1
         ENDIF

      END IF
      ! End of new code to contruct table containing melt curve information
      ! =========================================================================
      ! =========================================================================
C
C      LIST COLD TERMS IN EOS
C
      ACK(IT+71)=ZERO         !add melt terms
      CMLT7=ACK(IT+71)        !enable addition of melt terms for execution
      IF (NINT(ACK(IT+30)).EQ.2) GO TO 1210      !EOS type 2 only
      IF (NINT(ACK(IT+30)).LT.ZERO) GO TO 1210   !EOS type < 0
C
C      choose density step depending upon whether a solid-solid phase
C      transition is present or not.  ACK(IT+1) = eta1 of eq. 5.60
C
      IF (ACK(IT+1).LE.BIGNUM) THEN !write table
         WRITE(KLST,1410)       !write header for table
         S3=ACK(IT+11)          !reference density RHO0
C
C     compute density increment S2 for table
C
         IF (ACK(IT+1).EQ.BIGNUM) THEN
            S3=QCL50*S3         !no solid-solid phase transition
            S2=S3
         ELSE                   !solid-solid phase transition is present
            S2=(ACK(IT+2)*ACK(IT+19)-S3)/QCL51
         END IF
C  
         !ACK(IT+71)=-ONE        !do not add melt terms---GSC THIS LINE SEEMS TO CAUSE A PROBLEM
         CMLT7=-ONE             !disable addition of melt terms
         DO I=1,50              !loop over 50 densities
            BOOT=1.D-6          !very low temperature
            SQTS(1)=SQRT(BOOT)
            CALL ANEOS1 (BOOT,S3,S4,S5,S6,S7,S8,GAM,IT)
            S6=S3/ACK(IT+19)    !RHO/RHO00
            GIBB=S5+S4/S3       !Gibbs potential
C     
C     ZERO TEMPERATURE ISOTHERM
C     Outputs are: density, pressure, dP/dRHO,energy, eta, gibbs energy
C
            WRITE(KLST,1420) S3*DCONV,S4*PCONV,GAM*ECONV,S5*ECONV,S6,
     &           GIBB*ECONV
            S3=S3+S2            !increment density for next cycle
         END DO

      END IF
C
C       cold compression list for expanded edit
C
      IF(KEXPRT.NE.0) THEN
         SPS=ACK(IT+30)         !EOS type switch
         IF(NINT(ACK(IT+30)).GT.4) ACK(IT+30)=ACK(IT+30)-ONE   ! Convert EOS type 5 to type 4
         IF(NINT(ACK(IT+30)).GT.2) ACK(IT+30)=ACK(IT+30)-THREE
         S3=ACK(IT+11)*TEN      !start printout at ten times reference density
         DO I=1,101
            BOOT=1.D-6          !temperature at which EOS is evaluated
            SQTS(1)=SQRT(BOOT)
            CMLT7=-ONE          !set flag to skip addition of melt terms
C     
C     the following get redefined in ckeos only
C
            GAMMA=QCL56         !QCL56 = -1.11D22 (big negative number)
            PSI=QCL56           !these variables are transferred from ANEOS1
            THETA=QCL56         !through common /ANEDIS/
            CALL ANEOS1 (BOOT,S3,S4,S5,S6,S7,S8,GAM,IT)
            S6=S3/ACK(IT+11)    !S6 is now eta, not entropy
            GIBB=S5+S4/S3       !Gibbs potential
            IF(GAMMA.EQ.QCL56) THEN !This EOS did not compute GAMMA
               IF(I.EQ.1) WRITE(KLST,1410) !header for short table
C     
C     ZERO TEMPERATURE ISOTHERM
C     Outputs are: density, pressure, dP/dRHO,energy, eta, gibbs energy
C     
               WRITE(KLST,1420) S3*DCONV,S4*PCONV,GAM*ECONV,S5*ECONV,S6,
     &              GIBB*ECONV
            ELSE                !EOS type 0,1,3,4 with GAMMA
               IF(I.EQ.1) WRITE(KLST,1411) !header for long table
               PSI=BOOT/PSI
C     
C     ZERO TEMPERATURE ISOTHERM
C     Outputs are: density, pressure, dP/dRHO,energy, eta, gibbs energy 
C     gruneisen gamma, temperature at psi = 1, debye temperature
C     
               WRITE(KLST,1420) S3*DCONV,S4*PCONV,GAM*ECONV,S5*ECONV,S6,
     &              GAMMA,PSI*TCONV,THETA*TCONV,GIBB*ECONV
            END IF
            S3=S3*QCL52         !decrement density by .912 for next cycle
         END DO
         ACK(IT+30)=SPS
      END IF
 1210 IT=IT+99
      IZ=IZI+1
      ACK(IT+71)=ZERO   !add melt terms
      CMLT7=ACK(IT+71)  !enable addition of melt terms
      IF (THUG.LT.ZERO) THUG=DIN(12)   !use reference temp. for hugoniot
      IF (RHUG.LT.ZERO) RHUG=DIN(11)   !use reference density for hugoniot
C
C        compute hugoniot curve
C
      CALL ANHUG (MAT,RHUG,THUG,RHUGMX)
C
      IF((EMOVER.NE.ZERO).OR.(SMOVER.NE.ZERO)) THEN
C       restore energy & entropy shifts
        ACK(IT+64-99)=EMOVER
        ACK(IT+65-99)=SMOVER
        WRITE(KLST,1610) EMOVER, SMOVER
 1610   FORMAT(//,' Energy-Entropy shift employed for this material'
     &  ,/,' Delta Energy  =',1PE13.6
     &  ,/,' Delta Entropy =',E13.6,/)
      END IF
      IF(KEXPRT.NE.0) CALL ANPRTR   !list counter results for experts
      CALL ANZRTR  !zero all counters
 1220 CONTINUE  !end current definition and start loop for new material
      IF(IAMBAD.NE.0) THEN
        CALL ANMARK('ANEOS FATAL GENERATION ERROR')
        WRITE(KLST,1221) IAMBAD
 1221   FORMAT(' SEE PREVIOUS ERROR MESSAGE - IAMBAD=',I5)
        STOP
      END IF
C
C     END LOOP FOR EACH MATERIAL ***************************************
C
      CALL ANZRTR    !zero all counters
      IF((IZ.GT.30*MAXMAT).OR.
     1      (IT.GT.99*MAXMAT).OR.
     2      (IKPN.GT.100*MAXMAT)) GO TO 1230
      GO TO 2323  !return
 1230 CONTINUE
      CALL ANMARK('ARRAY OVERFLOW IN ANEOS')
      WRITE(KLST,1550) IZ,IT,IKPN,MAXMAT,30*MAXMAT,99*MAXMAT,100*MAXMAT
      STOP
C
C=======================================================================
C     WRITE RESTART DATA
 1240 CALL ANEOSS (IGK,ITAPE)
      GO TO 2323   !return
C
C=======================================================================
C     READ RESTART DATA
 1250 CALL ANEOSS (IGK,ITAPE)
      GO TO 2323   !return
C
C***********************************************************************
      ENTRY ANEXT2
C
C     set flag for extended setup print
C
      IEXPRT=1
      GO TO 2323   !return
 1260 FORMAT(' INPUT IN(',I2,') IS NOT ALLOWED',1PE13.6,
     &         ' RESET TO',E13.6)
 1270 FORMAT (//,35H  UNABLE TO INCLUDE MELT TRANSITION,5X,
     1 21HWILL CONTINUE WITHOUT)
 1280 FORMAT (/,34H  MELT TEMPERATURE INCREASED FROM ,1PD12.5,/)
 1290 FORMAT (//,7H  RHOL=,1PD13.6,4X,5HRHOS=,D13.6)
 1300 FORMAT (//,29H  HIGH TEMPERATURE MELT ERROR,1P3D15.6)
 1310 FORMAT (/,' IN(17) has been increased from',1PD12.5,
     1' to',D12.5,' for the melt transition',/)
 1320 FORMAT (//,28H  LOW TEMPERATURE MELT ERROR,1P3D15.6)
 1330 FORMAT (29H1 RECALCULATION OF EOS NUMBER,I5)
 1340 FORMAT ('  FAST ANELSM ITERATION FAILURE MAT=',I5,/,6H  RHO=,
     & 1PD14.7,3X,2HT=,D14.7,2X,7HKPHASE=,I5,2I7,2D15.7)
 1345 FORMAT (//' BOUNDS OF LIQUID/SOLID PHASE TRANSITION',/
     $ '   ENTHALPY OF FUSION (ERG/GM)...................',1PD12.5/,
     1 '   LOWER LIMIT OF METASTABLE LIQUID ACK(49)......',1PD12.5/,
     2 '   MELT TEMPERATURE (TRIPLE POINT) ACK(18).......',1PD12.5/,
     3 '   INTERMEDIATE TEMPERATURE ACK(52)..............',1PD12.5/,
     4 '   MAXIMUM TEMPERATURE ACK(48)...................',1PD12.5/,
     5 '   SLOPE OF LEFT  BOUNDARY, dT/dRHO ACK(50)......',1PD12.5/,
     6 '   SLOPE OF RIGHT BOUNDARY, dT/dRHO ACK(51)......',1PD12.5/,
     7 '   MINIMUM DENSITY ACK(23), GM/CC................',1PD12.5/,
     8 '   SOLID/MELT DENSITY AT Tm ACK(46), GM/CC.......',1PD12.5/,
     9 '   LIQUID/VAPOR DENSITY AT Tm ACK(47), GM/CC.....',1PD12.5/)
 1350 FORMAT (//' LIQUID/SOLID PHASE CURVE',//
     1 '       T         RLIQ       RSOLID      PLIQ       PSOLID',
     2 '      ELIQ        ESOLID       SLIQ       SOLID',
     3 '        GLIQ       GSOLID        #ITER'/
     4 '       K        kg/m**3     kg/m**3      GPa         GPa',
     5 '       J/kg         J/kg        J/kg-K     J/kg-K       J/kg',
     6 '        J/kg'/)
 1360 FORMAT(1P11D12.4,I10)
 1365 FORMAT (//'Constructing melt curve table for use',
     $     ' with solid-solid transition'//)
 1368 FORMAT( //' Second attempt at constructing melt table',
     $     ' with smaller density step size'//)
 1370 FORMAT (//' SUMMARY OF LIQUID/SOLID PHASE CURVE TABLE',//
     1 '       T         RLIQ       RSOLID      #ITER'/
     4 '       K        kg/m**3     kg/m**3  '/)
 1371 FORMAT (//'Constructing solid-solid phase transition table'//)
 1372 FORMAT (//' SUMMARY OF LPP/HPP PHASE CURVE TABLE',//
     1 '       T         RLPP       RHPP      #ITER'/
     4 '       K        kg/m**3     kg/m**3  '/)
 1375 FORMAT(1P3D12.4,I10,1P2D12.4)
 1380 FORMAT (
     1 ' THIS EOS IS NOT PROPERLY DEFINED AND SHOULD NOT BE USED')
 1382 FORMAT (/,' POSSIBILITY OF MULTIPLE CRITICAL POINTS')
 1384 FORMAT (/,' CRITICAL TEMPERATURE IS TOO HIGH',/,
     1 ' EOS MIGHT REQUIRE REDEFINITION',/)
 1390 FORMAT (/,'  MELT TEMPERATURE ERROR',/,1P7D13.4,I6)
 1400 FORMAT (1X,8D13.6)
 1410 FORMAT (//,'   ZERO-TEMPERATURE ISOTHERM',//
     1'        RHO          P         DPDR          E          ',
     2'ETA          GIBBS'/
     3'      kg/m**3       GPa	    (m/sec)**2     J/kg        ',
     4'             J/kg'/)
 1411 FORMAT (//,'   ZERO-TEMPERATURE ISOTHERM',//
     1'        RHO          P         DPDR          E          ',
     2'ETA       GAMMA      T(PSI=1)     THETA       GIBBS'/
     3'      kg/m**3       GPa	    (m/sec)**2     J/kg        ',
     4'                        K           K          J/kg'/)
 1420 FORMAT (2X,1P9D12.4)
 1440 FORMAT(//,'  Data for ANEOS number',I4,3X,'library number'
     1 ,I5,2X,'type',I3,5X,'ver ',A,//,2X,A,6X,
     2 'RHUG=',1PD12.5,/,58X,'THUG=',D12.5,/)
 1450 FORMAT(//,'  ERROR IN MATCHING OF MATERIAL AND EOS DATA',//,
     1 2X,'+ MATERIAL VALUES HAVE BEEN LOCATED',/,
     2 2X,'- MATERIAL VALUES HAVE NOT BEEN PROCESSED',/,
     3 2X,'CURRENT EOS NUMBER=',I8,//,2X,'MATERIALS',/,(5X,I8))
 1460 FORMAT (8D10.3)
 1470 FORMAT (3(5H  IN(,I2,2H)=,1PD16.9))
 1480 FORMAT (1X)
 1490 FORMAT (/,' THE IONIZATION POTENTIALS FOR Z=',I4,
     & ' ARE NOT IN TABLE')
 1510 FORMAT (3(4H  C(,I2,3H) =,1PD16.9))
 1520 FORMAT ('  Z(',I2,')=',F4.0,'  number fraction(',I2,')='
     1 ,1PD12.5,'  atoms/gram(',I2,')=',D12.5)
 1530 FORMAT(//,'  REFERENCE POINT CONDITIONS',//,
     1'  T=',1PD15.6,6X,'RHO=',D16.6 ,6X,'P=' ,D15.6,/,
     2'  E=',D15.6,  6X,'DPDT=',D15.6,6X,'CV=',D14.6,/,
     3'  S=',D15.6,  6X,'DPDR=',D15.6,6X,'B0=',D14.6,/,
     4'  CS=',D14.6, 6X,'GAMMA=',D14.6)
 1540 FORMAT (/,'  LOCKP(',I2,2H)=,I4,11H    LOCKPL(,I2,2H)=,I4,
     12X,'(two-phase array storage data)')
 1550 FORMAT ('1 ARRAY OVERFLOW IN ANEOS',/,' IZ=',I5,5H  IT=,I5,
     1  7H  IKPN=,I6,/,' LIMITS FOR',I5,' MATERIALS ARE',/,
     2  ' IZ=',I5,5H  IT=,I5,7H  IKPN=,I6)
 1560 FORMAT (5(F5.0,D10.3))
 1570 FORMAT(//,'  DATA ARRAY FOR MATERIAL',I4,/)
 2323 RETURN
      END
