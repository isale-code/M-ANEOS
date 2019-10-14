      PROGRAM ANEOSTEST
C
C      PROGRAM TO INITIALIZE ANEOS AND CALL IT ONCE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER STYLE
      PARAMETER (NDENS=0)
      PARAMETER (NTEMP=0)
      DIMENSION IZETL(21)
      COMMON /FILEOS/ KLST, KINP
      COMMON /ANGELX/ W,Y,ALFA,BETA,DADT,DADR,DBDT,ZZ,ET,
     & SN,CVN,EN,PN,PM,EM,SM,CVM,EC,PC,RHO0,RHO00,ZVIB,
     & SMLT,CVMLT,EMLT,PMLT,H1,H2
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
     &  ,PE,EE,SE,CVE,DPTE,DPRE
     &  ,NMATSX,IIZX
      COMMON /ANEDIS/ GAMMA,PSI,THETA
      PARAMETER (MATBUF=64)
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
      DIMENSION RHO(NDENS),TEMP(NTEMP),HISTDENS(100),HISTPRESS(100)
      DIMENSION HISTDPDR(100)
      DIMENSION TARR(2000)
      DIMENSION DARR(2000)
      DIMENSION PARR(2000,2000)
      DIMENSION SARR(2000,2000)
      DIMENSION EARR(2000,2000)
      DIMENSION CSARR(2000,2000)
      DIMENSION CVARR(2000,2000)
      DIMENSION ZKARR(2000,2000)
C
CSTS ADDS --help, --no_table, and file checks
      character(:), allocatable :: arg
      integer arglen, stat
      logical file_exists
      call get_command_argument(number=1, length=arglen)  ! Assume for simplicity success
      allocate (character(arglen) :: arg)
      call get_command_argument(number=1, value=arg, status=stat)
      IF (arg.EQ.'--help') THEN
           WRITE(*,9100)
           WRITE(*,9101)
           WRITE(*,9102)
           WRITE(*,9103)
           WRITE(*,9104)
           WRITE(*,9105)
           WRITE(*,9106)
           WRITE(*,9107)
           WRITE(*,9108)
           WRITE(*,9109)
           WRITE(*,9110)
           WRITE(*,9111)
           WRITE(*,9112)
           WRITE(*,9113)
           WRITE(*,9114)
           WRITE(*,9115)
           WRITE(*,9116)
           WRITE(*,9117)
           WRITE(*,9118)
           WRITE(*,9119)
           WRITE(*,9120)
           WRITE(*,9121)
           WRITE(*,9122)
           call EXIT(0)
      ENDIF
      INQUIRE(FILE="ANEOS.INPUT", EXIST=file_exists)
      IF (file_exists) THEN
C         WRITE(*,*)'ANEOS.INPUT exists'
      ELSE
         WRITE(*,*)'Missing ANEOS.INPUT file'
         WRITE(*,*)'See --help'
         call EXIT(0)
      ENDIF
      IF (arg.NE.'--no_table') THEN
        INQUIRE(FILE="tablegrid.txt", EXIST=file_exists)
        IF (file_exists) THEN
C         WRITE(*,*)'tablegrid.txt exists'
         ELSE
           WRITE(*,*)'Missing tablegrid.txt file'
           WRITE(*,*)'Use --no_table flag to run without file'
           WRITE(*,*)'See --help'
           call EXIT(0)
        ENDIF
      ENDIF
CSTS  END --help, --no_table checks
C
      OPEN(10,FILE='ANEOS.INPUT',STATUS='OLD')
      OPEN(12,FILE='ANEOS.OUTPUT')
C     corrected STSM 5/31/10=9
      TCONV=1.16045D4
C
C     DEFINE ARRAYS OF DENSITY, TEMPERATURE
C
      DO 10 I=1,NDENS
      RHO(I)=10.D0**(1-I)
C      RHO(I)=1.11D0    !REFERENCE DENSITY
  10  CONTINUE
C
      DO 20 J=1,NTEMP
      TEMP(J)=1.D3*DFLOAT(J)/TCONV
C      TEMP(J)=2.008991D-2
  20  CONTINUE
C
C     INITIALIZE ANALYTIC EOS
C
      KLST=12
      KINP=10
      M=0
      IZETL(1)=-1
      CALL ANEOS2 (1,1,0,IZETL)
C
C     ========================== STSM NEW EOS TABLE CONSTRUCTION ====================
      IF (arg.EQ.'--no_table') THEN
         call EXIT(0)
      ENDIF
C     READ IN TABLE GRID POINTS FROM FILE
C     FORMAT is a single column list of nden (int), ntemp (int), denarray, temparray. Arrays formatted as .6E in python
C     loops in this section begin at 50
C
      SESMATID=1.0
      DATE=1.0
      VERSION=1.0
      FMN=1.0
      FMW=1.0
      RHO0REF=1.0
      SESK0REF=1.0
      T0REF=1.0
      DCOUNT=1.0
      TCOUNT=1.0

      
      OPEN(13,FILE='tablegrid.txt',STATUS='OLD')
      READ(13, '(E12.6)') SESMATID
      READ(13, '(E12.6)') DATE
      READ(13, '(E12.6)') VERSION
      READ(13, '(E12.6)') FMN
      READ(13, '(E12.6)') FMW
      READ(13, '(E12.6)') RHO0REF
      READ(13, '(E12.6)') SESK0REF
      READ(13, '(E12.6)') T0REF
      READ(13, '(E12.6)') DCOUNT
      READ(13, '(E12.6)') TCOUNT
      DO 50 I=1,DCOUNT
      READ(13, '(E12.6)') DARR(I)
 50   CONTINUE
      DO 51 I=1,TCOUNT
      READ(13, '(E12.6)') TMP
      TARR(I) = TMP/TCONV
 51   CONTINUE
      CLOSE(13)
C
C     CALL ANEOSV TO FILL IN THE ARRAYS OF THERMODYNAMIC VARIABLES
C
C      PREF=1.D6      !REFERENCE PRESSURE =  1 ATMOSPHERE
C      TOL=1.D-7      !PRECISION OF CONVERGENCE TO REF PRESSURE
      
      DO 61 I=1,DCOUNT
      DO 60 J=1,TCOUNT
      CALL ANEOSV (1,TARR(J),DARR(I),1,P,E,S,CV,DPDT,DPDR,FKRO,CS,KPA,
     &  R2PL,R2PH,ZBAR)
      IF (ABS(P) .GE. 1E-20) THEN
         PARR(I,J) = P
      ELSE
         PARR(I,J) = 1.E-20
      ENDIF
      SARR(I,J) = S
      EARR(I,J) = E
      CSARR(I,J) = CS
      CVARR(I,J) = CV
      ZKARR(I,J) = KPA
 60   CONTINUE
 61   CONTINUE


C     WRITE SESAME TABLE TO FILE
      OPEN(14,FILE='NEW-SESAME-EXT.TXT')
C     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
C     could input matid, date, version with the grid
      NWDS=9
      SESNTABLES=2.0
      TABLE1 = 201.0
      TABLE2 = 301.0
C     5 entries in 201 table
      SESNWDS1=5.0
C     Number of entries in EXTENDED 301 table: 4 variables at each rho,T point: S, cs, cv, KPA
      SESNWDS2=2.+DCOUNT+TCOUNT+DCOUNT*TCOUNT*4.
C     HEADER SECTION
      WRITE(14, 8253) INT(SESMATID), NWDS
C     readf,unit,matid,date,date,version,ntables,format="(5(E16.8E2))"
C     readf,unit,t1,t2,nwds1,nwds,format="(4(E16.8E2))"
      WRITE(14, 8254) SESMATID, DATE, DATE, VERSION, SESNTABLES
      WRITE(14, 8255) TABLE1, TABLE2, SESNWDS1, SESNWDS2
C     201 SECTION
      WRITE(14, 8256) INT(TABLE1),INT(SESNWDS1)
C     readf,unit,fmn,fmw,rho0_ref,k0_ref,t0_ref,format="(5(E16.8E2))"
      WRITE(14, 8254) FMN, FMW, RHO0REF, SESK0REF, T0REF
      WRITE(14, 8256) INT(TABLE2),INT(SESNWDS2)
C     number of density points, number of temperature points in grid
      WRITE(14, 8251) DCOUNT, TCOUNT
      STYLE=2
C     density array g/cm3
      DO 70 K=1, DCOUNT
      WRITE(14,8250) DARR(K)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 70   CONTINUE
C     temperature array eV to K
      DO 71 J=1, TCOUNT
      WRITE(14,8250) TARR(J)*TCONV
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 71   CONTINUE
C     specific entropy in ergs/eV/g to MJ/K/kg
      DO 78 J=1, TCOUNT
      DO 79 K=1, DCOUNT
      WRITE(14,8250) SARR(K,J)*(1.D-10/TCONV)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 79   CONTINUE
 78   CONTINUE
C     sound speed in cm/s
      DO 80 J=1, TCOUNT
      DO 81 K=1, DCOUNT
      WRITE(14,8250) CSARR(K,J)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 81   CONTINUE
 80   CONTINUE
C     specific heat capacity in erg/eV/g to MJ/kg/K
      DO 82 J=1, TCOUNT
      DO 83 K=1, DCOUNT
      WRITE(14,8250) CVARR(K,J)*(1.D-10/TCONV)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 83   CONTINUE
 82   CONTINUE
C     KPA FLAG
      DO 84 J=1, TCOUNT
      DO 85 K=1, DCOUNT
      WRITE(14,8250) ZKARR(K,J)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
 85   CONTINUE
 84   CONTINUE
C
      CLOSE(14)

C     WRITE STANDARD (SHORT) SESAME FILE
C     WRITE SESAME TABLE TO FILE LIMITED TO P, E, HFE
      OPEN(14,FILE='NEW-SESAME-STD.TXT')
C     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
C     could input matid, date, version with the grid
      NWDS=9
      SESNTABLES=2.0
      TABLE1 = 201.0
      TABLE2 = 301.0
C     5 entries in 201 table
      SESNWDS1=5.0
C     Number of entries in STANDARD 301 table: 3 variables at each rho,T point: P, U, A
      SESNWDS2=2.+DCOUNT+TCOUNT+DCOUNT*TCOUNT*3.
C     HEADER SECTION
      WRITE(14, 8253) INT(SESMATID), NWDS
C     readf,unit,matid,date,date,version,ntables,format="(5(E16.8E2))"
C     readf,unit,t1,t2,nwds1,nwds,format="(4(E16.8E2))"
      WRITE(14, 8254) SESMATID, DATE, DATE, VERSION, SESNTABLES
      WRITE(14, 8255) TABLE1, TABLE2, SESNWDS1, SESNWDS2
C     201 SECTION
      WRITE(14, 8256) INT(TABLE1),INT(SESNWDS1)
C     readf,unit,fmn,fmw,rho0_ref,k0_ref,t0_ref,format="(5(E16.8E2))"
      WRITE(14, 8254) FMN, FMW, RHO0REF, SESK0REF, T0REF
      WRITE(14, 8256) INT(TABLE2),INT(SESNWDS2)
C     number of density points, number of temperature points in grid
      WRITE(14, 8251) DCOUNT, TCOUNT
      STYLE=2
C     density array g/cm3
      DO 100 K=1, DCOUNT
      WRITE(14,8250) DARR(K)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
100   CONTINUE
C     temperature array eV to K
      DO 101 J=1, TCOUNT
      WRITE(14,8250) TARR(J)*TCONV
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
101   CONTINUE
C     pressure array dynes/cm2 to GPa
      DO 102 J=1, TCOUNT
      DO 103 K=1, DCOUNT
C      IF (PARR(K,J) .LT. 1.D-20 .AND. PARR(K,J) .GT. 0.) 
C     &    PARR(K,J)=1.D-20
      WRITE(14,8250) PARR(K,J)*1.D-10
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
103   CONTINUE
102   CONTINUE
C     specific internal energy array in ergs/g to MJ/kg
      DO 104 J=1, TCOUNT
      DO 105 K=1, DCOUNT
      WRITE(14,8250) EARR(K,J)*1.D-10
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
105   CONTINUE
104   CONTINUE
C     Helmholtz free energy array in ergs/g to MJ/kg
      DO 106 J=1, TCOUNT
      DO 107 K=1, DCOUNT
      WRITE(14,8250) (EARR(K,J)-TARR(J)*SARR(K,J))*(1.D-10)
      STYLE=STYLE+1
      IF (mod(STYLE,5) .EQ. 0) WRITE(14,8200)
107   CONTINUE
106   CONTINUE
C
      CLOSE(14)
C
C      WRITE(*,9200)
C     =====================END STSM SECTION====================
C
CCCCCC========================================
1000  FORMAT(///'DENSITY = ',1PE15.5,'  TEMP(K) = ',1PE15.5//)
2000  FORMAT('  ANALYTIC COMPUTATIONAL RESULTS FROM ANEOS1'// 
     &'       PRESSURE = ',1PE15.5,'   AT RHO+DELTA RHO  = ',1PE15.5/
     &'       ENERGY   = ',1PE15.5,'   AT T+DELTA T      = ',1PE15.5/
     &'       ENTROPY  = ',1PE15.5,'   AT RHO+DELTA RHO  = ',1PE15.5/
     &'       HEAT CAP = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       DPDT     = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       DPDR     = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5)
3000  FORMAT(//'  TEST OF SECOND DERIVATIVES--THESE SHOULD BE EQUAL'// 
     &'       -DSDR*RHO**2 = ',1PE15.5/   
     &'        DPDT        = ',1PE15.5)
4000  FORMAT(//'  BOUND STATE AND EOS FACTORS'//
     &'       ZVIB FUNC = ',1PE15.5/
     &'       FACTOR Y  = ',1PE15.5/
     &'       UNBOUND W = ',1PE15.5/
     &'       RENORM Z  = ',1PE15.5/
     &'       PSI       = ',1PE15.5/
     &'       GAMMA     = ',1PE15.5/
     &'       THETA     = ',1PE15.5/)
5000  FORMAT(//'  INTERMEDIATE QUANTITIES FOR MOLECULAR COMPUTATION'//
     &'       ALPHA     = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       BETA      = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       DADT      = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       DADR      = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       DBDT      = ',1PE15.5,'   APPROX DERIVATIVE = ',1PE15.5/
     &'       H1        = ',1PE15.5/
     &'       H2        = ',1PE15.5)
6000  FORMAT(//' THERMODYNAMIC FUNCTIONS OF TEMPERATURE AT 1 BAR,'//
     & '        TEMP         DENSITY       PRESSURE       ENERGY',
     & '         ENTHALPY        ENTROPY        GIBBS       PHASE',
     & '     #ITER'/
     & '         K           kg/m**3         GPa            J/kg',
     & '           J/kg           J/kg-K         J/kg'/)
6100  FORMAT(//' THERMODYNAMIC FUNCTIONS OF TEMPERATURE AT 1 BAR,'//
     & '        TEMP         DENSITY       PRESSURE       ENERGY',
     & '         ENTROPY       PHASE',
     & '     #ITER'/
     & '         K           kg/m**3         GPa            J/kg',
     & '         J/kg-K '/)
7000  FORMAT(//' PRESSURE ITERATION FAILS TO CONVERGE AFTER ',
     & I5,' STEPS'//)
8000  FORMAT(7(1PE15.5),5X,I2,5X,I5)
8100  FORMAT(5(1PE15.5),5X,I2,5X,I5)
8200  FORMAT(3(1PE15.8))
8250  FORMAT(1(1PE16.8), $)
8251  FORMAT(2(1PE16.8), $)
8252  FORMAT(2(I16), $)
8253  FORMAT(' INDEX     MATID =',I7,'   NWDS = ',I8)
8254  FORMAT(5(1PE16.8E2))
8255  FORMAT(4(1PE16.8E2))
8256  FORMAT(' RECORD',5X,'TYPE =',I5,5X,'NWDS = ',I8)
8300  FORMAT(//' COLD PRESSURE DEPENDENCE ON DENSITY, CGS UNITS'//
     &'       DENSITY       PRESSURE       ENERGY         GIBBS',
     &'           DPDR           ETA'/)
8500  FORMAT(6(1PE15.5))
8700  FORMAT(/'ITERATION HISTORY:'/' STEP#    DENSITY      PRESSURE   ',
     &'   DPDR'/)
8900  FORMAT(I5,3(1PE15.5))
9000  FORMAT(//' INDIVIDUAL CONTRIBUTIONS TO THERMODYNAMIC FUNCTIONS:'//
     &'  QUANTITY         COLD          NUCLEAR      MOLECULAR',
     &'         MELT        ELECTRONIC'//
     &'   ENERGY   ',5(1PE15.5)/
     &'   PRESSURE ',5(1PE15.5)/
     &'   ENTROPY  ',15X,4(1PE15.5)/
     &'   HEAT CAP.',15X,4(1PE15.5))
 9001 FORMAT(//'EOS TABLE'//)
 9002 FORMAT(//'DENSITY'//)
 9003 FORMAT(//'TEMPERATURE'//)
 9004 FORMAT(//'PRESSURE'//)
 9005 FORMAT(//'ENERGY'//)
 9006 FORMAT(//'ENTROPY'//)
 9007 FORMAT(//'SOUND SPEED'//)
 9008 FORMAT(//'KPA FLAG'//)
 9009 FORMAT(//'END OF FILE'//)
CSTS --help
 9100 FORMAT('ANEOS INPUTS:')
 9101 FORMAT('   ANEOS.INPUT   = formatted material parameter file')
 9102 FORMAT('   tablegrid.txt = SESAME table parameters (optional)')
 9103 FORMAT('Optional command line argument:')
 9104 FORMAT('   --help')
 9105 FORMAT('   --no_table')
 9106 FORMAT('     Do not generate a SESAME table.')
 9107 FORMAT('tablegrid.txt format: single column of floats')
 9108 FORMAT('   SESMATID')
 9109 FORMAT('   DATE')
 9110 FORMAT('   VERSION')
 9111 FORMAT('   FORMULA TOTAL ATOMIC NUMBER')
 9112 FORMAT('   FORMULA WEIGHT')
 9113 FORMAT('   Reference density')
 9114 FORMAT('   Reference bulk modulus')
 9115 FORMAT('   Reference temperature')
 9116 FORMAT('   Number of density points in grid')
 9117 FORMAT('   Number of temperature poitns in grid')
 9118 FORMAT('   Density grid points')
 9119 FORMAT('   Temperature grid points')
 9120 FORMAT('Output files:')
 9121 FORMAT('   NEW-SESAME-STD.TXT  Standard 201 and 301 tables')
 9122 FORMAT('   NEW-SESAME-EXT.TXT  Extra variables 301-style table')
CSTS end --help
C
      END

