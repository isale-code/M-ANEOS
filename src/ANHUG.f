C
C
      SUBROUTINE ANHUG (M,RO,TO,RHUGMX)
C
C***********************************************************************
C
C     ANEOS PACKAGE     SIMPLE INLINE HUGONIOT CALCULATION
C
C************************************** 5/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE   /FILEOS/
      CHARACTER*7 WDF
      DIMENSION TS(48)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (QCO1=1.D-8)
      PARAMETER (TVLOW=1.D-6)
      DATA (TS(I),I=1,48)/.026D0,.0265D0,.0275D0,.0285D0,.03D0,.035D0,
     1 .04D0,.05D0,.06D0,.08D0,.1D0,.12D0,.14D0,.16D0,.18D0,.2D0,.25D0,
     2 .3D0,.35D0,.4D0,.45D0,.5D0,.55D0,.6D0,.65D0,.7D0,.75D0,.8D0,
     3 .85D0,.9D0,.95D0,1.D0,1.1D0,1.2D0,1.3D0,1.4D0,1.5D0,1.7D0,2.D0,
     4 2.5D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0/
      IF ((RO.LE.ZERO).OR.(TO.LE.ZERO)) RETURN
C
C    define some conversion factors for ease of readers:
C      tconv = converts temperature in eV to K
C      sconv = converts entropy in erg/gm-eV to J/kg-K
C      econv = converts energy in erg/gm to J/kg
C      dconv = converts density in gm/cc to Kg/m**3
C      pconv = converts pressure in dyne/cm to GPa
C      vconv = convertes velocity in cm/sec to km/sec
C
      TCONV = 1.160564D4
      SCONV = 8.617D-9
      ECONV = 1.0D-4
      DCONV = 1.0D3
      PCONV = 1.0D-10
      VCONV = 1.0D-5
C
      WRITE(KLST,40)
      CALL ANEOSV(1,TVLOW,RO,M,PC,EC,SC,CVC,PTC,PRC,D1C,D2C,KP,
     &   EXTRAA,EXTRAB,EXTRAC)
      CALL ANEOSV(1,TO,RO,M,PO,EO,SO,D1,D2,D3,D4,VO,KP,
     &   EXTRAA,EXTRAB,EXTRAC)
      D1=ZERO
      D2=ONE
      WRITE(KLST,60)
      WRITE(KLST,50) RO*DCONV,TO*TCONV,PO*PCONV,PC*PCONV,EO*ECONV,
     & SO*SCONV,VO*VCONV,D1,D2,0,'INITIAL'
      N=51
C
      DO 30 I=1,48
      T=TS(I)
      IF (T.LE.TO) GO TO 30
      IF (N.GT.50) R=RO
      R=MIN(R,RHUGMX)
      N=0
   10 CALL ANEOSV(1,T,R,M,P,E,S,CV,PT,PR,D1,D2,KP,EXTRAA,EXTRAB,EXTRAC)
      F=E-EO+HALF*(PO+P)*(RO-R)/(R*RO)
      DF=(P-T*PT)/R**2+HALF*PR*(RO-R)/(RO*R)-HALF*(PO+P)/R**2
      IF (DF.EQ.ZERO) GO TO 30
      DR=-F/DF
      IF (ABS(DR).LE.QCO1*R) GO TO 20
      D1=ONE
      IF (DR.LT.ZERO) D1=-ONE
      IF (ABS(DR).GT.HALF*R) DR=HALF*R*D1
      IF((R.GE.RHUGMX).AND.(DR.GT.ZERO)) GO TO 31
      R=MIN(R+DR,RHUGMX)
      N=N+1
      IF (N-50) 10,10,30
   20 IF(R.GE.RHUGMX) GO TO 31
      V=(P-PO)/(RO*(ONE-RO/R))
      V=SQRT(V)
      U=V*(ONE-RO/R)
      D1=R/RO
      IF (KP.EQ.1) THEN
         WDF='1 PHASE'
      ELSEIF (KP.EQ.2) THEN
         WDF='L/S+VAP'
      ELSEIF (KP.EQ.4) THEN
         WDF='SOLID  '
      ELSEIF (KP.EQ.5) THEN
         WDF='SOL+LIQ'
      ELSEIF (KP.EQ.6) THEN
         WDF='LIQUID '
      ELSE
        WDF='       '
      END IF
      CALL ANEOSV(1,TVLOW,R,M,PC,EC,SC,CVC,PTC,PRC,D1C,D2C,KP,
     & EXTRAA,EXTRAB,EXTRAC)
      WRITE(KLST,50) R*DCONV,T*TCONV,P*PCONV,PC*PCONV,E*ECONV,S*SCONV,
     & V*VCONV,U*VCONV,D1,N,WDF
   30 CONTINUE
C
   40 FORMAT (//,'   HUGONIOT')
   50 FORMAT (1P9D12.4,I3,2X,A7)
   60 FORMAT (/,'      RHO          T           P          PC        ',
     1'   E           S           V           U       RHO/RHOO',
     2'  #IT  STATE'/
     3'    kg/m**3        K          GPa        GPa          J/kg',
     4'      J/kg-K       km/sec      km/sec'/)
   31 CONTINUE
      END
