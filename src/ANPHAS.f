C
C
      SUBROUTINE ANPHAS (MAT,IT,IKPN)
C
C***********************************************************************
C
C     ANEOS PACKAGE     SET UP FOR LIQUID-VAPOR AND SOLID-VAPOR
C                       CALCULATION - DETERMINES CRITICAL POINT
C                       AND TWO-PHASE BOUNDARIES
C
C           INPUTS: MAT  = Material number
C                   IT   = Pointer to beginning of storage in ACK array
C                   IKPN = Pointer to beginning of storage in RSOL, etc arrays
C
C           OUTPUTS: TCT = Critical Temperature
C                    RCT = Critical Density
C                    RSOL = Solid/liquid density on phase boundary (table)
C                    RVAP = Vapor density on phase boundary (table)
C                    TTWO = Temperature on phase boundary (table)
C
C************************************** 8/87 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /CNES/ P1,E1,S1,G1,P2,E2,S2,G2,PSI1,PSI2
      COMMON /FNES/ CMLT7,ACK46
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     kexprt is flag for extended print
C     = 0, no print      =1, print
      COMMON /KEXARP/ KEXPRT
      COMMON /ANEDIS/ GAMMA,PSI,THETA
C     Terms needed for table generation (Added by GSC, 17/07/07)
      COMMON /TABLES/ RCRIT(MAXMAT),TCRIT(MAXMAT),RNORM(MAXMAT),
     1     TNORM(MAXMAT),RZERO(MAXMAT)
C
      SAVE   /FILEOS/,/ANES/,/ANE2/,/CNES/,/FNES/,/KEXARP/,/ANEDIS/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (S3S3=0.0001D0)   !minimum step for derivatives
      PARAMETER (QCM1=0.3D0)      !critical density/zero pressure density
      PARAMETER (QCM2=1.D-6)      !convergence precision
      PARAMETER (QCM3=0.1D0)
      PARAMETER (QCM4=1.D-3)
      PARAMETER (QCM5=0.95D0)
      PARAMETER (QCM6=0.005D0)
      PARAMETER (QCM7=1.D-4)
      PARAMETER (QCM8=0.15D0)
      PARAMETER (QCM9=0.25D0)
      PARAMETER (QCM10=0.99D0)
      PARAMETER (QCM11=0.015D0)   !temperature minimum for Maxwell construction (170K)
      PARAMETER (QCM12=1.D-100)
      PARAMETER (QCM13=1.D6)      !high vapor pressure limit, dyne/cm**2
      PARAMETER (KWAYP=2)         !number of first guesses available
      PARAMETER (ONEBAR=1.D6)     !one bar pressure in dyne/cm**2
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
C      skip this subroutine for EOS without detailed treatment
C      of liquid/vapor region. 
C
      IF (NINT(ACK(IT+30)).LE.2 .OR. NINT(ACK(IT+30)).EQ.5 ) RETURN 
C
C     if you want to see entire 2-phase, set next parameter to 1
C
      ISKIP=5                      !print out only every ISKIPth line
      IF(KEXPRT.NE.ZERO) ISKIP=1   !all lines for experts
C
C     FIND CRITICAL POINT
C
C     Method is to find density and temperature where both
C          dP/dRHO = d2P/d2RHO = 0
C
C    Uses Newtonian iteration on a mesh of 9 points from which
C    second and third derivatives are evaluated.
C    See Equations 5.47 and 5.48.
C
C       RCT = critical density
C       TCT = critical temperature
C
      KWAY=0
    5 NTY=0
      KLY=0
      CMLT7=ZERO   !turn on special liquid terms (Eq. 5.22)
   10 KLY=KLY+1
      IF(KLY.EQ.1) THEN          !first guess:
        RCT(MAT)=QCM1*ACK(IT+19) !critical density = (zero temp density)/3
        TCT(MAT)=TWO             !critical temperature = 2 eV
      ELSE                       !alternative guess:
        RCT(MAT)=ACK(IT+19)      !critical density = zero temp density
        TCT(MAT)=.05D0           !critical temperature = 0.05 eV
      END IF
      S3=0.001D0                 !fractional step size for mesh
   50 R1=S3*RCT(MAT)             !density step
      T1=S3*TCT(MAT)             !temperature step
      KK=-2
        DO 70 I=1,9  !seek solution in a square mesh of 9 points
        IF (3*((I-1)/3).NE.I-1) GO TO 60
        KK=KK+1
        KN=-2
   60   KN=KN+1
        T2=KK
        T2=TCT(MAT)*(ONE+HALF*S3*T2)
        R2=KN
        R2=RCT(MAT)*(ONE+HALF*S3*R2)
        SQTS(1)=SQRT(T2)
        IF(KWAY.NE.1) THEN   !Compute dP/dRho and P on grid points
          CALL ANEOS1 (T2,R2,P1,E1,S1,D1,D2,RSOL(IKPN+I),IT)
        ELSE                 !Alternative computation of dP/dRho
          CALL ANNDPR (T2,R2,P1,E1,S1,D1,D2,RSOL(IKPN+I),IT)
        END IF
c
c       RSOL(IKPN+I) now contains dP/dRho at point I, I=1,9
c
        IF (I.NE.5) GO TO 70 !save pressure,energy, entropy
          RSOL(IKPN+10)=P1   !at central point (5) of mesh
          RSOL(IKPN+11)=E1   !in memory locations above 9
          RSOL(IKPN+12)=S1
   70   CONTINUE
c
c    Evaluate density, temperature increments to seek critical point
c    See Equations 5.47 and 5.48
c
      D1=RSOL(IKPN+5)
      D2=(RSOL(IKPN+6)-RSOL(IKPN+4))/R1
      D3=(RSOL(IKPN+8)-RSOL(IKPN+2))/T1
      D4=FOUR*(RSOL(IKPN+6)-TWO*RSOL(IKPN+5)+RSOL(IKPN+4))/(R1**2)
      D5=(RSOL(IKPN+9)-RSOL(IKPN+7)-RSOL(IKPN+3)+RSOL(IKPN+1))/(R1*T1)
      DR2=D3*D4-D2*D5
      DR1=(D2*D2-D1*D4)/DR2
      DR2=(D1*D5-D2*D3)/DR2
c
      IF(KWAY.GE.KWAYP) THEN   !neither initial guess worked
        WRITE(KLST,103) NTY,   !write out current guess
     &  TCT(MAT),DR1,RCT(MAT),DR2,(RSOL(IKPN+I),I=1,12)
      END IF
  103 FORMAT(' CP',I5,1P8E13.6,/,8X,8E13.6)
c
      IF (ABS(DR1).GT.QCM2*TCT(MAT)) GO TO 80
      IF (ABS(DR2).LE.QCM2*RCT(MAT)) GO TO 130 !converged OK!
   80 IF (ABS(DR1).LE.QCM3*TCT(MAT)) GO TO 90
      DR1=SIGN(QCM3*TCT(MAT),DR1)
   90 IF (ABS(DR2).LE.QCM3*RCT(MAT)) GO TO 100
      DR2=SIGN(QCM3*RCT(MAT),DR2)
  100 CONTINUE
      RCT(MAT)=RCT(MAT)+DR2  !new guess at critical density
      TCT(MAT)=TCT(MAT)+DR1  !new guess at critical temperature
      IF (S3.EQ.S3S3) GO TO 110
      IF (ABS(DR1).GT.QCM4*TCT(MAT)) GO TO 110
      IF (ABS(DR2).GT.QCM4*RCT(MAT)) GO TO 110
      S3=S3S3    !use 10X finer mesh for derivatives as we zero in
  110 NTY=NTY+1  !increment loop counter: if lt 100, continue search
      IF (NTY-100) 50,10,120  !if eq 100, try second guess at cp
  120 IF (NTY-200) 50,140,140 !if ge 200, use last resort
  130 IF (RSOL(IKPN+10).GT.ZERO) THEN  !pressure is positive
        D4=1.001*TCT(MAT)
        SQTS(1)=SQRT(D4)
          DO 133 IKJ=1,9
          R1=0.05*REAL(IKJ-5)
          R1=RCT(MAT)*(ONE-R1)
          CALL ANEOS1 (D4,R1,P1,E1,S1,D1,D2,D3,IT)
          IF(D3.LT.ZERO) GO TO 140 !negative dP/dRHO: go to last resort
  133     CONTINUE
        GO TO 260
      END IF
C
C     LAST RESORT METHOD TO FIND CRITICAL POINT
C
C       Try to follow curve for dP/dRHO = 0, locate
C       first maximum temperature on this curve.  Search
C       is from highest possible density downward.
C       This may fail if multiple maxima exist!
C
  140 CONTINUE
      R1=ACK(IT+19)   !Zero temperature, pressure density RHO00
      D6=R1           !look for a solution below this maximum density
      D5=ZERO         !initialize variable that records previous step temperature
      NTY=200
  150 KLY=0
      T1=-9.999D0     !this starts the hunt at 0.00005 eV
      T2=10.D0        !don't look at temperatures above 10 eV
  160 D4=HALF*(T1+T2) !we will seek right temperature by bisection
      SQTS(1)=SQRT(D4)
C
C   Get dP/dRHO by previously successful method
C
      IF(KWAY.NE.1) THEN
        CALL ANEOS1 (D4,R1,P1,E1,S1,D1,D2,D3,IT)
      ELSE
        CALL ANNDPR (D4,R1,P1,E1,S1,D1,D2,D3,IT)
      END IF
C
      KLY=KLY+1                        !count iterations
      IF (KLY.GT.9900) GO TO 240       !no convergence.  Quit trying
      IF (T2-T1.LE.QCM2*D4) GO TO 210  !Temperature converged
      IF (D3) 170,210,180  !Test dP/dRHO
  170 T1=D4                !dP/dRHO negative: Raise minimum temperature
      GO TO 160            !and try again
  180 IF (T1) 220,220,190  !dP/dRHO positive: 
  190 T2=D4                !lower maximum temperature
      GO TO 160            !and try again
  210 CONTINUE             !come here when temperature converged
C
      IF(KWAY.GE.KWAYP) THEN  !out of guesses--write out current results
        WRITE(KLST,213) NTY,D4,R1,P1,E1,S1,D1,D2,D3
      END IF
  213 FORMAT(' CPB',I5,1P8E13.6)
C
      IF(D5.LT.D4) THEN !temperature at this step bigger than last
        D5=D4           !save temperature at this step
        D6=R1           !save density at this step
      END IF
      IF(D4.LT.0.8*D5) GO TO 230  !temperatue maximum is passed
C                                 !accept this as cp and leave
  220 NTY=NTY+1
      R2=MAX(QCM6*R1,QCM7)  !small density increment
      R1=R1-R2              !decrease upper density slightly
      IF (R1) 240,240,150   !if this is gt 0, go back to T search
  230 TCT(MAT)=D5           !record critical temperature
      RCT(MAT)=D6           !and critical density
      SQTS(1)=SQRT(D5)
C
C   Compute critical point Pressure, Energy and Entropy
C
      IF(KWAY.NE.1) THEN
        CALL ANEOS1 (TCT(MAT),RCT(MAT),RSOL(IKPN+10),RSOL(IKPN+11),
     &  RSOL(IKPN+12),D1,D2,D3,IT)
      ELSE
        CALL ANNDPR (TCT(MAT),RCT(MAT),RSOL(IKPN+10),RSOL(IKPN+11),
     7  RSOL(IKPN+12),D1,D2,D3,IT)
      END IF
C
      GO TO 260      !continue with two phase setup
  240 CONTINUE       !come here if critical pressure, density le zero
      IF(KWAY.LE.1) THEN
        KWAY=KWAY+1  !try a different initial guess, if you still have one
        GO TO 5
      END IF
      CALL ANWARN(0)
      WRITE(KLST,440) MAT         !iteration will not converge--warn user
  250 ACK(IT+30)=ACK(IT+30)-THREE !switch EOS type to one without liq/solid
      KEXPRT=1                    !turn on full debug outputs
      RETURN         !failed to find critical point. Give it up
c
c    Apparently successful convergence gets you here
c
  260 KN=IKPN+10 !pointer to critical pressure
      KK=KN+2    !next 2 locations are energy and entropy
c
c     Make sure critical pressure is positive.  If not, and alternative initial
c     guess has not yet been used, go back and try again
c
      IF((RSOL(KN).LE.ZERO).AND.(KWAY.LE.1)) GO TO 240
C
C     Write out successful critical point parameters to output file
C
      WRITE(KLST,450) MAT,RCT(MAT),TCT(MAT),(RSOL(I),I=KN,KK),NTY
C
C   extract other critical point data for later printing
C
      PCT=RSOL(KN)     !Critical Pressure
      ECT=RSOL(KN+1)   !Critical Energy
      SCT=RSOL(KN+2)   !Critical Entropy
      PSICT=PSI        !Interpolation ratio at CP
C
C     Store critical density and temperature, and reference density, 
C     for table generation (Added by GSC, 17/07/07)
C
      RCRIT(MAT) = RCT(MAT)*DCONV    ! Density at critical point
      TCRIT(MAT) = TCT(MAT)*TCONV    ! Temperature at critical point
      RNORM(MAT) = ACK(IT+11)*DCONV  ! Reference density
      TNORM(MAT) = ACK(IT+12)*TCONV  ! Reference temperature
      RZERO(MAT) = ACK(IT+19)*DCONV  ! Density at zero temperature and pressure
C
      ACK(IT+73)=TCT(MAT)  !store critical temperature in ACK
      ACK(IT+74)=RCT(MAT)  !store critical density in ACK
      IF((ZBAR.GT.HALF).AND.(ACK(IT+63).GE.ZERO)) THEN
        CALL ANWARN(0)
        WRITE(KLST,262)
  262   FORMAT(/,' Critical point is too high'
     1  ,' and EOS might not be useable.',/)
      END IF
      IF (RSOL(KN).LE.ZERO) GO TO 240
C
C    Critical density is negative.  Try last another method
C
      IF (TCT(MAT).GT.ACK(IT+18)) GO TO 270
C
C    note ACK(IT+18) is melt temperature of material pointer IT.
C    If you get here, critical temperature is larger than melt
C    temperature.  Issue warning of trouble, but go on
C
      CALL ANWARN(0)
      WRITE(KLST,460) ACK(IT+18)
      GO TO 250
C
C     FIND LIQUID-VAPOR PHASE BOUNDARIES
C
C      This segment constructs arrays RSOL(100), RVAP(100) of solid
C      and vapor densities on the phase curve, as a function of temperature
C      TTWO(100).  The first entry is the critical point, the last one is
C      extrapolated to zero temperature.
C
  270 KK=60   !number of steps between critical temp to melt temp
      KN=20   !number of steps between melt temp to zero
      IF (ACK(IT+18).GT.QCM8) KN=30
      IF (ACK(IT+18).GT.QCM9) KN=40  !use more steps for higher Tmelt
C
C    Storage note:  KK + KN must be le 100, since only 100 locations
C    are allocated in arrays RVAP, RSOL and TTWO.  If more is needed,
C    their dimensions must be increased.
C
      KLY=0
      RVAP(IKPN)=RCT(MAT)   !critical point--beginning of tables.
      RSOL(IKPN)=RCT(MAT)   !previous contents of RVAP and
      TTWO(IKPN)=TCT(MAT)   !RSOL are overwritten here.
      WRITE(KLST,470)       !header for output listing
C
C   start output of phase boundary with critical point data echo
C
      GCT=ECT+PCT/RCT(MAT)-TCT(MAT)*SCT  !Gibbs free energy at cp
      WRITE(KLST,490) TCT(MAT)*TCONV,RCT(MAT)*DCONV,RCT(MAT)*DCONV,
     & PCT*PCONV,PCT*PCONV,ECT*ECONV,ECT*ECONV,SCT*SCONV,SCT*SCONV,
     & GCT*ECONV,GCT*ECONV,PSICT,PSICT,0
C
      IK=IKPN+1           !next memory location in RVAP, etc. arrays
      D5=KK
      D5=(TCT(MAT)-ACK(IT+18))/D5  !step size from critical temp down to Tmelt
      D6=KN
      D6=ACK(IT+18)/D6             !step size from Tmelt to zero
      S3=QCM5
C
C    The phase curve tables are constructed in two parts, first from the
C    critical temperature down to the melt (or vaporization) temperature,
C    then from the melt temperature to zero.
C
      DO 392 JJJ=1,2
      D4=D5
      JJJJ=KK+10
      T=TCT(MAT)
      D1=ZERO
      IF (ACK46.GT.ZERO) D1=ACK(IT+47) !ACK(47) = liquid density at triple point
      IF (JJJ.EQ.1) GO TO 280  !compute phase curve from critical temp to Tmelt
      D4=D6                    !otherwise, compute phase curve from Tmelt to zero
      JJJJ=KN
      T=ACK(IT+18)
      IF (ACK46.LE.ZERO) GO TO 280
      JJJJ=KN+1
      T=QCM10*T+D4
      D1=ZERO
C
C   Determine temperature for next entry in tables
C
  280   DO 394 I=1,JJJJ
        IF (I.EQ.KK-9) D4=HALF*D4
        T=T-D4              !nominal temperature step
        IF (JJJ.EQ.1) THEN  !if this drops below Tmelt, then use Tmelt
          IF(T-HALF*D4.LT.ACK(IT+18)) T=ACK(IT+18)
        ELSE
          IF(T.GE.TTWO(IK-1)) GO TO 390
        END IF
        IF (T.LT.QCM11) GO TO 400  !skip Maxwell if temperature is below QCM11
        NTY=0
        IF (T.GT.S3*TCT(MAT)) NTY=-1
        R2=ZERO    !initial guess at vapor density
        R1=D1      !initial guess at liquid/solid density
        IF (ACK46.GT.ZERO.AND.T.EQ.ACK(IT+18)) R1=-R1
        IF (RVAP(IK-1).LE.QCM12) R2=-ONE
        IF (IK.GT.IKPN+1) GO TO 290
        IF (KLY.GE.12) GO TO 290
        NTY=-1
C
C       Call to Maxwell construction at given temperature to find
C       density of liquid and vapor on the phase curve
C
  290   CALL ANMAXW (T,R1,R2,IT,MAT,NTY)
C
C       If this succeeds OK, enter R1 R2 and T in tables at line 340.
C       If not, try to figure out what went wrong and inform the user.
C
        IF (NTY) 300,340,340
  300   KLY=KLY+1
        IF(KEXPRT.NE.0) THEN     !for expert users, write out error code
          WRITE(KLST,346) NTY,T  !and temperature where it occurred
        END IF
  346   FORMAT(' anphas liq-vap construction failure',I5,' T=',1PE13.6)
        IF (T.GT.S3*TCT(MAT)) KLY=0
        IF (IK.EQ.IKPN+1) GO TO 330
        IF (KLY-2) 310,320,320
  310   IF (T.EQ.ACK(IT+18)) GO TO 320
        IF (T.GT.S3*TCT(MAT)) GO TO 390
        CALL ANWARN(0)               !warning--no convergence at this point
        WRITE(KLST,410) T,R1,R2,NTY  !print out problem point
        GO TO 390                    !leave this point out and go on
  320   CONTINUE
        CALL ANWARN(1)     !fatal error
        WRITE(KLST,420)    !print error message
        GO TO 250          !change form of EOS to one without liq/solid phase
  330   IF (KLY-13) 390,310,320
C
C       Insert phase densities, temperature into tables
C
  340   RSOL(IK)=R1
        RVAP(IK)=R2
        TTWO(IK)=T
C
C       Check for other possible problems with this EOS
C
        IK=MIN(IK+1,100*MAXMAT)
        KLY=0
        IF (T-ACK(IT+18)) 370,360,350
  350   IF (R1.LE.ACK(IT+47)) GO TO 370   !for T above Tmelt,
        CALL ANWARN(0)                    !R1 gt density at triple point
        WRITE(KLST,430) T,R1              !print warning and location of trouble
        GO TO 370
  360   IF (R1.LT.ACK(IT+23)) THEN        !for T = Tmelt,
           CALL ANWARN(0)                 !solid density less than minimum preset
           WRITE(KLST,480) ACK(IT+23),R1  !warn the user to lower RHOmin.
        END IF
        IF (P2.GT.QCM13) THEN             !vapor pressure is more than preset limit,
           CALL ANWARN(0)                 !which is 1 bar at Tmelt
           WRITE(KLST,500)                !warn the user, but go on
        END IF
        GO TO 380
  370   IF (JJJ.EQ.2.AND.I.EQ.1.AND.ACK46.GT.ZERO) GO TO 380
        IF (ISKIP*(I/ISKIP).NE.I) GO TO 390
  380   CONTINUE
C
C    Check to see if pressure drops below 1 bar.  If so, interpolate to 1 bar
C    and write out interpolated data to approximate 1 bar thermodyanmic
C    parameters, for comparison to tabulated data.
C
        IF(.NOT.((JJJ.EQ.1).AND.(I.EQ.1))) THEN !don't do this for first entry
          IF((P1OLD.GT.ONEBAR).AND.(P1.LE.ONEBAR)) THEN !new pressure below 1 bar
            FRAC=(P1OLD-ONEBAR)/(P1OLD-P1)
            TBAR=TOLD-FRAC*(TOLD-T)
            R1BAR=R1OLD-FRAC*(R1OLD-R1)
            R2BAR=R2OLD-FRAC*(R2OLD-R2)
            P1BAR=P1OLD-FRAC*(P1OLD-P1)
            P2BAR=P2OLD-FRAC*(P2OLD-P2)
            E1BAR=E1OLD-FRAC*(E1OLD-E1)
            E2BAR=E2OLD-FRAC*(E2OLD-E2)
            S1BAR=S1OLD-FRAC*(S1OLD-S1)
            S2BAR=S2OLD-FRAC*(S2OLD-S2)
            G1BAR=G1OLD-FRAC*(G1OLD-G1)
            G2BAR=G2OLD-FRAC*(G2OLD-G2)
            PSI1BAR=PSI1OLD-FRAC*(PSI1OLD-PSI1)
            PSI2BAR=PSI2OLD-FRAC*(PSI2OLD-PSI2)
            WRITE(KLST,490) TBAR*TCONV,R1BAR*DCONV,R2BAR*DCONV,
     &      P1BAR*PCONV,P2BAR*PCONV,E1BAR*ECONV,E2BAR*ECONV,
     &      S1BAR*SCONV,S2BAR*SCONV,G1BAR*ECONV,G2BAR*ECONV,
     &      PSI1BAR,PSI2BAR,0
          ENDIF
        ENDIF           
C
C     Echo phase data to output file.  Note that data comes from last
C     call to ANMAXW, via labeled common /CNES/
C
        WRITE(KLST,490) T*TCONV,R1*DCONV,R2*DCONV,P1*PCONV,P2*PCONV,
     & E1*ECONV,E2*ECONV,S1*SCONV,S2*SCONV,G1*ECONV,G2*ECONV,PSI1,PSI2,
     & NTY
C
C     Save data just written to perform interpolation to 1 bar state
C
        TOLD=T
        R1OLD=R1
        R2OLD=R2
        P1OLD=P1
        P2OLD=P2
        E1OLD=E1
        E2OLD=E2
        S1OLD=S1
        S2OLD=S2
        G1OLD=G1
        G2OLD=G2
        PSI1OLD=PSI1
        PSI2OLD=PSI2
C
  390   IF(T.EQ.ACK(IT+18)) GO TO 392 !this prevents Tmelt from being entered twice
  394   CONTINUE
  392 CONTINUE
C
C   Come here if temperature drops below preset minimum QCM11.  This
C   is the last entry in the table
C
  400 RSOL(IK)=ACK(IT+19)   !Zero-temperature, pressure density
      RVAP(IK)=ZERO
      TTWO(IK)=ZERO
      IKPN=IK+1
C
  410 FORMAT ('  Iteration failure T=',1PE12.4,' Will leave point out.'
     1,/,'   R =',2E12.4,I5)
  420 FORMAT (/,'  FATAL ERROR - Will change form of eos.')
  430 FORMAT (
     1 '  Negative expansion coefficient in the liquid phase'
     2 ,/,2X,2HT=,1PD12.5,2X,4HRHO=,D12.5,5X,
     2 'Improper behavior will result')
  440 FORMAT (68H0 THE CRITICAL POINT ITERATION WILL NOT CONVERGE FOR MA
     1TERIAL NUMBER,I5,26H. WILL CHANGE FORM OF EOS.)
  450 FORMAT (//'  TWO-PHASE CALCULATION FOR MATERIAL ',I5,//,
     1 16H  CRITICAL POINT,/,6H  RHO=,1PD15.7,7X,2HT=,D15.7,9X,
     1 2HP=,D15.7,/,2X,2HE=,D15.7,9X,2HS=,D15.7,9X,4HNTY=,I5,/)
  460 FORMAT (26H0 THE MELTING TEMPERATURE(,1PD15.7,64H) IS GREATER THAN
     1 CRITICAL TEMPERATURE. WILL CHANGE FORM OF EOS.)
  470 FORMAT (/,'  TWO-PHASE BOUNDARIES',/
     1'       T         RHOLIQ        RHOVAP        PLIQ         PVAP',
     2'        ELIQ         EVAP         SLIQ         SVAP',
     3'        GLIQ         GVAP         PSILIQ      PSIVAP         ',
     4'NTY'/
     5'       K         kg/m**3       kg/m**3       GPa          GPa ',
     6'        J/kg         J/kg        J/kg-K       J/kg-K',
     7'       J/kg         J/kg'/)
  480 FORMAT (/,' The minimum solid density(',1PD12.5,
     1') is greater than the triple point density(',D12.5,').',/,
     2'  Improper solid behavior will result.',
     3' To correct use smaller value.',/)
  490 FORMAT (1P13D13.5,I7)
  500 FORMAT (//,' High vapor pressure',/)
      RETURN
      END
