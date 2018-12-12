C
C
      SUBROUTINE ANMAXW (T,RL,RV,L,MAT,IERR)
C
C***********************************************************************
C
C     ANEOS PACKAGE  LIQUID-VAPOR AND SOLID-VAPOR MAXWELL CONSTRUCTION
C
C************************************** 8/87 version slt ***************
C
C       INPUTS:     T = Temperature in eV
C                   L = Pointer to storage for material MAT in ACK array
C                 MAT = Material number
C                  RL = first guess at liquid/solid density on entry
C                       If RL = 0, code uses RHO00 for initial guess
C                       If RL < 0, liquid density is fixed at -RL,
C                          forcing code to adjust vapor density only
C                  RV = first guess at vapor density on entry
C                  If IERR lt 0 upon entry, convergence error message is
C                     suppressed.
C
C       OUTPUTS:   RL = Liquid/solid density in gm/cc
C                  RV = Vapor density in gm/cc
C                IERR = Error indicator:
C                     =  N >= 0, successful completion after N loops
C
C                     = -1, convergence failure for vapor density
C                     = -2, convergence failure for liquid density
C                      Either of these errors indicates that the code
C                      cannot find a density where dP/dRHO > 0!
C                      The EOS must be really pathological.
C
C                     = -3, Gibbs energy, pressure not equal after 500
C                          iterations
C                     = -4, Special low density treatment fails
C
C                The pressure, energy, entropy and gibbs energy
C                of each phase are available from the common /CNES/
C                upon exit from this subroutine
C    
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
C
C   This common transmits thermodynamic data for each phase
C
      COMMON /CNES/ P1,E1,S1,G1,P2,E2,S2,G2,PSI1,PSI2
      COMMON /ANEDIS/ GAMMA,PSI,THETA
C
C   CMLT7 is flag for addition of liquid terms in ANEOS1
C   ACK46 is density of solid at triple point
C
      COMMON /FNES/ CMLT7,ACK46
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
      SAVE   /FILEOS/,/ANES/,/CNES/,/FNES/,/ANEDIS/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (QCN1=1.5D0)
      PARAMETER (QCN2=1.D-3)
      PARAMETER (QCN3=1.D-100)
      PARAMETER (QCN4=0.99D0)    !step for vapor density
      PARAMETER (QCN5=1.005D0)   !step for liquid density
      PARAMETER (QCN6=1.D-6)     !convergence precision
      PARAMETER (QCN7=0.05D0)
      PARAMETER (QCN8=0.01D0)    !fractional deviation of Gibbs energy
      PARAMETER (QCN9=1.D4)
      PARAMETER (QCN10=0.1D0)
      PARAMETER (QCN11=1.D-9)
C
      SQTS(1)=SQRT(T)
      NP=0
      N=0
      CM7=CMLT7   !save value of CMLT7 flag on entry
      IF (IERR.LT.0) NP=1
      IF (RV.LT.ZERO) GO TO 230
      RLO=RCT(MAT)
      RVO=RCT(MAT)
C
C     use limiting EOS to estimate density of vapor phase
C     at temperature T as a first guess.  See Eq. 4.18 to 4.20
C
      RV=THREE*ACK(L+15)-ONE-ACK(L+10)/(ACK(L+27)*BOLTS*T)
      RV=ACK(L+25)**3*EXP(RV)/(ACK(L+13)*T)**QCN1
      IF(RV.LT.ZERO) RV=-RV  !this should never happen
      IF (RL) 10,20,30
   10 RL=-RL
      AL=ZERO
      GO TO 40
   20 RL=ACK(L+19)  !if input RL is zero, use RHO00
   30 AL=ONE
   40 RLM=RL
      DP2=QCN2*RL
      RV=MAX(QCN3,MIN(DP2,RV)) !RV not larger than 0.001 RL
C
C    Make repeated calls to ANEOS1 to find RV where dP/dRHO > 0
C    In one-phase regions dP/dRHO for fixed T must be positive,
C    but not necessarily in two phase region.
C
   50 IERR=0
   60 CALL ANEOS1 (T,RV,P2,E2,S2,D1,D2,DP2,L)
C
      IF (DP2.GT.ZERO) GO TO 80  !dP/dRHO is positive here
      RVO=RV         !save last estimate in RVO
      RV=QCN4*RV     !decrease vapor density estimate by 1%
      IF (IERR.GT.30) RV=HALF*RV !drop vapor density by half
      IERR=IERR+1
      IF (IERR-900) 60,60,70  
   70 IERR=-1        !we have looped 900 times!  Give up here
      GO TO 220
   80 G2=E2-T*S2+P2/RV   !Gibbs energy of vapor phase
      PSI2=PSI           !Extract interpolation raio of vapor
C
C   now make a similar search of the liquid region, again
C   seeking a density where dP/dRHO is positive.
C
      IERR=0
   90 IF (T.LT.ACK(L+18)) CMLT7=-ONE  !we are below triple point
C                                     !so do not add liquid terms
      CALL ANEOS1 (T,RL,P1,E1,S1,D1,D2,DP1,L)
      CMLT7=CM7                       !reset melt flag
      IF (DP1.GT.ZERO) GO TO 110      !dP/dRHO is positive here
      RLO=RL
      RL=QCN5*RL                      !increase liquid density by 0.5%
      IERR=IERR+1
      IF (IERR-900) 90,90,100         !we have looped 900 times!  Give up here
  100 IERR=-2
      GO TO 220
  110 G1=E1-T*S1+P1/RL    !Gibbs energy of liquid phase
      PSI1=PSI            !Extract interpolation ratio of liquid
C
C     compute the differences between the pressure and Gibbs energy
C     of the liquid and vapor phases.  Equations 5.50 and 5.51.
C     Use these to compute the density changes needed to bring them
C     into closer agreement, Equations 5.52 and 5.53
C
      SP=P1-P2
      SG=G1-G2
      DRL=AL*RL*(SP-RV*SG)/(DP1*(RV-RL))
      DRV=RV*(SP-RL*SG)/(DP2*(RV-RL))
C
      IF (ABS(DRL).GT.QCN6*RL) GO TO 120   !big step
      IF (ABS(DRV).LE.QCN6*RV) GO TO 200   !little step
  120 IF (N.GT.40) DRL=HALF*DRL
      IF (N.LT.60) GO TO 130
      DRL=QCN7*DRL
      IF (ABS(SP).GT.QCN8*(P1+P2+QCN9)) GO TO 130
      IF (ABS(SG).LE.QCN8*(ABS(G1)+ABS(G2))) GO TO 200
  130 SP=RL+DRL  !now SP is updated liquid density
      IF (SP.GT.RLO) GO TO 150
  140 DRL=HALF*DRL
      GO TO 130
  150 IF (DRL.GT.QCN10*RL) GO TO 140
  160 SG=RV+DRV  !now SG is updated vapor density
      IF (SG.GT.ZERO) GO TO 180
  170 DRV=HALF*DRV
      GO TO 160
  180 IF (SG.GE.RVO) GO TO 170
      RV=SG                !new guess for vapor density
      RL=SP                !new guess for liquid density
      N=N+1
      IF (N-500) 50,50,190 !if less than 500 loops, go back
  190 IERR=-3              !Gibbs energy, pressure do not agree
      GO TO 210
  200 IERR=N
      IF (RL.LE.RLM) RETURN !done--successful exit
      CALL ANWARN(0)        !liquid density is > RHO00!
      WRITE(KLST,310) T,RL  !possible negative expansion coefficient
      RETURN
  210 IF (RV.LT.QCN3) GO TO 230
  220 IF (NP.EQ.0) THEN
         CALL ANWARN(0)                !cannot get convergence
         WRITE(KLST,320) T,IERR,RV,RL  !write out location
      END IF
      RETURN
C
C     VAPOR DENSITY TOO SMALL TO CALCULATE
C     LIQUID-SOLID POINT AT P=0
C
  230 RLM=ACK(L+19)  !liquid density is RHO00
      RL=RLM
      N=0
      IF (T.LT.ACK(L+18)) CMLT7=-ONE  !below triple point
      P2=HALF*RL
      E2=RL
C
C     get thermodynamic properties of liquid
C
  240 CALL ANEOS1 (T,RL,P1,E1,S1,D1,D2,DP1,L)
      PSI1=PSI   !Extract interpolation ratio of liquid
C
      IF (P1.LT.ZERO.AND.N.LT.800) GO TO 250
      IF (ABS(P1).LE.QCN2) GO TO 300  !pressure less than 1 nannobar 
      IF (E2-P2.LE.QCN11)  GO TO 300
  250 IF (P1) 260,300,270   !try to fix negative pressure
  260 P2=RL
      GO TO 280
  270 E2=RL
  280 RL=HALF*(E2+P2)
      N=N+1
      IF (N-900) 240,240,290 !loop back again
  290 IERR=-4    !does not converge after 900 loops
      GO TO 220  !give up
C
  300 RV=QCN3    !successful completion comes here
C                !vapor density is given a very small value (1.D-100)
      CMLT7=CM7  !compute thermodynamic properties of vapor
      CALL ANEOS1 (T,RV,P2,E2,S2,D1,D2,DP1,L)
      PSI2=PSI   !Extract interpolation ratio of vapor
      GO TO 200  !we're done
C
  310 FORMAT(/,' Possible negative expansion coefficient  T='
     1,1PD12.5,5H RHO=,D12.5)
  320 FORMAT (42H  ANMAXW TWO-PHASE CONVERGENCE ERROR AT T=,1PD12.5,6H I
     1ERR=,I5,2D15.7)
      RETURN
      END
