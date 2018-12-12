C
C
      SUBROUTINE ANELSM (T,RHO,RL,RS,L,IERR)
C
C***********************************************************************
C
C     ANEOS PACKAGE  SOLID-LIQUID TWO-PHASE MAXWELL CONSTRUCTION (MELT)
C
C              INPUTS:   T    = TEMPERATURE IN EV
C                        RHO  = DENSITY IN GM/CC
C                        L    = MATERIAL NUMBER
C                        IERR = OPERATION CODE ON INPUT
C                             = 0 NORMAL VALUE ON RUNNING ENTRY
C                                --FATAL ERRORS CAUSE CODE TO PRINT ERROR
C                                  MSG AND THEN STOP
C                           IERR > 1 ARE FLAGS USED FOR DEBUGGING EOS:
C                             = 1 EXECUTION CONTINUES AFTER PRINTING ERROR MSG
C                             = 2 IERR SET TO NUMBER OF ITERATIONS ON EXIT
C                             = 3 DOES NOT HALT EXECUTION IF ZERO COMPRESSIBILITY
C                                 IS ENCOUNTERED (AND IERR IS SET TO -3)
C                        ACK(L+72) = OPERATION FLAG (CARRIED IN /ANES/)
C                             ACK(L+72) = ZERO, USE FAST ITERATION METHOD
C                             ACK(L+72) = ONE,  USE SLOW ITERATION METHOD
C                        ACK46 = DENSITY OF SOLID AT TRIPLE POINT (IN /FNES/)
C
C               OUTPUTS: RL  = DENSITY OF LIQUID AT T
C                        RS  = DENSITY OF SOLID AT T
C                        DATA ON MIXED PHASE IS TRANSMITTED IN COMMONS:
C                              /ANE2/  CRHS2P,CRLS2P,ZBAR,FKROS
C                              /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
C                        IERR = COMPLETION/ERROR CODE UPON RETURN:
C                             =  1 MATERIAL IS IN PURE SOLID REGION
C                             =  2 MATERIAL IS IN MIXED PHASE REGION
C                             =  3 MATERIAL IS IN PURE LIQUID REGION
C                        IERR < 0 ARE USUALLY FATAL ERROR CONDITIONS:
C                             = -1 DENOMINATOR IN EQ. 5.39 VANISHES
C                             = -2 TOO MANY ITERATIONS ( > 500)
C                             = -3 IERR = -3 IS NOT TREATED AS FATAL, RETURNS
C                                  TO CALLING PROGRAM, OTHERWISE SAME AS -4
C                             = -4 ZERO COMPRESSIBILITY OF SOLID OR LIQUID
C
C************************************** 8/87 version slt ***************
C************************************** 6/07 revised hjm ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FAST
      PARAMETER(MAXIT=500)     !maximum number of iterations
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10) !maximum number of materials
      PARAMETER (MAXHPP=3)  !maximum number of hpp phases per material
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /ANXTRA/ TMLT(1000*MAXMAT),RMLT(1000*MAXMAT)
     1 ,RSLD(1000*MAXMAT),THPP(MAXHPP,1000*MAXMAT)
     2 ,RLPP(MAXHPP,1000*MAXMAT),RHPP(MAXHPP,1000*MAXMAT)
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
      COMMON /FNES/ CMLT7,ACK46
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE   /FILEOS/,/ANES/,/ANE2/,/BNES/,/FNES/,/ANCNTR/
      CHARACTER*75 ERRMES
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (QCI1=0.05D0)  !large step size for iteration
      PARAMETER (QCI2=0.99D0)  !first guess at RL = QCI2*RHO
      PARAMETER (QCI3=1.01D0)  !first guess at RS = QCI2*RHO
C      PARAMETER (QCI4=0.02D0)  !small step size for iteration
      PARAMETER (QCI5=1.04D0)  !branch density in equations 5.44, 5.45
      PARAMETER (QCI6=1.D-7)   !criterion for convergence (delta rho)/rho
      PARAMETER (QCI7=0.1D0)   !iteration step for last-ditch iteration
      PARAMETER (QCI8=1.00001D0) !reset size for RS = QCI8* RHO if RS < RL
C
C        PROCESS INPUT FLAGS AND BRANCH TO APPROPRIATE AREAS OF PROGRAM
C
      QCI4=0.02D0               ! Small step size for iteration (first attempt)
      IERR_IN=IERR              ! Remember input value in case of second attempt
 10   NP=0                      ! flag for entry mode: 0 is normal, fast iteration
      FAST=.TRUE.               ! fast iteration if this remains true
      IF(ACK(L+72).NE.ZERO) FAST=.FALSE.  !nope, we're going slow
      IF(FAST.OR.(IERR.LE.0)) THEN        !fast iteration, or IERR = 0
        IF(T.GT.ACK(L+18)) THEN !T is > Tm; ACK(L+18) = melt temperature
          IF((RHO.LE.ACK(L+47)).OR.(T.GE.ACK(L+48))) GO TO 260 !exit as melt
          GOTO 90               !check melt region boundaries
        END IF                  !T is < = Tm
        IF (RHO.GE.ACK46)     GO TO 250  !RHO above triple point density, exit as solid
        IF (RHO.LT.ACK(L+23)) GO TO 260  !RHO less than RHOMIN, exit as melt
        IF (T.LT.ACK(L+49))   GO TO 250  !T below lower bound for melt, exit as solid
      ELSE                      !slow iteration and IERR > = 1
        IF( IERR.LE.3) NP=IERR
        IF((IERR.EQ.2).AND.(T.GT.ACK(L+18))) GOTO 90 !check melt region boundaries
      END IF
C
C       set initial estimates of RL, RS in case of T < Tm: this occurs only
C       for tensile pressure extension, P < 0.
C
      RL=ACK(L+47)        !set initial RL to liquid density at triple point
      RS=ACK46            !set initial RS to solid density at triple point
      GO TO 120           !branch to slow iteration
C
C       TEST WHETHER INPUT (T, RHO) LIES IN POSSIBLE MELT/SOLID REGION--
C           SKIP ITERATION IF IT DOES NOT
C
   90 RL=ACK(L+47)+(T-ACK(L+18))/ACK(L+50) !eq. 5.41 for boundary curve A
      IF (RHO.LE.RL) GO TO 260             !to left of boundary, exit as melt
      RS=RL*ACK46/ACK(L+47) ! RS=RL*(solid density/liquid density at triple pt.)
      IF (FAST) THEN
        IF (T.LE.ACK(L+52)) THEN     ! T in lower melt region
          SS46=ACK(L+18)+ACK(L+51)*(RHO-ACK(L+11))  !eq. 5.42 for boundary curve B
          IF (T.LT.SS46) GO TO 250   !exit as solid
        END IF
        IF (RHO.LE.RS) GO TO 120  !use slow iteration in this case
        GSX=QCI1      !use large iteration step GSX
        RL=QCI2*RHO   !initial guess at liquid density = 0.99*RHO
        RS=QCI3*RHO   !initial guess at solid density  = 1.01*RHO
        GOTO 130
      END IF
  120 GSX=QCI4        !use small iteration step GSX
C
C    MAIN ITERATION LOOP TO GET PROPERTIES OF COEXISTING MELT, SOLID PHASES
C    USING TWO-VARIABLE NEWTONIAN ITERATION METHOD TO EQUATE PRESSURE
C    AND GIBBS ENERGY IN BOTH PHASES.
C         SEE SECTION V-1 OF THOMPSON AND LAUSON (1972).
C
  130 SS46=ACK(L+46)  !density of solid at triple point, saved for later
C
      DO LOOP=1,MAXIT     !start of main iteration loop
C
        ACK46=ZERO  !set ACK46 to zero to prevent ANEOS1 from adding melt terms
C
C        call ANEOS1 to get pure solid thermodynamic state at given T, RS
C
        CALL ANEOS1 (T,RS,PS,ES,SS,CVS,DPDTS,DPDRS,L)
C
        LNKNTR(MATKTR,13)=LNKNTR(MATKTR,13)+1   !count access for link statistics
        IF (DPDRS.LE.ZERO) GO TO 300  !bombed because of infinite solid compressibility
        FKS=FKROS     !save solid opacity for later
        ACK46=SS46    !restore ACK46 so melt terms are now added
C
C        call ANEOS1 to get pure liquid thermodynamic state at given T, RL
C
        CALL ANEOS1 (T,RL,PL,EL,SL,CVL,DPDTL,DPDRL,L)
C
        LNKNTR(MATKTR,13)=LNKNTR(MATKTR,13)+1  !count access for link statistics
        IF (DPDRL.LE.ZERO) GO TO 300  !bombed because of infinite liquid compressibility
C
C        use solid, melt thermo. data to do two-variable newtonian search
C        for point at which pressure, gibbs energy of solid and melt match
C
        X1=PL-PS      !pressure difference in eq. 5.36 (actually with opposite sign)
        GL=EL-T*SL+PL/RL  !Gibbs energy of liquid
        GS=ES-T*SS+PS/RS  !Gibbs energy of solid
        X2=GL-GS      !Gibbs energy difference in eq. 5.37 (again, with flipped sign)
        DRL=DPDRS*(ONE/RS-ONE/RL)  !denominator of eq. 5.39
        IF (DRL.EQ.ZERO) GO TO 270 !...which had better not be zero!
        DRS=(X2-X1/RL)/DRL         !delta rho of solid in equation 5.39
        DRL=(DRS*DPDRS-X1)/DPDRL   !delta rho of liquid in equation 5.40
C
        IF (FAST) THEN !do the following tests only for the fast iteration
C
C       Implement equation 5.45 to make quick decision on whether we are
C       denfinitely in the solid region
C
          IF((RHO.GT.QCI5*RS).AND.(DRL.LT.ZERO).AND.(DRS.LT.ZERO))
     &        GOTO 250    !exit as solid
C
C       Implement equation 5.44 to make quick decision on whether we are
C       denfinitely in the melt region
C
          IF((QCI5*RHO.LT.RL).AND.(DRL.GT.ZERO).AND.(DRS.GT.ZERO))
     &        GOTO 260    !exit as melt
        END IF
C
        ADRS=ABS(DRS)
        ADRL=ABS(DRL)
        IF (ADRL.LE.QCI6*RL) EXIT !convergence achieved, exit iteration loop
C
        IF (LOOP.GE.(MAXIT-100)) THEN   !within 100 steps of maximum!
          DRS=QCI7*DRS                  !increase step size for last-ditch...
          DRL=QCI7*DRL                  !attempt to converge
          IF (LOOP.EQ.MAXIT) GO TO 280  !too many iterations, print error message
        END IF
C
        GL=GSX*RS     !biggest acceptable step size for solid
        IF (ADRS.GT.GL) DRS=SIGN(GL,DRS)  !limit solid step size if too big
        GL=GSX*RL     !biggest acceptable step size for melt
        IF (ADRL.GT.GL) DRL=SIGN(GL,DRL)  !limit melt step size if too big
C
        RS=RS+DRS    !update solid density
        RL=RL+DRL    !update liquid density
        IF (RS.LE.RL) RS=QCI8*RL  !liquid is denser than solid! try again
      END DO           !loop back for next iteration
C
C       WRAP UP ITERATION, MAKE FINAL DECISIONS ON PHASE
C
      IF (RHO.LE.RL) GO TO 260  !but we are definitely melt
      IF (RHO.GE.RS) GO TO 250  !au contraire, we are definitly solid
      IERR=2                    !actually, we are mixed
C
C      IN LIQUID-SOLID REGION HERE, SO COMPUTE THE DETAILED THERMODYNAMIC
C      PROPERTIES OF THE MELT/SOLID MIXTURE
C
      X2=RHO*(RS-RL)
      DRS=RS*(RHO-RL)/X2   !solid mass fraction
      DRL=RL*(RS-RHO)/X2   !liquid mass fraction
      DPDTM=(SL-SS)*((RS*RL)/(RS-RL))  !dP/dT by Clapeyron slope
      DRLDT=(DPDTM-DPDTL)/DPDRL  !dRL/dT in liquid phase
      DRSDT=(DPDTM-DPDTS)/DPDRS  !dRS/dT in solid phase
      X1=-RHO*(RL*(RHO-RL)*DRSDT+RS*(RS-RHO)*DRLDT)/X2**2
      EM=DRS*ES+DRL*EL          !energy of mixture
      SM=DRS*SS+DRL*SL          !entropy of mixture
      FKROS=DRS*FKS+DRL*FKROS   !rosseland mean opacity of mixture
      CVM=X1*(ES-EL)+DRS*(CVS+(PS-T*DPDTS)*DRSDT/RS**2)
     1    +DRL*(CVL+(PL-T*DPDTL)*DRLDT/RL**2) !heat capacity of mixture
      DPDRM=ZERO                !compressibility of *mixture* is infinite
      PM=DRS*PS+DRL*PL          !PS should equal PL, but just in case...
      CRHS2P=RS                 !save RS for transmission in common /ANE2/
      CRLS2P=RL                 !save RL for transmission in common /ANE2/
      GO TO 320
  250 IERR=1             !destination of all solid determinations
      GO TO 320
  260 IERR=3             !destination of all liquid determinations
      GO TO 320
  320 IF (NP.EQ.2) IERR=LOOP !when NP = 2, IERR is the number of iterations
      RETURN
C
C        PROCESS ERROR CONDITIONS BEFORE RETURNING TO CALLING PROGRAM
C
  270 IERR=-1    !fatal error, denominator of eq. 5.39 is zero
      GO TO 290
  280 IERR=-2    !fatal error, too many iterations (> MAXIT)
      GO TO 290
  300 IERR=-4    !fatal error, zero compressibility of solid or liquid
      IF (NP.EQ.3) IERR=-3   ! set error code for NP = 3
  290 CONTINUE
      ! If convergence fails with a large (small) density step size,
      ! reduce the step size by a factor of two and try again. . .
      IF (QCI4.GT.0.01D0) THEN 
         QCI4=0.01D0
         IERR=IERR_IN
         GOTO 10
      END IF
      IF (IERR.EQ.-3) RETURN !if NP = 3, no error message is printed
      CALL ANMARK('FATAL ERROR IN ANELSM  IERR=') !bad message to get
      WRITE(ERRMES,330) IERR,T,RHO,RS,RL,DRS,DRL  !echo conitions at failure 
      CALL ANMARK(ERRMES)
      IF (NP.EQ.0) THEN   !in execution mode
        CALL ANMARK('300 IN ANELSM')
        STOP   !we're dead, nothing more we can do during execution
      END IF
      RETURN   !return to calling program without stopping
  330 FORMAT (I3,1P6E12.5)
      END
