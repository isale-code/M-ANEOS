C
C
      SUBROUTINE ANEHPP (T,RHO,RLPP,RHPP,L,IERR,X1,X2)
C
C***********************************************************************
C
C     ANEOS PACKAGE  SOLID-SOLID TWO-PHASE MAXWELL CONSTRUCTION FOR
C                    HIGH PRESSURE PHASE TRANSFORMATION
C
C              INPUTS:   T    = TEMPERATURE IN EV
C                        RHO  = DENSITY IN GM/CC
C                        L    = MATERIAL NUMBER OF HPP
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
C               OUTPUTS: RLPP  = DENSITY OF LPP AT T
C                        RHPP  = DENSITY OF HPP AT T
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
!      COMMON /HPPTRANS/ THPP(1000*MAXMAT)
!     2 ,RLPP(1000*MAXMAT),RHPP(1000*MAXMAT)
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
      PARAMETER (QCI2=0.99D0)  !first guess at RLPP = QCI2*RHO
      PARAMETER (QCI3=1.01D0)  !first guess at RHPP = QCI2*RHO
      PARAMETER (QCI4=0.02D0)  !small step size for iteration
      PARAMETER (QCI5=1.04D0)  !branch density in equations 5.44, 5.45
      PARAMETER (QCI6=1.D-7)   !criterion for convergence (delta rho)/rho
      PARAMETER (QCI7=0.1D0)   !iteration step for last-ditch iteration
      PARAMETER (QCI8=1.00001D0) !reset size for RHPP = QCI8* RHO if RHPP < RLPP

  120 GSX=QCI4          !use small iteration step GSX
      ACK46=ZERO                !set ACK46 to something handed in from ANEOS2
C
C    MAIN ITERATION LOOP TO GET PROPERTIES OF COEXISTING LPP, HPP PHASES
C    USING TWO-VARIABLE NEWTONIAN ITERATION METHOD TO EQUATE PRESSURE
C    AND GIBBS ENERGY IN BOTH PHASES.
C         SEE SECTION V-1 OF THOMPSON AND LAUSON (1972).
C
      DO LOOP=1,MAXIT     !start of main iteration loop
C
C
C        call ANEOS1 to get LPP thermodynamic state at given T, RHPP
C
        CALL ANEOS1 (T,RHPP,PS,ES,SS,CVS,DPDTHPP,DPDRHPP,L) !L is for HPP
        ES=ES+1.2D10  !+1.D10
        SS=SS+9.07D9  !+17.7D9
C     
        LNKNTR(MATKTR,13)=LNKNTR(MATKTR,13)+1   !count access for link statistics
        IF (DPDRHPP.LE.ZERO) GO TO 300  !bombed because of infinite HPP compressibility
        FKS=FKROS     !save solid opacity for later
C
C        call ANEOS1 to get pure HPP thermodynamic state at given T, RLPP
C
        CALL ANEOS1 (T,RLPP,PL,EL,SL,CVL,DPDTLPP,DPDRLPP,L-99)  !L-1 is for LPP
C
        LNKNTR(MATKTR,13)=LNKNTR(MATKTR,13)+1  !count access for link statistics
        IF (DPDRLPP.LE.ZERO) GO TO 300  !bombed because of infinite LPP compressibility
C
C        use solid, melt thermo. data to do two-variable newtonian search
C        for point at which pressure, gibbs energy of solid and melt match
C
        X1=PL-PS      !pressure difference in eq. 5.36 (actually with opposite sign)
        GL=EL-T*SL+PL/RLPP  !Gibbs energy of low-pressure phase
        GS=ES-T*SS+PS/RHPP  !Gibbs energy of high-pressure phase
        X2=GL-GS            !Gibbs energy difference in eq. 5.37 (again, with flipped sign)
        DRLPP=DPDRHPP*(ONE/RHPP-ONE/RLPP)  !denominator of eq. 5.39
        IF (DRLPP.EQ.ZERO) GO TO 270 !...which had better not be zero!
        DRHPP=(X2-X1/RLPP)/DRLPP         !delta rho of solid in equation 5.39
        DRLPP=(DRHPP*DPDRHPP-X1)/DPDRLPP   !delta rho of liquid in equation 5.40
C
        ADRHPP=ABS(DRHPP)
        ADRLPP=ABS(DRLPP)
        IF (ADRLPP.LE.QCI6*RLPP) EXIT !convergence achieved, exit iteration loop
C
        IF (LOOP.GE.(MAXIT-100)) THEN   !within 100 steps of maximum!
          DRHPP=QCI7*DRHPP                  !increase step size for last-ditch...
          DRLPP=QCI7*DRLPP                  !attempt to converge
          IF (LOOP.EQ.MAXIT) GO TO 280  !too many iterations, print error message
        END IF
C
        GL=GSX*RHPP     !biggest acceptable step size for solid
        IF (ADRHPP.GT.GL) DRHPP=SIGN(GL,DRHPP)  !limit solid step size if too big
        GL=GSX*RLPP     !biggest acceptable step size for melt
        IF (ADRLPP.GT.GL) DRLPP=SIGN(GL,DRLPP)  !limit melt step size if too big
C
        RHPP=RHPP+DRHPP    !update solid density
        RLPP=RLPP+DRLPP    !update liquid density
        IF (RHPP.LE.RLPP) RHPP=QCI8*RLPP  !liquid is denser than solid! try again
      END DO           !loop back for next iteration
C
C       WRAP UP ITERATION, MAKE FINAL DECISIONS ON PHASE
C
      IF (RHO.LE.RLPP) GO TO 260  !but we are definitely LPP
      IF (RHO.GE.RHPP) GO TO 250  !au contraire, we are definitly HPP
      IERR=2                    !actually, we are mixed
C
C      IN LPP-HPP REGION HERE, SO COMPUTE THE DETAILED THERMODYNAMIC
C      PROPERTIES OF THE LPP/HPP MIXTURE
C
      X2=RHO*(RHPP-RLPP)
      DRHPP=RHPP*(RHO-RLPP)/X2   !solid mass fraction
      DRLPP=RLPP*(RHPP-RHO)/X2   !liquid mass fraction
      DPDTM=(SL-SS)*((RHPP*RLPP)/(RHPP-RLPP))  !dP/dT by Clapeyron slope
      DRLPPDT=(DPDTM-DPDTLPP)/DPDRLPP  !dRLPP/dT in liquid phase
      DRHPPDT=(DPDTM-DPDTHPP)/DPDRHPP  !dRHPP/dT in solid phase
      X1=-RHO*(RLPP*(RHO-RLPP)*DRHPPDT+RHPP*(RHPP-RHO)*DRLPPDT)/X2**2
      EM=DRHPP*ES+DRLPP*EL          !energy of mixture
      SM=DRHPP*SS+DRLPP*SL          !entropy of mixture
      FKROS=DRHPP*FKS+DRLPP*FKROS   !rosseland mean opacity of mixture
      CVM=X1*(ES-EL)+DRHPP*(CVS+(PS-T*DPDTHPP)*DRHPPDT/RHPP**2)
     1    +DRLPP*(CVL+(PL-T*DPDTLPP)*DRLPPDT/RLPP**2) !heat capacity of mixture
      DPDRM=ZERO                !compressibility of *mixture* is infinite
      PM=DRHPP*PS+DRLPP*PL          !PS should equal PL, but just in case...
      CRHS2P=RHPP                 !save RHPP for transmission in common /ANE2/
      CRLS2P=RLPP                 !save RLPP for transmission in common /ANE2/
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
  300 IERR=-4    !fatal error, zero compressibility of HPP or LPP
      IF (NP.EQ.3) THEN   !if NP = 3, no error message is printed
        IERR=-3           !set error code for NP = 3
        RETURN            !back to calling program
      END IF
  290 CALL ANMARK('FATAL ERROR IN ANEHPP  IERR=') !bad message to get
      WRITE(ERRMES,330) IERR,T,RHO,RHPP,RLPP,DRHPP,DRLPP  !echo conitions at failure 
      CALL ANMARK(ERRMES)
      IF (NP.EQ.0) THEN   !in execution mode
        CALL ANMARK('300 IN ANELSM')
        !STOP   !we're dead, nothing more we can do during execution
      END IF
      RETURN   !return to calling program without stopping
  330 FORMAT (I3,1P6E12.5)
      END
