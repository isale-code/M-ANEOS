C
C
	SUBROUTINE EOS(TEMP,TEMP0,P,RHO,E,M)
C
C     PROGRAM TO USE ANEOS TO MATCH INPUT ENERGY DENSITY AT CONSTANT DENSITY.
C     RHO AND E ARE NOT CHANGED, TEMP, P ARE OUTPUTS.  TEMP0 IS A TRIAL
C     VALUE THAT MUST BE DEFINED ON INPUT.
C
C        INPUTS:  TEMP0 = Initial estimate of temperature
C                 RHO   = Density
C                 E     = Internal energy to be matched by EOS
C                 M     = Material number
C
C       OUTPUTS:  TEMP  = Temperature (in eV) at given RHO, E
C                 P     = Pressure, P(RHO,E)
C
C       All units are cgs-eV
C
C     This routine generally uses Newton-Raphson iteration, which is fast
C     and ususally converges rapidly, because CV is positive-definite and E
C     is a monotonically increasing function of TEMP.  However, N-R can 
C     occasionally become trapped in an infinite loop, cycling around the
C     actual root.  This is especially likely where phase transitions are 
C     present.  In this case, this routine resorts to bisection to converge
C     on the actual root.
C
C     Note that this routine requires estimates of the maximum and minimum
C     temperatues likely to be encountered.  The user must take care that
C     these are reasonable limits for his/her particular problem!
C
C     The COST of this security is that this routine uses approximately TWICE
C     as many iterations as the simple Newton-Raphson EOS without bisection.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTEPS=100)
      LOGICAL TEST1,TEST2
      HALF=0.5D0
      ZERO=0.D0
      TWO=2.D0
      QDTOL=1.0D-6  !tolerance for temperature convergence
      TL=1.D-6      !minimum temperature for solution
      TH=1.0D2      !maximum temperature for solution
      DTOLD=HALF*(TH-TL)
      TEMP=TEMP0
      CALL ANEOSD(TEMP,RHO,P,ETRIAL,S,CV,DPDT,DPDR,FKROS,CS,KPA,M)
      DO IT=1,NSTEPS
        TEST1=((TEMP-TH)*CV-(ETRIAL-E))*((TEMP-TL)*CV-(ETRIAL-E))
     &         .GT.ZERO
        TEST2=DABS(TWO*(ETRIAL-E)).GT.DABS(DTOLD*CV)
        IF(TEST1.OR.TEST2) THEN        !Bisection step
          DTOLD=DT
          DT=HALF*(TH-TL)
          TEMP=TL+DT
        ELSE                           !Newton-Raphson step
          DTOLD=DT
          DT=(ETRIAL-E)/CV
          TEMP=TEMP-DT
        END IF
        IF(DABS(DT).LT.QDTOL) RETURN   !Successful completion
        CALL ANEOSD(TEMP,RHO,P,ETRIAL,S,CV,DPDT,DPDR,FKROS,CS,KPA,M)
        IF((ETRIAL-E).LT.ZERO) THEN
          TL=TEMP
        ELSE
          TH=TEMP
        ENDIF
      END DO
      STOP 'TOO MANY ITERATIONS--EOS'
      END
