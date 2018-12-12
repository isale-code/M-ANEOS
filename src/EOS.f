C
C
      SUBROUTINE EOS(TEMP,TEMP0,P,RHO,E,M)
C
C	PROGRAM TO USE ANEOS TO MATCH INPUT ENERGY DENSITY
C
C       This is the simplest possible algorithm, that sometimes
C       fails due to classic Newton-Raphson nonconverging loop,
C       especially near phase transitions.  If trouble is
C       encountered, use EOS_safe, a routine that backs the N-R
C       method up with bisection.
C
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
      IMPLICIT REAL*8 (A-H,O-Z)
      TEMP=TEMP0
      IT=0
 10   CONTINUE
      IT=IT+1
      IF(IT.GT.50) STOP 'TOO MANY ITERATIONS--EOS'
      CALL ANEOSD (TEMP,RHO,P,ETRIAL,S,CV,DPDT,DPDR,FKROS,CS,KPA,M)
      IF (ABS(E-ETRIAL) .LT. 1.0D-9*E) RETURN
      TEMP=TEMP+(E-ETRIAL)/CV
      GOTO 10
      END
