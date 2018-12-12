C
C
      SUBROUTINE ANNDPR(T,RHO,P,E,S,CV,DPT,DPR,IMAT)
C
C**********************************************************************
C
C     ANEOS PACKAGE - alternate call to aneos1 which calculates
C     numerical derivative for dP/dRHO - In some case analytic
C     derivate is not good because of numerics.
C     Parameters are the same as aneos1
C
C     DEL is 2X fractional step size
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DEL=0.001D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      CALL ANEOS1(T,RHO*(ONE+DEL),PCP,E,S,CV,DPT,DPR,IMAT)
      CALL ANEOS1(T,RHO*(ONE-DEL),PCN,E,S,CV,DPT,DPR,IMAT)
      CALL ANEOS1(T,RHO,P,E,S,CV,DPT,DPR,IMAT)
      DPR=(PCP-PCN)/(TWO*DEL*RHO)
      RETURN
      END
