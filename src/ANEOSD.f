C
C
      SUBROUTINE ANEOSD(TQQ,RHOQQ,PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,
     1                  FKROQQ,CSQQ,KPA,MATQQ)
C
C***********************************************************************
C
C     ANEOS PACKAGE - RUNNING ENTRY POINT - DOUBLE PRECISION ENTRY
C     See ANEOSV for argument list
C
C************************************** 9/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CALL ANEOSV (1,TQQ,RHOQQ,MATQQ,
     &  PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,FKROQQ,CSQQ,KPA,R2PLQQ,R2PHQQ,
     &  ZBARQQ)
      RETURN
      END
