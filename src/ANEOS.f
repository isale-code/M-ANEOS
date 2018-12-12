C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * *        DOUBLE PRECISION ANEOS PACKAGE             * * * * *
C * * * * * * * * * *         WARNING               * * * * * * * * * *
C *  ALL REAL VARIABLES (LOCAL, COMMON AND ARGUMENT) ARE TYPE DOUBLE  *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE ANEOS (T,RHO,P,E,S,CV,DPDT,DPDR,FKRO,CS,KPA,MAT)
C
C***********************************************************************
C
C     ANEOS PACKAGE - RUNNING ENTRY POINT - SINGLE PRECISION ENTRY
C     See ANEOSV for argument list
C
C************************************** 9/89 version slt ***************
C
C     ANEOS ARGUMENTS
      DOUBLE PRECISION TQQ,RHOQQ,PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,
     1                 FKROQQ,CSQQ,R2PLQQ,R2PHQQ,ZBARQQ
C
      TQQ=T
      RHOQQ=RHO
      CALL ANEOSV (1,TQQ,RHOQQ,MAT,
     &  PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,FKROQQ,CSQQ,KPA,R2PLQQ,R2PHQQ,
     &  ZBARQQ)
      P=PQQ
      E=EQQ
      S=SQQ
      CV=CVQQ
      DPDT=DPDTQQ
      DPDR=DPDRQQ
      FKRO=FKROQQ
      CS=CSQQ
      RETURN
      END
