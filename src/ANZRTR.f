C
C
      SUBROUTINE ANZRTR
C
C***********************************************************************
C
C     ZERO COUNTERS FOR ANEOS PACKAGE
C
C************************************** 5/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE  /ANCNTR/
      MATKTR=1
      LTKNTR=0
      DO 10 IE=1,KDIMTR
        DO 20 IM=1,MAXMAT+1
        LNKNTR(IM,IE)=0
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
