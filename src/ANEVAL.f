C
C
      SUBROUTINE ANEVAL(MAT,NUM,IWAY,VALUE,IERR)
C
C***********************************************************************
C
C     This routine added to allow external access to aneos common blocks
C     without having to change blocks in various programs. Values can be
C     returned or reset. Aneos2 or restart routines should be completed
C     before using this entry.
C
C     INPUTS:
C
C     MAT = aneos material number
C
C     NUM = number of variable
C
C     IWAY = 0, return value of variable
C          = 1, set new value of variable
C
C     VALUE = new value (only if IWAY=1)
C
C     OUTPUTS:
C
C     VALUE = value of variable (only if IWAY=0)
C
C     IERR = 0, no error
C          = 1, some type of error, value not set
C
C************************************** 8/87 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      SAVE   /ANES/
      IMAT=ABS(MAT)
C
C     check for allowed material
C
      IF((IMAT.EQ.0).OR.(IMAT.GT.MAXMAT)) THEN
        IF(IWAY.NE.0) VALUE=0.D0
        IERR=1
        GO TO 10
      END IF
      L=LOCSV(IMAT)
      IF(IWAY.EQ.0) THEN
        VALUE=ACK(L+NUM)
      ELSE
        ACK(L+NUM)=VALUE
      END IF
      IERR=0
   10 RETURN
      END
