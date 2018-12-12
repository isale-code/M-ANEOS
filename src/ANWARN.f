C
C
      SUBROUTINE ANWARN (IWAY)
C
C***********************************************************************
C
C     ANEOS PACKAGE  WRITE WARNING (IWAY=0)
C                         OR FATAL (IWAY=1) ERROR MESSAGE
C
C************************************* 10/88 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ANEST/ NOEOSU
      SAVE   /ANEST/
      CHARACTER*11 MESTYP
      CHARACTER*6  NESTYP
      IF(IWAY.EQ.0) THEN
         MESTYP='WARNING'
         NC=7
      ELSE
         MESTYP='FATAL error'
         NC=11
      END IF
      NESTYP=' '
      WRITE(NESTYP,'(I6)',ERR=10) NOEOSU
   10 CALL ANMARK(MESTYP(1:NC)
     & //' message from ANEOS package - eos number'//NESTYP)
      RETURN
      END
