C
C
      SUBROUTINE ANMARK (MESSAG)
C
C***********************************************************************
C
C     ANEOS PACKAGE  Write message to listing, log and terminal
C
C************************************* 10/88 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE   /FILEOS/
      CHARACTER*(*) MESSAG
      WRITE(KLST,10,ERR=20) MESSAG
   10 FORMAT(/,1X,A)
   20 CONTINUE
      RETURN
      END
