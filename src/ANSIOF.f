C
C
      SUBROUTINE ANSIOF (KOFILE,KIFILE)
C
C***********************************************************************
C
C     ANEOS PACKAGE  Set list and read files for package
C
C     Inputs:
C
C     KOFILE = Listing output file
C
C     KIFILE = Input file (used only during set up)
C
C     These files should be opened by the calling program
C
C************************************************ 10/88 slt ************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE   /FILEOS/
C
      KLST=KOFILE
      KINP=KIFILE
      RETURN
      END
