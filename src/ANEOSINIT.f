      SUBROUTINE ANEOSINIT(INPUT)
C**********************************************************************
C     INITIALIZE ANEOS FOR THE C WRAPPER
C**********************************************************************
      PARAMETER (NMATMAX=21)
      PARAMETER (NUMMAT=1)
      DIMENSION IZETL(NMATMAX)
      COMMON /FILEOS/ KLST, KINP

      CHARACTER*(*) INPUT

C     INITIALIZE INPUT AND OUTPUT FILES
      KLST=10
      KINP=12
      OPEN(KINP,FILE=INPUT,STATUS='OLD')
      OPEN(KLST,FILE='ANEOS.OUTPUT')

C     INITIALIZE ANEOS
      IZETL(1)=-1
      CALL ANEOS2 (1,NUMMAT,0,IZETL)


C     CLOSE INPUT AND OUTPUT FILES
      CLOSE(KINP)
      CLOSE(KLST)

      RETURN
      END
