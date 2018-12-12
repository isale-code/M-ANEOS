C
C
      SUBROUTINE ANDEBY(THETA,TK,D,D2,DS)
C
C***********************************************************************
C
C     ANEOS PACKAGE         EVALUATES DEBYE FUNCTIONS
C
C************************************** 8/87 version slt ***************
C
C     INPUTS:  TK    = Temperature of evaluation
C              THETA = Debye temperature
C
C     OUTPUTS: D  =  Debye integral, D(X)
C              D2 =  function ALOG(1-EXP(-X))
C              DS =  function 3*(D(X)-X/(EXP(X)-1))
C              where X = THETA/TK
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE /FILEOS/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (SIX=6.D0)
      PARAMETER (QCG1=0.1D0)
      PARAMETER (QCG2=0.375D0)
      PARAMETER (QCG3=0.05D0)
      PARAMETER (QCG4=5.952380953D-4)
      PARAMETER (QCG5=7.25D0)
      PARAMETER (QCG6=0.0946173D0)
      PARAMETER (QCG7=4.432582D0)
      PARAMETER (QCG8=85.07724D0)
      PARAMETER (QCG9=800.6087D0)
      PARAMETER (QCG10=3953.632D0)
      PARAMETER (QCG11=15.121491D0)
      PARAMETER (QCG12=143.155337D0)
      PARAMETER (QCG13=682.0012D0)
      PARAMETER (QCG14=3953.632D0)
      PARAMETER (QCG15=25.D0)
      PARAMETER (QCG16=6.493939402D0)
C
      X=THETA/TK
      EXX=EXP(-X)
      IF(X.GT.QCG1) GO TO 5
        X2=X*X
        D=ONE-QCG2*X+X2*(QCG3-QCG4*X2)
        GO TO 25
    5 IF(X-QCG5) 10,10,13
   10 D=((((QCG6*X-QCG7)*X+QCG8)*X-QCG9)*X+QCG10)
     1 /((((X+QCG11)*X+QCG12)*X+QCG13)*X+QCG14)
      GO TO 25
   13 N=QCG15/X
      D=ZERO
      IF (N.LE.0) GO TO 20
      D2=ONE
        DO 15 I=1,N
        DS=I
        D2=D2*EXX
        X3=DS*X
   15   D=D+D2*(SIX+X3*(SIX+X3*(THREE+X3)))/(DS**4)
   20 D=THREE*(QCG16-D)/(X**3)
   25 D2=LOG(ONE-EXX)
      DS=THREE*(D-X*(EXX/(ONE-EXX)))
      RETURN
      END
