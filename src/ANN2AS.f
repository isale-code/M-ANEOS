C
C
      SUBROUTINE ANN2AS (ETA,KTYPE,S1,AK1,AK2,AK3,AK4,GAMMAO,VALUE)
C
C***********************************************************************
C
C     Preform numerical integral for Th/Ph type -1 
C        (Low temperature solid)
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE   /FILEOS/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (TEN=10.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (KUM=200)
      IF(ABS(ETA).LE.1.E-7) THEN
        IF(KTYPE.EQ.0) THEN
          VALUE=TWO*S1/THREE
        ELSE
          VALUE=ONE
        END IF
        GO TO 20
      END IF
      NUM=MAX(KUM,KUM*INT(ABS(TEN*ETA)))
      DN=FLOAT(NUM)
      DN=ETA/DN
      VALUE=ZERO
      AL=ZERO
        DO 10 I=1,NUM
        ETAV=REAL(I)*DN
        IF(KTYPE.EQ.0) THEN
          AV=TWO*S1/(ONE-S1*ETAV)**3
        ELSE
          AV=AK1+TWO*AK2*ETAV+THREE*AK3*ETAV**2+FOUR*AK4*ETAV**3
        END IF
        AV=AV*EXP(-GAMMAO*ETAV)*ETAV**2
        VALUE=VALUE+HALF*(AL+AV)*DN
   10   AL=AV
      IF(KTYPE.EQ.0) THEN
        VALUE=VALUE*(ONE-S1*ETA)**3/(ETA**3)
      ELSE
        VALUE=THREE*VALUE/(AK1*ETA**3)
      END IF
   20 CONTINUE
c       write(klst,1313) eta,value
c1313   format(' int',1p2e15.7)
      RETURN
      END
