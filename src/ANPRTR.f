C
C
      SUBROUTINE ANPRTR
C
C***********************************************************************
C
C     LIST COUNTERS FOR ANEOS PACKAGE. IF CODE IS USED IN PARALLEL
C     PROCESSING MODE, THE COUNTERS MIGHT NOT BE CORRECT.
C
C************************************** 5/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
      SAVE   /FILEOS/
      PARAMETER (ZERO=0.D0)
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
C     aneos counters are as follows (for eos number M)
C     LNKNTR(M,1) = not used
C     LNKNTR(M,2) = calls to aneosv (once for each vector element)
C     LNKNTR(M,3) = total calls to aneos1
C     LNKNTR(M,4) = calls to anion1
C     LNKNTR(M,5) = calls to anion2
C     LNKNTR(M,6) = calls to antomf
C     LNKNTR(M,7) = calls to aneos1 from aneosv
C     LNKNTR(M,8) = calls to ane2ph from aneosv
C     LNKNTR(M,9) = calls to anelsm from aneosv
C     LNKNTR(M,10)= calls to anelsm from ane2ph
C     LNKNTR(M,11)= total calls to anelsm
C     LNKNTR(M,12)= calls to aneos1 from ane2ph
C     LNKNTR(M,13)= calls to aneos1 from anelsm
C     LNKNTR(M,14)= sum of iterations in anion2
C     LNKNTR(M,15)= calls to anueos
C     LTKNTR      = calls to eos (table routines)
C     MATKTR      = current value of M above
C
C                    | 2
C                    |
C                  ANEOSV
C                 /  |  \ \ ___
C              9 /   |8  \      \ 15
C               /    |    \      \
C              /   ANE2PH  \      ANUEOS (user defined)
C             /   /      \  \
C            /   / 10     \  \ 7
C           ANELSM      12 \  \
C                \          \  \
C                 \     13   \  \    3=7+12+13     11=9+10
C                  --------- ANEOS1
C                              |
C         --------------------------------------------------
C            |         |          | 4        | 5        | 6
C            |         |          |          |          |
C         ANDEBY     ANEI3      ANION1     ANION2     ANTOMF
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE  /ANCNTR/
      DIMENSION FNKNTR(MAXMAT+1)
      CHARACTER*20 ME
      CHARACTER*4 MAD(MAXMAT),MADA
      PARAMETER (MADA='ALL')
   10 FORMAT(1X,A,11I10)
   12 FORMAT(1X,A,1P11E10.3)
   20 FORMAT(//,' EOS package call numbers',/,2X,A,11(6X,A))
   22 FORMAT(' Calls to')
      IMX=1
        DO 100 IM=1,MAXMAT
        IF((LNKNTR(IM,1).GT.0).OR.(LNKNTR(IM,2).GT.0)) IMX=IM
        LNKNTR(IM,3)=LNKNTR(IM,7)+LNKNTR(IM,12)+LNKNTR(IM,13)
        LNKNTR(IM,11)=LNKNTR(IM,9)+LNKNTR(IM,10)
  100   CONTINUE
      DO 110 IE=1,KDIMTR
      LNKNTR(MAXMAT+1,IE)=0
        DO 120 IM=1,IMX
        LNKNTR(MAXMAT+1,IE)=LNKNTR(MAXMAT+1,IE)+LNKNTR(IM,IE)
  120   CONTINUE
  110 CONTINUE
      IF((LNKNTR(MAXMAT+1,2).EQ.0).AND.(LTKNTR.EQ.0)) THEN
        WRITE(KLST,114)
  114   FORMAT(//,' EOS package call numbers NULL',/)
        GO TO 140
      END IF
      IMX=MIN(IMX,10)
        DO 170 IEE=1,MIN(IMX,9)
        WRITE(MAD(IEE),'(I2,2X)') -IEE
  170   CONTINUE
      IF(IMX.EQ.10) MAD(10)='-10'
      ME='      ANEOS number ='
      WRITE(KLST,20) ME,(MAD(IE),IE=1,IMX),MADA
      WRITE(KLST,22)
      DO 130 IEE=1,KDIMTR
      IE=IEE+1
      IF(IEE.EQ.13) IE=15
      IF(IEE.EQ.14) IE=1
      IF(IEE.EQ.15) IE=14
      IF(IE.EQ.1) THEN
        ME='ANEOS1/ANEOSV'
          DO 126 IM=1,MAXMAT+1
          IF(LNKNTR(IM,2).GT.0) THEN
            FNKNTR(IM)=REAL(LNKNTR(IM,3))/REAL(LNKNTR(IM,2))
          ELSE
            FNKNTR(IM)=0.
          END IF
  126     CONTINUE
        IF(FNKNTR(MAXMAT+1).GT.ZERO) THEN
          WRITE(KLST,12) ME,(FNKNTR(IM),IM=1,IMX),FNKNTR(MAXMAT+1)
        END IF
        GO TO 130
      ELSEIF(IE.EQ.2) THEN
        ME='ANEOS'
      ELSEIF(IE.EQ.3) THEN
        ME='ANEOS1'
      ELSEIF(IE.EQ.4) THEN
        ME='ANION1'
      ELSEIF(IE.EQ.5) THEN
        ME='ANION2'
      ELSEIF(IE.EQ.6) THEN
        ME='ANTOMF'
      ELSEIF(IE.EQ.7) THEN
        ME='ANEOS1 from ANEOSV'
      ELSEIF(IE.EQ.8) THEN
        ME='ANE2PH from ANEOSV'
      ELSEIF(IE.EQ.9) THEN
        ME='ANELSM from ANEOSV'
      ELSEIF(IE.EQ.10) THEN
        ME='ANELSM from ANE2PH'
      ELSEIF(IE.EQ.11) THEN
        ME='ANELSM'
      ELSEIF(IE.EQ.12) THEN
        ME='ANEOS1 from ANE2PH'
      ELSEIF(IE.EQ.13) THEN
        ME='ANEOS1 from ANELSM'
      ELSEIF(IE.EQ.15) THEN
        ME='ANUEOS'
      ELSEIF(IE.EQ.14) THEN
        ME='ANION2 iter. count'
          DO 122 IM=1,MAXMAT+1
          IF(LNKNTR(IM,5).GT.0) THEN
            FNKNTR(IM)=REAL(LNKNTR(IM,14))/REAL(LNKNTR(IM,5))
          ELSE
            FNKNTR(IM)=0.
          END IF
  122     CONTINUE
        IF(FNKNTR(MAXMAT+1).GT.ZERO) THEN
          WRITE(KLST,12) ME,(FNKNTR(IM),IM=1,IMX),FNKNTR(MAXMAT+1)
        END IF
        GO TO 130
      END IF
      IF((LNKNTR(MAXMAT+1,IE).GT.0).OR.(IE.EQ.2)) THEN
        WRITE(KLST,10) ME,(LNKNTR(IM,IE),IM=1,IMX),LNKNTR(MAXMAT+1,IE)
      END IF
  130 CONTINUE
      WRITE(KLST,160) LTKNTR
  160 FORMAT(' Calls to eos tables =',I10)
  140 CONTINUE
      RETURN
      END
