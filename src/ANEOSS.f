C
C
      SUBROUTINE ANEOSS (IGK,ITAPE)
C
C***********************************************************************
C
C     PART OF SET UP FOR ANEOS PACKAGE
C
C     THIS ROUTINE WAS PART OF ANEOS2 - SEPARATED TO MINIMIZE STORAGE
C     WHEN ONLY RESET WRITE OR READ IS NEEDED - ALSO RAD TURN OFF.
C
C     INPUTS
C
C       IGK = 1
C           = 2   WRITE RESTART FILE ITAPE
C           = 3   READ RESTART FILE ITAPE
C           = 4   TURN RADIATION OPACITY CALCULATION OFF
C
C************************************** 8/87 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /ANZB/ ZB(92), ZRAT
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
      COMMON /MIXCON/  RHOMPH(MAXMAT),RHOMPL(MAXMAT),IPEOS
      COMMON /ANZLTF/ ZLTOMF
      SAVE   /FILEOS/,/MIXCON/,/ANES/,/ANZB/
      IF(IGK.NE.2) CALL ANZRTR
      ZRAT=0.000045D0
      IPSQTS=1
C
      IF(IGK.EQ.4) THEN
C       RADIATION OPACITY FLAG
        NOTRAD=0
        GO TO 130
C
      ELSEIF(IGK.EQ.2) THEN
C       write restart data
        WRITE (ITAPE) ACK,ZZS,COT,FNI,RCT,TCT,RSOL,RVAP,TTWO,SAVER
     1      ,BOLTS,EIP,LOCSV,LOCKP,LOCKPL,NOTRAD
C
      ELSE
C       read restart data
        READ (ITAPE,ERR=120,END=120)
     1       ACK,ZZS,COT,FNI,RCT,TCT,RSOL,RVAP,TTWO,SAVER
     2      ,BOLTS,EIP,LOCSV,LOCKP,LOCKPL,NOTRAD
        IPEOS=1
        ZLTOMF=-1.
      END IF
      GO TO 130
C
  120 CONTINUE
      CALL ANMARK('eof or error reading aneos restart')
      STOP
  130 CONTINUE
      RETURN
      END
