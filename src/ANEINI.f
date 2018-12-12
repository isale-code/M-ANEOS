C
C
      SUBROUTINE ANEINI
C
C***********************************************************************
C
C     This routine should only be used if aneos is used WITHOUT the
C     tabular EOS routines and under multitasking conditions.
C     A call to this routine should be made when starting a new task.
C     It sets a few variables in task common.
C
C************************************** 8/89 slt ***********************
C
C     multitasking note - In multitasking mode, the following commons
C     should be task common. Each time a new task is started,
C     routine ANEINI should be called.
C
C     Other task common save answers and return answers.
C
C     The following aneos common should be task common
C       /ANE2/ /BNES/ /FNES/ /ANZB/ /ANZLTF/ /ANEEL/ /ANESQT/
C
C     The following aneos common have call counters that will not
C     give an accurate count in multitask mode. This can be fixed
C       /ANCNTR/           with critical statements in the code
C                          if desired. Values are for print only.
C
C     The following task common should be included in the calling
C     routine.    /MIXCON/ /ANZLTF/
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /MIXCON/  RHOMPH(MAXMAT),RHOMPL(MAXMAT),IPEOS
      SAVE   /MIXCON/
      COMMON /ANZLTF/ ZLTOMF
      IPEOS=1
      ZLTOMF=-1.
      RETURN
      END
