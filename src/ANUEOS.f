C
C
      SUBROUTINE ANUEOS(
C         inputs
     &    TEMP,RHO,MATNUM,UARR,
C         outputs
     &    P,E,S,CV,DPDT,DPDR,FKROS,KP,R2PL,R2PH,ZB)
C
C***********************************************************************
C
C     User defined EOS running entry point
C     See ANUSET for set up of data for this routine.
C
C     WARNING - depending on the computer system
C     All variables are single or double precision.
C
C     Inputs: (all units cgs-ev)
C
C     TEMP   =  temperature
C     RHO    =  density
C     MATNUM =  material number
C     UARR   =  aneos data array for this material
C               (see ANUSET for array information)
C
C     Output:
C
C     P     =  pressure
C     E     =  specific energy
C     S     =  specific entropy
C     CV    =  specific heat capacity
C     DPDT  =  dP / dT
C     DPDR  =  dP / dRHO
C     FKROS =  Rosseland mean opacity (including conduction)
C     KP    =  phase indicator (see below)
C     R2PL  =  lower phase density for mixed-phase state
C              not used for single phase
C     R2PH  =  upper phase density for mixed-phase state
C              not used for single phase
C     ZB    =  ionization number
C
C     KP   =STATE INDICATOR      =1, 1p
C                                =2, 2p
C
C************************************** 5/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UARR(99)
C     klst is the output listing file
      COMMON /FILEOS/ KLST,KINP
      SAVE   /FILEOS/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
C
C     add your code here
C
C     example code for ideal gas ---------------------------------------
      CV=UARR(2)
      GAMMAN=UARR(1)-ONE
      E=CV*TEMP
      DPDR=GAMMAN*E
      DPDT=GAMMAN*CV*RHO
      P=DPDT*TEMP
      S=CV*(LOG(TEMP/UARR(12))-GAMMAN*LOG(RHO/UARR(11)))
      KP=1
      R2PL=ZERO
      R2PH=ZERO
      ZB=ZERO
      FKROS=1.E6
C     example code for ideal gas ---------------------------------------
C
      RETURN
      END
