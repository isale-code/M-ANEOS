C
C
      SUBROUTINE ANEOSV (NUMQQ,TQQ,RHOQQ,MATQQ,
     &  PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,FKROQQ,CSQQ,KPAQQ,R2PLQQ,R2PHQQ,
     &  ZBARQQ)
C
C***********************************************************************
C
C     ANEOS PACKAGE - VECTOR RUNNING ENTRY POINT
C
C     WARNING - depending on the computer system
C     All variables are single or double precision.
C
C     INPUTS (ALL UNITS CGS-EV)
C
C     NUMQQ =  dimensions of input and output arrays
C              (normally limited to 64)
C     TQQ   =  temperature array
C     RHOQQ =  density array
C     MATQQ =  material number array
C              skip any materials with matqq <= 0
C
C     Output arrays
C
C     PQQ   =  pressure
C     EQQ   =  specific energy
C     SQQ   =  specific entropy
C     CVQQ  =  specific heat capacity
C     DPDTQQ=  dP / dT
C     DPDRQQ=  dP / dRHO
C     FKROQQ=  Rosseland mean opacity (including conduction)
C     CSQQ  =  sound speed
C     KPAQQ =  phase indicator (see below)
C     R2PLQQ=  lower phase density for mixed-phase state
C              not used for single phase
C     R2PHQQ=  upper phase density for mixed-phase state
C              not used for single phase
C     ZBARQQ=  ionization number
C
C                                TABLE          ANEOS
C     KPAQQ=STATE INDICATOR      =1, 1p    =1, 1p    (eos without melt)
C                                =2, 2p lv =2, 2p liquid/solid plus vapor
C                                          =4, 1p solid  (eos with melt)
C                                          =5, 2p melt   (eos with melt)
C                                          =6, 1p liquid (eos with melt)
C                                =-1 bad value of temperature
C                                =-2 bad value of density
C                                =-3 bad value of material number
C
C                         Note: some versions stop for bad inputs
C                         instead of returning negative values of kpa.
C
C************************************** 6/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TQQ(NUMQQ),RHOQQ(NUMQQ),MATQQ(NUMQQ)
      DIMENSION PQQ(NUMQQ),EQQ(NUMQQ),SQQ(NUMQQ),CVQQ(NUMQQ)
      DIMENSION DPDTQQ(NUMQQ),DPDRQQ(NUMQQ),FKROQQ(NUMQQ),CSQQ(NUMQQ)
      DIMENSION KPAQQ(NUMQQ),R2PLQQ(NUMQQ),R2PHQQ(NUMQQ),ZBARQQ(NUMQQ)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /MIXCON/  RHOMPH(MAXMAT),RHOMPL(MAXMAT),IPEOS
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /MELTTRANS/ TMC(1000*MAXMAT),RSMC(1000*MAXMAT)
     1 ,RLMC(1000*MAXMAT),NUMTMC(MAXMAT),MELT_CURVE_TABLE
      COMMON /HPPTRANS/ THPP(1000*MAXMAT)
     2 ,RLPP(1000*MAXMAT),RHPP(1000*MAXMAT),NUMHPP(MAXMAT),NUMPHASE
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
      COMMON /FNES/ CMLT7,ACK46
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
C     next three commons are here for multitask only
      COMMON /ANZB/ ZB(92), ZRAT
C                    inputs
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
C        outputs
     &  ,PE,EE,SE,CVE,DPTE,DPRE
C        inputs
     &  ,NMATSX,IIZX
      COMMON /ANZLTF/ ZLTOMF
      COMMON /TPMEOS/ NPORFE(MAXMAT)
      SAVE   /FILEOS/,/MIXCON/,/ANES/,/ANE2/,/BNES/,/FNES/,/ANCNTR/,
     & /ANZB/,/ANEEL/,/TPMEOS/
      DIMENSION KOC(MATBUF), PHOFG(MATBUF), NMATSV(MATBUF)
      PARAMETER (ZERO=0.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (FIVE=5.D0)
      PARAMETER (ATHIRD=1.D0/3.D0)
      PARAMETER (QCB1=1.D-20)
      PARAMETER (QCB2=1.D-10)
      PARAMETER (QCB4=-100.00001D0)
      PARAMETER (QCB5=8.D0)
      PARAMETER (QCB6=15.D0)
      PARAMETER (QCB7=24.D0)
      PARAMETER (QCB8=0.2D0)
      PARAMETER (CSFLAG=-77.66D11)
      PARAMETER (TQQMN=0.001D0)
      PARAMETER (RQQMN=1.D-22)
      CHARACTER*50 MESS
      DATA L8BAD /0/

C
      DO LMAT=1,NUMQQ
         KOC(LMAT)=-1
         PHOFG(LMAT)=ZERO
         NMATSV(LMAT)=-77
         IF(MATQQ(LMAT).GT.0) THEN
            KOC(LMAT)=LOCSV(MATQQ(LMAT))
            CSQQ(LMAT)=CSFLAG
            NMATSV(LMAT)=NINT(ACK(KOC(LMAT)+30))
            PHOFG(LMAT)=ACK(KOC(LMAT)+22)
            KPAQQ(LMAT)=1
         END IF
      END DO
      KUSER=0
      KIDG=0
      KFSQT=0
      KKSOL=0
C
      DO 410 LMAT=1,NUMQQ
      IF(MATQQ(LMAT).LE.0) GO TO 410
      IF((TQQ(LMAT).LE.ZERO).OR.(RHOQQ(LMAT).LE.ZERO)) GO TO 480
      MATKTR=MATQQ(LMAT)
      LNKNTR(MATKTR,2)=LNKNTR(MATKTR,2)+1
C
C     check for special forms that are not processed in this loop
C
      IF(NMATSV(LMAT).LE.-1) THEN
        IF(NMATSV(LMAT).EQ.-1) THEN
          KKSOL=1
        ELSEIF(NMATSV(LMAT).EQ.-4) THEN
          KUSER=1
          LNKNTR(MATKTR,15)=LNKNTR(MATKTR,15)+1
        ELSE
          KIDG=1
        END IF
        GO TO 410
      END IF
      IF(KFSQT.EQ.0) THEN
        DO 337 KMAT=1,NUMQQ
        SQTS(KMAT)=MAX(QCB2,TQQ(KMAT))
        SQTS(KMAT)=SQRT(SQTS(KMAT))
  337   CONTINUE
        KFSQT=1
      END IF
      IPSQTS=LMAT
      LOC=LOCSV(MATKTR)
      CMLT7=ACK(LOC+71)
      ACK46=ACK(LOC+46)
      FKROS=1.D10
      CRHS2P=RHOQQ(LMAT)
      CRLS2P=RHOQQ(LMAT)
    5 IF(NMATSV(LMAT)-2) 30,80,6
C
C     TEST FOR POROUS MATERIAL (OLD MODEL)
C
    6 IF(NPORFE(MATKTR)) 7,10,7
    7 IF((RHOQQ(LMAT).GE.ACK(LOC+19)).OR.
     &     (TQQ(LMAT).GE.TCT(MATKTR))) GO TO 30
      GO TO 20
C
C     CHECK FOR LIQUID-VAPOR OR SOLID-VAPOR STATE
C
   10 IF((RHOQQ(LMAT).GE.ACK(LOC+47)).OR.          ! If density is greater than rho_lm
     &    (TQQ(LMAT).GE.TCT(MATKTR))) GO TO 30     ! or temperature is above critical point then not two-phase
      IF (TQQ(LMAT).GT.ACK(LOC+18)) GO TO 20       ! If temperature is greater than melt temp., check for mixed phase
      IF (RHOQQ(LMAT).GE.ACK(LOC+23)) GO TO 30     ! If density is greater than minimum density skip two-phase test
   20 CALL ANE2PH (TQQ(LMAT),RHOQQ(LMAT),MATKTR,PQQ(LMAT),EQQ(LMAT),
     &  SQQ(LMAT),CVQQ(LMAT),DPDTQQ(LMAT),DPDRQQ(LMAT),LOC,KPAQQ(LMAT))
      LNKNTR(MATKTR,8)=LNKNTR(MATKTR,8)+1
C
C     KPA=2 IF LIQUID-VAPOR OR SOLID-VAPOR STATE
      IF (KPAQQ(LMAT).EQ.2) GO TO 100
C
C     IS MELT TRANSITION INCLUDED
C
 30   CONTINUE

      ! New code for finding melt-transition state via tables
      IF (MELT_CURVE_TABLE.EQ.1) THEN

         ! Check for LIQUID-SOLID STATE
         CALL ANE2PH_MELT (TQQ(LMAT),RHOQQ(LMAT),MATKTR,PQQ(LMAT),
     $        EQQ(LMAT),SQQ(LMAT),CVQQ(LMAT),DPDTQQ(LMAT),
     $        DPDRQQ(LMAT),LOC,KPAQQ(LMAT))

         IF (KPAQQ(LMAT).EQ.5) GOTO 100  ! Mixed phase

         IF (KPAQQ(LMAT).EQ.4) GOTO 50  ! Solid phase

         IF (KPAQQ(LMAT).EQ.6) GOTO 90  ! Liquid phase
         
         GOTO 80  ! single phase (must be vapour?)

      ELSEIF (NUMPHASE.GT.1 .AND. MATKTR.EQ.1) THEN  ! New LPP/HPP method not yet compatible with new melt transition method.

         ! Check for LPP-HPP STATE
         CALL ANE2PH_LPPHPP (TQQ(LMAT),RHOQQ(LMAT),MATKTR,PQQ(LMAT),
     $        EQQ(LMAT),SQQ(LMAT),CVQQ(LMAT),DPDTQQ(LMAT),
     $        DPDRQQ(LMAT),LOC,KPAQQ(LMAT))
         
         GOTO 100   ! This routine defines the state explicitly

      END IF

      IF (ACK46) 80,80,40
   40 KPAQQ(LMAT)=0
C     FATAL ERROR FLAG SET TO STOP IN ANELSM
      CALL ANELSM (TQQ(LMAT),RHOQQ(LMAT),DPDTQQ(LMAT),DPDRQQ(LMAT),LOC,
     &   KPAQQ(LMAT))
      LNKNTR(MATKTR,9)=LNKNTR(MATKTR,9)+1
      IF (KPAQQ(LMAT)-2) 50,70,60
C
C     SOLID STATE (EOS WITH MELT)
C
   50 KPAQQ(LMAT)=4
      CMLT7=-ONE  ! Turn off melt terms
      CALL ANEOS1 (TQQ(LMAT),RHOQQ(LMAT),PQQ(LMAT),EQQ(LMAT),SQQ(LMAT),
     &  CVQQ(LMAT),DPDTQQ(LMAT),DPDRQQ(LMAT),LOC)
      LNKNTR(MATKTR,7)=LNKNTR(MATKTR,7)+1
      CMLT7=ZERO  ! Turn back on melt terms
      GO TO 100
C
C     LIQUID STATE (EOS WITH MELT)
C
   60 KPAQQ(LMAT)=6
      GO TO 90
C
C     LIQUID-SOLID STATE (EOS WITH MELT)
C
   70 KPAQQ(LMAT)=5
      PQQ(LMAT)=PM
      EQQ(LMAT)=EM
      SQQ(LMAT)=SM
      CVQQ(LMAT)=CVM
      DPDTQQ(LMAT)=DPDTM
      DPDRQQ(LMAT)=DPDRM
      GO TO 100
C
C     ONE-PHASE STATE (EOS WITHOUT MELT)
C
   80 KPAQQ(LMAT)=1
   90 CALL ANEOS1 (TQQ(LMAT),RHOQQ(LMAT),PQQ(LMAT),EQQ(LMAT),SQQ(LMAT),
     &   CVQQ(LMAT),DPDTQQ(LMAT),DPDRQQ(LMAT),LOC)
      LNKNTR(MATKTR,7)=LNKNTR(MATKTR,7)+1
      GO TO 100
C
C-----------------------------------------------------------------------
C     error in call parameters
C
  480 PQQ(LMAT)=ZERO
      EQQ(LMAT)=ZERO
      SQQ(LMAT)=ZERO
      CVQQ(LMAT)=ONE
      PHOFG(LMAT)=ZERO
      DPDTQQ(LMAT)=ZERO
      DPDRQQ(LMAT)=ZERO
      FKROQQ(LMAT)=ONE
      R2PHQQ(LMAT)=ZERO
      R2PLQQ(LMAT)=ZERO
      ZBARQQ(LMAT)=ZERO
      CSQQ(LMAT)=ZERO
      IF(TQQ(LMAT).LE.ZERO) THEN
        KPAQQ(LMAT)=-1
      ELSE
        KPAQQ(LMAT)=-2
      END IF
      IF(L8BAD.EQ.0) THEN
        L8BAD=1
        CALL ANMARK
     1  ('ZERO OR NEGATIVE DENSITY OR TEMPERATURE IN ANEOS')
        MESS=
     1   'MATERIAL12345  T=1234567890123   RHO=1234567890123'
        WRITE(MESS(9:13),'(I5)',ERR=694) MATQQ(LMAT)
        WRITE(MESS(18:30),'(1PE13.6)',ERR=694) TQQ(LMAT)
        WRITE(MESS(38:50),'(1PE13.6)',ERR=694) RHOQQ(LMAT)
  694   CONTINUE
        CALL ANMARK(MESS)
      END IF
C     stop code for bad input to eos package
      STOP
C
C----------------------------------------------
C     set arrays for return
C
  100 FKROQQ(LMAT)=FKROS
      R2PHQQ(LMAT)=CRHS2P
      R2PLQQ(LMAT)=CRLS2P
      ZBARQQ(LMAT)=ZBAR
C
  410 CONTINUE
C
      IF(KIDG.EQ.1) THEN
C
C       IDEAL GAS EOS
C
        DO 422 LMAT=1,NUMQQ
        IF(NMATSV(LMAT).EQ.-2) THEN
          CVQQ(LMAT)=ACK(KOC(LMAT)+4)
          EQQ(LMAT)=CVQQ(LMAT)*TQQ(LMAT)
          DPDRQQ(LMAT)=ACK(KOC(LMAT)+5)*TQQ(LMAT)
          DPDTQQ(LMAT)=ACK(KOC(LMAT)+5)*RHOQQ(LMAT)
          PQQ(LMAT)=DPDTQQ(LMAT)*TQQ(LMAT)
          SQQ(LMAT)=CVQQ(LMAT)*(LOG(TQQ(LMAT))-ACK(KOC(LMAT)+15)*
     &              LOG(RHOQQ(LMAT))-ACK(KOC(LMAT)+6))
          FKROQQ(LMAT)=QCB8+ACK(KOC(LMAT)+20)/(ACK(KOC(LMAT)+25)+
     &              TQQ(LMAT))**3
          R2PHQQ(LMAT)=ZERO
          R2PLQQ(LMAT)=ZERO
          ZBARQQ(LMAT)=ZERO
        END IF
  422   CONTINUE
        CRHS2P=ZERO
        CRLS2P=ZERO
        ZBAR=ZERO
      END IF
C
      IF(KUSER.EQ.1) THEN
C
C       USER DEFINED EOS
C
        DO 424 LMAT=1,NUMQQ
        IF(NMATSV(LMAT).EQ.-4) THEN
          CALL ANUEOS(
C         inputs
     &    TQQ(LMAT),RHOQQ(LMAT),MATQQ(LMAT),ACK(KOC(LMAT)+1),
C         outputs
     &    PQQ(LMAT),EQQ(LMAT),SQQ(LMAT),CVQQ(LMAT),DPDTQQ(LMAT),
     &    DPDRQQ(LMAT),FKROQQ(LMAT),KPAQQ(LMAT),
     &    R2PLQQ(LMAT),R2PHQQ(LMAT),ZBARQQ(LMAT))
C
        END IF
  424   CONTINUE
        CRHS2P=R2PHQQ(1)
        CRLS2P=R2PLQQ(1)
        ZBAR=ZBARQQ(1)
        FKROS=FKROQQ(1)
        RHOMPH(IPEOS)=CRHS2P
        RHOMPL(IPEOS)=CRLS2P
        GO TO 428
      END IF
      IF(KKSOL.EQ.1) THEN
C
C     type -1 (Low temperature solid)
C
        CALL ANN1VL (NUMQQ,TQQ,RHOQQ,KOC,NMATSV,
     &  PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,FKROQQ,KPAQQ,R2PLQQ,R2PHQQ,
     &  ZBARQQ)
        CRHS2P=R2PHQQ(1)
        CRLS2P=R2PLQQ(1)
        ZBAR=ZBARQQ(1)
        FKROS=FKROQQ(1)
        RHOMPH(IPEOS)=CRHS2P
        RHOMPL(IPEOS)=CRLS2P
      END IF
C
C     phonon conduction term
C
      IF(NOTRAD.NE.0) THEN
        DO 320 LMAT=1,NUMQQ
        IF(PHOFG(LMAT).GT.ZERO) THEN
          PHOFG(LMAT)=ACK(KOC(LMAT)+22)*
     &     TQQ(LMAT)**(THREE-ACK(KOC(LMAT)+41))/RHOQQ(LMAT)
          FKROQQ(LMAT)=FKROQQ(LMAT)*PHOFG(LMAT)/
     &     (FKROQQ(LMAT)+PHOFG(LMAT))
        END IF
  320   CONTINUE
      END IF
C
C     old entry variables
C
      IF(NUMQQ.EQ.1) THEN
        RHOMPH(IPEOS)=CRHS2P
        RHOMPL(IPEOS)=CRLS2P
        FKROS=FKROQQ(1)
      END IF
C
C     sound speed, energy-entropy shifts
C
  428   DO 430 LMAT=1,NUMQQ
        IF(CSQQ(LMAT).EQ.CSFLAG) THEN
          CSQQ(LMAT)=SQRT(MAX(QCB1,(DPDRQQ(LMAT)
     &    +(TQQ(LMAT)*DPDTQQ(LMAT)**2)/(CVQQ(LMAT)*RHOQQ(LMAT)**2))))
          EQQ(LMAT)=EQQ(LMAT)+ACK(KOC(LMAT)+64)
          SQQ(LMAT)=SQQ(LMAT)+ACK(KOC(LMAT)+65)
        END IF
  430   CONTINUE
C
      RETURN
      END
