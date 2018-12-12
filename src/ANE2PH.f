C
C
      SUBROUTINE ANE2PH (T,R,MAT,P,E,S,CV,DPDT,DPDR,LOC,KPA)
C
C***********************************************************************
C
C     ANEOS PACKAGE      EVALUATES THERMODYNAMIC FUNCTIONS
C                        IN THE LIQUID-VAPOR AND SOLID-VAPOR REGIONS
C
C************************************* 10/88 version slt ***************
C
C      This routine uses stored tables of density and temperature of the
C      two phase region to evaluate the quantities of each. 
C        
C      Tables are:  RSOL = solid/melt density on phase curve
C                   RVAP = vapor density on phase curve
C                   TTWO = temperature on phase curve
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
      COMMON /FNES/ CMLT7,ACK46
      COMMON /MIX/FM1,FM2
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE /FILEOS/,/ANES/,/ANE2/,/BNES/,/FNES/,/ANCNTR/,/MIX/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (ATHIRD=1.D0/3.D0)
      PARAMETER (SLP=1.03D0)
      PARAMETER (QCH1=0.995D0)
      PARAMETER (QCH2=1.D0-QCH1)
C
C     find temperature interval in stored 2 phase tables
C
      K1=LOCKP(MAT)   ! Pointer to first point (highest T) in table for material
      KJ=LOCKPL(MAT)  ! Pointer to last point (lowest T) in table for material
      ! Bi-sector splitting search for temperature T in table TTWO
   10   KK=(K1+KJ)/2
        IF(K1.NE.KK) THEN
          IF(T.GT.TTWO(KK)) THEN
            KJ=KK
          ELSE
            K1=KK
          END IF
          GO TO 10
        END IF
      TL=TTWO(KJ)  ! Table temperature below
      TU=TTWO(KK)  ! Table temperature above
      IF (KK.GT.LOCKP(MAT)) GO TO 30  ! Upper table point is not the critical point

      ! Interpolation in the region just below the critical point
      X1=((TU-T)/(TU-TL))**ATHIRD  ! Normalized temperature within table cell
      R1=RSOL(KK)+(RSOL(KJ)-RSOL(KK))*X1  ! Interpolate density of solid
      IF (R.GE.SLP*R1) GO TO 40 ! Density is above incipient vaporization (one phase)
      R2=RVAP(KK)-(RVAP(KK)-RVAP(KJ))*X1  ! Interpolate density of vapor
      IF (R.LE.R2) GO TO 40     ! Density is below complete vaporization (one phase)

      ! Definitely a two-phase mixture, compute the density at each phase boundary
      R1P=(RSOL(KK)-RSOL(KJ))*X1/(THREE*(TU-T))
      R2P=(RVAP(KK)-RVAP(KJ))*X1/(THREE*(TU-T))
      GO TO 50  

      ! Interpolation away from the critical point
   30 DT=TU-TL
      R1=((T-TL)*RSOL(KK)+(TU-T)*RSOL(KJ))/DT  !Interpolate R on incipient vapor curve at T
      IF (R.GE.SLP*R1) GO TO 40
      R2=((T-TL)*RVAP(KK)+(TU-T)*RVAP(KJ))/DT  !Interpolate R on complete vapor curve at T
      IF (R.GT.R2) GO TO 50

      ! Not in 2 phase region so exit
   40 KPA=1
      GO TO 180

      ! Definitely a two-phase region so now compute the properties. . .
   50 KPA=2
      
      ! get gas phase properties, density R2 (complete vapor)
      CALL ANEOS1 (T,R2,P2,E2,S2,CV2,DPDT2,DPDR2,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FK2=FKROS
      CRHS2P=ZBAR
      IF (R.LE.R1) GO TO 60  ! Density less than incipient vaporisation

      ! We get here if density is between R1 (incipient vapor) and SLP*R1
      ! In this case compute the pressure
      CALL ANEOS1 (T,R,P,E,S,CV,DPDT,DPDR,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      IF (P.GE.P2) GO TO 40  ! If pressure at R is greater than complete vapor pressure, single phase & exit
      ! Compute the pressure at R1 (incipient vapor)
      CALL ANEOS1 (T,R1,P1,E1,S1,CV1,DPDT1,DPDR1,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      P=P2
      DPDR=ZERO
      DPDT=(S2-S1)*(R1*R2)/(R1-R2)
      R1=R  ! Two phase with R as the incipient vapor density
      GO TO 170 !Exit

      ! see if we are also in melt region
   60 IF (ACK46) 110,110,70          ! Check density of solid at triple point is postive
   70 IF (TU-ACK(LOC+18)) 80,90,110  ! Above, at, or below the melt temperature

      ! Temperature is below triple point temperature
   80 CMLT7=-ONE  ! Turn off special liquid terms
      CALL ANEOS1 (T,R1,P1,E1,S1,CV1,DPDT1,DPDR1,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FK1=FKROS
      CMLT7=ZERO  ! Turn back on special liquid terms
      GO TO 120

      ! Temperature is at triple point temperature
   90 K9=0
C     FATAL FLAG SET TO STOP IN ANELSM
      CALL ANELSM(T,R1,X1,X2,LOC,K9)
      LNKNTR(MATKTR,10)=LNKNTR(MATKTR,10)+1
      IF(K9-2) 80,100,110
  100 P1=PM
      E1=EM
      S1=SM
      CV1=CVM
      DPDT1=DPDTM
      DPDR1=DPDRM
      FK1=FKROS
      GO TO 120

      ! Temperature is above the triple point temperature
      ! get liquid phase properties, density R1
  110 CALL ANEOS1 (T,R1,P1,E1,S1,CV1,DPDT1,DPDR1,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FK1=FKROS

      ! Get the mixed-phase state
      ! Get derived parameters
  120 X3=R1-R2                   !use lever rule and input density R
      X1=(R1-R)/X3               !along with gas density R2, liguid density R1
      X2=(R-R2)/X3               !to get proportion of liquid, gas in mixture
      FM1=R1*X2/R                !mass mixing ratio of liquid
      FM2=R2*X1/R                !mass mixing ratio of gas
      E=FM1*E1+FM2*E2            !internal energy of mixture
      S=FM1*S1+FM2*S2            !entropy of mixture
      FKROS=FM1*FK1+FM2*FK2      !Rosseland mean opacity of mixture
      ZBAR=FM1*CRHS2P+FM2*ZBAR
      IF (P1.LE.P2) GO TO 130
      DPDR=QCH1*R1+QCH2*R2       !DPDR of mixture
      IF (R.LE.DPDR) GO TO 130
      X4=R1-DPDR
      P=(P1*(R-DPDR)+P2*(R1-R))/X4
      DPDR=(P1-P2)/X4            !DPDR of mixture
      GO TO 140
  130 P=P2
      DPDR=ZERO                  !DPDR of mixture
 140  DPDT=(S2-S1)*R1*R2/X3     !DPDT of mixture
      IF (KK.NE.LOCKP(MAT)) THEN
        X4=(RVAP(KK)-RVAP(KJ))/DT
        X5=(RSOL(KK)-RSOL(KJ))/DT
      ELSE
        X4=R2P
        X5=R1P
      END IF
      X3=-(R1*X1*X4+R2*X2*X5)/(R*X3)
      X1=CV1+(P1-T*DPDT1)*X5/R1**2
      X2=CV2+(P2-T*DPDT2)*X4/R2**2   
      CV=X3*(E1-E2)+FM1*X1+FM2*X2    !Heat capacity of mixture
  170 CRHS2P=R1
      CRLS2P=R2


 180  CONTINUE

      RETURN
      END

C
C
      SUBROUTINE ANE2PH_MELT (T,RHO,MAT,P,E,S,CV,DPDT,DPDR,LOC,KPA)
C
C***********************************************************************
C
C     ANEOS PACKAGE      EVALUATES THERMODYNAMIC FUNCTIONS
C                        IN THE SOLID-LIQUID REGION
C
C************************************* 10/88 version slt ***************
C
C      This routine uses stored tables of density and temperature of the
C      two phase region to evaluate the quantities of each. 
C        
C      Tables are:  RSMC = solid density on solidus
C                   RLMC = liquid density on liquidus
C                   TMC  = temperature on solidus/liquidus
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /MELTTRANS/ TMC(1000*MAXMAT),RSMC(1000*MAXMAT)
     1 ,RLMC(1000*MAXMAT),NUMTMC(MAXMAT),MELT_CURVE_TABLE
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
      COMMON /FNES/ CMLT7,ACK46
      COMMON /MIX/FM1,FM2
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE /FILEOS/,/ANES/,/ANE2/,/BNES/,/FNES/,/ANCNTR/,/MIX/
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (ATHIRD=1.D0/3.D0)
      PARAMETER (SLP=1.03D0)
      PARAMETER (QCH1=0.995D0)
      PARAMETER (QCH2=1.D0-QCH1)
      
      ! Simple tests to identify phase (duplicated from ANELSM)
      IF (T.LE.ACK(LOC+18)) THEN ! T =< T_m0
         IF (RHO.GE.ACK46) THEN ! If density is above RHO_SM  (and T =< Tm0) then SOLID
            KPA = 4
            RETURN
         ELSEIF (RHO.LT.ACK(LOC+23)) THEN ! If density is below RHO_MIN (and T =< Tm0) then LIQUID(?)
            KPA = 6
            RETURN
         ELSEIF (T.LT.ACK(LOC+49)) THEN ! If temperature is less than minimum on solidus, then SOLID
            KPA = 4
            RETURN
         END IF
      ELSE  ! T > T_m0
         IF((RHO.LE.ACK(LOC+47)).OR.(T.GE.ACK(LOC+48))) THEN ! If RHO < RHO_LM or T > T_max => LIQUID
            KPA = 6
            RETURN
         END IF
      END IF

      ! Now check to see if within the LIQUID-SOLID two-phase region
      K1=1+1000*(MAT-1)        ! Index for lowest temperature in table
      KJ=K1+NUMTMC(MAT)-1      ! Index for highest temperature in table

      IF (T.LE.TMC(K1)) THEN         ! Temperature is below lower limit of table => SOLID
         KPA = 4
         RETURN
      ELSEIF (T.GE.TMC(KJ)) THEN     ! Temperature is above upper limit of table => LIQUID
         KPA = 6
         RETURN
      ELSEIF (RHO.LE.RLMC(K1)) THEN  ! Density is less than lower limit of liquidus => LIQUID 
         KPA = 6
         RETURN
      ELSEIF (RHO.GE.RSMC(KJ)) THEN  ! Density is more than upper limit of solidus (and T<Tmax) => SOLID 
         KPA = 4
         RETURN
      END IF

      ! Bi-sector splitting search for temperature T in table TMC
      KK=(K1+KJ)/2           ! Mid-point
      DO WHILE (K1.NE.KK)
         IF(T.GT.TMC(KK)) THEN
            K1=KK
         ELSE
            KJ=KK
         END IF
         KK=(K1+KJ)/2           ! Mid-point
      END DO
      TL=TMC(KK)  ! Table temperature below
      TU=TMC(KJ)  ! Table temperature above      

      ! Now the interpolation
      DT=TU-TL
      RS=((T-TL)*RSMC(KJ)+(TU-T)*RSMC(KK))/DT  !Interpolate R on solidus at T
      RL=((T-TL)*RLMC(KJ)+(TU-T)*RLMC(KK))/DT  !Interpolate R on liquidus at T

      ! Single phase solid/liquid
      IF (RHO.GE.RS) THEN         !SOLID
         KPA=4
         RETURN
      ELSEIF (RHO.LE.RL) THEN     !LIQUID
         KPA=6
         RETURN
      ENDIF

      ! Material is in the two-phase melt region
      KPA=5

      ! Get the solid-phase state
      CMLT7=-ONE                ! Turn off special liquid terms
      CALL ANEOS1 (T,RS,PS,ES,SS,CVS,DPDTS,DPDRS,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FKS=FKROS
      CMLT7=ZERO                ! Turn back on special liquid terms      

      ! Get the liquid-phase state
      CALL ANEOS1 (T,RL,PL,EL,SL,CVL,DPDTL,DPDRL,LOC)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FKL=FKROS

      ! Compute mixed-phase state here and store in BNES named common
  310 X2=RHO*(RS-RL)
      DRS=RS*(RHO-RL)/X2
      DRL=RL*(RS-RHO)/X2
      DPDT=(SL-SS)*((RS*RL)/(RS-RL)) ! DPDT of mixture
      DRLDT=(DPDTM-DPDTL)/DPDRL
      DRSDT=(DPDTM-DPDTS)/DPDRS
      X1=-RHO*(RL*(RHO-RL)*DRSDT+RS*(RS-RHO)*DRLDT)/X2**2
      E=DRS*ES+DRL*EL  ! Energy of mixture
      S=DRS*SS+DRL*SL  ! Entropy of mixture
      FKROS=DRS*FKS+DRL*FKL 
      CV=X1*(ES-EL)+DRS*(CVS+(PS-T*DPDTS)*DRSDT/RS**2)
     1    +DRL*(CVL+(PL-T*DPDTL)*DRLDT/RL**2) ! Heat capacity of mixture
      DPDR=ZERO        ! DPDR of mixture
      P=DRS*PS+DRL*PL  ! Pressure of mixture
      CRHS2P=RS
      CRLS2P=RL

      RETURN
      END

      SUBROUTINE ANE2PH_LPPHPP (T,RHO,MAT,P,E,S,CV,DPDT,DPDR,LOC,KPA)
C
C***********************************************************************
C
C     ANEOS PACKAGE      EVALUATES THERMODYNAMIC FUNCTIONS
C                        IN THE LPP-HPP REGIONS
C
C************************************* 10/88 version slt ***************
C
C      This routine uses stored tables of density and temperature of the
C      two phase region to evaluate the quantities of each. 
C        
C      Tables are:  RLPP = density on the low-pressure phase boundary
C                   RHPP = density on the high-pressure phase boundary
C                   THPP = temperature on solid-solid phase boundary
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /HPPTRANS/ THPP(1000*MAXMAT)
     2 ,RLPP(1000*MAXMAT),RHPP(1000*MAXMAT),NUMHPP(MAXMAT),NUMPHASE
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /BNES/ PM,EM,SM,CVM,DPDTM,DPDRM
      COMMON /FNES/ CMLT7,ACK46
      COMMON /MIX/FM1,FM2
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE /FILEOS/,/ANES/,/ANE2/,/BNES/,/FNES/,/ANCNTR/,/MIX/
      PARAMETER (ZERO=0.D0)
      
      ! Define the low(est) pressure phase material and ACK location
      MAT_LPP = MAT
      LOC_LPP = LOC     ! LPP is the current material
      KPA_LPP = 41

 100  K1=1000*(MAT_LPP-1)+1     ! Index for lowest temperature in table
      KJ=NUMHPP(MAT_LPP)      ! Index for highest temperature in table

      ! Define the ACK location of the high-pressure phase
      LOC_HPP=LOC_LPP+99  ! HPP is the next material

      !write(*,*)K1,KJ,MAT_LPP
      !write(*,*)THPP(K1),THPP(KJ)
      !write(*,*)RLPP(K1),RHPP(KJ)
      
      IF (T.LT.THPP(K1)) THEN   ! Temperature is below lower limit of table => LPP
         KPA=KPA_LPP
         CALL ANEOS1 (T,RHO,P,E,S,CV,DPDT,DPDR,LOC_LPP)
         LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
         RETURN
      ENDIF
c$$$      ELSEIF (T.GT.THPP(KJ)) THEN   ! Temperature is above upper limit of table => LPP
c$$$         KPA=46
c$$$         CALL ANEOS1 (T,RHO,P,E,S,CV,DPDT,DPDR,LOC_HPP)
c$$$         LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
c$$$         RETURN
c$$$      END IF
c$$$c$$$      ELSEIF (RHO.LE.RLPP(K1)) THEN  ! Density is less than lower limit of LPP => LPP 
c$$$         KPA = 6
c$$$         RETURN
c$$$      ELSEIF (RHO.GE.RHPP(KJ)) THEN  ! Density is more than upper limit of HPP (and T<Tmax) => HPP 
c$$$         KPA = 4
c$$$         RETURN
c$$$      END IF

      ! Bi-sector splitting search for temperature T in table THPP
      KK=(K1+KJ)/2              ! Mid-point
      DO WHILE (K1.NE.KK)
         IF(T.GT.THPP(KK)) THEN
            K1=KK
         ELSE
            KJ=KK
         END IF
         KK=(K1+KJ)/2           ! Mid-point
      END DO
      TL=THPP(KK)  ! Table temperature below
      TU=THPP(KJ)  ! Table temperature above      

      ! Now the interpolation
      DT=TU-TL
      RS=((T-TL)*RHPP(KJ)+(TU-T)*RHPP(KK))/DT  !Interpolate R on HPP at T
      RL=((T-TL)*RLPP(KJ)+(TU-T)*RLPP(KK))/DT  !Interpolate R on LPP at T

      !write(*,*)KK,KJ,TL,T,TU,RL,RHO,RS
      !read(*,*)
      !write(KLST,*)RS,RL

      ! Single LPP/HPP 
      IF (RHO.GE.RS) THEN         !HPP
         ! Phase must be one of the high-pressure phases
         ! If there are more than one high-pressure phases consider
         ! the next phase boundary until the highest phase is reached
         IF (MAT_LPP .LT. NUMPHASE-1) THEN
            MAT_LPP = MAT_LPP + 1
            LOC_LPP = LOC_LPP + 99
            KPA_LPP = KPA_LPP + 3
            GOTO 100
         END IF
         ! If phase is the highest available, use this.
         KPA=KPA_LPP+2
         CALL ANEOS1 (T,RHO,P,E,S,CV,DPDT,DPDR,LOC_HPP)
         E=E+1.2D10  !+1.D10
         S=S+9.07D9  !+17.7D9
         LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
         !write(*,*)"Leaving ANE2PH_LPPHPP as HPP"
         RETURN
      ELSEIF (RHO.LE.RL) THEN     !LPP
         KPA=KPA_LPP
         CALL ANEOS1 (T,RHO,P,E,S,CV,DPDT,DPDR,LOC_LPP)
         LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
         !write(*,*)"Leaving ANE2PH_LPPHPP as LPP"
         RETURN
      ENDIF

      ! Material is in the mixed-phase region
      KPA=KPA_LPP+1

      ! Get the solid-phase state
      CALL ANEOS1 (T,RS,PS,ES,SS,CVS,DPDTS,DPDRS,LOC_HPP)
      ES=ES+1.2D10  !+1.D10
      SS=SS+9.07D9  !+17.7D9
      !SS=SS+1.D10/T
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FKS=FKROS

      ! Get the liquid-phase state
      CALL ANEOS1 (T,RL,PL,EL,SL,CVL,DPDTL,DPDRL,LOC_LPP)
      LNKNTR(MATKTR,12)=LNKNTR(MATKTR,12)+1
      FKL=FKROS

      ! Compute mixed-phase state here and store in BNES named common
  310 X2=RHO*(RS-RL)
      DRS=RS*(RHO-RL)/X2
      DRL=RL*(RS-RHO)/X2
      DPDT=(SL-SS)*((RS*RL)/(RS-RL)) ! DPDT of mixture
      DRLDT=(DPDTM-DPDTL)/DPDRL
      DRSDT=(DPDTM-DPDTS)/DPDRS
      X1=-RHO*(RL*(RHO-RL)*DRSDT+RS*(RS-RHO)*DRLDT)/X2**2
      E=DRS*ES+DRL*EL  ! Energy of mixture
      S=DRS*SS+DRL*SL  ! Entropy of mixture
      FKROS=DRS*FKS+DRL*FKL 
      CV=X1*(ES-EL)+DRS*(CVS+(PS-T*DPDTS)*DRSDT/RS**2)
     1    +DRL*(CVL+(PL-T*DPDTL)*DRLDT/RL**2) ! Heat capacity of mixture
      DPDR=ZERO        ! DPDR of mixture
      P=DRS*PS+DRL*PL  ! Pressure of mixture
      CRHS2P=RS
      CRLS2P=RL
      
      !write(*,*)"Leaving ANE2PH_LPPHPP as LPP/HPP mix"

      RETURN
      END

