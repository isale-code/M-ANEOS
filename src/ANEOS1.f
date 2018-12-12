C
C
      SUBROUTINE ANEOS1 (T,RHO,P,E,S,CV,DPDT,DPDR,L)
C
C**********************************************************************
C
C     ANEOS PACKAGE     NUCLEAR AND COLD COMPONENTS
C
C     INPUTS (ALL UNITS CGS-EV)
C
C     T    =  TEMPERATURE
C     RHO  =  DENSITY
C     L    =  MATERIAL POINTER
C
C     For speed reasons there is one input variable in /ANESQT/
C     SQTS(ipsqts)=sqrt(t). This should be set in calling routine.
C
C     OUTPUTS (ALSO SEE /ANE2/)
C
C     P    =  PRESSURE
C     E    =  SPECIFIC ENERGY
C     S    =  SPECIFIC ENTROPY
C     CV   =  SPECIFIC HEAT CAPACITY
C     DPDT =  DP / DT
C     DPDR =  DP / DRHO
C
C************************************** 4/30/90 version slt ************
C********** 10/7/99 MODIFIED FOR MOLECULAR CLUSTERS BY HJM *************
C********** 4/27/05 MODIFIED FOR VIBRATIONAL DOF BY HJM ****************
C********** 5/01/06 MODIFIED FOR TRIATOMIC MOLECULES BY HJM ************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
C                    inputs
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
     &  ,PE,EE,SE,CVE,DPTE,DPRE
     &  ,NMATSX,IIZX
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      COMMON /FNES/ CMLT7,ACK46
      COMMON /ANZB/ ZB(92), ZRAT
C     vector buffer size
      PARAMETER (MATBUF=64)
C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
C     following aneos1 variables are for aneos2 edit only
      COMMON /ANEDIS/ GAMMA,PSIZ,THETA
C
C     the following common is used by hjm for debugging only--
C      it is not part of the usual ANEOS package
C
      COMMON /ANGELX/ W,Y,ALFA,BETA,DADT,DADR,DBDT,ZZ,ET,
     & SN,CVN,EN,PN,PM,EM,SM,CVM,EC,PC,RHO0,RHO00,ZVIB,
     & SMLT,CVMLT,EMLT,PMLT,H1,H2
C
      SAVE   /FILEOS/,/ANES/,/ANE2/,/FNES/,/ANCNTR/,
     & /ANZB/,/ANEEL/,/ANEDIS/
      DIMENSION XXA(3)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (ATHIRD=1.D0/3.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (FIVE=5.D0)
      PARAMETER (THALF=3.D0/2.D0)
      PARAMETER (FLGJWL=0.0252525D0)
C      PARAMETER (QCC1=1.D-10)    !SUPPLIED VALUE DOES NOT WORK WELL ON MAC
      PARAMETER (QCC1=1.D-6)      !low density cutoff for solid terms
      PARAMETER (QCC2=0.2D0)
      PARAMETER (QCC3=30.D0)
      PARAMETER (QCC4=20.D0)
      PARAMETER (QCC5=1.D-5)      !min temperature for gas treatment, units of theta
      PARAMETER (QCC6=1.D5)       !psi value at switch to pure gas treatment
      PARAMETER (QCC7=2.5D0)
      PARAMETER (QCC8=0.07D0)     !maximum temperature for Saha computation in eV
      PARAMETER (QCC9=0.5D11)
      PARAMETER (QCC10=0.4D0)
      PARAMETER (QCC11=6.18D7)
      PARAMETER (QCC12=1.41421356D0)
      PARAMETER (QCC13=0.34657359D0)
      PARAMETER (QCC14=416.D0)
      PARAMETER (QCC15=1.D8)      !Rosseland mean opacity for type 0 and 3
      PARAMETER (QCC16=1.D0)
CHJM
      PARAMETER (QCC17=100.D0)        !limit on size of exponential
      PARAMETER (QCC18=27.D0/4.D0)    !27/4, constant for triatomic function
      LOGICAL MOLEC
CHJM
C
      RHO0=ACK(L+11)    !Reference density at reference temperature
      RHO00=ACK(L+19)   !Reference density at zero temperature
      NTYPE=NINT(ACK(L+30))
C
C     test for mixed-phase jwl (high explosive EOS)
C
      IF (ACK(L+25).EQ.FLGJWL) GO TO 806
C
      FT=BOLTS*ACK(L+27)  !Boltzman k * N0
CHJM
CRMC    Be careful to initialize molecular quantities--thanks Robin!
CHJM                                                    11/4/01
      PM=ZERO
      PN=ZERO
      PMLT=ZERO
      PE=ZERO
      EN=ZERO
      EM=ZERO
      EMLT=ZERO
      EE=ZERO
      PC=ZERO
      EC=ZERO
      SN=ZERO
      SM=ZERO
      SE=ZERO
      SMLT=ZERO
      CVN=ZERO
      CVM=ZERO
      CVMLT=ZERO
      CVE=ZERO
CHJM
CHJM  compute auxiliary quantities to describe molecular clusters
CHJM
      ZZ=ONE
      ATOM=ACK(L+86)    !number of atoms in atomic cluster
      MOLEC=(ATOM.EQ.TWO.OR.ATOM.EQ.THREE)  !flag for atomic cluster addition
      IF(MOLEC) THEN                 !only for diatomic or triatomic
        ET=ACK(L+80)/T               !Binding energy/T
        TVT=ACK(L+85)/T              !Vibrational energy level/T
        VDOF=ACK(L+84)               !# of vibrational DOF
        TPWR=ACK(L+81)          !ACK(81) = temperature exponent in Y
        RPWR=ACK(L+87)          !ACK(87) = density exponent in Y
C
C     construct functions for vibrational energy contribution,
C     making sure that exponentials do not get too big
C
        ZVIB=ONE
        DZVIB=ZERO
        D2ZVIB=ZERO
        IF(VDOF.NE.ZERO) THEN
          IF(ET.GT.QCC17) THEN
            EXEB=ONE
            EXXEB=ZERO
          ELSE
            EXEB=ONE-EXP(-ET)
            EXXEB=EXP(-ET)
          ENDIF
          IF(TVT.GT.QCC17) THEN
            EXTV=ONE
            EXXTV=ZERO
          ELSE
            EXTV=ONE-EXP(-TVT)
            EXXTV=EXP(-TVT)
          ENDIF
          ZVIB=(EXEB/EXTV)**VDOF !vibrational partition function
          DZVIB=VDOF*(ET*EXXEB/EXEB-TVT*EXXTV/EXTV) !first derivative in T
          D2ZVIB=VDOF*(ET*(EXXEB-ONE+ET)*EXXEB/EXEB**2-
     &           TVT*(EXXTV-ONE+TVT)*EXXTV/EXTV**2) !second derivative
        ENDIF
C
C      construct function Y, except for exponential part, which requires
C      careful treatment
C
        ARG=ACK(L+79)*(T**TPWR)/(ZVIB*RHO**RPWR ) !ACK(79) = CXX
C
C     construct the crucial Y and W functions,  Y is the
C     inverse of the equilibrium function, W is the fraction
C     of *unbound* molecules
C
        IF(ET.GT.QCC17) THEN
           Y=ZERO
        ELSE
           Y=ARG*EXP(-ET)
        ENDIF
C
C     W for diatomic clusters
C
        IF(ATOM.EQ.TWO) THEN
           AFAC=TWO             !constant in Y - W equation = 2
           SQY=DSQRT(Y)
           IF(Y.LT.1.D2) THEN   !small Y formula
              W=SQY*(DSQRT(TWO+Y)-SQY)
           ELSE
              Z=ONE/Y           !good approximation for big Y
              W=ONE-(HALF-(HALF-(0.625D0-(0.875D0-1.3125*Z)
     &             *Z)*Z)*Z)*Z
           ENDIF
        ENDIF
C
C       W for triatomic clusters
C
        IF(ATOM.EQ.THREE) THEN
           AFAC=QCC18           !constant in Y - W equation = 27/4
           CUBERTY=Y**ATHIRD
           IF(Y.LT.1.D2) THEN
              G=DSQRT(ONE+Y)
              W=THALF*CUBERTY*((G+ONE)**ATHIRD-(G-ONE)**ATHIRD)
           ELSE
              Z=ONE/Y           !good approximation for big Y
              W=ONE+(-1.48148148D-1+(6.58436214D-2+(-3.9018442D-2+
     &             (2.6494D-2-1.948246D-2*Z)*Z)*Z)*Z)*Z
           ENDIF
        ENDIF
C     
C     construct the PSI renormalization function ZZ and its
C     derivatives for the computation of pressure, entropy, etc.
C     The following equations are good for any number of atoms.
C     
        EFAC=TWO*ATHIRD/ATOM
        IF(ET.GT.QCC17) THEN    !very small Y limit
           ZZ=(AFAC*ARG*EXP(RPWR*(ONE-W)))**EFAC
        ELSE                    !normal Y
           ZZ=(W**(TWO*ATHIRD))*EXP(EFAC*(RPWR*(ONE-W)+ET))
        ENDIF
        DYDR=-RPWR
        DYDT=TPWR+ET+DZVIB
        DYDT2=-ET+D2ZVIB
        H1=EFAC*(ONE-W)
        H2=-W/(ONE+RPWR*(ONE-W))
        ALFA=-THALF*H1*DYDR
        BETA=H1*DYDT-EFAC*ET
        DBDT=H1*(H2*DYDT**2+DYDT2)+EFAC*ET
        DADT=-THALF*H1*H2*DYDT*DYDR
        DADR=-THALF*H1*H2*DYDR**2
      ENDIF
CHJM
      IF (NTYPE.EQ.2) THEN
C
C       type 2:  Ideal gas only with electronic terms
C
         DPDR=FT*T
         E=THALF*DPDR
         P=DPDR*RHO
         CV=THALF*FT
         DPDT=RHO*FT
         S=FT*(THALF*LOG(T/ACK(L+13))-LOG(RHO)+QCC7)
         GO TO 170              !now add electronic components
C
      ELSEIF (RHO.LE.QCC1) THEN
C
C       other types at density le 1.e-6
C         --treat as a perfect gas with no cold contribution
C
         PN=RHO*FT*T
         EN=ACK(L+10)+THALF*FT*T
         CVN=THALF*FT
         DPDT=RHO*FT
         DPDR=FT*T
         SN=FT*(THALF*LOG(T/ACK(L+13))-LOG(RHO)+QCC7)
         PSIZ=1.D100            !Report really big PSI for debugging purposes
CHJM
CHJM      modified for molecular clusters
CHJM
         IF(MOLEC) THEN
            PM=-PN*ALFA
            EM=-THALF*FT*T*BETA
            SM=-THALF*FT*(LOG(ZZ)+BETA)
            CVM=-CVN*(BETA+DBDT)
            DPDT=DPDT*(ONE-ALFA-DADT)
            DPDR=DPDR*(ONE-ALFA-DADR)
         ENDIF
C
C       Update thermodynamic quantities
C
         P=PN+PM
         E=EN+EM
         S=SN+SM
         CV=CVN+CVM
CHJM
         GO TO 130              !now add liquid phase components, if requested
      END IF
C
C     all other conditions for types 0,1,3,4
C
      FTT=FT*T                  ! N0 k T
      RHOO=ACK(L+11)            !reference density
      RHOOO=ACK(L+19)           !zero temperature, pressure density
C
C     determine gamma, gamma derivative, theta
C     and interpolation function
C
      IF (RHO.LE.RHOO) THEN     !below reference density
         X2=RHO*ACK(L+16)
         X3=ACK(L+61)+ONE
         GAMMA=RHO*(ACK(L+17)+X2)+X3
         GAMP=ACK(L+17)+TWO*X2
         THETA=ACK(L+14)*EXP(RHO*(ACK(L+17)+HALF*X2))*RHO**X3
      ELSE                      !above reference density
         X3=RHOO/RHO
         X4=ONE-X3
         X6=ONE-ACK(L+60)
         X2=ACK(L+60)*ACK(L+15)
         X1=ACK(L+24)-X2
         GAMMA=X6*X3*ACK(L+15)+X2+X1*X4**2
         GAMP=X3*(-X6*ACK(L+15)+TWO*X1*X4)/RHO
         THETA=ACK(L+25)*EXP(X6*X4*ACK(L+15)
     &        -HALF*X1*(THREE-X3*(FOUR-X3)))*(RHO/RHOO)**ACK(L+24)
      END IF
C     
      X1=RHO**ATHIRD
      X2=RHO/RHOOO
      X3=X2**ATHIRD
      X4=X2/X3
      X6=ONE/X3
      PSI=ACK(L+13)*T*(X1/THETA)**2 !original definition of PSI
CHJM
C      modification for molecular clusters
C
      PSIZ=PSI*ZZ
CHJM
C
C     cold compression terms
C
      IF (X2.LE.ONE) THEN       !density is less than cold reference density
C
C     computation for expanded states using modified morse potential
C
         IF(ACK(L+82).EQ.ZERO) THEN
            X5=ONE-X6
            X7=EXP(ACK(L+5)*X5)
            X8=EXP(ACK(L+6)*X5)
            PC=ACK(L+4)*(X7-X8)*X4
            DPDR=PC/(THALF*RHO)+ACK(L+4)*(ACK(L+5)*X7-ACK(L+6)*X8)
     &           /(THREE*X4*RHOOO)
            EC=THREE*ACK(L+4)*
     &           ((X7-ONE)/ACK(L+5)-(X8-ONE)/ACK(L+6))/RHOOO
CHJM
C       alternative expanded eos using lennard-jones potential
C
         ELSE
            PC=ACK(L+4)*(X2**ACK(L+5)-X2**ACK(L+6))
            X7=ACK(L+5)-ONE
            X8=ACK(L+6)-ONE
            EC=ACK(L+4)*((X2**X7-ONE)/X7-(X2**X8-ONE)/X8)/RHOOO
            DPDR=ACK(L+4)*(ACK(L+5)*X2**X7-ACK(L+6)*X2**X8)/RHOOO
         ENDIF
CHJM
C
C        modify pressure to change shape of cold compression curve
C        and move critical point, if requested
C
         IF ((ACK(L+53).NE.ZERO).AND.(X2.LT.ACK(L+54))) THEN
            X3=X2/ACK(L+54)
            X4=ONE-X3
            X5=X4**2
            X6=ACK(L+53)*X2*X5/(FIVE*RHOOO)
            EC=EC-X6*X5
            PC=PC+ACK(L+53)*(X3-QCC2)*X4*X5*X2**2
            DPDR=DPDR-X6*(X3*(QCC3*X3-QCC4)+TWO)
         ENDIF
C     
      ELSE
C
C     computation for compressed states
C       --density is above cold reference density
C
         X8=ACK(L+33)*X6
         X5=EXP(-X8)
         X7=X5*ACK(L+32)
         IF ((ACK(L+1).LT.QCC16).OR.(X2.LE.ACK(L+1))) THEN
C     
C         no solid-solid phase change or lower density solid phase
C
            PC=X2*X4*X7-(ACK(L+34)+ACK(L+35)*X3+ACK(L+36)*X4)
            DPDR=(X7*X3*(FIVE*X3+ACK(L+33))-X6*(ACK(L+35)*X6
     &           +TWO*ACK(L+36)))/(THREE*RHOOO)
            CALL ANEI3 (X8,X5,X9) !third exponential integral
            EC=(THREE*ACK(L+32)*X4*X9+(ACK(L+34)+THALF*ACK(L+35)*X3
     &           +THREE*ACK(L+36)*X4)/X2-ACK(L+37))/RHOOO

         ELSEIF (X2.LE.ACK(L+2)) THEN
C     
C         in mixed solid-solid phase change region
C
            PC=ACK(L+7)
            DPDR=ZERO
            EC=ACK(L+8)+PC*(X2-ACK(L+1))/(RHOOO*X2*ACK(L+1))
            CRLS2P=RHOOO*ACK(L+1)
            CRHS2P=RHOOO*ACK(L+2)
         ELSE
C
C         higher density solid phase in solid-solid phase change
C
            PC=X2*X4*X7-(ACK(L+38)+ACK(L+39)*X3+ACK(L+40)*X4)
            DPDR=(X7*X3*(FIVE*X3+ACK(L+33))-X6*(ACK(L+39)*X6
     &           +TWO*ACK(L+40)))/(THREE*RHOOO)
            CALL ANEI3 (X8,X5,X9) !third exponential integral
            EC=ACK(L+9)+(THREE*ACK(L+32)*X4*X9+(ACK(L+38)
     &           +THALF*ACK(L+39)*X3+THREE*ACK(L+40)*X4)/X2)/RHOOO
         END IF
      END IF
C
C     thermal terms
C
      IF ((T.GE.QCC5*ABS(THETA)).AND.(PSIZ.GT.QCC6)) THEN
C
C       temperature is above interpolation region
C         --add perfect gas terms to cold eos and branch to
C           addition of electronic terms.
C
         PN=RHO*FTT
         EN=THALF*FTT
         SN=FT*(THALF*LOG(T/ACK(L+13))-LOG(RHO)+QCC7)
         CVN=THALF*FT
         DPDT=RHO*FT
         DPDR=DPDR+FTT
CHJM
CHJM      modified for molecular clusters
CHJM
         IF(MOLEC) THEN
            PM=-ALFA*RHO*FTT
            EM=-THALF*FT*T*BETA
            SM=-THALF*FT*(LOG(ZZ)+BETA)
            CVM=-CVN*(BETA+DBDT)
            DPDT=DPDT*(ONE-ALFA-DADT)
            DPDR=DPDR*(ONE-ALFA-DADR)
         ENDIF
CHJM
C
      ELSE
C
C    Here we are in usual liquid/solid-vapor interpolation region
C
C        Use general form of PSI interpolation of Section IV-3, supplemeted
C        by 1990 SLT update.  Input variable 32, stored in ACK(L+62), is
C        used to define variable BB that multiplies PSI to adjust how
C        low to high temperature interpolation is performed.
C           (variable 32) = 1.0 forces PSI = 0, no transition to vapor
C                         = 0.0 original interpolation
C       values between 0. and 1. delay onset of vapor.  Note that BB=0 is
C       singular because of division by zero!
C
         BB=ONE-ACK(L+62)
         IF(BB.EQ.ZERO) BB=ONE  !protection from singularity
         PSIB=PSIZ**BB
         X1=THREE*GAMMA
         X2=ONE/(ONE+PSIB)
         X3=PSIB*X2
         X4=TWO+PSIB
         X5=THREE*GAMMA+PSIB
         X6=ONE-THREE*GAMMA
C
         IF (THETA.GE.ZERO) THEN
C
C         normal approximation of debye functions
C
            PN=RHO*FTT*X2*X5
            EN=THALF*FTT*X2*X4
            SN=FT*(FOUR-THREE*LOG(THETA/T)+THALF*(LOG(X2)/BB-X3))
            CVN=EN*(ONE-BB*X3/X4)/T
            DPDT=PN*(ONE+BB*X3*X6/X5)/T
            DPDR=DPDR+PN*(ONE+BB*X3*X6**2/(THALF*X5))/RHO
     &           +THREE*RHO*GAMP*FTT*X2
         ELSE
C
C         full treatment of debye functions
C
            CALL ANDEBY(-THETA,T,DFUNC,D2,DS)
C
            PN=RHO*FTT*(X1*DFUNC+X6*X3)
            EN=FTT*(THREE*DFUNC-THALF*X3)
            SN=FT*(FOUR*DFUNC-THREE*D2+THALF*(LOG(X2)/BB-X3))
            CVN=FT*(THREE*DS-THALF*BB*X3*X2)+EN/T
            DPDT=RHO*FT*(X1*DS+X6*BB*X3*X2)+PN/T
            DPDR=PN/RHO+FTT*(THREE*RHO*(DFUNC-X3)*GAMP
     &           -X1*DS*GAMMA+BB*PSIB*(X6*X2)**2/THALF)+DPDR
         END IF
CHJM
CHJM      modified for molecular clusters
CHJM
         IF(MOLEC) THEN
            PM=-RHO*FTT*X3*ALFA
            EM=-THALF*FTT*X3*BETA
            SM=-THALF*FT*X3*BETA
            CVM=-THALF*FT*X3*(BETA*(ONE+BB*(TWO+BETA)*X2)+DBDT)
            DPDT=DPDT-RHO*FT*X3*(ALFA+BB*(ALFA*(ONE+BETA)-BETA*X6)*X2
     $           +DADT)
            DPDR=DPDR-FTT*X3*(ALFA*(ONE+TWO*BB*ATHIRD*(TWO*X6-ALFA)*X2)
     $           +DADR)
         END IF
      END IF
C
C      Update thermodynamic variables for nuclear and molecular contributions
C
      E=EC+EN+EM
      P=PC+PN+PM
      S=SN+SM
      CV=CVN+CVM
C
CHJM
C
C     Add liquid correction terms if appropriate.
C     This implements Equation 5.22
C
C      CMLT7 is a flag for addition of these terms:
C        CMLT7 < 0 Skip addition
C        CMLT7 = 0 add the terms
C        CMLT7 > 0 Skip addition if T is too low
C
 130  IF (ACK46.LE.ZERO) GO TO 170 !add electronic components
      IF (CMLT7) 170,160,150
 150  IF (T.GE.ACK(L+18)) GO TO 160 !ACK(L=18) is melt temperature
      IF ((RHO.GE.ACK46).OR.(T.LT.ACK(L+49))) GO TO 170
C
C       Addition of liquid correction terms, see section V-1
C
C           Definitions:
C            C43,C44,C45 = ACK(L+43), ACK(L+44), ACK(L+45)
C            alpha, beta, gamma = ACK(L+57), ACK(L+58), ACK(L+59)
C
 160  DO I=1,3 !construct three terms of liquid free energy, eq. 5.22
         XXA(I)=RHO**ACK(I+L+56)*ACK(I+L+42)
      END DO
      X6=XXA(1)/(TWO*SQTS(IPSQTS))  !derivative of Fm wrt temperature
      XXA(1)=SQTS(IPSQTS)*XXA(1)    !dFm/dT in eq. 5.22
      SMLT=-X6                      !liquid addition to entropy
      DPDT=DPDT+ACK(L+57)*X6*RHO    !equation 5.29
      CVMLT=HALF*X6                 !equation 5.28
      EMLT=-HALF*XXA(1)             !equation 5.26, half of first term
      PMLT=ZERO
      DO I=1,3
         EMLT=EMLT+XXA(I)       !equation 5.26
         PMLT=PMLT+ACK(I+L+56)*XXA(I)*RHO !equation 5.25
         DPDR=DPDR+ACK(I+L+56)*(ACK(I+L+56)+ONE)*XXA(I) !equation 5.30
      END DO
C        Again, update basic thermodynamic quantities for melt addition
      E=E+EMLT
      P=P+PMLT
      S=S+SMLT
      CV=CV+CVMLT
C
      GO TO 170                 !now add electronic components
C
C-----------------------------------------------------------------------
C     modification of JWL High Explosive form for very high densities
C     and mixed-phase construction
C
  806 X3=ONE/RHO
      S=ACK(L+32)*EXP(-ACK(L+5)*X3)
      X1=ACK(L+33)*EXP(-ACK(L+6)*X3)
      X2=ACK(L+15)
      P=RHO**X2
      CV=ACK(L+4)
      E=S/ACK(L+5)+X1/ACK(L+6)-ACK(L+35)*P+ACK(L+3)+CV*T
      DPDR=((ACK(L+5)*S+ACK(L+6)*X1)*X3)*X3
     1     -X2*(X2+ONE)*ACK(L+35)*P+X2*CV*T
      DPDT=X2*CV*RHO
      P=S+X1-X2*ACK(L+35)*P*RHO+DPDT*T
      S=CV*(ONE+LOG(T/ACK(L+12))-X2*LOG(RHO/ACK(L+11)))
      IF (RHO.GT.ACK(L+53)) THEN
        X1=RHO**ATHIRD
        P=(ACK(L+24)*X1+ACK(L+54))*RHO*X1+DPDT*T
        DPDR=(FIVE*ACK(L+24)*X1+FOUR*ACK(L+54))*X1/THREE+ACK(L+15)*CV*T
        E=ACK(L+14)+THREE*(HALF*ACK(L+24)*X1+ACK(L+54))*X1+CV*T
      END IF
      IF (NTYPE.EQ.3) GO TO 270
      GO TO 175    !now add electronic components
C
C-----------------------------------------------------------------------
C     electronic and opacity terms
C
  170 IF ((NTYPE.LE.0).OR.(NTYPE.EQ.3)) GO TO 270
  175 CONTINUE
      IF (ACK(L+63).LT.ZERO) GO TO 180  !Flag for detailed T-F computation
C
C     Saha
C
      IF (T.LE.QCC8) THEN     !temp is below 812 K--don't bother with Saha
        ZBAR=ZERO
        GO TO 215
      END IF
      TEVX=T
      RHOX=RHO
      T32X=T*SQTS(IPSQTS)
      ZBARM=ACK(L+26)
      FNX=ACK(L+27)
      NMATSX=NINT(ACK(L+28))
      IF (NMATSX.EQ.1) THEN
C
C       one-element Saha
C
        CALL ANION1
        LNKNTR(MATKTR,4)=LNKNTR(MATKTR,4)+1
        IF (ZBAR.LE.ZERO) GO TO 215
        ZBAR2=ZBAR**2
      ELSE
C
C       multi-element Saha
C
        IIZX=NINT(ACK(L+31))
        CALL ANION2
        LNKNTR(MATKTR,5)=LNKNTR(MATKTR,5)+1
        IF (ZBAR.LE.ZERO) GO TO 215
        ZBAR2=ZERO
          DO 210 I=1,NMATSX
  210     ZBAR2=ZBAR2+COT(I+IIZX-1)*ZB(I)**2
      END IF
      GO TO 182
C
C     DETAILED TREATMENT BY THOMAS-FERMI TABULAR EOS
C
  180 TEVX=T
      RHOX=RHO
      ZBARM=ACK(L+26)
      ABARX=ACK(L+29)
      CALL ANTOMF
      LNKNTR(MATKTR,6)=LNKNTR(MATKTR,6)+1
      ZBAR2=ZBAR**2
      T32X=T*SQTS(IPSQTS)
C
C   Update thermodynamic quantities for electronic components
C
  182 CONTINUE
      P=P+PE
      E=E+EE
      S=S+SE
      CV=CV+CVE
      DPDT=DPDT+DPTE
      DPDR=DPDR+DPRE
      IF (NOTRAD.EQ.0) RETURN
C
C    Calculate Rosseland Mean opacity
C
       FKROS=
     1  (QCC9*RHO*ZBAR*ZBAR2/(ACK(L+29)*T32X*T*T)+QCC10*ACK(L+26))
     2   / ACK(L+29)
      GO TO 220
C
  215 IF (NOTRAD.EQ.0) RETURN
      FKROS=QCC10*ACK(L+26)/ACK(L+29)
C
C-----------------------------------------------------------------------
C     electronic conduction term
C
  220 X4=QCC11*SQTS(IPSQTS)/((RHO*ACK(L+27))**ATHIRD)
      IF (X4.LE.QCC12) THEN
        X4=QCC13
      ELSE
        X4=LOG(X4)
      END IF
      X4=QCC14*MAX(ZBAR,ACK(L+42))*X4*SQTS(IPSQTS)/RHO
      FKROS=FKROS*X4/(X4+FKROS)
      RETURN
C
  270 FKROS=QCC15
      ZBAR=ZERO
C
      RETURN
      END
