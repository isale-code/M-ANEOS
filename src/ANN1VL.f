C
C
      SUBROUTINE ANN1VL (NUMQQ,TQQ,RHOQQ,KOC,NMATSV,
     &  PQQ,EQQ,SQQ,CVQQ,DPDTQQ,DPDRQQ,FKROQQ,KPAQQ,R2PLQQ,R2PHQQ,
     &  ZBARQQ)
C
C***********************************************************************
C
C     ANEOS PACKAGE - TYPE -1
C
C     NUMQQ =  dimensions of input and output arrays
C              (normally limited to 64)
C     TQQ   =  temperature array
C     RHOQQ =  density array
C     KOC   =  array set in aneosv (pointer to material data)
C     NMATSV=  array set in aneosv (material number)
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
C     KPAQQ =  scratch array (set to 1 before exit)
C     R2PLQQ=  lower phase density for mixed-phase state
C              not used for single phase
C     R2PHQQ=  upper phase density for mixed-phase state
C              not used for single phase
C     ZBARQQ=  ionization number
C
C************************************** 11/89 version slt **************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TQQ(NUMQQ),RHOQQ(NUMQQ)
      DIMENSION PQQ(NUMQQ),EQQ(NUMQQ),SQQ(NUMQQ),CVQQ(NUMQQ)
      DIMENSION DPDTQQ(NUMQQ),DPDRQQ(NUMQQ),FKROQQ(NUMQQ)
      DIMENSION KPAQQ(NUMQQ),R2PLQQ(NUMQQ),R2PHQQ(NUMQQ),ZBARQQ(NUMQQ)
      DIMENSION KOC(NUMQQ), NMATSV(NUMQQ)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      SAVE   /FILEOS/,/ANES/
      PARAMETER (ZERO=0.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (SIX=6.D0)
      PARAMETER (OV3=1.D0/3.D0)
      PARAMETER (OV5=1.D0/5.D0)
C
C     determine form to be used if this eos is desired
C
        DO 10 LMAT=1,NUMQQ
        IF(NMATSV(LMAT).EQ.-1) THEN
          IF(RHOQQ(LMAT).LE.ACK(KOC(LMAT)+48)) THEN
            KPAQQ(LMAT)=-102
          ELSEIF(RHOQQ(LMAT).GE.ACK(KOC(LMAT)+37)) THEN
            KPAQQ(LMAT)=-101
          ELSEIF(ACK(KOC(LMAT)+33).GT.ZERO) THEN
            KPAQQ(LMAT)=-100
          ELSE
            KPAQQ(LMAT)=-103
          END IF
        END IF
   10   CONTINUE
C
C     determine cold terms for Pc, Ec and dPc/dRHO
C
C     Us-Up form in allowed density range
C
        DO 20 LMAT=1,NUMQQ
        IF(KPAQQ(LMAT).EQ.-100) THEN
C         eta stored in eqq
          EQQ(LMAT)=MAX(ACK(KOC(LMAT)+42),MIN(ACK(KOC(LMAT)+43),
     &    ONE-ACK(KOC(LMAT)+11)/RHOQQ(LMAT)))
C         delta eta for table
          R2PLQQ(LMAT)=(ACK(KOC(LMAT)+43)-ACK(KOC(LMAT)+42))/
     &    ACK(KOC(LMAT)+34)
C         find interval in stored table
          KPAQQ(LMAT)=INT(
     &    (EQQ(LMAT)-ACK(KOC(LMAT)+42))/R2PLQQ(LMAT))
C         interpolate in table
          R2PHQQ(LMAT)=EQQ(LMAT)-KPAQQ(LMAT)*R2PLQQ(LMAT)
     &    -ACK(KOC(LMAT)+42)
          R2PLQQ(LMAT)=
     &    (ACK(KOC(LMAT)+KPAQQ(LMAT)+75)*(R2PLQQ(LMAT)-R2PHQQ(LMAT))+
     &    ACK(KOC(LMAT)+KPAQQ(LMAT)+76)*R2PHQQ(LMAT))/R2PLQQ(LMAT)
C         cvqq is Th interpolated
          CVQQ(LMAT)=EXP(ACK(KOC(LMAT)+15)*EQQ(LMAT))*
     &    (ACK(KOC(LMAT)+12)+HALF*ACK(KOC(LMAT)+32)**2*
     &    (EQQ(LMAT)/(ONE-ACK(KOC(LMAT)+33)*EQQ(LMAT)))**3
     &    *R2PLQQ(LMAT)/ACK(KOC(LMAT)+4))
C         Ph
          PQQ(LMAT)=ACK(KOC(LMAT)+21)*EQQ(LMAT)/
     &    (ONE-ACK(KOC(LMAT)+33)*EQQ(LMAT))**2
C            write(klst,9191) lmat,rhoqq(lmat),eqq(lmat),
C    &       tqq(lmat),cvqq(lmat),pqq(lmat),r2plqq(lmat)
C    &       ,kpaqq(lmat)
C9191        format(' ann1vl',i4,1p6e13.6,i6)
C         dPc/dRHO
          DPDRQQ(LMAT)=(ACK(KOC(LMAT)+21)*
     &    (ONE+ACK(KOC(LMAT)+33)*EQQ(LMAT))*
     &    (ONE-HALF*ACK(KOC(LMAT)+15)*EQQ(LMAT))
     &    /(ONE-ACK(KOC(LMAT)+33)*EQQ(LMAT))**3-ACK(KOC(LMAT)+15)*
     &    (ACK(KOC(LMAT)+5)*CVQQ(LMAT)-HALF*PQQ(LMAT)))
     &    *ACK(KOC(LMAT)+11)/RHOQQ(LMAT)**2
C         Ec
          EQQ(LMAT)=PQQ(LMAT)*EQQ(LMAT)/(TWO*ACK(KOC(LMAT)+11))
     &    +ACK(KOC(LMAT)+3)-ACK(KOC(LMAT)+4)*CVQQ(LMAT)
C         Pc
          PQQ(LMAT)=PQQ(LMAT)-ACK(KOC(LMAT)+5)*CVQQ(LMAT)
        END IF
   20   CONTINUE
C
C     high density form
C
        DO 30 LMAT=1,NUMQQ
        IF(KPAQQ(LMAT).EQ.-101) THEN
C         dPc/dRHO
          DPDRQQ(LMAT)=TWO*ACK(KOC(LMAT)+39)*RHOQQ(LMAT)
C         Pc
          PQQ(LMAT)=ACK(KOC(LMAT)+38)+ACK(KOC(LMAT)+39)*RHOQQ(LMAT)**2
C         Ec
          EQQ(LMAT)=(ACK(KOC(LMAT)+38)/(RHOQQ(LMAT)*ACK(KOC(LMAT)+37))
     &    +ACK(KOC(LMAT)+39))*(RHOQQ(LMAT)-ACK(KOC(LMAT)+37))
     &    +ACK(KOC(LMAT)+40)
        END IF
   30   CONTINUE
C
C     low density form
C
        DO 40 LMAT=1,NUMQQ
        IF(KPAQQ(LMAT).EQ.-102) THEN
C         dPc/dRHO
          DPDRQQ(LMAT)=-ACK(KOC(LMAT)+50)/RHOQQ(LMAT)**2
C         Pc
          PQQ(LMAT)=ACK(KOC(LMAT)+49)+ACK(KOC(LMAT)+50)/RHOQQ(LMAT)
C         Ec
          CVQQ(LMAT)=ACK(KOC(LMAT)+48)/RHOQQ(LMAT)
          EQQ(LMAT)=( ACK(KOC(LMAT)+49)*(ONE-CVQQ(LMAT)) +
     &      HALF*ACK(KOC(LMAT)+50)*(ONE-CVQQ(LMAT)**2)/ACK(KOC(LMAT)+48)
     &      )/ACK(KOC(LMAT)+48)
     &       +ACK(KOC(LMAT)+51)
        END IF
   40   CONTINUE
C
C     power series form
C
        DO 50 LMAT=1,NUMQQ
        IF(KPAQQ(LMAT).EQ.-103) THEN
C         eta stored in eqq
          EQQ(LMAT)=MAX(ACK(KOC(LMAT)+42),MIN(ACK(KOC(LMAT)+43),
     &    ONE-ACK(KOC(LMAT)+11)/RHOQQ(LMAT)))
C         delta eta for table
          R2PLQQ(LMAT)=(ACK(KOC(LMAT)+43)-ACK(KOC(LMAT)+42))/
     &    ACK(KOC(LMAT)+34)
C         find interval in stored table
          KPAQQ(LMAT)=INT(
     &    (EQQ(LMAT)-ACK(KOC(LMAT)+42))/R2PLQQ(LMAT))
C         interpolate in table
          R2PHQQ(LMAT)=EQQ(LMAT)-KPAQQ(LMAT)*R2PLQQ(LMAT)
     &    -ACK(KOC(LMAT)+42)
          R2PLQQ(LMAT)=
     &    (ACK(KOC(LMAT)+KPAQQ(LMAT)+75)*(R2PLQQ(LMAT)-R2PHQQ(LMAT))+
     &    ACK(KOC(LMAT)+KPAQQ(LMAT)+76)*R2PHQQ(LMAT))/R2PLQQ(LMAT)
C         cvqq is Th interpolated
          CVQQ(LMAT)=EXP(ACK(KOC(LMAT)+15)*EQQ(LMAT))*
     &    (ACK(KOC(LMAT)+12)+ACK(KOC(LMAT)+32)**2*ACK(KOC(LMAT)+44)*
     &    EQQ(LMAT)**3*R2PLQQ(LMAT)/(SIX*ACK(KOC(LMAT)+4)))
C         Ph
          DPDRQQ(LMAT)=ACK(KOC(LMAT)+21)*
     &    (ONE+EQQ(LMAT)*(ACK(KOC(LMAT)+44)+EQQ(LMAT)*(ACK(KOC(LMAT)+45)
     &    +EQQ(LMAT)*(ACK(KOC(LMAT)+53)+EQQ(LMAT)*ACK(KOC(LMAT)+54)))))
          PQQ(LMAT)=DPDRQQ(LMAT)*EQQ(LMAT)
C            write(klst,9193) lmat,rhoqq(lmat),eqq(lmat),
C    &       tqq(lmat),cvqq(lmat),pqq(lmat),r2plqq(lmat)
C    &       ,kpaqq(lmat)
C9193        format(' ANN1VL',i4,1p6e13.6,i6)
C         dPc/dRHO
          DPDRQQ(LMAT)=(DPDRQQ(LMAT)-ACK(KOC(LMAT)+15)*ACK(KOC(LMAT)+5)
     &    *CVQQ(LMAT)+EQQ(LMAT)*(ONE-HALF*ACK(KOC(LMAT)+15)*EQQ(LMAT))
     &    *ACK(KOC(LMAT)+21)*(ACK(KOC(LMAT)+44)+EQQ(LMAT)*
     &    (TWO*ACK(KOC(LMAT)+45)+EQQ(LMAT)*
     &    (THREE*ACK(KOC(LMAT)+53)+EQQ(LMAT)*FOUR*ACK(KOC(LMAT)+54)))))
     &    *ACK(KOC(LMAT)+11)/RHOQQ(LMAT)**2
C         Ec
          EQQ(LMAT)=PQQ(LMAT)*EQQ(LMAT)/(TWO*ACK(KOC(LMAT)+11))
     &    +ACK(KOC(LMAT)+3)-ACK(KOC(LMAT)+4)*CVQQ(LMAT)
C         Pc
          PQQ(LMAT)=PQQ(LMAT)-ACK(KOC(LMAT)+5)*CVQQ(LMAT)
        END IF
   50   CONTINUE
C
C     now have cold terms for Pc, Ec and dPc/dRHO
C     add in thermal terms
C
        DO 70 LMAT=1,NUMQQ
        IF(NMATSV(LMAT).EQ.-1) THEN
          CVQQ(LMAT)=ACK(KOC(LMAT)+4)
          DPDTQQ(LMAT)=ACK(KOC(LMAT)+5)
          EQQ(LMAT)=CVQQ(LMAT)*TQQ(LMAT)+EQQ(LMAT)
          PQQ(LMAT)=DPDTQQ(LMAT)*TQQ(LMAT)+PQQ(LMAT)
          SQQ(LMAT)=ACK(KOC(LMAT)+5)/RHOQQ(LMAT)
     &      +CVQQ(LMAT)*(LOG(TQQ(LMAT)/ACK(KOC(LMAT)+12))+ONE)
          FKROQQ(LMAT)=1.E11
          R2PHQQ(LMAT)=ZERO
          R2PLQQ(LMAT)=ZERO
          ZBARQQ(LMAT)=ZERO
          KPAQQ(LMAT)=1
        END IF
   70   CONTINUE
C
C     cut out for P < 0 for solid tensions
C
        DO 110 LMAT=1,NUMQQ
        IF((NMATSV(LMAT).EQ.-1).AND.
     &      (PQQ(LMAT).LT.ZERO).AND.
     &        ((SQQ(LMAT).GT.ACK(KOC(LMAT)+96)).OR.
     &         (RHOQQ(LMAT).LT.ACK(KOC(LMAT)+23)))) THEN
C          phase density
           R2PHQQ(LMAT)=-ACK(KOC(LMAT)+50)/
     &                  (ACK(KOC(LMAT)+49)+ACK(KOC(LMAT)+5)*TQQ(LMAT))
           DPDRQQ(LMAT)=ACK(KOC(LMAT)+48)/R2PHQQ(LMAT)
           EQQ(LMAT)=( ACK(KOC(LMAT)+49)*(ONE-DPDRQQ(LMAT)) +
     &       HALF*ACK(KOC(LMAT)+50)*(ONE-DPDRQQ(LMAT)**2)
     &       /ACK(KOC(LMAT)+48) )/ACK(KOC(LMAT)+48)
     &       +ACK(KOC(LMAT)+51) + CVQQ(LMAT)*TQQ(LMAT)
           SQQ(LMAT)=ACK(KOC(LMAT)+5)*
     &              (ONE/R2PHQQ(LMAT)-ONE/RHOQQ(LMAT)) + SQQ(LMAT)
           CVQQ(LMAT)=CVQQ(LMAT)-TQQ(LMAT)*(ACK(KOC(LMAT)+5)**2)
     &                                     /ACK(KOC(LMAT)+50)
           R2PLQQ(LMAT)=1.2345E-33
           DPDTQQ(LMAT)=1.2345E-22
           DPDRQQ(LMAT)=ZERO
           PQQ(LMAT)=ZERO
           KPAQQ(LMAT)=2
        END IF
  110   CONTINUE
C
      RETURN
      END
