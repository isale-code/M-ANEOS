      SUBROUTINE ANMTAB

      

      IF (ACK(IT+46).LE.ZERO) GO TO 1190
      DFN1=ACK(IT+72)
      I1=0
      DFN2=ZERO
      PCP1=ZERO
      ACK(IT+71)=ZERO
      CMLT7=ACK(IT+71)
      ACK(IT+72)=ONE
      WRITE(KLST,1350)       !write header for melt curve output
      BOOT=ACK(IT+49)        !define initial temperature for melt curve
      DTBELOW=(ACK(IT+18)-ACK(IT+49))/99.D0  ! 99 linear steps between Tmin and Tm
      DTABOVE=(ACK(IT+48)-ACK(IT+18))/900.D0 ! 900 linear steps between Tm and Tmax
      DO WHILE (BOOT.LT.ACK(IT+48)  ! Loop over temperature until upper bound
         I2=2
         S7=1.D10
         SQTS(1)=SQRT(BOOT)
         CALL ANELSM (BOOT,S7,S1,S2,IT,I2) !get density of solid, liquid
         IF (I2.GE.0) THEN      !no error from ANELSM
                                !Store results in table   
            WRITE(KLST,1360) BOOT*TCONV,S1*DCONV,S2*DCONV,S3*PCONV,
     $           S3S*PCONV,S4*ECONV,S4S*ECONV,S5*SCONV,S5S*SCONV,
     $           S9*ECONV,S9S*ECONV,I2
         ELSE
            I1=1                ! Error
         END IF

         DFN2=DFN2+ONE          ! Counter for points in table

         ! Increment the temperature differently above and below the triple point
         IF (DFN2-100. < 0.) THEN  ! Temp below the melt temperature
            BOOT=BOOT+DTBELOW
         ELSEIF (DFN2-100. > 0.) THEN  ! Temperature above the melt temperature
            BOOT=BOOT+DTABOVE
         ELSE
            BOOT=ACK(IT+18)  ! The melt temperature
         END IF

      END DO
      
      ! Report error if EoS is not properly defined
      IF (I1.EQ.1) THEN
         CALL ANWARN(0)
         WRITE(KLST,1380)
      END IF

      RETURN
      END
