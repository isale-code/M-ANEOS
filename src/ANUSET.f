C
C
      SUBROUTINE ANUSET(
C         inputs
     &    NEOS, URINP,
C         outputs
     &    UARR)
C
C***********************************************************************
C
C     User defined EOS setup entry point. See ANUEOS for running
C     entry point.
C
C     WARNING - depending on the computer system
C     All variables are single or double precision.
C
C     Inputs: (all units cgs-ev)
C
C     NEOS   =  User input eos number (negative)
C
C     URINP  =  input array
C
C     Output:
C
C     UARR   =  aneos data array for this material
C
C     The input array is the aneos input record - either 24 or 40
C     variables are defined depending on aneos input record format
C     that you selected. Of these the following must be defined:
C
C     urinp(1) = number of elements in this material
C     urinp(2) = -4.  (to get here)
C     urinp(3) = reference density for this material
C     urinp(4) = reference temperature for this material
C     The others can be anything that you want to input to your model.
C
C     The other input array is of length 99. Certain location in this
C     array must not be redefined by your model. Any of the remaining
C     locations can be used by your model. Note that this array is
C     written to the restart file so that you don't have to do anything
C     about restarts. * indicates values that you MUST set.
C
C     UARR location       value
C         11        reference density
C         12        reference temperature
C     *   18        melt temperature
C         22        zero
C         23        internal value
C         26        average atomic number
C         27        number of atoms per gram
C         28        number of elements
C         29        average atomic weight
C         30        internal switch
C         31        internal pointer
C     *   47        liquid density at melt point for solid form
C         64        internal value
C         65        internal value
C         71        internal value
C         72        internal value
C         73        internal value
C         74        internal value
C
C     EOS values must be defined for all positive temperatures and
C     densities. Returned values must be realistic - negative dP/dRHO
C     is not allowed.
C
C************************************** 4/90 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION URINP(48), UARR(99)
C     klst is the output listing file
      COMMON /FILEOS/ KLST,KINP
      SAVE   /FILEOS/
      PARAMETER (ZERO=0.D0)
C
C     add your code here
C
C     example code for ideal gas ---------------------------------------
C     input assumed
C     aneos #elements -4  rho0   temp0   gamma   Cv
C     aneos 0
C     aneos 0
C     aneos element records
C
      UARR(1)=URINP(5)
      UARR(2)=URINP(6)
      UARR(18)=ZERO
C
C     example code for ideal gas ---------------------------------------
C
      RETURN
      END
