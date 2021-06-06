/*
 * The header file for the ANEOS wrapper.
 */
#ifndef LIBANEOS_HINCLUDED
#define LIBANEOS_HINCLUDED

#define FALSE 0
#define TRUE 1

/* Convert from K to eV (1 K = 1 eV in Joule / kbolz in J) */
#define KELVIN_IN_EV 1.16054e4

/**
 * FORTRAN subroutine "ANEOSINIT" to init ANEOS
 */
void aneosinit_(
	const char*,  /// materials file
	int           /// filename string length
);

/**
 * FORTRAN subroutine "ANEOS" for p( rho, T, mat )
 *
 * This methods takes an array of items to achieve multiple calculation with
 * the same call, the number of elements being given in the first parameter.
 */
void aneosv_(
	const int*,     /// no of array elements
	const double*,  /// T (input)
	const double*,  /// rho (input)
	const int*,     /// material (input)
	double*,        /// p
	double*,        /// E
	double*,        /// S
	double*,        /// c_v
	double*,        /// dp/dT
	double*,        /// dp/dRho
	double*,        /// fkro (Rosseland mean opacity)
	double*,        /// cs
	int*,           /// phase
	double*,        /// lower  phase density
	double*,        /// higher phase density
	double*         /// ionization number
);

/**
 * Initialize the ANEOS library.
 *
 * This is a wrapper arround aneosinit Fortran subroutine.
 */
void initaneos(char *matFilename);

/**
 * Compute the ANEOS equation of state for given ( T, rho, mat ) values.
 *
 * This is a wrapper around ANEOSV Fortran subroutine.
 * All units are CGS-eV.
 *
 * @param T Temperature (input)
 * @param rho Density (input)
 * @param mat Material number, as defined in the input file (input)
 * @param p Pressure (output)
 * @param u Specific internal energy (output)
 * @param S Specific entropy (output)
 * @param cv Specific heat capacity at constant volume (output)
 * @param dpdt Temperature derivative of the pressure (output)
 * @param dpdrho Density derivative of the pressure (output)
 * @param fkros Rossland mean opacity (output)
 * @param cs Speed of sound (output)
 * @param kpa Phase (output)
 * @param rhoL Density of the lower phase, for multi-phase states (output)
 * @param rhoH Density of the higher phase, for multi-phase states (output)
 * @param ion Ionisation number (output)
 */
void callaneos (
	const double T, const double rho, const int mat,
	double* p, double* u, double* S, double* cv,
	double* dpdt, double* dpdrho, double* fkros, double* cs,
	int* kpa, double* rhoL, double* rhoH, double* ion
);

/*
 * The same function as callaneos() but temperature has to be provided in K not eV.
 */
void callaneos_cgs (
	const double T, const double rho, const int mat,
	double* p, double* u, double* S, double* cv,
	double* dpdt, double* dpdrho, double* fkros, double* cs,
	int* kpa, double* rhoL, double* rhoH, double* ion
);
#endif

