/*
 * This is a C wrapper to the Fortran functions of the ANEOS equation of state.
 * It is based on a C++ code that was provided by F. Benitez.
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "aneos.h"

/*
 * Define a global variable that stores if ANEOS was properly initialized.
 */
static int bANEOSInitialized = FALSE;

/**
 * Initialize the ANEOS library.
 *
 * This is a wrapper arround the ANEOSINIT Fortran subroutine.
 */
void initaneos(char *matFilename) {
	if ( bANEOSInitialized ) {
		return;
	}

	bANEOSInitialized = TRUE;
	fprintf(stderr,"Initializing ANEOS with input file %s.\n", matFilename);

	// Use aneos.input for ANEOS
    assert(matFilename != NULL);
    assert(strlen(matFilename) > 0);
    
    // At some point it would be smart to check, if the desired file exists.
	aneosinit_( matFilename, strlen(matFilename) );
}

/**
 * Compute the ANEOS equation of state for given ( T, rho, mat ) values.
 *
 * This is a wrapper around ANEOSV Fortran subroutine.
 * The library will be initialized is not already done, so there is no
 * requirement to call initaneos() before this function.
 * All units are CGS-eV.
 *
 * @param T         Temperature (input)
 * @param rho       Density (input)
 * @param mat       Material number, as defined in the input file (input)
 * @param p         Pressure (output)
 * @param u         Specific internal energy (output)
 * @param S         Specific entropy (output)
 * @param cv        Specific heat capacity at constant volume (output)
 * @param dpdt      Temperature derivative of the pressure (output)
 * @param dpdrho    Density derivative of the pressure (output)
 * @param fkros     Rossland mean opacity (output)
 * @param cs        Speed of sound (output)
 * @param kpa       Phase (output)
 * @param rhoL      Density of the lower phase, for multi-phase states (output)
 * @param rhoH      Density of the higher phase, for multi-phase states (output)
 * @param ion       Ionisation number (output)
 */
void callaneos (
	const double T, const double rho, const int mat,
	double* p, double* u, double* S, double* cv,
	double* dpdt, double* dpdrho, double* fkros, double* cs,
	int* kpa, double* rhoL, double* rhoH, double* ion
) {
	const int n = 1;

//	initaneos();
    // Make sure that the library was initalized.
    assert(bANEOSInitialized);

	aneosv_(
		&n, &T, &rho, &mat, p, u, S, cv, dpdt, dpdrho, fkros,
		cs, kpa, rhoL, rhoH, ion
	);
}

/**
 * Compute the ANEOS equation of state for given ( T, rho, mat ) values in CGS units.
 *
 * This function basically converts the input variables from cgs and Kelvin to cgs and
 * eV then calls the EOS and converts the results back to cgs/K.
 *
 * @param T         Temperature [K] (input)
 * @param rho       Density [g/cc] (input)
 * @param mat       Material number, as defined in the input file (input)
 * @param p         Pressure [erg/cc] (output)
 * @param u         Specific internal energy [erg/g] (output)
 * @param S         Specific entropy [erg/g/K] (output)
 * @param cv        Specific heat capacity at constant volume [erg/g/K] (output)
 * @param dpdt      Temperature derivative of the pressure (output)
 * @param dpdrho    Density derivative of the pressure (output)
 * @param fkros     Rossland mean opacity (output)
 * @param cs        Speed of sound [cm/s] (output)
 * @param kpa       Phase (output)
 * @param rhoL      Density of the lower phase, for multi-phase states (output)
 * @param rhoH      Density of the higher phase, for multi-phase states (output)
 * @param ion       Ionisation number (output)
 */
void callaneos_cgs (
	const double T, const double rho, const int mat,
	double* p, double* u, double* S, double* cv,
	double* dpdt, double* dpdrho, double* fkros, double* cs,
	int* kpa, double* rhoL, double* rhoH, double* ion
) {
    double T_in_eV;
	const int n = 1;

    // Make sure that the library was initalized.
    assert(bANEOSInitialized);

    // Convert from K to eV
    T_in_eV = T/KELVIN_IN_EV;

	aneosv_(
		&n, &T_in_eV, &rho, &mat, p, u, S, cv, dpdt, dpdrho, fkros,
		cs, kpa, rhoL, rhoH, ion
	);

    // Convert the output to K
    *S /= KELVIN_IN_EV;
    *cv /= KELVIN_IN_EV;
}
