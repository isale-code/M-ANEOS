/*
 * Calculate all thermodynamic variables from ANEOS.
 *
 * Author:   Christian Reinhardt
 * Created:  29.07.2020
 * Modified:  
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "aneos.h"

int main(int argc, char **argv) {
    double rho;
    double T;
    int iMat;
    char matFilename[256] = "ANEOS.INPUT";
    // Uncomment below to use M-ANEOS
    //char matFilename[256] = "maneos.in";
    double p;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;

    if (argc != 4) {
        fprintf(stderr, "Usage: aneoscall <rho> <T> <iMat>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);
    iMat = atoi(argv[3]);

    assert(rho > 0.0);
    assert(T > 0.0);
    assert(iMat >= 0);

    fprintf(stderr, "ANEOS: Initializing material...\n");
    initaneos(matFilename);

    callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH,
                  &ion);

    printf("Input:\n");
    printf("rho = %15.7E\n", rho);
    printf("T   = %15.7E\n", T);
    printf("iMat= %i\n", iMat);
    printf("\n");

    printf("Output:\n");
    printf("p      = %15.7E\n", p);
    printf("u      = %15.7E\n", u);
    printf("s      = %15.7E\n", s);
    printf("cv     = %15.7E\n", cv);
    printf("dPdT   = %15.7E\n", dPdT);
    printf("dPdrho = %15.7E\n", dPdrho);
    printf("fkros  = %15.7E\n", fkros);
    printf("cs     = %15.7E\n", cs);
    printf("iPhase = %i\n", iPhase);
    printf("rhoL   = %15.7E\n", rhoL);
    printf("rhoH   = %15.7E\n", rhoH);
    printf("ion    = %15.7E\n", ion);
    printf("\n");

    return 0;
}
