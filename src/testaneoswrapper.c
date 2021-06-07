/*
 * Generate an EOS table and compare to the result of the Fortran code.
 *
 * Author:   Christian Reinhardt
 * Created:  06.06.2021
 * Modified:  
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "aneos.h"

#define SUCCESS 0
#define FAIL -1

int ReadTableGrid(char *chFile, double *rho, double *T, int *pnRho, int *pnT) {
    double sesameid;
    double date;
    double version;
    double fmn;
    double fmw;
    double rho0;
    double K0;
    double T0;
    int nRho;
    int nT;
    FILE *fp;
    int iRet;
    int i;

    if (rho != NULL) free(rho);
    if (T != NULL) free(T);

    fp = fopen(chFile, "r");
    if (fp == NULL) return FAIL;

    /* Read the file. */
    iRet = fscanf(fp, "%lf", &sesameid);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &date);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &version);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &fmn);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &fmw);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &rho0);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &K0);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &T0);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%d", &nRho);
    if (iRet != 1) return FAIL;
 
    iRet = fscanf(fp, "%d", &nT);
    if (iRet != 1) return FAIL;

    iRet = fscanf(fp, "%lf", &K0);
    if (iRet != 1) return FAIL;

    /* Allocate memory and read rho and T. */
    rho = calloc(nRho, sizeof(double));
    assert(rho != NULL);

    T = calloc(nT, sizeof(double));
    assert(rho != NULL);

    for (i=0; i<nRho; i++) {
        iRet = fscanf(fp, "%lf", &rho[i]);
        if (iRet != 1) return FAIL;
    }

    for (i=0; i<nT; i++) {
        iRet = fscanf(fp, "%lf", &T[i]);
        if (iRet != 1) return FAIL;
    }

    fclose(fp);

    *pnRho = nRho;
    *pnT = nT;

    return SUCCESS;
}

int main(int argc, char **argv) {
    char matFilename[256] = "ANEOS.INPUT";
    /* Variables for aneoscall(). */
    double rho;
    double T;
    int iMat = 1;
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
    /* EOS table. */
    int nRho;
    int nT;
    double *rhoAxis;
    double *TAxis;
    int i;

    /* Test if reading the file tablegrid.txt works. */
    if (!ReadTableGrid("tablegrid.txt", rhoAxis, TAxis, &nRho, &nT)) {
        fprintf(stderr, "Failed reading tablegrid.txt.\n");
        exit(1);
    }
    
    fprintf(stderr, "nRho= %i nT= %i\n", nRho, nT);
    
    exit(1);

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
