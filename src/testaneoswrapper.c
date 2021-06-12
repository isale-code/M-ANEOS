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

struct ANEOSTable {
    int iMat;
    double SesameId;
    double Date;
    double Version;
    double fmn;
    double fmw;
    double rho0;
    double K0;
    double T0;
    int nRho;
    int nT;
    double *rho;
    double *T;
    double **P;
    double **u;
    double **s;
    double **cs;
    double **cv;
    int **phase;
};

/*
 * Allocate memory and read rho and T from a file.
 */
struct ANEOSTable *ANEOSTableInit(char *chFile) {
    struct ANEOSTable *Table;
    double tmp;
    FILE *fp;
    int iRet;
    int i;

    Table = (struct ANEOSTable *) calloc(1, sizeof(struct ANEOSTable));
    assert(Table != NULL);

    fp = fopen(chFile, "r");
    assert(fp != NULL);

    /* Read the file. */
    iRet = fscanf(fp, "%lf", &Table->SesameId);
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->Date);
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->Version);
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->fmn);
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->fmw);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->rho0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->K0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &Table->T0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);
    
    iRet = fscanf(fp, "%lf", &tmp);
    assert(iRet == 1);
    Table->nRho = (int) tmp;
    assert(Table->nRho > 0);

 
    iRet = fscanf(fp, "%lf", &tmp);
    assert(iRet == 1);
    Table->nT = (int) tmp;
    assert(Table->nT > 0);

    /* Allocate memory and read rho and T. */
    Table->rho = calloc(Table->nRho, sizeof(double));
    assert(Table->rho != NULL);

    Table->T = calloc(Table->nT, sizeof(double));
    assert(Table->rho != NULL);

    for (i=0; i<Table->nRho; i++) {
        iRet = fscanf(fp, "%lf", &Table->rho[i]);
        assert(iRet == 1);
    }

    for (i=0; i<Table->nT; i++) {
        iRet = fscanf(fp, "%lf", &Table->T[i]);
        assert(iRet == 1);
    }

    fclose(fp);

    /* Allocate memory for p, u, s, cs and cv. */
    Table->P = (double **) calloc(Table->nRho, sizeof(double *));
    Table->u = (double **) calloc(Table->nRho, sizeof(double *));
    Table->s = (double **) calloc(Table->nRho, sizeof(double *));
    Table->cs = (double **) calloc(Table->nRho, sizeof(double *));
    Table->cv = (double **) calloc(Table->nRho, sizeof(double *));
    Table->phase = (int **) calloc(Table->nRho, sizeof(int *));

    for (i=0; i<Table->nRho; i++) {
        Table->P[i] = (double *) calloc(Table->nT, sizeof(double));
        Table->u[i] = (double *) calloc(Table->nT, sizeof(double));
        Table->s[i] = (double *) calloc(Table->nT, sizeof(double));
        Table->cs[i] = (double *) calloc(Table->nT, sizeof(double));
        Table->cv[i] = (double *) calloc(Table->nT, sizeof(double));
        Table->phase[i] = (int *) calloc(Table->nT, sizeof(int));
    }

    return Table;
}

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
    double tmp;
    FILE *fp;
    int iRet;
    int i;

    if (rho != NULL) free(rho);
    if (T != NULL) free(T);

    fp = fopen(chFile, "r");
    //if (fp == NULL) return FAIL;
    assert(fp != NULL);

    /* Read the file. */
    iRet = fscanf(fp, "%lf", &sesameid);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &date);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &version);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &fmn);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &fmw);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &rho0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &K0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);

    iRet = fscanf(fp, "%lf", &T0);
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);
    
    iRet = fscanf(fp, "%lf", &tmp);
    nRho = (int) tmp;
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);
    assert(nRho > 0);

 
    iRet = fscanf(fp, "%lf", &tmp);
    nT = (int) tmp;
    //if (iRet != 1) return FAIL;
    assert(iRet == 1);
    assert(nT > 0);

    /* Allocate memory and read rho and T. */
    rho = calloc(nRho, sizeof(double));
    assert(rho != NULL);

    T = calloc(nT, sizeof(double));
    assert(rho != NULL);

    for (i=0; i<nRho; i++) {
        iRet = fscanf(fp, "%lf", &rho[i]);
        //if (iRet != 1) return FAIL;
        assert(iRet == 1);
    }

    for (i=0; i<nT; i++) {
        iRet = fscanf(fp, "%lf", &T[i]);
        //if (iRet != 1) return FAIL;
        assert(iRet == 1);
    }

    fclose(fp);

    *pnRho = nRho;
    *pnT = nT;

#if 0
    /* Write the data to tablegrid.new.*/
    fp = fopen("tablegrid.new", "w");
    assert(fp != NULL);

    /* Write the data to the file. */
    fprintf(fp, "%11.6e\n", sesameid);
    fprintf(fp, "%11.6e\n", date);
    fprintf(fp, "%11.6e\n", version);
    fprintf(fp, "%11.6e\n", fmn);
    fprintf(fp, "%11.6e\n", fmw);
    fprintf(fp, "%11.6e\n", rho0);
    fprintf(fp, "%11.6e\n", K0);
    fprintf(fp, "%11.6e\n", T0);
    fprintf(fp, "%11.6e\n", (double) nRho);
    fprintf(fp, "%11.6e\n", (double) nT);

    for (i=0; i<nRho; i++) {
        fprintf(fp, "%11.6e\n", rho[i]);
    }

    for (i=0; i<nT; i++) {
        fprintf(fp, "%11.6e\n", T[i]);
    }

    fclose(fp);
#endif
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
    struct ANEOSTable *Table;
    int nwds;
    int nTables;
    int Table1;
    int Table2;
    int nwds1;
    int nwds2;
    int nValues;
    FILE *fp;
    /* Convert from eV to K. */
    double Tconv = 1.16045e4;
    int i, j;

#if 0
    /* Reading the file tablegrid.txt. */
    if (ReadTableGrid("tablegrid.txt", rhoAxis, TAxis, &nRho, &nT)) {
        fprintf(stderr, "Failed reading tablegrid.txt.\n");
        exit(1);
    }
    
    fprintf(stderr, "nRho= %i nT= %i\n", nRho, nT);
#endif

    /* Initialize ANEOS. */
    fprintf(stderr, "ANEOS: Initializing material.\n");
    initaneos(matFilename);
    fprintf(stderr, "ANEOS: Done.\n");

    /*
     * Generate the EOS table.
     */
    Table = ANEOSTableInit("tablegrid.txt");
    assert(Table != NULL);

    fprintf(stderr, "Open input file.\n");
    fprintf(stderr, "nRho= %i nT= %i\n", Table->nRho, Table->nT);

    fprintf(stderr, "rho_min= %16.8E rho_max= %16.8E\n", Table->rho[0], Table->rho[Table->nRho-1]);
    fprintf(stderr, "T_min= %16.8E T_max= %16.8E\n", Table->T[0], Table->T[Table->nT-1]);

//#if 0 
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {
            rho = Table->rho[i];
            T = Table->T[j];

            callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL,
                          &rhoH, &ion);

            /* Limit the pressure. */
            if (fabs(p) < 1e-20) p = 1e-20;

            /* Convert pressure from cgs to GPa. */
            p *= 1e-10;

            /* Convert specific internal energy from erg/g to MJ/kg. */
            u *= 1e-10;

            /* Convert specific entropy from ergs/g/K to MJ/kg/K. */
            s *= 1e-10;

            /* Convert specific heat capacity from erg/g/K to MJ/kg/K. */
            cv *= 1e-10;

            Table->P[i][j] = p;
            Table->u[i][j] = u;
            Table->s[i][j] = s;
            Table->cs[i][j] = cs;
            Table->cv[i][j] = cv;
            Table->phase[i][j] = iPhase;
        }
    }
//#endif

#if 0
    i = 0;
    j = 0;

    fprintf(stderr, "i= %i, j= %i, T= %16.8E K = %16.8E eV, rho= %16.8E\n", i, j, Table->T[j],
            Table->T[j]/Tconv, Table->rho[i]);

    fprintf(stderr, "T= %16.8E eV (normal), T= %16.8E eV (define) frac= %16.8E\n", Table->T[j]/Tconv, Table->T[j]/KELVIN_IN_EV, Table->T[j]/Tconv/Table->T[j]/KELVIN_IN_EV);

    /* Call ANEOS in the usual units. */
    fprintf(stderr, "Call ANEOS: T= %16.8E eV, rho= %16.8E\n", Table->T[j]/Tconv, Table->rho[i]);

    callaneos(Table->T[j]/Tconv, Table->rho[i], iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs,
                  &iPhase, &rhoL, &rhoH, &ion);
    
    fprintf(stderr, "ANEOS units: cgs and eV\n");
    fprintf(stderr, "P = %16.8E [erg/cc]\n", p);
    fprintf(stderr, "u = %16.8E [erg/g]\n", u);
    fprintf(stderr, "s = %16.8E [erg/g/eV] \n", s);
    fprintf(stderr, "cs = %16.8E [cm/s]\n", cs);
    fprintf(stderr, "cv = %16.8E [erg/g/eV]\n", cv);
    fprintf(stderr, "phase = %i\n", iPhase);
    fprintf(stderr, "\n"); 

    /* Call ANEOS in cgs and K. */
    fprintf(stderr, "Call ANEOS: T= %16.8E K, rho= %16.8E\n", Table->T[j], Table->rho[i]);
    
    callaneos_cgs(Table->T[j], Table->rho[i], iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs,
                  &iPhase, &rhoL, &rhoH, &ion);
      
    fprintf(stderr, "ANEOS units: cgs and K\n");
    fprintf(stderr, "P = %16.8E [erg/cc]\n", p);
    fprintf(stderr, "u = %16.8E [erg/g]\n", u);
    fprintf(stderr, "s = %16.8E [erg/g/K] \n", s);
    fprintf(stderr, "cs = %16.8E [cm/s]\n", cs);
    fprintf(stderr, "cv = %16.8E [erg/g/K]\n", cv);
    fprintf(stderr, "phase = %i\n", iPhase);
    fprintf(stderr, "\n"); 

    /* Convert pressure from cgs to GPa. */
    p *= 1e-10;
    /* Convert specific internal energy from erg/g to MJ/kg. */
    u *= 1e-10;
    /* Convert specific entropy from ergs/g/K to MJ/kg/K. */
    s *= 1e-10;
    /* Convert specific heat capacity from erg/g/K to MJ/kg/K. */
    cv *= 1e-10;

    fprintf(stderr, "ANEOS units: SESAME table\n");
    fprintf(stderr, "P = %16.8E [GPa]\n", p);
    fprintf(stderr, "u = %16.8E [MJ/kg]\n", u);
    fprintf(stderr, "s = %16.8E [MJ/kg/K] \n", s);
    fprintf(stderr, "cs = %16.8E [cm/s]\n", cs);
    fprintf(stderr, "cv = %16.8E [MJ/kg/K]\n", cv);
    fprintf(stderr, "phase = %i\n", iPhase);
    fprintf(stderr, "\n"); 

    exit(1);
#endif
#if 0
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {

            callaneos_cgs(Table->T[j], Table->rho[i], iMat, &Table->P[i][j], &Table->u[i][j],
                    &Table->s[i][j], &Table->cv[i][j], &dPdT, &dPdrho, &fkros, &cs,
                    &Table->phase[i][j], &rhoL, &rhoH, &ion);

#if 0
            /* Convert specific entropy from ergs/g/K to MJ/kg/K. */
            s *= 1e-10;
            callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL,
                    &rhoH, &ion);
#endif
            /* Limit the pressure. */
            if (fabs(p) < 1e-20) p = 1e-20;

            /* Convert pressure from cgs to GPa. */
            Table->P[i][j] *= 1e-10;

            /* Convert specific internal energy from erg/g to MJ/kg. */
            Table->u[i][j] *= 1e-10;

            /* Convert specific entropy from ergs/g/K to MJ/kg/K. */
            Table->s[i][j] *= 1e-10;

            /* Convert specific heat capacity from erg/g/K to MJ/kg/K. */
            Table->cv[i][j] *= 1e-10;
        }
    }
#endif


    /*
     * Write the SESAME table.
     */
    fp = fopen("NEW-SESAME-EXT.NEW", "w");
    assert(fp != NULL);

    nwds = 9;
    nTables = 2;
    Table1 = 201;
    Table2 = 301;
    nwds1 = 5;
    nwds2 = 2+Table->nRho+Table->nT+Table->nRho*Table->nT*4;

    /* Write the header. */
    fprintf(fp, " INDEX     MATID =%7i   NWDS = %8i\n", (int) Table->SesameId, nwds);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E%16.8E\n", Table->SesameId, Table->Date, Table->Date,
                                                    Table->Version, (double) nTables);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E\n", (double) Table1, (double) Table2, (double) nwds1, (double) nwds2);

    /* Table 1. */
    fprintf(fp, " RECORD     TYPE =%5i     NWDS = %8i\n", Table1, nwds1);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E%16.8E\n", Table->fmn, Table->fmw, Table->rho0, Table->K0,
                                                    Table->T0);
    /* Table 2. */
    fprintf(fp, " RECORD     TYPE =%5i     NWDS = %8i\n", Table2, nwds2);

    /* Number of grid points in rho and T. */
    fprintf(fp, "%16.8E%16.8E", (double) Table->nRho, (double) Table->nT);

    /* Print rho and T axis. */
    nValues = 2;

    /* Density */
    for (i=0; i<Table->nRho; i++) {
        fprintf(fp, "%16.8E", Table->rho[i]);
        nValues++;

        /* Print a new line after every 5 grid points. */
        if ((nValues % 5) == 0) fprintf(fp, "\n");
    }
    
    /* Temperature */
    for (i=0; i<Table->nT; i++) {
        fprintf(fp, "%16.8E", Table->T[i]);
        nValues++;

        if ((nValues % 5) == 0) fprintf(fp, "\n");
    }

    /* Specific entropy */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {

            fprintf(fp, "%16.8E", Table->s[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Sound speed */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", Table->cs[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Specific heat capacity */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", Table->cv[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Phase */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", (double) Table->phase[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    fclose(fp);
    /*
     * Write the standard (short) SESAME table.
     */
    fp = fopen("NEW-SESAME-STD.NEW", "w");
    assert(fp != NULL);

    nwds = 9;
    nTables = 2;
    Table1 = 201;
    Table2 = 301;
    nwds1 = 5;
    nwds2 = 2+Table->nRho+Table->nT+Table->nRho*Table->nT*3;

    /* Write the header. */
    fprintf(fp, " INDEX     MATID = %7i   NWDS = %8i\n", (int) Table->SesameId, nwds);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E%16.8E\n", Table->SesameId, Table->Date, Table->Date,
                                                    Table->Version, (double) nTables);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E\n", (double) Table1, (double) Table2, (double) nwds1, (double) nwds2);

    /* Table 1. */
    fprintf(fp, " RECORD     TYPE =%5i     NWDS = %8i\n", Table1, nwds1);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E%16.8E\n", Table->fmn, Table->fmw, Table->rho0, Table->K0,
                                                    Table->T0);
    /* Table 2. */
    fprintf(fp, " RECORD     TYPE =%5i     NWDS = %8i\n", Table2, nwds2);

    /* Number of grid points in rho and T. */
    fprintf(fp, "%16.8E%16.8E", (double) Table->nRho, (double) Table->nT);

    /* Print rho and T axis. */
    nValues = 2;

    /* Density */
    for (i=0; i<Table->nRho; i++) {
        fprintf(fp, "%16.8E", Table->rho[i]);
        nValues++;

        /* Print a new line after every 5 grid points. */
        if ((nValues % 5) == 0) fprintf(fp, "\n");
    }

    /* Temperature */
    for (i=0; i<Table->nT; i++) {
        fprintf(fp, "%16.8E", Table->T[i]);
        nValues++;

        if ((nValues % 5) == 0) fprintf(fp, "\n");
    }

    /* Pressure */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", Table->P[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Specific internal energy */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", Table->u[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Helmholtz free energy */
    for (j=0; j<Table->nT; j++) {
        for (i=0; i<Table->nRho; i++) {
            fprintf(fp, "%16.8E", Table->u[i][j]-Table->T[j]*Table->s[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

#if 0
    /* Specific entropy */
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {
            fprintf(fp, "%16.8E", Table->s[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Sound speed */
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {
            fprintf(fp, "%16.8E", Table->cs[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Specific heat capacity */
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {
            fprintf(fp, "%16.8E", Table->cv[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }

    /* Phase */
    for (i=0; i<Table->nRho; i++) {
        for (j=0; j<Table->nT; j++) {
            fprintf(fp, "%15.8E", (double) Table->phase[i][j]);
            nValues++;

            if ((nValues % 5) == 0) fprintf(fp, "\n");
        }
    }
#endif

    fclose(fp);

    exit(1);
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
