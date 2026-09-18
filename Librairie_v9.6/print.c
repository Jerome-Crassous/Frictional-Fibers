#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//---------------------------------------------------------------------
void PrintPositionsLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int index) {

	unsigned  int		i, iFiber;
	char                fileName[256];
	FILE				* filePtr;

	//----------------------------------------------
	sprintf(fileName, "%sposition_%d_%d.txt", pathName, theParameterPtr->line, index);
	filePtr = fopen(fileName, "w");
	if (filePtr == NULL) {
		printf("### opening %s failed\n", fileName);
		return;
	}
	fprintf(filePtr, "iter fiber i xx yy zz thetax thetay thetaz m1barx m1bary m1barz\n");

	//----------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			if (theFiberPtr[iFiber].status != STATUS_VIRTUAL) {
				fprintf(filePtr, "%d %d %d %e %e %e %e %e %e %e %e %e\n",
					theParameterPtr->iter, iFiber, i,
					theFiberPtr[iFiber].xt[i][0], theFiberPtr[iFiber].xt[i][1], theFiberPtr[iFiber].xt[i][2],
					theFiberPtr[iFiber].thetat[i][0], theFiberPtr[iFiber].thetat[i][1], theFiberPtr[iFiber].thetat[i][2],
					theFiberPtr[iFiber].m1[i][0], theFiberPtr[iFiber].m1[i][1], theFiberPtr[iFiber].m1[i][2]);
			}
		}
	}
	fflush(filePtr);
	fclose(filePtr);
}

//---------------------------------------------------------------------
void PrintForcesLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index) {

    unsigned  int		i, iFiber;
    char                fileName[256];
    FILE* filePtr;

    //----------------------------------------------
    sprintf(fileName, "force_%d.dat", index);
    filePtr = fopen(fileName, "w");
    fprintf(filePtr, "FIter Ffiber Fi Fx Fy Fz Mz\n");

    //----------------------------------------------
    for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
        for (i = 0; i < theFiberPtr[iFiber].n; i++) {
            fprintf(filePtr, "%d %d %d %e %e %e %e\n", theParameterPtr->iter, iFiber, i,
                theFiberPtr[iFiber].f[i][0], theFiberPtr[iFiber].f[i][1], theFiberPtr[iFiber].f[i][2],theFiberPtr[iFiber].moment[i][2]);
        }
    }
    fflush(filePtr);
    fclose(filePtr);
}

//---------------------------------------------------------------------
void Print_leLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index) {

    int					i, iFiber;
    char                fileName[256];
    FILE* filePtr;

    //----------------------------------------------
    sprintf(fileName, "le_%d.dat", index);
    filePtr = fopen(fileName, "w");
    fprintf(filePtr, "le_%d\n", index);

    //-----------------------------------------
    for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
        if (theFiberPtr[iFiber].status == STATUS_FREE) {
            for (i = 0; i < theFiberPtr[iFiber].n - 1; i++) {
                fprintf(filePtr, "%e \n", theFiberPtr[iFiber].lt[i] - theFiberPtr[iFiber].l0);
            }
        }
    }
    fflush(filePtr);
    fclose(filePtr);
}

//---------------------------------------------------------------------
void Print_thetaLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index) {

	int					i, iFiber;
	char                fileName[256];
	FILE* filePtr;

	//----------------------------------------------
	sprintf(fileName, "theta_%d.dat", index);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "ii_%d theta_%d\n", index, index);

	//-----------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			for (i = 0; i < theFiberPtr[iFiber].n - 1; i++) {
				fprintf(filePtr, "%d %e \n", i, theFiberPtr[iFiber].thetat[i][2]);
			}
		}
	}
	fflush(filePtr);
	fclose(filePtr);
}

//---------------------------------------------------------------------
double PrintBendingEnergyPerFiber(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		i, iFiber;
	double				aux = 0, kappa, energy = 0;
	static FILE* filePtr;
	char					fileName[256];


	//-------------------------------------------------------------------
	sprintf(fileName, "C:\\data\\bending_energy_%d_%d.txt", theParameterPtr->line, theParameterPtr->iter);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "iFiber_eb eb\n");

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			aux = 0;
			for (i = 1; i < theFiberPtr[iFiber].n - 1; i++) {
				kappa = KappaLib(theFiberPtr[iFiber], i, 1) - theFiberPtr[iFiber].kappa1_bar[i];
				aux += kappa * kappa;
				kappa = KappaLib(theFiberPtr[iFiber], i, 2) - theFiberPtr[iFiber].kappa2_bar[i];
				aux += kappa * kappa;
			}
			aux *= theFiberPtr[iFiber].bending;
			fprintf(filePtr, "%d %e\n", iFiber, aux / 2);
		}
		energy += aux;
	}
	fclose(filePtr);
	return energy / 2;
}