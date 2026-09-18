#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//-------------------------------------------------------
void SaveConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n) {

	unsigned int		iFiber, i, j;
	FILE* filePtr;
	char				fileName[256];

	//---------------------------------------------
	if (n < 0) sprintf(fileName, "config.dat");
	else sprintf(fileName, "%sconfig_%d_%d.tmp", pathName, theParameterPtr->line, n);

	filePtr = fopen(fileName, "wb");
	if (filePtr == NULL) {
		printf("### opening %s failed\n", fileName);
		return;
	}

	//----------------------------------------------------------------
	fwrite(theParameterPtr, sizeof(parameter), 1, filePtr);

	//--------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {

		fwrite(&(theFiberPtr[iFiber].n), sizeof(int), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].status), sizeof(int), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].radius), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].masse), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].l0), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].k0), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].bending), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].j), sizeof(double), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].c), sizeof(double), 1, filePtr);

		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				fwrite(&(theFiberPtr[iFiber].xt[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].xtm[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].thetat[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].thetatm[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].e[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].f[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].moment[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].m1_bar[i][j]), sizeof(double), 1, filePtr);
				fwrite(&(theFiberPtr[iFiber].m1[i][j]), sizeof(double), 1, filePtr);
			}
			fwrite(&(theFiberPtr[iFiber].lt[i]), sizeof(double), 1, filePtr);
			fwrite(&(theFiberPtr[iFiber].ltm[i]), sizeof(double), 1, filePtr);;
			fwrite(&(theFiberPtr[iFiber].kappa1_bar[i]), sizeof(double), 1, filePtr);
			fwrite(&(theFiberPtr[iFiber].kappa2_bar[i]), sizeof(double), 1, filePtr);
		}
	}
	fclose(filePtr);
}

//-------------------------------------------------------
void ReadConfigLib(fiber** theFiberPtr, parameter* theParameterPtr, char* pathName, int n) {

	unsigned int				iFiber, i, j, nSegment;
	static int					inited = 0;	
	char						fileName[256];
	FILE						* filePtr;

	//---------------------------------------------
	if (n < 0) sprintf(fileName, "config.dat");
	else sprintf(fileName, "%sconfig_%d_%d.tmp", pathName, theParameterPtr->line, n);

	filePtr = fopen(fileName, "rb");
	if (filePtr == NULL) {
		printf("### ReadConfigLib opening %s failed\n", fileName);
		return;
	}
	else
		printf("ReadConfigLib opening %s ok\n", fileName);

	//----------------------------------------------------------------
	fread(theParameterPtr, sizeof(parameter), 1, filePtr);

	//-----------------------------------------------------------------
	if (*theFiberPtr == NULL)			//the fiber does not exist
		*theFiberPtr = (fiber*)malloc((theParameterPtr->nFiber) * sizeof(fiber));

	//--------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {

		//--------------------------------------------------
		fread(&nSegment, sizeof(int), 1, filePtr);
		if (inited == 0)
			AllocateOneFiberLib(*theFiberPtr, iFiber, nSegment);

		fread(&((*theFiberPtr)[iFiber].status), sizeof(int), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].radius), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].masse), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].l0), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].k0), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].bending), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].j), sizeof(double), 1, filePtr);
		fread(&((*theFiberPtr)[iFiber].c), sizeof(double), 1, filePtr);

		for (i = 0; i < (*theFiberPtr)[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				fread(&((*theFiberPtr)[iFiber].xt[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].xtm[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].thetat[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].thetatm[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].e[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].f[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].moment[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].m1_bar[i][j]), sizeof(double), 1, filePtr);
				fread(&((*theFiberPtr)[iFiber].m1[i][j]), sizeof(double), 1, filePtr);
			}
			fread(&((*theFiberPtr)[iFiber].lt[i]), sizeof(double), 1, filePtr);
			fread(&((*theFiberPtr)[iFiber].ltm[i]), sizeof(double), 1, filePtr);;
			fread(&((*theFiberPtr)[iFiber].kappa1_bar[i]), sizeof(double), 1, filePtr);
			fread(&((*theFiberPtr)[iFiber].kappa2_bar[i]), sizeof(double), 1, filePtr);
		}
	}
	fclose(filePtr);
	inited = 1;
}

//-------------------------------------------------------
void SaveReducedConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n) {

	unsigned int	iFiber, i, j;
	FILE			* filePtr;
	char			fileName[256];
	float			aux;
	
	//---------------------------------------------
	if (n < 0) sprintf(fileName, "reduced_Config.tmp");
	else sprintf(fileName, "%sreduced_Config_%d_%d.tmp", pathName,theParameterPtr->line, n);

	filePtr = fopen(fileName, "wb");
	if (filePtr == NULL) {
		printf("### opening %s failed\n", fileName);
		return;
	}

	//----------------------------------------------------------------
	fwrite(theParameterPtr, sizeof(parameter), 1, filePtr);

	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {

		fwrite(&(theFiberPtr[iFiber].n), sizeof(int), 1, filePtr);
		fwrite(&(theFiberPtr[iFiber].status), sizeof(int), 1, filePtr);
		aux = (float)theFiberPtr[iFiber].radius;
		fwrite(&aux, sizeof(float), 1, filePtr);

		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				aux = (float)theFiberPtr[iFiber].xt[i][j];
				fwrite(&aux, sizeof(float), 1, filePtr);
				aux = (float)theFiberPtr[iFiber].e[i][j];
				fwrite(&aux, sizeof(float), 1, filePtr);
				aux = (float)theFiberPtr[iFiber].thetat[i][j];
				fwrite(&aux, sizeof(float), 1, filePtr);
				aux = (float)theFiberPtr[iFiber].m1[i][j];
				fwrite(&aux, sizeof(float), 1, filePtr);
			}
		}
	}
	fflush(filePtr);
	fclose(filePtr);
}

//-------------------------------------------------------
int ReadReducedConfigLib(fiber** theFiberPtr, parameter* theParameterPtr, char* pathName, int n) {

	unsigned int	iFiber, i, j, iterOld,nSegment ;
	int				ok;
	FILE			* filePtr;
	float			aux;
	static int		inited = 0;
	char			fileName[100];

	//---------------------------------------------
	if (n < 0) sprintf(fileName, "reduced_Config.tmp");
	else sprintf(fileName, "%sreduced_Config_%d_%d.tmp", pathName,theParameterPtr->line, n);

	filePtr = fopen(fileName, "rb");
	if (filePtr == NULL) {
		printf("### opening %s failed\n", fileName);
		return 0;
	}

	//----------------------------------------------------------------
	iterOld = theParameterPtr->iter;
	ok = fread(theParameterPtr, sizeof(parameter), 1, filePtr);

	if (theParameterPtr->iter == iterOld) {
		fclose(filePtr);
		return 1;
	}

	//------------------- allocation du tableau de fibre
	if (*theFiberPtr == NULL) {
		*theFiberPtr = (fiber*)malloc(theParameterPtr->nFiber * sizeof(fiber));
		if (theFiberPtr != NULL)
			printf("### allocation theFiberPtr ok in ReadReducedConfigLib\n");
		else {
			printf("### allocation theFiberPtr failed\n");
			exit(1253);
		}
	}

	//--------------------------------------------
	for (iFiber=0; iFiber< theParameterPtr->nFiber; iFiber++) {

		ok = fread(&nSegment, sizeof(unsigned int), 1, filePtr);

		if (inited == 0)
			AllocateOneFiberLib(*theFiberPtr, iFiber, nSegment);

		fread(&((*theFiberPtr)[iFiber].status), sizeof(int), 1, filePtr);
		fread(&aux, sizeof(float), 1, filePtr);
		(*theFiberPtr)[iFiber].radius = (double) aux;

		for (i = 0; i < (*theFiberPtr)[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				fread(&aux, sizeof(float), 1, filePtr);
				(*theFiberPtr)[iFiber].xt[i][j] = (double) aux;
				fread(&aux, sizeof(float), 1, filePtr);
				(*theFiberPtr)[iFiber].e[i][j] = (double)aux;
				fread(&aux, sizeof(float), 1, filePtr);
				(*theFiberPtr)[iFiber].thetat[i][j] = (double) aux;
				fread(&aux, sizeof(float), 1, filePtr);
				(*theFiberPtr)[iFiber].m1[i][j] = (double) aux;
			}
		}
	//	CalculateLeLib(&((*theFiberPtr)[iFiber]));
	}
	fclose(filePtr);
	inited = 1;
	return 2;
}

//-------------------------------------------------------
void PrintReducedConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n) {

	unsigned int	iFiber, i;
	FILE			* filePtr;
	char			fileName[256];

	//---------------------------------------------
	sprintf(fileName, "check_Reduced_Config.txt");

	filePtr = fopen(fileName, "w");
	if (filePtr == NULL) {
		printf("### opening %s failed\n", fileName);
		return;
	}

	//----------------------------------------------------------------
	fprintf(filePtr, "//--------- parameters\n");
	fprintf(filePtr, "theParameterPtr->iter %d\n", theParameterPtr->iter);
	fprintf(filePtr, "theParameterPtr->line %d\n", theParameterPtr->line);
	fprintf(filePtr, "theParameterPtr->periodic %d\n", theParameterPtr->periodic);
	fprintf(filePtr, "theParameterPtr->nFiber %d\n", theParameterPtr->nFiber);
	fprintf(filePtr, "theParameterPtr->nSegment %d\n", theParameterPtr->nSegment);
	fprintf(filePtr, "theParameterPtr->nContactMax %d\n", theParameterPtr->nContactMax);
	
	//----------------------------------------------------------------
	fprintf(filePtr,"//--------- fibres\n");
	fprintf(filePtr, "iFiber n ii xx yy zz ex ey ez thetax thetay thetaz m1x m1y m1z\n");

	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			fprintf(filePtr, "%d %d %d ", iFiber, theFiberPtr[iFiber].n, i);
			fprintf(filePtr, "%.2e %.2e %.2e -- %.2e %.2e %.2e -- %.2e %.2e %.2e -- %.2e %.2e %.2e\n",
				theFiberPtr[iFiber].xt[i][0], theFiberPtr[iFiber].xt[i][1], theFiberPtr[iFiber].xt[i][2],
				theFiberPtr[iFiber].e[i][0], theFiberPtr[iFiber].e[i][1], theFiberPtr[iFiber].e[i][2],
				theFiberPtr[iFiber].thetat[i][0], theFiberPtr[iFiber].thetat[i][1], theFiberPtr[iFiber].thetat[i][2],
				theFiberPtr[iFiber].m1[i][0], theFiberPtr[iFiber].m1[i][1], theFiberPtr[iFiber].m1[i][2]);
		}
	}

	//----------------------------------------------------------------
	fflush(filePtr);
	fclose(filePtr);
}