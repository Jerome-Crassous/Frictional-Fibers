#include	"main.h"

extern fiber* gFiber;
extern parameter* gParameterPtr;

double shift[3] = { 0.,0.,0 };
extern double * rPovRay, * gPovRay, * bPovRay;

void PrintOnePovRayScene(FILE* filePtr);

//-----------------------------------------------------------------
void Print_PovRay(int line, int iter){

	FILE			*filePtr;
	char			fileName[256];
	int				i, iLoaded, iFiber;
	double			pointOnAxis[3] = { 0.,0. ,gParameterPtr->lx / gParameterPtr->alpha };
	double			direction[3] = { 0.,1.,0. };
	double			angle = gParameterPtr->alpha;
	double			dx, dy;

	printf("input Print_PovRay\n");

	//---------------------------------------------
	if((line<0)||(iter<0))
		sprintf(fileName, "C:\\data\\PovRay\\fiber_povray.txt");
	else
		sprintf(fileName, "C:\\data\\PovRay\\fiber_povray_%d_%d.txt", line,iter);

	filePtr = fopen(fileName, "w");
	if (!filePtr)
		printf("Opening %s failed in PrintPovRay\n", fileName);
	else
		printf("Print_PovRay file ok\n");

	PrintOnePovRayScene(filePtr);
	
/*	for (dx = -2;dx < 3;dx++) {
		for (dy = -2;dy < 3;dy++) {

			for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {
				for (i = 0; i < gFiber[iFiber].n; i++) {
					Rotate_Point_Around_Axis(gFiber[iFiber].xt[i],
						pointOnAxis, direction, dx * angle);
					gFiber[iFiber].xt[i][1] += dy * gParameterPtr->ly;
				}
			}
			PrintOnePovRayScene(filePtr);
			for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {
				for (i = 0; i < gFiber[iFiber].n; i++) {
					gFiber[iFiber].xt[i][1] -= dy * gParameterPtr->ly;	
					Rotate_Point_Around_Axis(gFiber[iFiber].xt[i],
						pointOnAxis, direction, -dx * angle);
				}
			}
		}
	}*/

	fflush(filePtr);
	fclose(filePtr);
	exit(5);
}

//---------------------------------------------
void PrintOnePovRayScene(FILE* filePtr) {

	int				iFiber, i;

	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {

		if (gFiber[iFiber].status != STATUS_VIRTUAL) {
			for (i = 0; i < gFiber[iFiber].n - 1; i++) {
				fprintf(filePtr, "cylinder{\n<%f,%f,%f>\n<%f,%f,%f>\n %f\n pigment { hsl2rgb(<%f,%f,%f>) }\n }\n",
					gFiber[iFiber].xt[i][0], gFiber[iFiber].xt[i][1], gFiber[iFiber].xt[i][2],
					gFiber[iFiber].xt[i + 1][0], gFiber[iFiber].xt[i + 1][1], gFiber[iFiber].xt[i + 1][2],
					gFiber[iFiber].radius, rPovRay[iFiber], gPovRay[iFiber], bPovRay[iFiber]);
			}


			for (i = 0; i < gFiber[iFiber].n; i++) {
				fprintf(filePtr, "sphere{\n<%f,%f,%f>\n %f\n pigment  {hsl2rgb(<%f,%f,%f>) }\n }\n",
					gFiber[iFiber].xt[i][0], gFiber[iFiber].xt[i][1], gFiber[iFiber].xt[i][2],
					gFiber[iFiber].radius, rPovRay[iFiber], gPovRay[iFiber], bPovRay[iFiber]);
			}
		}
	}
}