#include "main.h"

int							gIndexMovie, gStartIndexMovie, gStep, gSleep, gSave, gPrintContinuousPovray;
fiber                       * gFiber;
parameter                   *gParameterPtr;
char                        *g_PathName;
double						* r, * g, * b, *rPovRay, *gPovRay, *bPovRay;

void SetColors(void);

//----------------------------------------------------------------------
int main(int argc, char** argv) {
    
    FILE* filePtr;
    
    //-------------- init 
    glutInit(&argc, argv);
    InitGlut();
    gParameterPtr = (parameter*)malloc(sizeof(parameter));
    g_PathName = (char*)malloc(256 * sizeof(char));
    gParameterPtr->iter = -1;
    gPrintContinuousPovray = 0;
	gSave = 0;
	gStep = 1;
	gSleep = 10;
    SetColors();

    //------------- lecture movie si present
    filePtr = fopen("movie.txt", "r");
    if (filePtr == NULL)
        gIndexMovie = -1;
    else {
        gIndexMovie = 0;
        fscanf(filePtr, "%s", g_PathName);
        fscanf(filePtr, "%d", &(gParameterPtr->line));
        fscanf(filePtr, "%d", &(gStartIndexMovie));
        fscanf(filePtr, "%d", &(gStep));
		fscanf(filePtr, "%d", &(gSleep));
		fclose(filePtr); 
        printf("g_PathName = %s\n", g_PathName);
        printf("gParameterPtr->line = %u\n", gParameterPtr->line);
        printf("gStartMovie = %d\n", gStartIndexMovie);
        printf("gParameterPtr->gStep = %u\n", gStep);
        printf("gParameterPtr->gSleep = %u\n", gSleep);
        gIndexMovie = gStartIndexMovie;
    }
    printf("gIndexMovie = %d\n", gIndexMovie);

    DoOneIteration();

     //-------------- main loop
    glutIdleFunc(DoOneIteration);
    glutMainLoop();

    return 0;
}

//---------------------------------------
void SetColors(void) {

	int				iFiber;

	r = (double*)malloc(NCOLOR * sizeof(double));
	g = (double*)malloc(NCOLOR * sizeof(double));
	b = (double*)malloc(NCOLOR * sizeof(double));
	rPovRay = (double*)malloc(NCOLOR * sizeof(double));
	gPovRay = (double*)malloc(NCOLOR * sizeof(double));
	bPovRay = (double*)malloc(NCOLOR * sizeof(double));

	for (iFiber = 0; iFiber < NCOLOR; iFiber++) {
		rPovRay[iFiber] = (((double)iFiber) * 0.33 * 7. / 8.);
		rPovRay[iFiber] -= floor(rPovRay[iFiber]);
		gPovRay[iFiber] = 1.;
		bPovRay[iFiber] = 0.5 * MyRand();

		r[iFiber] = MyRand();
		g[iFiber] = MyRand();
		b[iFiber] = MyRand();
	}
}