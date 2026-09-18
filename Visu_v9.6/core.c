#include	"main.h"

extern int							gIndexMovie, gStartIndexMovie,  gStep, gSleep, gSave;
extern char							*g_PathName;
extern fiber						*gFiber;
extern parameter					*gParameterPtr;
extern int gPrintContinuousPovray;

//-----------------------------------------------------------------
void DoOneIteration(void) {

	int				ok;

	//-----------------------------------
	Sleep(gSleep);

	if (gIndexMovie < 0) {
		ReadReducedConfigLib(&gFiber, gParameterPtr, NULL, -1);
		PrintReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
		if (gPrintContinuousPovray == 1)
			Print_PovRay(gParameterPtr->line, gParameterPtr->iter);
	}
	else {
		ok=ReadReducedConfigLib(&gFiber, gParameterPtr, g_PathName, gIndexMovie);

		if (gSave == 1)
			PrintPositionsLib(gFiber, gParameterPtr, g_PathName, gIndexMovie);
		if (gPrintContinuousPovray == 1)
			Print_PovRay(gParameterPtr->line, gParameterPtr->iter);

		if (ok == 0)
			gIndexMovie = gStartIndexMovie;
		else
			gIndexMovie += gStep;
	}

	glutPostRedisplay();
}