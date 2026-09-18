#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define NCOLOR				30000

//---------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#include <winbase.h>
#include <math.h>
//#include "GL/glew.h"
#include "GL/freeglut.h"
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//---------------- core.c
void			DoOneIteration(void);

//----------------- glut.c
void			InitGlut(void);
void			SetColor(void);
void			Display(void);
void			PrintGlutText(void);
void			mouse(int button, int state, int x, int y);
void			keyboard(unsigned char key, int x, int y);
void			reshape(int w, int h);


void Print_PovRay(int line, int iter);