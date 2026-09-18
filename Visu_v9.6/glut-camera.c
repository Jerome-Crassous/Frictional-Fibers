#include "main.h"

#define SIZE				600
#define COLOR_NORMAL		0
#define COLOR_GRAY			1

extern fiber				* gFiber;
extern parameter			* gParameterPtr;
extern int gPrintContinuousPovray, gSave;
extern double* r, * g, * b;


double						distance,roll,pitch,heading,randomColor;
int							drawFixed,j,p;

GLubyte mipmapImage128[128][128][3]; 
GLubyte mipmapImage64[64][64][3];
GLubyte mipmapImage32[32][32][3];
GLubyte mipmapImage16[16][16][3];
GLubyte mipmapImage8[8][8][3];
GLubyte mipmapImage4[4][4][3];
GLubyte mipmapImage2[2][2][3];
GLubyte mipmapImage1[1][1][3];
GLuint tex;
GLUquadric* sphere, *cylinder;

void MakeImages(void);
void MakeImages_m1m2(void);
void MakeTextures(void);
void DrawOneScene(int flagColor);

//-----------------------------------------------------------------
void InitGlut(void)  {

	int			ok;
	
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(SIZE, SIZE);
	glutInitWindowPosition(10, 10);
	printf("before glutCreateWindows\n");
	ok=glutCreateWindow("");
	printf("glutCreateWindows returns ok %d\n", ok);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_FLAT);
	glutReshapeFunc(reshape);
	glutDisplayFunc(Display);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(DoOneIteration);

	MakeTextures();
	sphere = gluNewQuadric();
	cylinder = gluNewQuadric();

	gSave = 0;
	drawFixed = 1;

	roll=0, pitch=-80, heading=0, distance=700;
	roll = 0, pitch = 90., heading = 0, distance = 120;
	roll = 0, pitch = 0., heading = 0, distance = 150;
	roll = 0, pitch = 0., heading = 0, distance = 20;
//	roll = 0, pitch = 0., heading = 0, distance = 450;
}

//-----------------------------------------------------------------
void keyboard(unsigned char key, int x, int y){
	if (key == 'a') roll += 5.;
	if (key == 'z') roll -= 5.;
	if (key == 'q') pitch += 5.;
	if (key == 's') pitch -= 5.;
	if (key == 'w') heading += 5.;
	if (key == 'x') heading -= 5.;
	if (key == 'f') drawFixed = 1 - drawFixed;
	if (key == 't') MakeTextures();
	if (key == 'p') {
		gPrintContinuousPovray = 1- gPrintContinuousPovray;
		printf("gPrintContinuousPovray = %d\n", gPrintContinuousPovray);
	}
	if (key == 'o') {
		gSave = 1 - gSave;
		printf("gSave = %d\n", gSave);
	}
	if (key == 'r') randomColor = 1 - randomColor;
	printf("roll %f - pitch %f - heading %f - distance %f \n", roll, pitch, heading, distance);
}

//-----------------------------------------------------------------
void mouse(int button, int state, int x, int y) {
	if (button == 3)
		distance += 1.;
	else if (button == 4)
		distance -= 1.;
	if (distance < 1)
		distance = 1;
	printf("roll %f - pitch %f - heading %f - distance %f \n", roll, pitch, heading, distance);
}

//-----------------------------------------------------------------
void reshape(int w, int h){
	//double l = 80;
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50, 1, 0.1, 5000.);	
//	glOrtho(-l, l, -l, l, -200, 200);
}

//-----------------------------------------------------------------
void Display(void) {

	int			ix=0, iy=0, iz=0;
	double		lz;

	//--------------------------------
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//--------- rotation de la camera
	glTranslated(0.0, 0.0, -distance);
	glRotated(heading, 0, 0, 1);
	glRotated(pitch, 1, 0, 0);
	glRotated(roll, 0, 0, 1);

	//--------- origine + axes
	glLineWidth(3.0);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0), glVertex3f(0, 0, 0), glVertex3f(100, 0, 0);
	glColor3f(0, 1, 0), glVertex3f(0, 0, 0), glVertex3f(0, 100, 0);
	glColor3f(0, 0, 1), glVertex3f(0, 0, 0), glVertex3f(0, 0, 100);
	glEnd();
	glLineWidth(1.0);

	//-------------------- draws fibers
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tex);
	gluQuadricTexture(cylinder, GL_TRUE);
	gluQuadricDrawStyle(cylinder, GLU_FILL);
	gluQuadricTexture(sphere, GL_TRUE);
	gluQuadricDrawStyle(sphere, GLU_FILL);
	glMatrixMode(GL_MODELVIEW);

	//---------------------------------------	
	if (gParameterPtr->periodicVisu == PERIODIC_VISU_NO) {
		DrawOneScene(COLOR_NORMAL);
	}
	else if (gParameterPtr->periodicVisu == PERIODIC_VISU_XY){
		for (ix = -1; ix <= 1; ix++) {
			for (iy = -1; iy <= 1; iy++){
				glPushMatrix();
				glTranslatef(((double)ix) * gParameterPtr->lx,((double)iy) * gParameterPtr->ly,0.);
				DrawOneScene(((ix==0)&&(iy==0)&&(iz==0)) ? COLOR_NORMAL : COLOR_GRAY);
				glPopMatrix();
			}
		}
	}
	else if (gParameterPtr->periodicVisu == PERIODIC_VISU_ALPHA) {
		for (ix = -2; ix <= 2; ix++) {
			for (iy = -2; iy <= 2; iy++) {
				glPushMatrix();
				lz = gParameterPtr->lx / gParameterPtr->alphaX;
				glTranslatef(0., 0, lz);		//------ rotation autour de l'axe = 3 etapes
				glRotatef(((double)ix)* gParameterPtr->alphaX * 180. / PI, 0, 1., 0);
				glTranslatef(0., 0, -lz);
				glTranslatef(0., ((double)iy) * gParameterPtr->ly, 0.);	//--- translation simple y
				DrawOneScene(ix*iy*iz);
				glPopMatrix();
			}
		}
	}

	//----------------------------------------------------
	glDisable(GL_TEXTURE_2D);		
	PrintGlutText();
	glFlush();
	glutSwapBuffers();
}

//-----------------------------------------------------------------
void DrawOneScene(int flagColor) {

	int			i, iLoaded, iFiber;
	unsigned int	k;
	double		l0, m1[3], m2[3], m3[3], theta, c, s, u[3];

	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {

		if (gFiber[iFiber].status != STATUS_VIRTUAL) {

			//---------- couleur fibre
			/*if (flagColor == COLOR_NORMAL) {
				glColor3f(theColorArray[iFiber].r, theColorArray[iFiber].g, theColorArray[iFiber].b);
				printf("%d %e %e %e\n", iFiber,
					theColorArray[iFiber].r,
					theColorArray[iFiber].g,
					theColorArray[iFiber].b);
			}
			else
				glColor3f(0.8, 0.8, 0.8);*/
	//		glColor3d(MyRand(), MyRand(), MyRand());
			glColor3d(r[iFiber], g[iFiber], b[iFiber]);

			//--------- tracé fibre
			for (i = 0; i < gFiber[iFiber].n; i++) {
				
				glPushMatrix();

				//------------------ translation du cylindre
				glTranslatef(gFiber[iFiber].xt[i][0], gFiber[iFiber].xt[i][1], gFiber[iFiber].xt[i][2]);

				//------------------ rotation du cylindre
				for (k = 0; k < 3; k++) {
					if ((i == gFiber[iFiber].n - 1) && (gFiber[iFiber].n !=1)) 
						iLoaded = i - 1;
					else iLoaded = i;				
					m1[k] = gFiber[iFiber].m1[iLoaded][k];
					m3[k] = gFiber[iFiber].e[iLoaded][k];
				}
				CrossProduct(m3, m1, m2);

				c = (m1[0] + m2[1] + m3[2] - 1.)/2.;

				if (c >= 1) theta = 0.;
				else if (c <= -1) theta = 180.;
				else theta = 180. / PI * acos(c);
				s = sin(theta*PI/180.);

				if (fabs(s)> 0.001) {
					u[0] = (m2[2] - m3[1]) / 2./ s;
					u[1] = (m3[0] - m1[2]) / 2./ s;
					u[2] = (m1[1] - m2[0]) / 2./ s;
				}
				else if (c <  -0.999) {
					u[0] = sqrt((m1[0] + 1) / 2.);
					u[1] = sqrt((m2[1] + 1) / 2.);
					u[2] = sqrt((m3[2] + 1) / 2.);
				}
				glRotatef(theta, u[0], u[1], u[2]);

				//----------- dessine le cylindre ou cylindre + speher
				if (gFiber[iFiber].n != 1) {
					if (i < gFiber[iFiber].n - 1) {
						l0 = Distance(gFiber[iFiber].xt[i], gFiber[iFiber].xt[i + 1]);
						gluCylinder(cylinder, gFiber[iFiber].radius, gFiber[iFiber].radius, l0, 8, 8);
					}
					if ((i == 0) || (i == gFiber[iFiber].n - 1))
						gluSphere(sphere, gFiber[iFiber].radius, 16, 8);
				}
				else {
					gluSphere(sphere, gFiber[iFiber].radius, 16, 8);
				}

				//--------- restaure l'orientation
				glPopMatrix();
			}
		}
	}
}

//-----------------------------------------------------------------
void PrintGlutText(void){

	char		string[256];

	sprintf(string, "gIter=%d", gParameterPtr->iter);
	glutSetWindowTitle(string);
}

//-----------------------------------------------------------------
void MakeImages(void) {

	int i, j,color;


	for (i = 0; i < 128; i++) {
		for (j = 0; j < 128; j++) {
			for (color = 0; color < 3; color++) {
				if ((i < 32) || (i > 96)) mipmapImage128[i][j][color] = 0;
			/*	else if ((i < 16+64) || (i > 112+64)) mipmapImage128[i][j][color] = 0;
				else if ((i < 16 + 2*64) || (i > 112 + 2*64)) mipmapImage128[i][j][color] = 0;
				else if ((i < 16 + 3*64) || (i > 112 + 3*64)) mipmapImage128[i][j][color] = 0;*/
				else mipmapImage128[i][j][color] = 255;
			}
		}
	}

	for (i = 0; i < 64; i++) {
		for (j = 0; j < 64; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage64[i][j][color] =
					(mipmapImage128[2 * i][2 * j][color] + mipmapImage128[2 * i][2 * j + 1][color]
						+ mipmapImage128[2 * i + 1][2 * j][color] + mipmapImage128[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 32; i++) {
		for (j = 0; j < 32; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage32[i][j][color] =
					(mipmapImage64[2 * i][2 * j][color] + mipmapImage64[2 * i][2 * j + 1][color]
						+ mipmapImage64[2 * i + 1][2 * j][color] + mipmapImage64[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}
	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage16[i][j][color] =
					(mipmapImage32[2 * i][2 * j][color] + mipmapImage32[2 * i][2 * j + 1][color]
						+ mipmapImage32[2 * i + 1][2 * j][color] + mipmapImage32[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage8[i][j][color] =
					(mipmapImage16[2 * i][2 * j][color] + mipmapImage16[2 * i][2 * j + 1][color]
						+ mipmapImage16[2 * i + 1][2 * j][color] + mipmapImage16[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage4[i][j][color] =
					(mipmapImage8[2 * i][2 * j][color] + mipmapImage8[2 * i][2 * j + 1][color]
						+ mipmapImage8[2 * i + 1][2 * j][color] + mipmapImage8[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage2[i][j][color] =
					(mipmapImage4[2 * i][2 * j][color] + mipmapImage4[2 * i][2 * j + 1][color]
						+ mipmapImage4[2 * i + 1][2 * j][color] + mipmapImage4[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (color = 0; color < 3; color++) {
		mipmapImage1[i][j][color] =
			(mipmapImage2[2 * i][2 * j][color] + mipmapImage2[2 * i][2 * j + 1][color]
				+ mipmapImage2[2 * i + 1][2 * j][color] + mipmapImage2[2 * i + 1][2 * j + 1][color]) / 4.;
	}
}

//-----------------------------------------------------------------
void MakeImages_m1m2(void) {

	int i, j, color;

	for (i = 0; i < 128; i++) {
		for (j = 0; j < 128; j++) {
			for (color = 0; color < 3; color++) {
				if ((i < 16) || (i > 112)) mipmapImage128[i][j][color] = 0;			//debut fin
			//	else if ((j > 124) || (j < 4)) mipmapImage128[i][j][color] = 0;	//marque m1
			//	else if ((j > 28) && (j < 36)) mipmapImage128[i][j][color] = 0;	//marque m1
			//	else if ((j > 60) && (j < 68)) mipmapImage128[i][j][color] = 0;	//marque m1
				else if ((j > 92) && (j < 104)) mipmapImage128[i][j][color] = 0;	//marque m1
				else mipmapImage128[i][j][color] = 255;
			}
		}
	}

	for (i = 0; i < 64; i++) {
		for (j = 0; j < 64; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage64[i][j][color] =
					(mipmapImage128[2 * i][2 * j][color] + mipmapImage128[2 * i][2 * j + 1][color]
						+ mipmapImage128[2 * i + 1][2 * j][color] + mipmapImage128[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 32; i++) {
		for (j = 0; j < 32; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage32[i][j][color] =
					(mipmapImage64[2 * i][2 * j][color] + mipmapImage64[2 * i][2 * j + 1][color]
						+ mipmapImage64[2 * i + 1][2 * j][color] + mipmapImage64[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage16[i][j][color] =
					(mipmapImage32[2 * i][2 * j][color] + mipmapImage32[2 * i][2 * j + 1][color]
						+ mipmapImage32[2 * i + 1][2 * j][color] + mipmapImage32[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage8[i][j][color] =
					(mipmapImage16[2 * i][2 * j][color] + mipmapImage16[2 * i][2 * j + 1][color]
						+ mipmapImage16[2 * i + 1][2 * j][color] + mipmapImage16[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage4[i][j][color] =
					(mipmapImage8[2 * i][2 * j][color] + mipmapImage8[2 * i][2 * j + 1][color]
						+ mipmapImage8[2 * i + 1][2 * j][color] + mipmapImage8[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			for (color = 0; color < 3; color++) {
				mipmapImage2[i][j][color] =
					(mipmapImage4[2 * i][2 * j][color] + mipmapImage4[2 * i][2 * j + 1][color]
						+ mipmapImage4[2 * i + 1][2 * j][color] + mipmapImage4[2 * i + 1][2 * j + 1][color]) / 4.;
			}
		}
	}

	for (color = 0; color < 3; color++) {
		mipmapImage1[i][j][color] =
			(mipmapImage2[2 * i][2 * j][color] + mipmapImage2[2 * i][2 * j + 1][color]
				+ mipmapImage2[2 * i + 1][2 * j][color] + mipmapImage2[2 * i + 1][2 * j + 1][color]) / 4.;
	}
}

//-----------------------------------------------------------------
void MakeTextures(void){
	static int flag = 0;

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glShadeModel(GL_FLAT);

	if (flag == 1)
		MakeImages();
	else
		MakeImages_m1m2();
	flag = 1 - flag;

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 128, 128, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage128[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 1, 3, 64, 64, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage64[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 2, 3, 32, 32, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage32[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 3, 3, 16, 16, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage16[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 4, 3, 8, 8, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage8[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 5, 3, 4, 4, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage4[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 6, 3, 2, 2, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage2[0][0][0]);
	glTexImage2D(GL_TEXTURE_2D, 7, 3, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, &mipmapImage1[0][0][0]);/**/
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_LOD, 0);
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_LOD, max_lod);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glEnable(GL_TEXTURE_2D);
}