//-----------------------------------------------------
__attribute__((always_inline))  uint NSegment(__global uint *n , uint fiber, uint i) {

	uint				fiberLoc=0, aux=0;

	for (fiberLoc = 0; fiberLoc < fiber; fiberLoc++)
		aux += n[fiberLoc];
	return aux + i;
}

//-----------------------------------------------------
void __attribute__((always_inline)) ComputeShift(uint periodic, uint iShift, float dx, float dy, float dz, float3* shift) {
	
	if (periodic == PERIODIC_NO) {
		if (iShift == 0)
			(*shift)[0] = 0., (*shift)[1] = 0., (*shift)[2] = 0.;
		else
			printf("#### error in ComputeShift PERIODIC_NO %u %u\n", periodic, iShift);
		return;
	}
	else if (periodic == PERIODIC_XY) {
		if (iShift == 0) { (*shift)[0] = 0., (*shift)[1] = 0., (*shift)[2] = 0.; }
		else if (iShift == 1) { (*shift)[0] = dx, (*shift)[1] = 0., (*shift)[2] = 0.; }
		else if (iShift == 2) { (*shift)[0] = dx, (*shift)[1] = dy, (*shift)[2] = 0.; }
		else if (iShift == 3) { (*shift)[0] = 0., (*shift)[1] = dy, (*shift)[2] = 0.; }
		else if (iShift == 4) { (*shift)[0] = -dx, (*shift)[1] = dy, (*shift)[2] = 0.; }
		else if (iShift == 5) { (*shift)[0] = -dx, (*shift)[1] = 0., (*shift)[2] = 0.; }
		else if (iShift == 6) { (*shift)[0] = -dx, (*shift)[1] = -dy, (*shift)[2] = 0.; }
		else if (iShift == 7) { (*shift)[0] = 0., (*shift)[1] = -dy, (*shift)[2] = 0.; }
		else if (iShift == 8) { (*shift)[0] = dx, (*shift)[1] = -dy, (*shift)[2] = 0.; }
		else printf("#### error in ComputeShift PERIODIC_XY %u %u\n", periodic, iShift);
		return;
	}
}

//-----------------------------------------------------
double3 rotate_point_around_axis(double3 P, double3 A, double3 u, double alpha){
	//P le point, A un point sur l'axe, u le directeur normalisé, alpha l'angle
	double3 p = P - A;
    double cosA = cos(alpha);
	double sinA = sin(alpha);
	double3 p_rot = p * cosA + cross(u, p) * sinA + u * dot(u, p) * (1.0f - cosA);
    return A + p_rot;
}

//-----------------------------------------------------
inline float3 parallel_transport(float3 x, float3 a, float3 b) {
	//----- parallel transport of x (a is transported in b) 
	//----- suppose a \simeq b et ||a||=||b||=1
	return x + cross(cross(a, b), x);
}

//-----------------------------------------------------
inline float3 rotate(float3 m1_bar, float3 e, float theta) {
	//----- rotate m1_bar of an angle theta along e
	//----- suppose ||e||=1, et e.m1_bar=0
	return cos(theta) * m1_bar + sin(theta) * cross(e, m1_bar);
}