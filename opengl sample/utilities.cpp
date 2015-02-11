#include "utilities.h"
#include <stack>
#include <map>
#include <algorithm>
#include <string>
#include <fstream>
#include <gdiplus.h>

using namespace Gdiplus;
using namespace std;


double PerformanceCounter::gTickFrequency = 0;

float MTH_Round(float value) {
	int intVal = (int)value;
	if( value - intVal >= 0.5f ) {
		return (float)(intVal+1);
	} else if( value - intVal <= -0.5f ) {
		return (float)(intVal-1);
	} else {
		return (float)intVal;
	}
}

double MTH_Round(double value) {
	int intVal = (int)value;
	if( value - intVal >= .5 ) {
		return (double)(intVal+1);
	} else if( value - intVal <= -0.5 ) {
		return (double)(intVal-1);
	} else {
		return (double)intVal;
	}
}


double RND_DoubleInRange(double min, double max) {
	return (rand() / (double)RAND_MAX) * (max-min) + min;
}
int RND_IntInRange(int min, int max) {
	return (rand() % (max-min+1)) + min;
}

TransformMatrix RND_MakeTransform(double minx, double miny, double minz, double maxx, double maxy, double maxz) {
	WorldPt3 x = RND_GenerateNormalizedVector();

	WorldPt3 y = RND_GenerateNormalizedVector();

	WorldPt3 z = x.Cross(y).Normal();
	y = z.Cross(x).Normal();

	TransformMatrix mat;
	mat.R0() = x;
	mat.R1() = y;
	mat.R2() = z;

	mat.v1.xOff = RND_DoubleInRange(minx, maxx);
	mat.v1.yOff = RND_DoubleInRange(miny, maxy);
	mat.v1.zOff = RND_DoubleInRange(minz, maxz);

	return mat;
}



void VEC3_Equate(double * vec3, double * output_vec3) {
	output_vec3[0] = vec3[0];
	output_vec3[1] = vec3[1];
	output_vec3[2] = vec3[2];
}
void VEC3_Negate(double * vec3) {
	vec3[0] = -vec3[0];
	vec3[1] = -vec3[1];
	vec3[2] = -vec3[2];
}
double VEC3_Dot(double * first_vec3, double * second_vec3) {
	return first_vec3[0]*second_vec3[0] + first_vec3[1]*second_vec3[1] + first_vec3[2]*second_vec3[2];
}


double VEC3_Magnitude(double * vec3) {
	return sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2]);
}
void VEC3_Normalize(double * vec3, double * normalized_vec3) {
	double magnitude = VEC3_Magnitude(vec3);
	if( magnitude == 0 ) {
		magnitude = 1;
	}

	normalized_vec3[0] = vec3[0] / magnitude;
	normalized_vec3[1] = vec3[1] / magnitude;
	normalized_vec3[2] = vec3[2] / magnitude;
}

void VEC3_Normalize(double * vec3) {
	VEC3_Normalize(vec3, vec3);
}

void VEC3_Cross(double * a_vec3, double * b_vec3, double * crossed_vec3) {
	crossed_vec3[0] = a_vec3[1]*b_vec3[2] - a_vec3[2]*b_vec3[1];
	crossed_vec3[1] = -(a_vec3[0]*b_vec3[2] - a_vec3[2]*b_vec3[0]);
	crossed_vec3[2] = a_vec3[0]*b_vec3[1] - a_vec3[1]*b_vec3[0];
}

void VEC3_Difference(double * a_vec3, double * b_vec3, double * difference_vec3) {
	difference_vec3[0] = b_vec3[0] - a_vec3[0];
	difference_vec3[1] = b_vec3[1] - a_vec3[1];
	difference_vec3[2] = b_vec3[2] - a_vec3[2];
}

void VEC3_Swap(double * a_vec3, double * b_vec3) {
	double vec3[3];
	vec3[0] = a_vec3[0];
	vec3[1] = a_vec3[1];
	vec3[2] = a_vec3[2];

	a_vec3[0] = b_vec3[0];
	a_vec3[1] = b_vec3[1];
	a_vec3[2] = b_vec3[2];

	b_vec3[0] = vec3[0];
	b_vec3[1] = vec3[1];
	b_vec3[2] = vec3[2];
}

void VEC3_Add(double * first_vec3, double * second_vec3, double * output_vec3) {
	output_vec3[0] = first_vec3[0] + second_vec3[0];
	output_vec3[1] = first_vec3[1] + second_vec3[1];
	output_vec3[2] = first_vec3[2] + second_vec3[2];
}

void VEC3_Subtract(double * first_vec3, double * second_vec3, double * output_vec3) {
	output_vec3[0] = first_vec3[0] - second_vec3[0];
	output_vec3[1] = first_vec3[1] - second_vec3[1];
	output_vec3[2] = first_vec3[2] - second_vec3[2];
}

void VEC3_Multiply(double * first_vec3, double * second_vec3, double * output_vec3) {
	output_vec3[0] = first_vec3[0] * second_vec3[0];
	output_vec3[1] = first_vec3[1] * second_vec3[1];
	output_vec3[2] = first_vec3[2] * second_vec3[2];
}

void VEC3_Multiply(double * vec3, double scalar, double * output_vec3) {
	output_vec3[0] = vec3[0]*scalar;
	output_vec3[1] = vec3[1]*scalar;
	output_vec3[2] = vec3[2]*scalar;
}
double VEC3_AngleBetween(double * first_vec3, double * second_vec3) {
	double mag1 = VEC3_Magnitude(first_vec3);

	double mag2 = VEC3_Magnitude(second_vec3);

	double dot = VEC3_Dot(first_vec3, second_vec3);

	double angleBetween = acos(max(min(dot / (mag1*mag2), 1.0), -1.0) );

	return angleBetween;
}
void MAT4_ExtractCol3(double * mat4, int col, double * vec3) {
    vec3[0] = mat4[col*4+0];
    vec3[1] = mat4[col*4+1];
    vec3[2] = mat4[col*4+2];
}

void MAT4_ExtractRow3(double * mat4, int row, double * vec3) {
	vec3[0] = mat4[row];
	vec3[1] = mat4[row+4];
	vec3[2] = mat4[row+8];
}
void MAT4_ExtractFrustumPoints(float * inVPInvertMat, float outFrustumPoints[3*8])
// frustum points are coordinates in world space in this order: NearBottomLeft, NTopL, NTRight, NBR, FarBL, FTL, FTR, FBR
{
	float ndcFrustumPts[] = {
		-1,-1,-1, 1,  -1, 1,-1, 1,  1, 1,-1, 1,  1,-1,-1,1, // near
		-1,-1, 1, 1,  -1, 1, 1, 1,  1, 1, 1, 1,  1,-1, 1, 1 // far
	};

	for( int i=0; i<8; i++ ) {
		float p[4];
		MAT4_VEC4_Multiply(inVPInvertMat, &ndcFrustumPts[i*4], p);
		VEC3_Divide(p, p[3], &outFrustumPoints[i*3]);
	}

}

void MAT4_SetCol3(int col, double * vec3, double * mat4) {
	mat4[col*4+0] = vec3[0];
	mat4[col*4+1] = vec3[1];
	mat4[col*4+2] = vec3[2];
}

void MAT4_LoadIdentity(double * mat4) {
	mat4[0] = 1;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 1;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = 1;
	mat4[11] = 0;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = 0;
	mat4[15] = 1;
}

void MAT4_MatrixVectorMultiply(double * mat4, double * column_vec3, double * row_vec3) {
	row_vec3[0] = mat4[ 0]*column_vec3[ 0] + mat4[ 4]*column_vec3[ 1] + mat4[ 8]*column_vec3[ 2] + mat4[12]*1;
	row_vec3[1] = mat4[ 1]*column_vec3[ 0] + mat4[ 5]*column_vec3[ 1] + mat4[ 9]*column_vec3[ 2] + mat4[13]*1;
	row_vec3[2] = mat4[ 2]*column_vec3[ 0] + mat4[ 6]*column_vec3[ 1] + mat4[10]*column_vec3[ 2] + mat4[14]*1;
}
void MAT4_VectorMatrixMultiply(double * row_vec3, double * mat4, double * column_vec3) {
	column_vec3[0] = row_vec3[ 0]*mat4[ 0] + row_vec3[ 1]*mat4[ 1] + row_vec3[ 2]*mat4[ 2] + 1*mat4[ 3];
	column_vec3[1] = row_vec3[ 0]*mat4[ 4] + row_vec3[ 1]*mat4[ 5] + row_vec3[ 2]*mat4[ 6] + 1*mat4[ 7];
	column_vec3[2] = row_vec3[ 0]*mat4[ 8] + row_vec3[ 1]*mat4[ 9] + row_vec3[ 2]*mat4[10] + 1*mat4[11];
}
void MAT4_VEC4_MatrixVectorMultiply(double * mat4, double * column_vec4, double * out_vec4) {
	out_vec4[0] = mat4[ 0]*column_vec4[ 0] + mat4[ 4]*column_vec4[ 1] + mat4[ 8]*column_vec4[ 2] + mat4[12]*column_vec4[ 3];
	out_vec4[1] = mat4[ 1]*column_vec4[ 0] + mat4[ 5]*column_vec4[ 1] + mat4[ 9]*column_vec4[ 2] + mat4[13]*column_vec4[ 3];
	out_vec4[2] = mat4[ 2]*column_vec4[ 0] + mat4[ 6]*column_vec4[ 1] + mat4[10]*column_vec4[ 2] + mat4[14]*column_vec4[ 3];
	out_vec4[3] = mat4[ 3]*column_vec4[ 0] + mat4[ 7]*column_vec4[ 1] + mat4[11]*column_vec4[ 2] + mat4[15]*column_vec4[ 3];
}

void MAT4_Multiply(double * first_mat4, double * second_mat4, double * output_mat4) {
	double result[16];

	result[ 0] = first_mat4[ 0]*second_mat4[ 0] + first_mat4[ 4]*second_mat4[ 1] + first_mat4[ 8]*second_mat4[ 2] + first_mat4[12]*second_mat4[ 3];
	result[ 1] = first_mat4[ 1]*second_mat4[ 0] + first_mat4[ 5]*second_mat4[ 1] + first_mat4[ 9]*second_mat4[ 2] + first_mat4[13]*second_mat4[ 3];
	result[ 2] = first_mat4[ 2]*second_mat4[ 0] + first_mat4[ 6]*second_mat4[ 1] + first_mat4[10]*second_mat4[ 2] + first_mat4[14]*second_mat4[ 3];
	result[ 3] = first_mat4[ 3]*second_mat4[ 0] + first_mat4[ 7]*second_mat4[ 1] + first_mat4[11]*second_mat4[ 2] + first_mat4[15]*second_mat4[ 3];

	result[ 4] = first_mat4[ 0]*second_mat4[ 4] + first_mat4[ 4]*second_mat4[ 5] + first_mat4[ 8]*second_mat4[ 6] + first_mat4[12]*second_mat4[ 7];
	result[ 5] = first_mat4[ 1]*second_mat4[ 4] + first_mat4[ 5]*second_mat4[ 5] + first_mat4[ 9]*second_mat4[ 6] + first_mat4[13]*second_mat4[ 7];
	result[ 6] = first_mat4[ 2]*second_mat4[ 4] + first_mat4[ 6]*second_mat4[ 5] + first_mat4[10]*second_mat4[ 6] + first_mat4[14]*second_mat4[ 7];
	result[ 7] = first_mat4[ 3]*second_mat4[ 4] + first_mat4[ 7]*second_mat4[ 5] + first_mat4[11]*second_mat4[ 6] + first_mat4[15]*second_mat4[ 7];

	result[ 8] = first_mat4[ 0]*second_mat4[ 8] + first_mat4[ 4]*second_mat4[ 9] + first_mat4[ 8]*second_mat4[10] + first_mat4[12]*second_mat4[11];
	result[ 9] = first_mat4[ 1]*second_mat4[ 8] + first_mat4[ 5]*second_mat4[ 9] + first_mat4[ 9]*second_mat4[10] + first_mat4[13]*second_mat4[11];
	result[10] = first_mat4[ 2]*second_mat4[ 8] + first_mat4[ 6]*second_mat4[ 9] + first_mat4[10]*second_mat4[10] + first_mat4[14]*second_mat4[11];
	result[11] = first_mat4[ 3]*second_mat4[ 8] + first_mat4[ 7]*second_mat4[ 9] + first_mat4[11]*second_mat4[10] + first_mat4[15]*second_mat4[11];

	result[12] = first_mat4[ 0]*second_mat4[12] + first_mat4[ 4]*second_mat4[13] + first_mat4[ 8]*second_mat4[14] + first_mat4[12]*second_mat4[15];
	result[13] = first_mat4[ 1]*second_mat4[12] + first_mat4[ 5]*second_mat4[13] + first_mat4[ 9]*second_mat4[14] + first_mat4[13]*second_mat4[15];
	result[14] = first_mat4[ 2]*second_mat4[12] + first_mat4[ 6]*second_mat4[13] + first_mat4[10]*second_mat4[14] + first_mat4[14]*second_mat4[15];
	result[15] = first_mat4[ 3]*second_mat4[12] + first_mat4[ 7]*second_mat4[13] + first_mat4[11]*second_mat4[14] + first_mat4[15]*second_mat4[15];

	memcpy(output_mat4, result, sizeof(double)*16);
}

void MAT4_TransposeUpper3(double * mat4, double * transpose_mat4) {
	transpose_mat4[0] = mat4[0];
	transpose_mat4[1] = mat4[4];
	transpose_mat4[2] = mat4[8];
	transpose_mat4[3] = mat4[3];
	transpose_mat4[4] = mat4[1];
	transpose_mat4[5] = mat4[5];
	transpose_mat4[6] = mat4[9];
	transpose_mat4[7] = mat4[7];
	transpose_mat4[8] = mat4[2];
	transpose_mat4[9] = mat4[6];
	transpose_mat4[10] = mat4[10];
	transpose_mat4[11] = mat4[11];
	transpose_mat4[12] = mat4[3];
	transpose_mat4[13] = mat4[7];
	transpose_mat4[14] = mat4[11];
	transpose_mat4[15] = mat4[15];
}

void MAT4_Transpose(double * mat4) {
	swap(mat4[1], mat4[5]);
	swap(mat4[2], mat4[8]);
	swap(mat4[3], mat4[12]);
	swap(mat4[6], mat4[9]);
	swap(mat4[7], mat4[13]);
	swap(mat4[11], mat4[14]);
}

void MAT4_SetEqual(double * mat4, double * output_mat4) {
	for( int i=0; i<16; i++ ) {
		output_mat4[i] = mat4[i];
	}
}

void MAT4_InvertRotation(double * rotation_mat4, double * invert_mat4) {
	MAT4_TransposeUpper3(rotation_mat4, invert_mat4);

	double t_vec3[3];
	MAT4_ExtractCol3(rotation_mat4, 3, t_vec3);
	VEC3_Negate(t_vec3);

	double x_vec3[3];
	MAT4_ExtractCol3(rotation_mat4, 0, x_vec3);
	double y_vec3[3];
	MAT4_ExtractCol3(rotation_mat4, 1, y_vec3);
	double z_vec3[3];
	MAT4_ExtractCol3(rotation_mat4, 2, z_vec3);
	
	invert_mat4[12] = VEC3_Dot(t_vec3, x_vec3);
	invert_mat4[13] = VEC3_Dot(t_vec3, y_vec3);
	invert_mat4[14] = VEC3_Dot(t_vec3, z_vec3);
}


double MTH_DegToRad(double deg) {
    return deg*3.141592654/180.0;    
}

double MTH_RadToDeg(double rad) {
    return rad*180.0/3.141592654;
}

double MTH_Sign(double inValue) {
	return (inValue >= 0) ? 1.0 : -1.0;
}

TransformMatrix MTH_MakeArbitraryRotation(WorldPt3 & axis_vec3, double angle_rad) {
	double c = cos(angle_rad);
	double s = sin(angle_rad);
	double t = 1-c;
	double & x = axis_vec3.x;
	double & y = axis_vec3.y;
	double & z = axis_vec3.z;

	TransformMatrix result;
	result.v1.a00 = t*x*x + c;
	result.v1.a10 = t*x*y - s*z;
	result.v1.a20 = t*x*y + s*y;
	result.v1.a01 = t*x*y + s*z;
	result.v1.a11 = t*y*y + c;
	result.v1.a21 = t*y*z - s*x;
	result.v1.a02 = t*x*z - s*y;
	result.v1.a12 = t*y*z + s*x;
	result.v1.a22 = t*z*z + c;
	result.v1.xOff = 0;
	result.v1.yOff = 0;
	result.v1.zOff = 0;
	return result;
}

void MAT4_MakeZRotation(double deg, double * zRotation_mat4) {
	double rad = MTH_DegToRad(deg);

	double c = cos(rad);
	double s = sin(rad);
	zRotation_mat4[0] = c;
	zRotation_mat4[1] = s;
	zRotation_mat4[2] = 0;
	zRotation_mat4[3] = 0;
	zRotation_mat4[4] = -s;
	zRotation_mat4[5] = c;
	zRotation_mat4[6] = 0;
	zRotation_mat4[7] = 0;
	zRotation_mat4[8] = 0;
	zRotation_mat4[9] = 0;
	zRotation_mat4[10] = 1;
	zRotation_mat4[11] = 0;
	zRotation_mat4[12] = 0;
	zRotation_mat4[13] = 0;
	zRotation_mat4[14] = 0;
	zRotation_mat4[15] = 1;

}

void MAT4_MakeXRotation(double deg, double * xRotation_mat4) {
	double rad = MTH_DegToRad(deg);

	double c = cos(rad);
	double s = sin(rad);
	xRotation_mat4[0] = 1;
	xRotation_mat4[1] = 0;
	xRotation_mat4[2] = 0;
	xRotation_mat4[3] = 0;
	xRotation_mat4[4] = 0;
	xRotation_mat4[5] = c;
	xRotation_mat4[6] = s;
	xRotation_mat4[7] = 0;
	xRotation_mat4[8] = 0;
	xRotation_mat4[9] = -s;
	xRotation_mat4[10] = c;
	xRotation_mat4[11] = 0;
	xRotation_mat4[12] = 0;
	xRotation_mat4[13] = 0;
	xRotation_mat4[14] = 0;
	xRotation_mat4[15] = 1;

}

void MAT4_MakeArbitraryRotation(double * axis_vec3, double angle_rad, double * rotation_mat4) {
	double c = cos(angle_rad);
	double s = sin(angle_rad);
	double t = 1-c;
	double & x = axis_vec3[0];
	double & y = axis_vec3[1];
	double & z = axis_vec3[2];

	rotation_mat4[0] = t*x*x + c;
	rotation_mat4[1] = t*x*y - s*z;
	rotation_mat4[2] = t*x*y + s*y;
	rotation_mat4[3] = 0;
	rotation_mat4[4] = t*x*y + s*z;
	rotation_mat4[5] = t*y*y + c;
	rotation_mat4[6] = t*y*z - s*x;
	rotation_mat4[7] = 0;
	rotation_mat4[8] = t*x*z - s*y;
	rotation_mat4[9] = t*y*z + s*x;
	rotation_mat4[10] = t*z*z + c;
	rotation_mat4[11] = 0;
	rotation_mat4[12] = 0;
	rotation_mat4[13] = 0;
	rotation_mat4[14] = 0;
	rotation_mat4[15] = 1;
}


void MAT4_MakeArbitraryRotation(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * rotation_mat4) {
	float c = cosf(inAngleRad);
	float s = sinf(inAngleRad);
	float t = 1-c;
	float & x = inAxisX;
	float & y = inAxisY;
	float & z = inAxisZ;

	rotation_mat4[0] = t*x*x + c;
	rotation_mat4[1] = t*x*y - s*z;
	rotation_mat4[2] = t*x*y + s*y;
	rotation_mat4[3] = 0;
	rotation_mat4[4] = t*x*y + s*z;
	rotation_mat4[5] = t*y*y + c;
	rotation_mat4[6] = t*y*z - s*x;
	rotation_mat4[7] = 0;
	rotation_mat4[8] = t*x*z - s*y;
	rotation_mat4[9] = t*y*z + s*x;
	rotation_mat4[10] = t*z*z + c;
	rotation_mat4[11] = 0;
	rotation_mat4[12] = 0;
	rotation_mat4[13] = 0;
	rotation_mat4[14] = 0;
	rotation_mat4[15] = 1;
}

void MAT4_PostRotate(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * inOutMat) {
    float mat[16];
    MAT4_MakeArbitraryRotation(inAxisX, inAxisY, inAxisZ, inAngleRad, mat);
    
    MAT4_Multiply(inOutMat, mat, inOutMat);
    
}

void MAT4_PreRotate(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * inOutMat) {
    float mat[16];
    MAT4_MakeArbitraryRotation(inAxisX, inAxisY, inAxisZ, inAngleRad, mat);
    
    MAT4_Multiply(mat, inOutMat, inOutMat);
    
}

void MAT4_MakeTranslation(float inX, float inY, float inZ, float * outMat) {
    outMat[0] = 1;
    outMat[1] = 0;
    outMat[2] = 0;
    outMat[3] = 0;
    outMat[4] = 0;
    outMat[5] = 1;
    outMat[6] = 0;
    outMat[7] = 0;
    outMat[8] = 0;
    outMat[9] = 0;
    outMat[10] = 1;
    outMat[11] = 0;
    outMat[12] = inX;
    outMat[13] = inY;
    outMat[14] = inZ;
    outMat[15] = 1;
}

void MAT4_PostTranslate(float inX, float inY, float inZ, float * inOutMat) 
// inOutMat = inOutMat * translateMat
{
    float mat[16];
    MAT4_MakeTranslation(inX, inY, inZ, mat);
    MAT4_Multiply(inOutMat, mat, inOutMat);
}

void MAT4_PreTranslate(float inX, float inY, float inZ, float * inOutMat) 
// inOutMat = inOutMat * translateMat
{
    float mat[16];
    MAT4_MakeTranslation(inX, inY, inZ, mat);
    MAT4_Multiply(mat, inOutMat, inOutMat);
}

void MAT4_MakeScale(float inScaleX, float inScaleY, float inScaleZ, float * outMat) {
    outMat[0] = inScaleX;
    outMat[1] = 0;
    outMat[2] = 0;
    outMat[3] = 0;
    outMat[4] = 0;
    outMat[5] = inScaleY;
    outMat[6] = 0;
    outMat[7] = 0;
    outMat[8] = 0;
    outMat[9] = 0;
    outMat[10] = inScaleZ;
    outMat[11] = 0;
    outMat[12] = 0;
    outMat[13] = 0;
    outMat[14] = 0;
    outMat[15] = 1;
}

void MAT4_PreScale(float inScaleX, float inScaleY, float inScaleZ, float * inOutMat) {
    float mat[16];
    MAT4_MakeScale(inScaleX, inScaleY, inScaleZ, mat);
    MAT4_Multiply(mat, inOutMat, inOutMat);
}

void MAT4_PostScale(float inScaleX, float inScaleY, float inScaleZ, float * inOutMat) {
    float mat[16];
    MAT4_MakeScale(inScaleX, inScaleY, inScaleZ, mat);
    MAT4_Multiply(inOutMat, mat, inOutMat);
}

void MAT4_MakeInversePerspective(double fovy_rad, double aspect, double zNear, double zFar, double * mat4) {
	double f = 1.0 / tan(fovy_rad / 2);
	double dp = zNear - zFar;
	mat4[0] = aspect/f;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 1/f;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = 0;
	mat4[11] = dp/(2*zFar*zNear);
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = -1;
	mat4[15] =(zFar+zNear)/(2*zNear*zFar);
}

void MAT4_MakePerspective(double fovy_rad, double aspect, double zNear, double zFar, double * mat4) {
	double f = 1.0 / tan(fovy_rad / 2);
	mat4[0] = f/aspect;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = f;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = (zFar+zNear)/(zNear-zFar);
	mat4[11] = -1;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = (2*zFar*zNear)/(zNear-zFar);
	mat4[15] = 0;
}

void MAT4_MakeOrtho(double left, double right, double bottom, double top, double zNear, double zFar, double * mat4) {
	double tx = -(right+left) / (right-left);
	double ty = -(top+bottom) / (top-bottom);
	double tz = -(zFar+zNear) / (zFar-zNear);

	mat4[0] = 2.0 / (right-left);
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 2.0 / (top-bottom);
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = -2.0 / (zFar-zNear);
	mat4[11] = 0;
	mat4[12] = tx;
	mat4[13] = ty;
	mat4[14] = tz;
	mat4[15] = 1;
}

void MAT4_MakeFrustum(double left, double right, double bottom, double top, double zNear, double zFar, double * mat4) {
	double A = (right+left) / (right-left);
	double B = (top+bottom) / (top-bottom);
	double C = -(zFar+zNear) / (zFar-zNear);
	double D = -(2*zFar*zNear) / (zFar-zNear);

	mat4[0] = (2*zNear) / (right-left);
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = (2*zNear) / (top-bottom);
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = A;
	mat4[9] = B;
	mat4[10] = C;
	mat4[11] = -1;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = D;
	mat4[15] = 0;
}

void MAT4_LookAt(float eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ, double * viewMatrix) {
	double * m = viewMatrix;

	double eye[] = { eyeX, eyeY, eyeZ };
	double center[] = { centerX, centerY, centerZ };
	double up[] = { upX, upY, upZ };
	VEC3_Normalize(up);

	double f[3];
	VEC3_Subtract(center, eye, f);
	VEC3_Normalize(f);

	double s[3];
	VEC3_Cross(f, up, s);
	VEC3_Normalize(s);

	double u[3];
	VEC3_Cross(s, f, u);

	double mat[16] = {
		s[0], u[0], -f[0], 0,
		s[1], u[1], -f[1], 0,
		s[2], u[2], -f[2], 0,
		0, 0, 0, 1 };

	double t[16];
	MAT4_LoadIdentity(t);
	MAT4_Translate(-eyeX, -eyeY, -eyeZ, t);

	MAT4_Multiply(mat, t, m);	
}

void MAT4_Translate(double inX, double inY, double inZ, double * inOutMat4) {
	inOutMat4[12] += inX;
	inOutMat4[13] += inY;
	inOutMat4[14] += inZ;
}

void FPS_CAM_MoveForward(double amount, bool fixedZ, double * camera_mat4) {
	double forward_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 2, forward_vec3);
	VEC3_Negate(forward_vec3);

	double offset_vec3[3];
	if( fixedZ ) {
		//VEC3_Equal(
		//VEC3_Dot(forward_vec3, 
	} else {
		VEC3_Multiply(forward_vec3, amount, offset_vec3);

	}
	double pos_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 3, pos_vec3);

	VEC3_Add(pos_vec3, offset_vec3, pos_vec3);

	MAT4_SetCol3(3, pos_vec3, camera_mat4);
}

void FPS_CAM_StrafeRight(double amount, double * camera_mat4) {
	double right_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 0, right_vec3);

	double offset_vec3[3];
	VEC3_Multiply(right_vec3, amount, offset_vec3);

	double pos_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 3, pos_vec3);

	VEC3_Add(pos_vec3, offset_vec3, pos_vec3);

	MAT4_SetCol3(3, pos_vec3, camera_mat4);
}

void FPS_CAM_LookRight(double amount_deg, double * camera_mat4) {
    double zRotation_mat4[16];
	MAT4_MakeZRotation(-amount_deg, zRotation_mat4);

	double pos_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 3, pos_vec3);

	double mat4[16];
	MAT4_SetEqual(camera_mat4, mat4);
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = 0;

	MAT4_Multiply(zRotation_mat4, mat4, camera_mat4);
	camera_mat4[12] = pos_vec3[0];
	camera_mat4[13] = pos_vec3[1];
	camera_mat4[14] = pos_vec3[2];
}

void FPS_CAM_LookUp(double amount_deg, double * camera_mat4) {
	
	if( amount_deg > 0 ) {
		double up_vec3[3] = { 0, 0, 1 };

		double forward_vec3[3];
		MAT4_ExtractCol3(camera_mat4, 2, forward_vec3);
		VEC3_Negate(forward_vec3);

		double dot = VEC3_Dot(up_vec3, forward_vec3);

		double angleBetween_rad = acos(dot);
		double angleBetween_deg = MTH_RadToDeg(angleBetween_rad);

		if( angleBetween_deg < 10 ) {
			amount_deg = angleBetween_deg-10;
		}
	} else if( amount_deg < 0 ) {
		double down_vec3[3] = { 0, 0, -1 };

		double forward_vec3[3];
		MAT4_ExtractCol3(camera_mat4, 2, forward_vec3);
		VEC3_Negate(forward_vec3);

		double dot = VEC3_Dot(down_vec3, forward_vec3);

		double angleBetween_rad = acos(dot);
		double angleBetween_deg = MTH_RadToDeg(angleBetween_rad);

		if( angleBetween_deg < 10 ) {
			amount_deg = 10-angleBetween_deg;
		}

	}

    double xRotation_mat4[16];
	MAT4_MakeXRotation(amount_deg, xRotation_mat4);

	double pos_vec3[3];
	MAT4_ExtractCol3(camera_mat4, 3, pos_vec3);

	double mat4[16];
	MAT4_SetEqual(camera_mat4, mat4);
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = 0;

	MAT4_Multiply(mat4, xRotation_mat4, camera_mat4);
	camera_mat4[12] = pos_vec3[0];
	camera_mat4[13] = pos_vec3[1];
	camera_mat4[14] = pos_vec3[2];
}



void RND_GenerateNormalizedVector(double * vec3) {
	vec3[0] = RND_DoubleInRange(0.0, 1.0);
	vec3[1] = RND_DoubleInRange(0.0, 1-vec3[0]);
	vec3[2] = sqrt(1 - vec3[0]*vec3[0] - vec3[1]*vec3[1]);
}

WorldPt3 RND_GenerateNormalizedVector() {
	WorldPt3 vec;
	RND_GenerateNormalizedVector((double *)&vec);
	return vec;
}

bool RaySphereIntersection(WorldPt3 & sphereCenter, double radius, WorldPt3 & rayOrigin, WorldPt3 & rayDir, /*out*/ double * t, /*out*/ WorldPt3 * intersection) {
	WorldPt3 s = sphereCenter;
	WorldPt3 o = rayOrigin;
	WorldPt3 d = rayDir.Normal();
	double r2 = radius*radius;



	WorldPt3 dst = o - s;
	double B = dst.Dot(d);
	double C = dst.Dot(dst) - r2;
	double D = B*B - C;
	if( D <= 0 ) return false;

	double sqrtD = sqrt(D);
	double tt = -B - sqrtD;
	double tt2 = -B + sqrtD;
	if( tt > 0 && tt2 > 0) {
		tt = min(tt, tt2);
	} else if( tt2 > 0 ) {
		tt = tt2;
	}
	
	if( t != 0 ) {
		*t = tt;
	}
	if( intersection != 0 ) {
		*intersection = rayOrigin + d * tt;
	}

	return true;
}



bool RayCylinderIntersection(double radius, double * minZ, double * maxZ, WorldPt3 & rayOrigin, WorldPt3 & rayDir, double * t, WorldPt3 * intersection, bool findNegativeIntersections) { // cylinder base on xy plane extruded height in +z

	double a = rayDir.x*rayDir.x + rayDir.y*rayDir.y;
	double b = 2*(rayOrigin.x*rayDir.x + rayOrigin.y*rayDir.y);
	double c = rayOrigin.x*rayOrigin.x + rayOrigin.y*rayOrigin.y - radius*radius;

	vector<double> solution = MTH_QuadraticEquation(a, b, c);

	bool intersect = false;
	double bestT = 0;

	double tt;

	for( size_t i=0; i<solution.size(); i++ ) {
		tt = solution[i];
		WorldPt3 p = rayOrigin + tt*rayDir;

		if( tt>0 || findNegativeIntersections) {
			if( minZ && p.z < *minZ ) {
				continue;
			}
			if( maxZ && p.z > *maxZ ) {
				continue;
			}
			if( !intersect || tt<bestT ) {
				intersect = true;
				bestT = tt;
			}
		}
	}


	WorldPt3 intersectionPt;
	
	if( minZ ) {
		if( RayPlaneIntersection(WorldPt3(0, 0, *minZ), WorldPt3(0, 0, 1), rayOrigin, rayDir, &tt, &intersectionPt) ) {
			if( intersectionPt.MagnitudeSquared() < radius*radius ) {
				if( !intersect || tt<bestT ) {
					intersect = true;
					bestT = tt;
				}
			}
		}
	}

	if( maxZ ) {
		if( RayPlaneIntersection(WorldPt3(0, 0, *maxZ), WorldPt3(0, 0, 1), rayOrigin, rayDir, &tt, &intersectionPt) ) {
			if( WorldPt3(intersectionPt.x, intersectionPt.y, 0).MagnitudeSquared() < radius*radius ) {
				if( !intersect || tt<bestT ) {
					intersect = true;
					bestT = tt;
				}
			}
		}
	}

	if( intersect ) { // if ray intersects the cylinder
		if( t ) { // if client wants to know t value at intersection
			*t = bestT;
		}

		if( intersection ) { // if client wants to know intersection point
			*intersection = rayOrigin + bestT*rayDir;
		}

		return true;
	}

	return false;

}

void MAT4_VWViewMatToGLViewMat(TransformMatrix & inMat, double inCrsPersDistance, float * outMat) {
	float * m = outMat;
	m[ 0] = (float)inMat.v1.a00;
	m[ 1] = (float)inMat.v1.a01;
	m[ 2] = (float)inMat.v1.a02;
	m[ 3] = 0;
	m[ 4] = (float)inMat.v1.a10;
	m[ 5] = (float)inMat.v1.a11;
	m[ 6] = (float)inMat.v1.a12;
	m[ 7] = 0;
	m[ 8] = (float)inMat.v1.a20;
	m[ 9] = (float)inMat.v1.a21;
	m[10] = (float)inMat.v1.a22;
	m[11] = 0;
	m[12] = (float)inMat.v1.xOff;
	m[13] = (float)inMat.v1.yOff;
	m[14] = (float)inMat.v1.zOff;
	m[15] = 1;

	// in VW the world eye point is located 'inCrsPersDistance' back from what the matrix says
	m[14] -= (float)inCrsPersDistance;
}
void MAT4_VWToGLMatrix(TransformMatrix & inMat, float * outMat) {
	float * m = outMat;
	m[ 0] = (float)inMat.v1.a00;
	m[ 1] = (float)inMat.v1.a01;
	m[ 2] = (float)inMat.v1.a02;
	m[ 3] = 0;
	m[ 4] = (float)inMat.v1.a10;
	m[ 5] = (float)inMat.v1.a11;
	m[ 6] = (float)inMat.v1.a12;
	m[ 7] = 0;
	m[ 8] = (float)inMat.v1.a20;
	m[ 9] = (float)inMat.v1.a21;
	m[10] = (float)inMat.v1.a22;
	m[11] = 0;
	m[12] = (float)inMat.v1.xOff;
	m[13] = (float)inMat.v1.yOff;
	m[14] = (float)inMat.v1.zOff;
	m[15] = 1;
}

bool GEO_LineLineIntersection(WorldPt & inP1, WorldPt & inV1, WorldPt & inP2, WorldPt & inV2, bool inSegmentIntersection, bool & outParallel, double & outT1, double & outT2, WorldPt & outIntersection)
//   Determine the intersection point of two line segments
//   Return FALSE if the lines don't intersect
// original code: http://paulbourke.net/geometry/lineline2d/
{
	outParallel = false;
	outT1 = 0;
	outT2 = 0;

	WorldPt & p1 = inP1;
	WorldPt & v1 = inV1;

	WorldPt & p2 = inP2;
	WorldPt & v2 = inV2;

	double dx = p1.x-p2.x;
	double dy = p1.y-p2.y;

	double denom  = v2.y * v1.x - v2.x * v1.y;
	double numera = v2.x * dy - v2.y * dx;
	double numerb = v1.x * dy - v1.y * dx;

	double EPS = 1e-4;

	// Are the line coincident?
	if (abs(numera) < EPS && abs(numerb) < EPS && abs(denom) < EPS) {
		outIntersection = (p1 + p2) * 0.5;
		outParallel = true;
		return true;
	}

	// Are the line parallel
	if (abs(denom) < EPS) {
		outIntersection = WorldPt(0, 0);
		outParallel = true;
		return true;
	}

	// Is the intersection along the the segments
	outT1 = numera / denom;
	outT2 = numerb / denom;
	if( inSegmentIntersection && (outT1 < 0 || outT1 > 1 || outT2 < 0 || outT2 > 1) ) {
		outIntersection = WorldPt(0, 0);
		return false;
	}

	outIntersection = p1 + inV1 * outT1;

	return true;
}
void GEO_GenerateSphere(int inPrimitiveType, float * inOrigin, float inRadius, int inSamples, vector<float> & outVertices, vector<float> & outNormals, vector<float> & outTexCoords, vector<int> & outIndices, int * outNumTriangles) {
	
	if( inPrimitiveType == GL_TRIANGLES ) { // with indices
		int kSamples = inSamples;
		int kNumUniquePoints = (kSamples-2)*(kSamples) + 2; // 2 poles
		int numTriangles = (kSamples-3)*(2*kSamples) + 2*kSamples;
		outVertices.resize(kNumUniquePoints*3); // xyz
		outNormals.resize(kNumUniquePoints*3); // xyz
		outTexCoords.resize(kNumUniquePoints*2); // uv
		outIndices.resize(numTriangles*3); // indices 3 per triangle

		float radius = inRadius;
		float * v = &outVertices[0];
		float * n = &outNormals[0];
		float * t = &outTexCoords[0];
		int * i = &outIndices[0];
	
		if( outNumTriangles ) {
			*outNumTriangles = numTriangles;
		}

		for( int phiI=0; phiI<kSamples; phiI++ ) {
			float phi = 3.141592654f * phiI / (float)(kSamples-1);
			float z = inOrigin[2] + radius * cosf(phi);

			if( phiI == 0 ) {
				*v++ = inOrigin[0];
				*v++ = inOrigin[1];
				*v++ = inOrigin[2] + radius;
				*n++ = 0;
				*n++ = 0;
				*n++ = 1;
				*t++ = 0.5f;
				*t++ = 0;
			} else if( phiI == kSamples-1 ) {
				*v++ = inOrigin[0];
				*v++ = inOrigin[1];
				*v++ = inOrigin[2]-radius;
				*n++ = 0;
				*n++ = 0;
				*n++ = -1;
				*t++ = 0.5f;
				*t++ = 1;

				for( int thetaI=0; thetaI<kSamples; thetaI++ ) {
					*i++ = kNumUniquePoints-1;
					*i++ = 1+(phiI-2)*kSamples + (thetaI+1)%kSamples;
					*i++ = 1+(phiI-2)*kSamples + thetaI;
				}
			} else {

				for( int thetaI=0; thetaI<kSamples; thetaI++ ) {
					float theta = 2.0f*3.141592654f * thetaI / (float)kSamples;

					float x = radius * cosf(theta) * sinf(phi);
					float y = radius * sinf(theta) * sinf(phi);
					float magnitude = sqrtf(x*x + y*y + z*z);
			
					*v++ = inOrigin[0] + x;
					*v++ = inOrigin[1] + y;
					*v++ = inOrigin[2] + z;

					*n++ = x / magnitude;
					*n++ = y / magnitude;
					*n++ = z / magnitude;

					*t++ = thetaI / (float)kSamples;
					*t++ = phiI / (float)kSamples;

					if( phiI==1 ) {
						*i++ = 1+thetaI;
						*i++ = 1+(thetaI+1)%kSamples;
						*i++ = 0;
					} else {
						int v0 = 1+(phiI-1)*kSamples + thetaI;
						int v1 = 1+(phiI-1)*kSamples + (thetaI+1)%kSamples;
						int v2 = 1+(phiI-2)*kSamples + thetaI;
						int v3 = 1+(phiI-2)*kSamples + (thetaI+1)%kSamples;

						*i++ = v0;
						*i++ = v1;
						*i++ = v2;

						*i++ = v2;
						*i++ = v1;
						*i++ = v3;
					}
				}
			}
		}


	} else if( inPrimitiveType == GL_TRIANGLE_STRIP ) {
		const int Size = (int)((3.141592654 / .2 + 1)*(2*3.141592654 / .2 + 1)) * 6;
		outVertices.resize(Size);
		outNormals.resize(Size);

		if( outNumTriangles ) {
			*outNumTriangles = Size-2;
		}

		int i = 0;
		float radius = inRadius;

		for( float phi = 0.2f; phi <= 3.141592654f; phi += 0.2f ) {
			for( float theta = 0.0f; theta <= 2.0f*3.141592654f; theta += 0.2f ) {
				float px = radius * cosf(theta) * sinf(phi-0.2f);
				float py = radius * sinf(theta) * sinf(phi-0.2f);
				float pz = radius * cosf(phi-0.2f);
				float pmagnitude = sqrtf(px*px + py*py + pz*pz);
			
				outNormals[i] = px/pmagnitude;
				outVertices[i++] = inOrigin[0]+px;
				outNormals[i] = py/pmagnitude;
				outVertices[i++] = inOrigin[1]+py;
				outNormals[i] = pz/pmagnitude;
				outVertices[i++] = inOrigin[2]+pz;

				float x = radius * cosf(theta) * sinf(phi);
				float y = radius * sinf(theta) * sinf(phi);
				float z = radius * cosf(phi);
				float magnitude = sqrtf(x*x + y*y + z*z);

				outNormals[i] = x/magnitude;
				outVertices[i++] = inOrigin[0]+x;
				outNormals[i] = y/magnitude;
				outVertices[i++] = inOrigin[1]+y;
				outNormals[i] = z/magnitude;
				outVertices[i++] = inOrigin[2]+z;
			}
		}
	}
}


void MAT4_VWToGLMatrix(TransformMatrix & mat, double * glMatrix) {
	double * m = glMatrix;
	m[ 0] = mat.v1.a00;
	m[ 1] = mat.v1.a01;
	m[ 2] = mat.v1.a02;
	m[ 3] = 0;
	m[ 4] = mat.v1.a10;
	m[ 5] = mat.v1.a11;
	m[ 6] = mat.v1.a12;
	m[ 7] = 0;
	m[ 8] = mat.v1.a20;
	m[ 9] = mat.v1.a21;
	m[10] = mat.v1.a22;
	m[11] = 0;
	m[12] = mat.v1.xOff;
	m[13] = mat.v1.yOff;
	m[14] = mat.v1.zOff;
	m[15] = 1;
}
void MAT4_GLToVWMatrix(double * glMatrix, TransformMatrix & mat) {
	double * m = glMatrix;

	mat.v1.a00 = m[0];
	mat.v1.a01 = m[1];
	mat.v1.a02 = m[2];
	mat.v1.a10 = m[4];
	mat.v1.a11 = m[5];
	mat.v1.a12 = m[6];
	mat.v1.a20 = m[8];
	mat.v1.a21 = m[9];
	mat.v1.a22 = m[10];
	mat.v1.xOff = m[12];
	mat.v1.yOff = m[13];
	mat.v1.zOff = m[14];
}


void GL_Frustum(double left, double right, double bottom, double top, double zNear, double zFar, double * projectionMatrix) {
	double * m = projectionMatrix;
	
	double a = (right+left) / (right-left);
	double b = (top+bottom) / (top-bottom);
	double c = (zFar+zNear) / (zFar-zNear);
	double d = (2*zFar*zNear) / (zFar-zNear);

	m[ 0] = 2*zNear / (right-left);
	m[ 1] = 0;
	m[ 2] = 0;
	m[ 3] = 0;
	m[ 4] = 0;
	m[ 5] = 2*zNear / (top-bottom);
	m[ 6] = 0;
	m[ 7] = 0;
	m[ 8] = a;
	m[ 9] = b;
	m[10] = c;
	m[11] = -1;
	m[12] = 0;
	m[13] = 0;
	m[14] = d;
	m[15] = 0;
}

void GLU_LookAt(double eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ, double * viewMatrix) {
	double * m = viewMatrix;

	WorldPt3 eye(eyeX, eyeY, eyeZ);
	WorldPt3 center(centerX, centerY, centerZ);
	WorldPt3 up = WorldPt3(upX, upY, upZ).Normal();

	WorldPt3 f = (center - eye).Normal();
	WorldPt3 s = CrossProduct(f, up);
	WorldPt3 u = CrossProduct(s, f);

	TransformMatrix mat;
	mat.v1.a00 = s.x;
	mat.v1.a01 = s.y;
	mat.v1.a02 = s.z;
	mat.v1.a10 = u.x;
	mat.v1.a11 = u.y;
	mat.v1.a12 = u.z;
	mat.v1.a20 = -f.x;
	mat.v1.a21 = -f.y;
	mat.v1.a22 = -f.z;
	mat.v1.xOff = 0;
	mat.v1.yOff = 0;
	mat.v1.zOff = 0;

	TransformMatrix t;
	TranslateMatrix(t, -eyeX, -eyeY, -eyeZ);

	mat = MatrixMultiply(mat, t);

	m[ 0] = mat.v1.a00;
	m[ 1] = mat.v1.a10;
	m[ 2] = mat.v1.a20;
	m[ 3] = 0;
	m[ 4] = mat.v1.a01;
	m[ 5] = mat.v1.a11;
	m[ 6] = mat.v1.a21;
	m[ 7] = 0;
	m[ 8] = mat.v1.a02;
	m[ 9] = mat.v1.a12;
	m[10] = mat.v1.a22;
	m[11] = 0;
	m[12] = mat.v1.xOff;
	m[13] = mat.v1.yOff;
	m[14] = mat.v1.zOff;
	m[15] = 1;
}



void MAT4_VEC4_Multiply(double * mat4, double * column_vec4, double * row_vec4) {
	double result[4];

	result[0] = mat4[ 0]*column_vec4[ 0] + mat4[ 4]*column_vec4[ 1] + mat4[ 8]*column_vec4[ 2] + mat4[12]*column_vec4[ 3];
	result[1] = mat4[ 1]*column_vec4[ 0] + mat4[ 5]*column_vec4[ 1] + mat4[ 9]*column_vec4[ 2] + mat4[13]*column_vec4[ 3];
	result[2] = mat4[ 2]*column_vec4[ 0] + mat4[ 6]*column_vec4[ 1] + mat4[10]*column_vec4[ 2] + mat4[14]*column_vec4[ 3];
	result[3] = mat4[ 3]*column_vec4[ 0] + mat4[ 7]*column_vec4[ 1] + mat4[11]*column_vec4[ 2] + mat4[15]*column_vec4[ 3];

	memcpy(row_vec4, result, 4*sizeof(double));
}



void MAT4_RotationFromSphericalCoords(double theta, double phi, TransformMatrix & cameraMat) {

	double ct = cos(theta);
	double st = sin(theta);
	double cp = cos(phi);
	double sp = sin(phi);

	double xx = ct*sp;
	double yy = st*sp;
	double zz = cp;

	WorldPt3 back(xx, yy, zz);
	WorldPt3 up(0, 0, 1);
	WorldPt3 right = CrossProduct(up, back).Normal();
	up = CrossProduct(back, right).Normal();

	
	cameraMat.v1.a00 = right.x;
	cameraMat.v1.a01 = right.y;
	cameraMat.v1.a02 = right.z;
	cameraMat.v1.a10 = up.x;
	cameraMat.v1.a11 = up.y;
	cameraMat.v1.a12 = up.z;
	cameraMat.v1.a20 = back.x;
	cameraMat.v1.a21 = back.y;
	cameraMat.v1.a22 = back.z;
}

vector<double> MTH_QuadraticEquation(double a, double b, double c) { // ax^2 + bx + c = 0; return contains real solutions
	vector<double> solutions;

	double two_a = 2*a;
	if( abs(two_a) < 1e-7 ) { // if no solutions
		return solutions;
	}

	double root_squared = b*b - 4*a*c;
	if( root_squared < 0 ) { // if no real solutions
		return solutions;
	}

	double root = sqrt(root_squared);

	if( sqrt(root) < 1e-7 ) { // if single solution
		solutions.push_back(-b / two_a);

	} else { // if two solutions
		solutions.push_back((-b + root) / two_a);
		solutions.push_back((-b - root) / two_a);
	}

	return solutions;
}


bool MAT4_Invert(const double * mat4, double * out4) {
	/* Compute inverse of 4x4 transformation matrix.
	Code contributed by Jacques Leroy jle@star.be
	return - true for success, false for failure (singular matrix)
	
	modified from http://webcvs.freedesktop.org/mesa/Mesa/src/glu/mesa/project.c?revision=1.4&view=markup
	*/

/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { double *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(mat4,r,c) (mat4)[(c)*4+(r)]

	if( mat4 == out4 ) { // if input and output point to same memory
		return false;
	}

   double wtmp[4][8];
   double m0, m1, m2, m3, s;
   double *r0, *r1, *r2, *r3;

   r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

   r0[0] = MAT(mat4, 0, 0), r0[1] = MAT(mat4, 0, 1),
      r0[2] = MAT(mat4, 0, 2), r0[3] = MAT(mat4, 0, 3),
      r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,
      r1[0] = MAT(mat4, 1, 0), r1[1] = MAT(mat4, 1, 1),
      r1[2] = MAT(mat4, 1, 2), r1[3] = MAT(mat4, 1, 3),
      r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,
      r2[0] = MAT(mat4, 2, 0), r2[1] = MAT(mat4, 2, 1),
      r2[2] = MAT(mat4, 2, 2), r2[3] = MAT(mat4, 2, 3),
      r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,
      r3[0] = MAT(mat4, 3, 0), r3[1] = MAT(mat4, 3, 1),
      r3[2] = MAT(mat4, 3, 2), r3[3] = MAT(mat4, 3, 3),
      r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

   /* choose pivot - or die */
   if (fabs(r3[0]) > fabs(r2[0]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[0]) > fabs(r1[0]))
      SWAP_ROWS(r2, r1);
   if (fabs(r1[0]) > fabs(r0[0]))
      SWAP_ROWS(r1, r0);
   if (0.0 == r0[0])
      return false;

   /* eliminate first variable     */
   m1 = r1[0] / r0[0];
   m2 = r2[0] / r0[0];
   m3 = r3[0] / r0[0];
   s = r0[1];
   r1[1] -= m1 * s;
   r2[1] -= m2 * s;
   r3[1] -= m3 * s;
   s = r0[2];
   r1[2] -= m1 * s;
   r2[2] -= m2 * s;
   r3[2] -= m3 * s;
   s = r0[3];
   r1[3] -= m1 * s;
   r2[3] -= m2 * s;
   r3[3] -= m3 * s;
   s = r0[4];
   if (s != 0.0) {
      r1[4] -= m1 * s;
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r0[5];
   if (s != 0.0) {
      r1[5] -= m1 * s;
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r0[6];
   if (s != 0.0) {
      r1[6] -= m1 * s;
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r0[7];
   if (s != 0.0) {
      r1[7] -= m1 * s;
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[1]) > fabs(r2[1]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[1]) > fabs(r1[1]))
      SWAP_ROWS(r2, r1);
   if (0.0 == r1[1])
      return GL_FALSE;

   /* eliminate second variable */
   m2 = r2[1] / r1[1];
   m3 = r3[1] / r1[1];
   r2[2] -= m2 * r1[2];
   r3[2] -= m3 * r1[2];
   r2[3] -= m2 * r1[3];
   r3[3] -= m3 * r1[3];
   s = r1[4];
   if (0.0 != s) {
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r1[5];
   if (0.0 != s) {
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r1[6];
   if (0.0 != s) {
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r1[7];
   if (0.0 != s) {
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[2]) > fabs(r2[2]))
      SWAP_ROWS(r3, r2);
   if (0.0 == r2[2])
      return GL_FALSE;

   /* eliminate third variable */
   m3 = r3[2] / r2[2];
   r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
      r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

   /* last check */
   if (0.0 == r3[3])
      return GL_FALSE;

   s = 1.0 / r3[3];		/* now back substitute row 3 */
   r3[4] *= s;
   r3[5] *= s;
   r3[6] *= s;
   r3[7] *= s;

   m2 = r2[3];			/* now back substitute row 2 */
   s = 1.0 / r2[2];
   r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
      r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
   m1 = r1[3];
   r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
      r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
   m0 = r0[3];
   r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
      r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

   m1 = r1[2];			/* now back substitute row 1 */
   s = 1.0 / r1[1];
   r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
      r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
   m0 = r0[2];
   r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
      r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

   m0 = r0[1];			/* now back substitute row 0 */
   s = 1.0 / r0[0];
   r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
      r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

   MAT(out4, 0, 0) = r0[4];
   MAT(out4, 0, 1) = r0[5], MAT(out4, 0, 2) = r0[6];
   MAT(out4, 0, 3) = r0[7], MAT(out4, 1, 0) = r1[4];
   MAT(out4, 1, 1) = r1[5], MAT(out4, 1, 2) = r1[6];
   MAT(out4, 1, 3) = r1[7], MAT(out4, 2, 0) = r2[4];
   MAT(out4, 2, 1) = r2[5], MAT(out4, 2, 2) = r2[6];
   MAT(out4, 2, 3) = r2[7], MAT(out4, 3, 0) = r3[4];
   MAT(out4, 3, 1) = r3[5], MAT(out4, 3, 2) = r3[6];
   MAT(out4, 3, 3) = r3[7];

   return true;

#undef MAT
#undef SWAP_ROWS
}

bool GL_DrawString(const char * text, unsigned int textLen, int x, int y, int columnWidth, double * textWidth, double * textHeight) 
/* Draws left-justified text to the back buffer.  Uses the font of the current window's dc.  Assumes the 
		modelview matrix is active and that lighting is disabled.

		creates display lists 1000 to 1254 for fonts first time this function is called.
	text - the text to draw
	textLen - length of 'text'
	x - x coordinate of where to draw 'text'; 0 is the left side of the viewport
	y - y coordinate of where to darw 'text'; 0 is the top of the viewport
	columnWidth - the width from 'x' to the right side of an imaginary vertical line that is the word wrap break.
	return - true on succes; false on error
*/
{
	if( textWidth ) {
		*textWidth = 0;
	}

	if( textHeight ) {
		*textHeight = 0;
	}

#if GS_WIN && WIDGET_DEMO_APP

	glGetError(); // clear error flag

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	


	// assume gl_modelview
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, viewport[2]-1, 0, viewport[3]-1, -1, 1);
	

	GLuint displayListID = 1000;

	static TEXTMETRIC tm;
	static int widths[255];
	if( !glIsList(displayListID) ) {
		HDC dc = CreateDC(L"DISPLAY", 0, 0, 0);

		HFONT font = CreateFont(15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, L"Tahoma");
		SelectObject(dc, font);
	
		ABC abcs[256];
		GetCharABCWidths(dc, 0, 255, abcs);

		for( int i=0; i<255; i++ ) {
			widths[i] = abcs[i].abcA + abcs[i].abcB + abcs[i].abcC;
		}
	
		wglUseFontBitmaps(dc, 0, 255, displayListID);

		GetTextMetrics(dc, &tm);

		DeleteDC(dc);
	}
	glListBase(displayListID);

	int dy = tm.tmHeight;

	if( textHeight ) { // if client wants text height
		*textHeight = tm.tmHeight;
	}

	y += tm.tmAscent;
	
	int left = x;
	int right = x+columnWidth;


	/*
	glPushAttrib(GL_CURRENT_BIT);
	glBegin(GL_QUADS);
	glColor4d(1, 1, 1, .5);
	glVertex2i(x, y-tm.tmAscent);//- tm.tmAscent - 2*tm.tmDescent);
	glVertex2i(x, y+tm.tmDescent);//y+tm.tmh2*tm.tmDescent);
	glColor4d(1, 1, 1, .1);
	glVertex2i(x+300, y+tm.tmDescent);//y+2*tm.tmDescent );
	glVertex2i(x+300, y-tm.tmAscent);//y-tm.tmAscent - 2*tm.tmDescent);
	glEnd();
	glPopAttrib();
*/

	unsigned int currentI = 0;
	do {
		glRasterPos2d(0, 0);
		glBitmap(0, 0, 0, 0, (float)x, (float)(viewport[3]-y), 0);

		unsigned int chI = currentI;
		do {
			char ch = text[chI];
			if( ch == ' ' || ch=='\t' || ch=='\n' || ch=='\r' ) {
				break;
			}
			chI++;
		} while( chI < textLen );

		int width = 0;
		for( unsigned int i=currentI; i<chI; i++ ) {
			char ch = text[i];
			width += widths[ch];
		}

		bool newLine = false;

		if( x+width > right ) {
			newLine = true;
		} else {
			x += width;
		}

		glCallLists(chI - currentI, GL_UNSIGNED_BYTE, &text[currentI]);

		currentI = chI;

		do {
			char ch = text[chI];

			if( !(ch == ' ' || ch=='\t' || ch=='\n' || ch=='\r') ) {
				break;
			} else if( ch=='\n' || ch=='\r' ) { 
				newLine = true;
			}
			
			chI++;

		} while( chI < textLen );

		width = 0;
		if( !newLine ) {
			for( unsigned int i=currentI; i<chI; i++ ) {
				char ch = text[i];
				width += widths[ch];
			}

			if( x+width > right ) {
				newLine = true;
			}
		}

		if( newLine ) {
			if( textWidth ) {
				*textWidth = right-left;
			}
			y += dy;
			x = left;
		} else {
			x += width;
		}

		if( textWidth ) {
			*textWidth = max(*textWidth, (double)(x-left));
		}

		currentI = chI;
		

	} while( currentI < textLen );


	//glDeleteLists(displayListID, 255);

	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

#else
	TDrawPad * drawPad = gpRegenerationManager->GetOGLDrawPad();
	if( !drawPad ) {													// if error
		return false;
	}

	drawPad->TextMoveTo(x, y);
	drawPad->DrawTextString(text, textLen);
	if( textWidth ) {
		*textWidth = drawPad->GetTextWidth(text, textLen);
	}

	if( textHeight ) {
		DPPen pen;
		drawPad->GetPen(pen);
		*textHeight = pen.size;
	}

#endif
	return true;
}





void StringToWString(string & str, wstring & wstr) {
	wstr.resize(str.size());
	for( size_t i=0; i<str.size(); i++ ) {
		wstr[i] = str[i];
	}
}

void GL_DrawDiamond(const WorldPt & center, double width, double height) {
	glBegin(GL_LINE_LOOP);
	WorldPt p1 = center + WorldPt(width, 0);
	WorldPt p2 = center + WorldPt(0, height);
	WorldPt p3 = center + WorldPt(-width, 0);
	WorldPt p4 = center + WorldPt(0, -height);
	glVertex2dv((GLdouble *)&p1);
	glVertex2dv((GLdouble *)&p2);
	glVertex2dv((GLdouble *)&p3);
	glVertex2dv((GLdouble *)&p4);
	glEnd();
}
void GL_FillDiamond(const WorldPt & center, double width, double height) {
	glBegin(GL_POLYGON);
	WorldPt p1 = center + WorldPt(width, 0);
	WorldPt p2 = center + WorldPt(0, height);
	WorldPt p3 = center + WorldPt(-width, 0);
	WorldPt p4 = center + WorldPt(0, -height);
	glVertex2dv((GLdouble *)&p1);
	glVertex2dv((GLdouble *)&p2);
	glVertex2dv((GLdouble *)&p3);
	glVertex2dv((GLdouble *)&p4);
	glEnd();
}
void GL_DrawCircle(double radius, const WorldPt3 & normal, int segments) {

	TransformMatrix transform;
	transform.R2() = normal;
	WorldPt3 xAxis(1, 0, 0);

	if( ApproxEqualVectors(WorldPt3(1, 0, 0), normal) ) {
		xAxis = WorldPt3(0, 1, 0);
	}	
	transform.R1() = xAxis.Cross(normal).Normal();
	transform.R0() = transform.R1().Cross(normal).Normal();

	glBegin(GL_LINE_LOOP);
	for( int i=0; i<segments; i++ ) {
		double angle = i / (double)segments * 2*PI;
		WorldPt3 p(radius*cos(angle), radius*sin(angle), 0);
		p = VectorTransformN(p, transform);

		glVertex3dv((GLdouble *)&p);
	}
	glEnd();

}

void GL_DrawRect(const WorldPt3 & inOrigin, const WorldPt3 & inNormal, const WorldRect & inRect) {

	TransformMatrix transform;
	transform.R2() = inNormal;
	WorldPt3 xAxis(1, 0, 0);

	if( DoublesAreNearlyEqual(1, abs(xAxis.Dot(inNormal))) ) {
		xAxis = WorldPt3(0, 1, 0);
	}	
	transform.R1() = xAxis.Cross(inNormal).Normal();
	transform.R0() = transform.R1().Cross(inNormal).Normal();
	transform.P() = inOrigin;

	glBegin(GL_LINE_LOOP);

	WorldPt3 p0(inRect.TopLeft(), 0);
	p0 = PointTransformN(p0, transform);
	WorldPt3 p1(inRect.BotLeft(), 0);
	p1 = PointTransformN(p1, transform);
	WorldPt3 p2(inRect.BotRight(), 0);
	p2 = PointTransformN(p2, transform);
	WorldPt3 p3(inRect.TopRight(), 0);
	p3 = PointTransformN(p3, transform);

	glVertex3dv((GLdouble *)&p0);
	glVertex3dv((GLdouble *)&p1);
	glVertex3dv((GLdouble *)&p2);
	glVertex3dv((GLdouble *)&p3);

	glEnd();
}

void GL_Circle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments) {

	TransformMatrix transform;
	transform.R2() = normal;
	WorldPt3 xAxis(1, 0, 0);

	if( ApproxEqualVectors(WorldPt3(1, 0, 0), normal) ) {
		xAxis = WorldPt3(0, 1, 0);
	}	
	transform.R1() = xAxis.Cross(normal).Normal();
	transform.R0() = transform.R1().Cross(normal).Normal();

	for( int i=0; i<segments; i++ ) {
		double angle = i / (double)segments * 2*PI;
		WorldPt3 p(radius*cos(angle), radius*sin(angle), 0);
		p = center + VectorTransformN(p, transform);


		glVertex3dv((GLdouble *)&p);
	}

}

void GL_DrawCircle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments) {
	glBegin(GL_LINE_LOOP);
	GL_Circle(center, radius, normal, segments);
	glEnd();
}

void GL_FillCircle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments) {
	glBegin(GL_POLYGON);
	GL_Circle(center, radius, normal, segments);
	glEnd();
}

void GL_FillDisc(const WorldPt & center, double r1, double r2, int segments) {
	if( segments < 3 || r1 <= 0 || r2 <= 0 ) {							// erroneous conditions
		return;
	}

	glBegin(GL_TRIANGLE_STRIP);

	for( int i=0; i<=segments; i++ ) {									// for each segment
		double t = -(i/(double)segments) * 2*PI;
		WorldPt p1 = center + WorldPt(r1*cos(t), r1*sin(t));
		WorldPt p2 = center + WorldPt(r2*cos(t), r2*sin(t));
		
		glVertex2d(p1.x, p1.y);
		glVertex2d(p2.x, p2.y);
	}

	glEnd();
}

void GL_FillCircle(const WorldPt & center, double radius, int segments) {
	if( segments < 3 || radius <= 0 ) {				// erroneous conditions
		return;
	}

	glBegin(GL_POLYGON);

	for( int i=0; i<segments; i++ ) {				// for each segment
		double t = -(i/(double)segments) * 2*PI;
		WorldPt p = center + WorldPt(radius*cos(t), radius*sin(t));
		
		glVertex2d(p.x, p.y);
	}

	glEnd();
}
void GL_DrawCircle(const WorldPt & center, double radius, int segments) {
	if( segments < 3 || radius <= 0 ) {				// erroneous conditions
		return;
	}

	glBegin(GL_LINE_LOOP);

	for( int i=0; i<segments; i++ ) {				// for each segment
		double t = -(i/(double)segments) * 2*PI;
		WorldPt p = center + WorldPt(radius*cos(t), radius*sin(t));
		
		glVertex2d(p.x, p.y);
	}

	glEnd();
}
void GL_DrawDepthCuedLine(const WorldPt3 & origin, const WorldPt3 & dir, const WorldPt3 & localViewer, double maxAlpha, double minAlpha, double maxDist, double maxWidth, int segments) {
	double maxA = maxAlpha;
	double minA = minAlpha;

	double currColor[4];
	glGetDoublev(GL_CURRENT_COLOR, currColor);
	glPushAttrib(GL_LINE_BIT);

	for( int i=0; i<segments; i++ ) {				// for each line segment of depth cue'd line
		double t1 = i / (double)segments;
		double t2 = (i+1) / (double)segments;

		WorldPt3 p0 = origin + dir*t1;
		WorldPt3 p1 = origin + dir*t2;
		
		double d = p1.DistanceFrom(localViewer);
		d = maxDist / d;
		d = (d-1); // d is in 0 to 1
		
		float lineWidth = (float)(d*maxWidth+1);
		glLineWidth(lineWidth);

		double a = d*(maxA-minA)+minA;
		currColor[3] = a;
		
		glColor4dv(currColor);

		glBegin(GL_LINES);
		glVertex3dv((GLdouble *)&p0);
		glVertex3dv((GLdouble *)&p1);
		glEnd();
	}

	glPopAttrib();
}



WorldCube TransformAxesAlignedCube(WorldCube & worldCube, TransformMatrix & transform) {
	WorldCube cube;
	cube.Unite(PointTransformN(worldCube.PointXYZ(), transform));
	cube.Unite(PointTransformN(worldCube.Pointxyz(), transform));
	cube.Unite(PointTransformN(worldCube.PointXyz(), transform));
	cube.Unite(PointTransformN(worldCube.PointxYz(), transform));
	cube.Unite(PointTransformN(worldCube.PointxyZ(), transform));
	cube.Unite(PointTransformN(worldCube.PointxYZ(), transform));
	cube.Unite(PointTransformN(worldCube.PointXyZ(), transform));
	cube.Unite(PointTransformN(worldCube.PointXYz(), transform));

	return cube;
}

WorldCube InverseTransformAxesAlignedCube(WorldCube & worldCube, TransformMatrix & transform) {
	WorldCube cube;
	cube.Unite(InversePointTransformN(worldCube.PointXYZ(), transform));
	cube.Unite(InversePointTransformN(worldCube.Pointxyz(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointXyz(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointxYz(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointxyZ(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointxYZ(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointXyZ(), transform));
	cube.Unite(InversePointTransformN(worldCube.PointXYz(), transform));

	return cube;
}

float VEC3_DistanceBetween(float * inV1, float * inV2) {
	float diff[3];
	VEC3_Subtract(inV1, inV2, diff);
	return VEC3_Magnitude(diff);
}

float VEC3_DistanceBetween_Sq(float * inV1, float * inV2) {
	float diff[3];
	VEC3_Subtract(inV1, inV2, diff);
	return VEC3_Magnitude_Sq(diff);
}

void VEC3_Set(float x, float y, float z, float * outVec) {
	outVec[0] = x;
	outVec[1] = y;
	outVec[2] = z;
}
void VEC4_Set(float x, float y, float z, float w, float * outVec) {
	outVec[0] = x;
	outVec[1] = y;
	outVec[2] = z;
	outVec[3] = w;
}

void VEC3_Set(float * outVec, float x, float y, float z) {
	outVec[0] = x;
	outVec[1] = y;
	outVec[2] = z;
}
void VEC4_Set(float * outVec, float x, float y, float z, float w) {
	outVec[0] = x;
	outVec[1] = y;
	outVec[2] = z;
	outVec[3] = w;
}


void VEC4_Equate(float * inVec, float * outVec) {
	outVec[0] = inVec[0];
	outVec[1] = inVec[1];
	outVec[2] = inVec[2];
	outVec[3] = inVec[3];
}
void VEC3_Equate(float * vec3, float * output_vec3) {
	output_vec3[0] = vec3[0];
	output_vec3[1] = vec3[1];
	output_vec3[2] = vec3[2];
}
void VEC3_Negate(float * vec3) {
	vec3[0] = -vec3[0];
	vec3[1] = -vec3[1];
	vec3[2] = -vec3[2];
}
float VEC3_Dot(float * first_vec3, float * second_vec3) {
	return first_vec3[0]*second_vec3[0] + first_vec3[1]*second_vec3[1] + first_vec3[2]*second_vec3[2];
}

float VEC3_Magnitude_Sq(float * vec3) {
	return vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2];
}

float VEC3_Magnitude(float * vec3) {
	return sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2]);
}
void VEC3_Normalize(float * vec3, float * normalized_vec3) {
	float magnitude = VEC3_Magnitude(vec3);
	if( magnitude == 0 ) {
		magnitude = 1;
	}

	normalized_vec3[0] = vec3[0] / magnitude;
	normalized_vec3[1] = vec3[1] / magnitude;
	normalized_vec3[2] = vec3[2] / magnitude;
}

void VEC3_Normalize(float * vec3) {
	VEC3_Normalize(vec3, vec3);
}

void VEC3_Cross(float * a_vec3, float * b_vec3, float * crossed_vec3) {
	crossed_vec3[0] = a_vec3[1]*b_vec3[2] - a_vec3[2]*b_vec3[1];
	crossed_vec3[1] = -(a_vec3[0]*b_vec3[2] - a_vec3[2]*b_vec3[0]);
	crossed_vec3[2] = a_vec3[0]*b_vec3[1] - a_vec3[1]*b_vec3[0];
}

void VEC3_Difference(float * a_vec3, float * b_vec3, float * difference_vec3) {
	difference_vec3[0] = b_vec3[0] - a_vec3[0];
	difference_vec3[1] = b_vec3[1] - a_vec3[1];
	difference_vec3[2] = b_vec3[2] - a_vec3[2];
}

void VEC3_Swap(float * a_vec3, float * b_vec3) {
	float vec3[3];
	vec3[0] = a_vec3[0];
	vec3[1] = a_vec3[1];
	vec3[2] = a_vec3[2];

	a_vec3[0] = b_vec3[0];
	a_vec3[1] = b_vec3[1];
	a_vec3[2] = b_vec3[2];

	b_vec3[0] = vec3[0];
	b_vec3[1] = vec3[1];
	b_vec3[2] = vec3[2];
}

void VEC3_Add(float * first_vec3, float * second_vec3, float * output_vec3) {
	output_vec3[0] = first_vec3[0] + second_vec3[0];
	output_vec3[1] = first_vec3[1] + second_vec3[1];
	output_vec3[2] = first_vec3[2] + second_vec3[2];
}

void VEC3_Subtract(float * first_vec3, float * second_vec3, float * output_vec3) {
	output_vec3[0] = first_vec3[0] - second_vec3[0];
	output_vec3[1] = first_vec3[1] - second_vec3[1];
	output_vec3[2] = first_vec3[2] - second_vec3[2];
}

void VEC3_Multiply(float * first_vec3, float * second_vec3, float * output_vec3) {
	output_vec3[0] = first_vec3[0] * second_vec3[0];
	output_vec3[1] = first_vec3[1] * second_vec3[1];
	output_vec3[2] = first_vec3[2] * second_vec3[2];
}

void VEC3_Multiply(float * vec3, float scalar, float * output_vec3) {
	output_vec3[0] = vec3[0]*scalar;
	output_vec3[1] = vec3[1]*scalar;
	output_vec3[2] = vec3[2]*scalar;
}
void VEC4_Multiply(float * inVec4, float scalar, float * outVec) {
	outVec[0] = inVec4[0]*scalar;
	outVec[1] = inVec4[1]*scalar;
	outVec[2] = inVec4[2]*scalar;
	outVec[3] = inVec4[3]*scalar;
}

void VEC3_Divide(float * inVec3, float scalar, float * outVec) {
	outVec[0] = inVec3[0]/scalar;
	outVec[1] = inVec3[1]/scalar;
	outVec[2] = inVec3[2]/scalar;
}
void VEC4_Divide(float * inVec4, float scalar, float * outVec) {
	outVec[0] = inVec4[0]/scalar;
	outVec[1] = inVec4[1]/scalar;
	outVec[2] = inVec4[2]/scalar;
	outVec[3] = inVec4[3]/scalar;
}

float VEC3_AngleBetween(float * first_vec3, float * second_vec3) {
	float mag1 = VEC3_Magnitude(first_vec3);

	float mag2 = VEC3_Magnitude(second_vec3);

	float dot = VEC3_Dot(first_vec3, second_vec3);

	float angleBetween = acos(max(min(dot / (mag1*mag2), 1.0f), -1.0f) );

	return angleBetween;
}
void VEC3_MAT4_Multiply(float * row_vec3, float * mat4, float * column_vec3) {
	column_vec3[0] = row_vec3[ 0]*mat4[ 0] + row_vec3[ 1]*mat4[ 1] + row_vec3[ 2]*mat4[ 2] + 1*mat4[ 3];
	column_vec3[1] = row_vec3[ 0]*mat4[ 4] + row_vec3[ 1]*mat4[ 5] + row_vec3[ 2]*mat4[ 6] + 1*mat4[ 7];
	column_vec3[2] = row_vec3[ 0]*mat4[ 8] + row_vec3[ 1]*mat4[ 9] + row_vec3[ 2]*mat4[10] + 1*mat4[11];
}

void MAT4_Set(float in0, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float in8, float in9, float in10, float in11, float in12, float in13, float in14, float in15, float * outMat) {
	outMat[0] = in0;
	outMat[1] = in1;
	outMat[2] = in2;
	outMat[3] = in3;
	outMat[4] = in4;
	outMat[5] = in5;
	outMat[6] = in6;
	outMat[7] = in7;
	outMat[8] = in8;
	outMat[9] = in9;
	outMat[10] = in10;
	outMat[11] = in11;
	outMat[12] = in12;
	outMat[13] = in13;
	outMat[14] = in14;
	outMat[15] = in15;
}

void MAT4_ExtractCol3(float * mat4, int col, float * vec3) {
    vec3[0] = mat4[col*4+0];
    vec3[1] = mat4[col*4+1];
    vec3[2] = mat4[col*4+2];
}

void MAT4_ExtractRow3(float * mat4, int row, float * vec3) {
	vec3[0] = mat4[row];
	vec3[1] = mat4[row+4];
	vec3[2] = mat4[row+8];
}

void MAT4_Equate(float * inMat, float * outMat) {
	memcpy(outMat, inMat, sizeof(float)*16);
}

void MAT4_SetCol3(int col, float * vec3, float * mat4) {
	mat4[col*4+0] = vec3[0];
	mat4[col*4+1] = vec3[1];
	mat4[col*4+2] = vec3[2];
}

void MAT4_LoadIdentity(float * mat4) {
	mat4[0] = 1;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 1;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = 1;
	mat4[11] = 0;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = 0;
	mat4[15] = 1;
}

void MAT3_VEC3_Multiply(float * mat4, float * column_vec3, float * row_vec3) {
	float result[3];

	result[0] = mat4[ 0]*column_vec3[ 0] + mat4[ 4]*column_vec3[ 1] + mat4[ 8]*column_vec3[ 2];
	result[1] = mat4[ 1]*column_vec3[ 0] + mat4[ 5]*column_vec3[ 1] + mat4[ 9]*column_vec3[ 2];
	result[2] = mat4[ 2]*column_vec3[ 0] + mat4[ 6]*column_vec3[ 1] + mat4[10]*column_vec3[ 2];

	row_vec3[0] = result[0];
	row_vec3[1] = result[1];
	row_vec3[2] = result[2];
}
void MAT4_VEC3_Multiply(float * mat4, float * column_vec3, float * row_vec3) {
	row_vec3[0] = mat4[ 0]*column_vec3[ 0] + mat4[ 4]*column_vec3[ 1] + mat4[ 8]*column_vec3[ 2] + mat4[12]*1;
	row_vec3[1] = mat4[ 1]*column_vec3[ 0] + mat4[ 5]*column_vec3[ 1] + mat4[ 9]*column_vec3[ 2] + mat4[13]*1;
	row_vec3[2] = mat4[ 2]*column_vec3[ 0] + mat4[ 6]*column_vec3[ 1] + mat4[10]*column_vec3[ 2] + mat4[14]*1;
}
void MAT4_VEC4_MatrixVectorMultiply(float * mat4, float * column_vec4, float * out_vec4) {
	out_vec4[0] = mat4[ 0]*column_vec4[ 0] + mat4[ 4]*column_vec4[ 1] + mat4[ 8]*column_vec4[ 2] + mat4[12]*column_vec4[ 3];
	out_vec4[1] = mat4[ 1]*column_vec4[ 0] + mat4[ 5]*column_vec4[ 1] + mat4[ 9]*column_vec4[ 2] + mat4[13]*column_vec4[ 3];
	out_vec4[2] = mat4[ 2]*column_vec4[ 0] + mat4[ 6]*column_vec4[ 1] + mat4[10]*column_vec4[ 2] + mat4[14]*column_vec4[ 3];
	out_vec4[3] = mat4[ 3]*column_vec4[ 0] + mat4[ 7]*column_vec4[ 1] + mat4[11]*column_vec4[ 2] + mat4[15]*column_vec4[ 3];
}

void MAT4_Multiply(float * first_mat4, float * second_mat4, float * output_mat4) {
	float result[16];

	result[ 0] = first_mat4[ 0]*second_mat4[ 0] + first_mat4[ 4]*second_mat4[ 1] + first_mat4[ 8]*second_mat4[ 2] + first_mat4[12]*second_mat4[ 3];
	result[ 1] = first_mat4[ 1]*second_mat4[ 0] + first_mat4[ 5]*second_mat4[ 1] + first_mat4[ 9]*second_mat4[ 2] + first_mat4[13]*second_mat4[ 3];
	result[ 2] = first_mat4[ 2]*second_mat4[ 0] + first_mat4[ 6]*second_mat4[ 1] + first_mat4[10]*second_mat4[ 2] + first_mat4[14]*second_mat4[ 3];
	result[ 3] = first_mat4[ 3]*second_mat4[ 0] + first_mat4[ 7]*second_mat4[ 1] + first_mat4[11]*second_mat4[ 2] + first_mat4[15]*second_mat4[ 3];

	result[ 4] = first_mat4[ 0]*second_mat4[ 4] + first_mat4[ 4]*second_mat4[ 5] + first_mat4[ 8]*second_mat4[ 6] + first_mat4[12]*second_mat4[ 7];
	result[ 5] = first_mat4[ 1]*second_mat4[ 4] + first_mat4[ 5]*second_mat4[ 5] + first_mat4[ 9]*second_mat4[ 6] + first_mat4[13]*second_mat4[ 7];
	result[ 6] = first_mat4[ 2]*second_mat4[ 4] + first_mat4[ 6]*second_mat4[ 5] + first_mat4[10]*second_mat4[ 6] + first_mat4[14]*second_mat4[ 7];
	result[ 7] = first_mat4[ 3]*second_mat4[ 4] + first_mat4[ 7]*second_mat4[ 5] + first_mat4[11]*second_mat4[ 6] + first_mat4[15]*second_mat4[ 7];

	result[ 8] = first_mat4[ 0]*second_mat4[ 8] + first_mat4[ 4]*second_mat4[ 9] + first_mat4[ 8]*second_mat4[10] + first_mat4[12]*second_mat4[11];
	result[ 9] = first_mat4[ 1]*second_mat4[ 8] + first_mat4[ 5]*second_mat4[ 9] + first_mat4[ 9]*second_mat4[10] + first_mat4[13]*second_mat4[11];
	result[10] = first_mat4[ 2]*second_mat4[ 8] + first_mat4[ 6]*second_mat4[ 9] + first_mat4[10]*second_mat4[10] + first_mat4[14]*second_mat4[11];
	result[11] = first_mat4[ 3]*second_mat4[ 8] + first_mat4[ 7]*second_mat4[ 9] + first_mat4[11]*second_mat4[10] + first_mat4[15]*second_mat4[11];

	result[12] = first_mat4[ 0]*second_mat4[12] + first_mat4[ 4]*second_mat4[13] + first_mat4[ 8]*second_mat4[14] + first_mat4[12]*second_mat4[15];
	result[13] = first_mat4[ 1]*second_mat4[12] + first_mat4[ 5]*second_mat4[13] + first_mat4[ 9]*second_mat4[14] + first_mat4[13]*second_mat4[15];
	result[14] = first_mat4[ 2]*second_mat4[12] + first_mat4[ 6]*second_mat4[13] + first_mat4[10]*second_mat4[14] + first_mat4[14]*second_mat4[15];
	result[15] = first_mat4[ 3]*second_mat4[12] + first_mat4[ 7]*second_mat4[13] + first_mat4[11]*second_mat4[14] + first_mat4[15]*second_mat4[15];

	memcpy(output_mat4, result, sizeof(float)*16);
}

void MAT4_TransposeUpper3(float * mat4, float * transpose_mat4) {
	transpose_mat4[0] = mat4[0];
	transpose_mat4[1] = mat4[4];
	transpose_mat4[2] = mat4[8];
	transpose_mat4[3] = mat4[3];
	transpose_mat4[4] = mat4[1];
	transpose_mat4[5] = mat4[5];
	transpose_mat4[6] = mat4[9];
	transpose_mat4[7] = mat4[7];
	transpose_mat4[8] = mat4[2];
	transpose_mat4[9] = mat4[6];
	transpose_mat4[10] = mat4[10];
	transpose_mat4[11] = mat4[11];
	transpose_mat4[12] = mat4[3];
	transpose_mat4[13] = mat4[7];
	transpose_mat4[14] = mat4[11];
	transpose_mat4[15] = mat4[15];
}

void MAT4_Transpose(float * mat4) {
	swap(mat4[1], mat4[5]);
	swap(mat4[2], mat4[8]);
	swap(mat4[3], mat4[12]);
	swap(mat4[6], mat4[9]);
	swap(mat4[7], mat4[13]);
	swap(mat4[11], mat4[14]);
}

void MAT4_SetEqual(float * mat4, float * output_mat4) {
	for( int i=0; i<16; i++ ) {
		output_mat4[i] = mat4[i];
	}
}

void MAT4_InvertQuick(float * inMat, float * outMat) {
	memcpy(outMat, inMat, sizeof(float)*16);

	float offset[3];
	MAT4_ExtractCol3(inMat, 3, offset);
	VEC3_Negate(offset);

	outMat[12] = VEC3_Dot(offset, inMat+0);
	outMat[13] = VEC3_Dot(offset, inMat+4);
	outMat[14] = VEC3_Dot(offset, inMat+8);
	swap(outMat[1], outMat[4]);
	swap(outMat[2], outMat[8]);
	swap(outMat[6], outMat[9]);
}

float MTH_DegToRad(float deg) {
    return deg*3.141592654f/180.0f;
}

float MTH_RadToDeg(float rad) {
    return rad*180.0f/3.141592654f;
}

void MAT4_MakeZRotation(float deg, float * zRotation_mat4) {
	float rad = MTH_DegToRad(deg);

	float c = cos(rad);
	float s = sin(rad);
	zRotation_mat4[0] = c;
	zRotation_mat4[1] = s;
	zRotation_mat4[2] = 0;
	zRotation_mat4[3] = 0;
	zRotation_mat4[4] = -s;
	zRotation_mat4[5] = c;
	zRotation_mat4[6] = 0;
	zRotation_mat4[7] = 0;
	zRotation_mat4[8] = 0;
	zRotation_mat4[9] = 0;
	zRotation_mat4[10] = 1;
	zRotation_mat4[11] = 0;
	zRotation_mat4[12] = 0;
	zRotation_mat4[13] = 0;
	zRotation_mat4[14] = 0;
	zRotation_mat4[15] = 1;

}

void MAT4_MakeXRotation(float deg, float * xRotation_mat4) {
	float rad = MTH_DegToRad(deg);

	float c = cos(rad);
	float s = sin(rad);
	xRotation_mat4[0] = 1;
	xRotation_mat4[1] = 0;
	xRotation_mat4[2] = 0;
	xRotation_mat4[3] = 0;
	xRotation_mat4[4] = 0;
	xRotation_mat4[5] = c;
	xRotation_mat4[6] = s;
	xRotation_mat4[7] = 0;
	xRotation_mat4[8] = 0;
	xRotation_mat4[9] = -s;
	xRotation_mat4[10] = c;
	xRotation_mat4[11] = 0;
	xRotation_mat4[12] = 0;
	xRotation_mat4[13] = 0;
	xRotation_mat4[14] = 0;
	xRotation_mat4[15] = 1;

}
void MAT4_MakeYRotation(float deg, float * xRotation_mat4) {
	float rad = MTH_DegToRad(deg);

	float c = cos(rad);
	float s = sin(rad);
	xRotation_mat4[0] = c;
	xRotation_mat4[1] = 0;
	xRotation_mat4[2] = s;
	xRotation_mat4[3] = 0;
	xRotation_mat4[4] = 0;
	xRotation_mat4[5] = 1;
	xRotation_mat4[6] = 0;
	xRotation_mat4[7] = 0;
	xRotation_mat4[8] = -s;
	xRotation_mat4[9] = 0;
	xRotation_mat4[10] = c;
	xRotation_mat4[11] = 0;
	xRotation_mat4[12] = 0;
	xRotation_mat4[13] = 0;
	xRotation_mat4[14] = 0;
	xRotation_mat4[15] = 1;

}
void MAT4_MakeArbitraryRotation(float * axis_vec3, float angle_rad, float * rotation_mat4) {
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	float t = 1-c;
	float & x = axis_vec3[0];
	float & y = axis_vec3[1];
	float & z = axis_vec3[2];

	rotation_mat4[0] = t*x*x + c;
	rotation_mat4[1] = t*x*y - s*z;
	rotation_mat4[2] = t*x*y + s*y;
	rotation_mat4[3] = 0;
	rotation_mat4[4] = t*x*y + s*z;
	rotation_mat4[5] = t*y*y + c;
	rotation_mat4[6] = t*y*z - s*x;
	rotation_mat4[7] = 0;
	rotation_mat4[8] = t*x*z - s*y;
	rotation_mat4[9] = t*y*z + s*x;
	rotation_mat4[10] = t*z*z + c;
	rotation_mat4[11] = 0;
	rotation_mat4[12] = 0;
	rotation_mat4[13] = 0;
	rotation_mat4[14] = 0;
	rotation_mat4[15] = 1;
}

void MAT4_MakeInversePerspective(float fovy_rad, float aspect, float zNear, float zFar, float * mat4) {
	float f = 1.0f / tan(fovy_rad / 2);
	float dp = zNear - zFar;
	mat4[0] = aspect/f;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 1/f;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = 0;
	mat4[11] = dp/(2*zFar*zNear);
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = -1;
	mat4[15] =(zFar+zNear)/(2*zNear*zFar);
}

void MAT4_MakePerspective(float fovy_rad, float aspect, float zNear, float zFar, float * mat4) {
	float f = 1.0f / tan(fovy_rad / 2);
	mat4[0] = f/aspect;
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = f;
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = (zFar+zNear)/(zNear-zFar);
	mat4[11] = -1;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = (2*zFar*zNear)/(zNear-zFar);
	mat4[15] = 0;
}

void MAT4_MakeOrtho(float left, float right, float bottom, float top, float zNear, float zFar, float * mat4) {
	float tx = -(right+left) / (right-left);
	float ty = -(top+bottom) / (top-bottom);
	float tz = -(zFar+zNear) / (zFar-zNear);

	mat4[0] = 2.0f / (right-left);
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = 2.0f / (top-bottom);
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = 0;
	mat4[9] = 0;
	mat4[10] = -2.0f / (zFar-zNear);
	mat4[11] = 0;
	mat4[12] = tx;
	mat4[13] = ty;
	mat4[14] = tz;
	mat4[15] = 1;
}

void MAT4_MakeFrustum(float left, float right, float bottom, float top, float zNear, float zFar, float * mat4) {
	float A = (right+left) / (right-left);
	float B = (top+bottom) / (top-bottom);
	float C = -(zFar+zNear) / (zFar-zNear);
	float D = -(2*zFar*zNear) / (zFar-zNear);

	mat4[0] = (2*zNear) / (right-left);
	mat4[1] = 0;
	mat4[2] = 0;
	mat4[3] = 0;
	mat4[4] = 0;
	mat4[5] = (2*zNear) / (top-bottom);
	mat4[6] = 0;
	mat4[7] = 0;
	mat4[8] = A;
	mat4[9] = B;
	mat4[10] = C;
	mat4[11] = -1;
	mat4[12] = 0;
	mat4[13] = 0;
	mat4[14] = D;
	mat4[15] = 0;
}

void MAT4_LookAt(float eyeX, float eyeY, float eyeZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ, float * viewMatrix) {
	float * m = viewMatrix;

	float eye[] = { eyeX, eyeY, eyeZ };
	float center[] = { centerX, centerY, centerZ };
	float up[] = { upX, upY, upZ };
	VEC3_Normalize(up);

	float f[3];
	VEC3_Subtract(center, eye, f);
	VEC3_Normalize(f);

	float s[3];
	VEC3_Cross(f, up, s);
	VEC3_Normalize(s);

	float u[3];
	VEC3_Cross(s, f, u);

	float mat[16] = {
		s[0], u[0], -f[0], 0,
		s[1], u[1], -f[1], 0,
		s[2], u[2], -f[2], 0,
		0, 0, 0, 1 };

	float t[16];
	MAT4_LoadIdentity(t);
	MAT4_Translate(-eyeX, -eyeY, -eyeZ, t);

	MAT4_Multiply(mat, t, m);	
}


void MAT4_Translate(float inX, float inY, float inZ, float * inOutMat4) {
	inOutMat4[12] += inX;
	inOutMat4[13] += inY;
	inOutMat4[14] += inZ;
}

void MAT4_SetUpper3x3Rowwise(float * outMat, float * inR0, float * inR1, float * inR2) {
	outMat[0] = inR0[0];
	outMat[4] = inR0[1];
	outMat[8] = inR0[2];

	outMat[1] = inR1[0];
	outMat[5] = inR1[1];
	outMat[9] = inR1[2];

	outMat[2] = inR2[0];
	outMat[6] = inR2[1];
	outMat[10] = inR2[2];
}
void MAT4_SetUpper3x3(float * outMat, float * inR0, float * inR1, float * inR2) {
	outMat[0] = inR0[0];
	outMat[1] = inR0[1];
	outMat[2] = inR0[2];

	outMat[4] = inR1[0];
	outMat[5] = inR1[1];
	outMat[6] = inR1[2];

	outMat[8] = inR2[0];
	outMat[9] = inR2[1];
	outMat[10] = inR2[2];
}

void MAT4_AxisAngleToMatrix(float * inAxis, float inAngleDeg, float * outMat) 
/* AxisAngleToMatrix( ) - Generates a 4x4 rotation matrix representing rotation about 'axis' 
		'axisDeg' degrees.
	axis - axis of rotation.  Does not need to be normalized.  Should have magnitude > 0.
	angleDeg - angle of rotation about 'axis' in degrees.  Any degree is valid.
	return - a 4x4 rotation matrix representing rotation about 'axis' 'axisDeg' degrees.
*/		
{
	MAT4_LoadIdentity(outMat);

	if( VEC3_Magnitude_Sq(inAxis) < 1e-3 ) {
		return;
	}

	float axis[3];
	VEC3_Equate(inAxis, axis);
	VEC3_Normalize(axis);

	float angleRad = MTH_DegToRad(inAngleDeg);
	float c = cosf(angleRad);
	float s = sinf(angleRad);
	float t = 1-c;

	// variable rename
	float & x = axis[0];
	float & y = axis[1];
	float & z = axis[2];

	// equation from http://www.gamedev.net/reference/articles/article1199.asp
	outMat[ 0] = t*x*x + c;
	outMat[ 1] = t*x*y + s*z;
	outMat[ 2] = t*x*z - s*y;
	outMat[ 4] = t*x*y - s*z;
	outMat[ 5] = t*y*y + c;
	outMat[ 6] = t*y*z + s*x;
	outMat[ 8] = t*x*z + s*y;
	outMat[ 9] = t*y*z - s*x;
	outMat[10] = t*z*z + c;
}

void QT_Equate(float * inQ, float * outQ) {
	memcpy(outQ, inQ, sizeof(float)*4);
}

void QT_Set(float inX, float inY, float inZ, float inW, float * outQ) {
	outQ[0] = inX;
	outQ[1] = inY;
	outQ[2] = inZ;
	outQ[3] = inW;
}

void QT_MAT4_MakeQuaternion(float * inMat, float * outQ) {
	float * m = inMat;
	
	float & x = outQ[0];
	float & y = outQ[1];
	float & z = outQ[2];
	float & w = outQ[3];

	float trace = m[0] + m[5] + m[10] + 1;
	float s = 0;
	// if trace is positive then perform instant calculation
	if( trace > 0 ) {
		s = 0.5f / sqrt(trace);
		w = 0.25f / s;
		x = ( m[9] - m[6] ) * s;
		y = ( m[2] - m[8] ) * s;
		z = ( m[4] - m[1] ) * s;
	}
	else {
		int majorColumn = 0;
		if( m[5] > m[majorColumn*5] ) majorColumn = 1;
		if( m[10] > m[majorColumn*5] ) majorColumn = 2;

		switch( majorColumn ) {
			case 0:
				s = sqrt( 1.0f + m[0] - m[5] - m[10] ) * 2.0f;
				x = 0.25f * s;
				y = (m[4] + m[1] ) / s;
				z = (m[2] + m[8] ) / s;
				w = (m[9] - m[6] ) / s;
				break;

			case 1:
				s = sqrt( 1.0f + m[5] - m[0] - m[10] ) * 2.0f;
				x = (m[4] + m[1] ) / s;
				y = 0.25f * s;
				z = (m[9] + m[6] ) / s;
				w = (m[2] - m[8] ) / s;
				break;

			case 2:
				s = sqrt( 1.0f + m[10] - m[0] - m[5] ) * 2.0f;
				x = (m[2] + m[8] ) / s;
				y = (m[9] + m[6] ) / s;
				z = 0.25f * s;
				w = (m[4] - m[1] ) / s;
				break;
		}
	}
}

void QT_MAT4_MakeMatrix(float * inQ, float * outMat) {
	float * m = outMat;
	
	float & x = inQ[0];
	float & y = inQ[1];
	float & z = inQ[2];
	float & w = inQ[3];

	float xx = x * x;
	float xy = x * y;
	float xz = x * z;
	float xw = x * w;

	float yy = y * y;
	float yz = y * z;
	float yw = y * w;

	float zz = z * z;
	float zw = z * w;

	m[0] = 1 - 2 * ( yy + zz );
	m[1] = 2 * ( xy - zw );
	m[2] = 2 * ( xz + yw );
	m[3] = 0;

	m[4] = 2 * ( xy + zw );
	m[5] = 1 - 2 * ( xx + zz );
	m[6] = 2 * ( yz - xw );
	m[7] = 0;

	m[8] = 2 * ( xz - yw );
	m[9] = 2 * ( yz + xw );
	m[10] = 1 - 2 * ( xx + yy );
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;
}

void QT_GetAxisAngle(float * inQ, float * outAxis, float outAngle /*rad*/) {
	outAngle = 2.0f*acosf(inQ[3]);
	float sinHalfAngle = sinf(outAngle/2.0f);
	VEC3_Set(inQ[0] / sinHalfAngle, inQ[1] / sinHalfAngle, inQ[2] / sinHalfAngle, outAxis);
}
void QT_SLERP(float * inQ0, float * inQ1, float inT, float * outQ) {

	float * q0 = inQ0;
	float * q1 = inQ1;

	float cosom = q0[0]*q1[0] + q0[1]*q1[1] + q0[2]*q1[2] + q0[3]*q1[3];
	float delta = 1e-10f;

	float adjustedEnd[4];
	if( cosom < 0.0f ) {
		QT_Set(-q1[0], -q1[1], -q1[2], -q1[3], adjustedEnd);
		cosom *= -1.0f;
	}
	else {
		QT_Equate(q1, adjustedEnd);
	}

	float scaleBegin, scaleEnd;

	if( 1.0f-cosom > delta ) { // normal slerp
		float omega = acosf(cosom);
		float sinom = sinf(omega);
		scaleBegin = sinf((1.0f-inT)*omega) / sinom;
		scaleEnd = sinf(inT * omega) / sinom;
	} else { // linear interp because of divide by zero
		scaleBegin = 1.0f - inT;
		scaleEnd = inT;
	}

	QT_Set(
		scaleBegin*q0[0] + scaleEnd*adjustedEnd[0],
		scaleBegin*q0[1] + scaleEnd*adjustedEnd[1],
		scaleBegin*q0[2] + scaleEnd*adjustedEnd[2],
		scaleBegin*q0[3] + scaleEnd*adjustedEnd[3],
		outQ);
}


void GLU_LookAt(float eyeX, float eyeY, float eyeZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ, float * viewMatrix) {
	float * m = viewMatrix;

	float eye[] = { eyeX, eyeY, eyeZ };
	float center[] = { centerX, centerY, centerZ };
	float up[] = { upX, upY, upZ };
	VEC3_Normalize(up);

	float f[3];
	VEC3_Subtract(center, eye, f);
	VEC3_Normalize(f);

	float s[3];
	VEC3_Cross(f, up, s);

	float u[3];
	VEC3_Cross(s, f, u);

	float mat[16] = {
		s[0], u[0], -f[0], 0,
		s[1], u[1], -f[1], 0,
		s[2], u[2], -f[2], 0,
		0, 0, 0, 1 };

	float t[16];
	MAT4_LoadIdentity(t);
	MAT4_Translate(-eyeX, -eyeY, -eyeZ, t);

	MAT4_Multiply(mat, t, m);	
}


void MAT4_VEC4_Multiply(float * mat4, float * column_vec4, float * row_vec4) {
	float result[4];

	result[0] = mat4[ 0]*column_vec4[ 0] + mat4[ 4]*column_vec4[ 1] + mat4[ 8]*column_vec4[ 2] + mat4[12]*column_vec4[ 3];
	result[1] = mat4[ 1]*column_vec4[ 0] + mat4[ 5]*column_vec4[ 1] + mat4[ 9]*column_vec4[ 2] + mat4[13]*column_vec4[ 3];
	result[2] = mat4[ 2]*column_vec4[ 0] + mat4[ 6]*column_vec4[ 1] + mat4[10]*column_vec4[ 2] + mat4[14]*column_vec4[ 3];
	result[3] = mat4[ 3]*column_vec4[ 0] + mat4[ 7]*column_vec4[ 1] + mat4[11]*column_vec4[ 2] + mat4[15]*column_vec4[ 3];

	memcpy(row_vec4, result, 4*sizeof(float));
}

bool MAT4_Invert(const float * mat4, float * outMat) {
	/* Compute inverse of 4x4 transformation matrix.
	Code contributed by Jacques Leroy jle@star.be
	return - true for success, false for failure (singular matrix)
	
	modified from http://webcvs.freedesktop.org/mesa/Mesa/src/glu/mesa/project.c?revision=1.4&view=markup
	*/

/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { float *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(mat4,r,c) (mat4)[(c)*4+(r)]

	float * out4 = outMat;
	float inv[16];
	if( mat4 == out4 ) { // if input and output point to same memory
		out4 = inv;
	}

   float wtmp[4][8];
   float m0, m1, m2, m3, s;
   float *r0, *r1, *r2, *r3;

   r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

   r0[0] = MAT(mat4, 0, 0), r0[1] = MAT(mat4, 0, 1),
      r0[2] = MAT(mat4, 0, 2), r0[3] = MAT(mat4, 0, 3),
      r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,
      r1[0] = MAT(mat4, 1, 0), r1[1] = MAT(mat4, 1, 1),
      r1[2] = MAT(mat4, 1, 2), r1[3] = MAT(mat4, 1, 3),
      r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,
      r2[0] = MAT(mat4, 2, 0), r2[1] = MAT(mat4, 2, 1),
      r2[2] = MAT(mat4, 2, 2), r2[3] = MAT(mat4, 2, 3),
      r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,
      r3[0] = MAT(mat4, 3, 0), r3[1] = MAT(mat4, 3, 1),
      r3[2] = MAT(mat4, 3, 2), r3[3] = MAT(mat4, 3, 3),
      r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

   /* choose pivot - or die */
   if (fabs(r3[0]) > fabs(r2[0]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[0]) > fabs(r1[0]))
      SWAP_ROWS(r2, r1);
   if (fabs(r1[0]) > fabs(r0[0]))
      SWAP_ROWS(r1, r0);
   if (0.0 == r0[0])
      return false;

   /* eliminate first variable     */
   m1 = r1[0] / r0[0];
   m2 = r2[0] / r0[0];
   m3 = r3[0] / r0[0];
   s = r0[1];
   r1[1] -= m1 * s;
   r2[1] -= m2 * s;
   r3[1] -= m3 * s;
   s = r0[2];
   r1[2] -= m1 * s;
   r2[2] -= m2 * s;
   r3[2] -= m3 * s;
   s = r0[3];
   r1[3] -= m1 * s;
   r2[3] -= m2 * s;
   r3[3] -= m3 * s;
   s = r0[4];
   if (s != 0.0) {
      r1[4] -= m1 * s;
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r0[5];
   if (s != 0.0) {
      r1[5] -= m1 * s;
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r0[6];
   if (s != 0.0) {
      r1[6] -= m1 * s;
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r0[7];
   if (s != 0.0) {
      r1[7] -= m1 * s;
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[1]) > fabs(r2[1]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[1]) > fabs(r1[1]))
      SWAP_ROWS(r2, r1);
   if (0.0 == r1[1])
      return GL_FALSE;

   /* eliminate second variable */
   m2 = r2[1] / r1[1];
   m3 = r3[1] / r1[1];
   r2[2] -= m2 * r1[2];
   r3[2] -= m3 * r1[2];
   r2[3] -= m2 * r1[3];
   r3[3] -= m3 * r1[3];
   s = r1[4];
   if (0.0 != s) {
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r1[5];
   if (0.0 != s) {
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r1[6];
   if (0.0 != s) {
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r1[7];
   if (0.0 != s) {
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[2]) > fabs(r2[2]))
      SWAP_ROWS(r3, r2);
   if (0.0 == r2[2])
      return GL_FALSE;

   /* eliminate third variable */
   m3 = r3[2] / r2[2];
   r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
      r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

   /* last check */
   if (0.0 == r3[3])
      return GL_FALSE;

   s = 1.0f / r3[3];		/* now back substitute row 3 */
   r3[4] *= s;
   r3[5] *= s;
   r3[6] *= s;
   r3[7] *= s;

   m2 = r2[3];			/* now back substitute row 2 */
   s = 1.0f / r2[2];
   r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
      r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
   m1 = r1[3];
   r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
      r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
   m0 = r0[3];
   r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
      r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

   m1 = r1[2];			/* now back substitute row 1 */
   s = 1.0f / r1[1];
   r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
      r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
   m0 = r0[2];
   r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
      r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

   m0 = r0[1];			/* now back substitute row 0 */
   s = 1.0f / r0[0];
   r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
      r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

   MAT(out4, 0, 0) = r0[4];
   MAT(out4, 0, 1) = r0[5], MAT(out4, 0, 2) = r0[6];
   MAT(out4, 0, 3) = r0[7], MAT(out4, 1, 0) = r1[4];
   MAT(out4, 1, 1) = r1[5], MAT(out4, 1, 2) = r1[6];
   MAT(out4, 1, 3) = r1[7], MAT(out4, 2, 0) = r2[4];
   MAT(out4, 2, 1) = r2[5], MAT(out4, 2, 2) = r2[6];
   MAT(out4, 2, 3) = r2[7], MAT(out4, 3, 0) = r3[4];
   MAT(out4, 3, 1) = r3[5], MAT(out4, 3, 2) = r3[6];
   MAT(out4, 3, 3) = r3[7];

	if( mat4 == outMat ) { // if input and output point to same memory
		memcpy(outMat, out4, sizeof(float)*16);
	}

   return true;

#undef MAT
#undef SWAP_ROWS
}
void MAT4_Transpose(float * inMat4, float * outMat4) {
	outMat4[ 0] = inMat4[ 0];
	outMat4[ 1] = inMat4[ 4];
	outMat4[ 2] = inMat4[ 8];
	outMat4[ 3] = inMat4[12];
	outMat4[ 4] = inMat4[ 1];
	outMat4[ 5] = inMat4[ 5];
	outMat4[ 6] = inMat4[ 9];
	outMat4[ 7] = inMat4[13];
	outMat4[ 8] = inMat4[ 2];
	outMat4[ 9] = inMat4[ 6];
	outMat4[10] = inMat4[10];
	outMat4[11] = inMat4[14];
	outMat4[12] = inMat4[ 3];
	outMat4[13] = inMat4[ 7];
	outMat4[14] = inMat4[11];
	outMat4[15] = inMat4[15];
}
void MAT4_ExtractUpper3x3(float * inMat4, float * outMat3) {
	outMat3[ 0] = inMat4[ 0];
	outMat3[ 1] = inMat4[ 1];
	outMat3[ 2] = inMat4[ 2];
	outMat3[ 3] = inMat4[ 4];
	outMat3[ 4] = inMat4[ 5];
	outMat3[ 5] = inMat4[ 6];
	outMat3[ 6] = inMat4[ 8];
	outMat3[ 7] = inMat4[ 9];
	outMat3[ 8] = inMat4[10];
}

float MTH_GetFraction(float inF) {
	return inF - (int)inF;
}

int ExtractFrustumCorners(float * inProjectionMat, float * inModelviewMat, float * outCorners) {
	if( !inProjectionMat || !inModelviewMat || !outCorners ) {
		return 1;
	}

	float MPMat[16];
	MAT4_Multiply(inProjectionMat, inModelviewMat, MPMat);
	MAT4_Invert(MPMat, MPMat);

	float corners[] = { 
		-1, -1, -1,
		 1, -1, -1,
		 1,  1, -1,
		-1,  1, -1,
		-1, -1,  1,
		 1, -1,  1,
		 1,  1,  1,
		-1,  1,  1
	};

	MAT4_VEC4_Multiply(MPMat, corners+ 0, corners+ 0);
	MAT4_VEC4_Multiply(MPMat, corners+ 3, corners+ 3);
	MAT4_VEC4_Multiply(MPMat, corners+ 6, corners+ 6);
	MAT4_VEC4_Multiply(MPMat, corners+ 9, corners+ 9);
	MAT4_VEC4_Multiply(MPMat, corners+12, corners+12);
	MAT4_VEC4_Multiply(MPMat, corners+15, corners+15);
	MAT4_VEC4_Multiply(MPMat, corners+18, corners+18);
	MAT4_VEC4_Multiply(MPMat, corners+21, corners+21);

	memcpy(outCorners, corners, sizeof(float)*24);

	return 0;
}

void FPSCAM_MoveForward(float amount, bool fixedZ, float * inOutViewMat) {
	float offset[4];

	if( fixedZ ) {
	
		offset[0] = inOutViewMat[2];
		offset[1] = inOutViewMat[6];
		offset[2] = 0;
		offset[3] = 0;
		MAT4_VEC4_Multiply(inOutViewMat, offset, offset);
		VEC3_Normalize(offset);
		VEC3_Multiply(offset, amount, offset);

	} else {
		VEC3_Set(offset, 0, 0, -amount);
	}
	inOutViewMat[12] += offset[0];
	inOutViewMat[13] += offset[1];
	inOutViewMat[14] += offset[2];
	
}

void FPSCAM_MoveUp(float amount, bool fixedZ, float * inOutViewMat) {
	float offset[3];
	
	if( fixedZ ) {
		MAT4_ExtractCol3(inOutViewMat, 2, offset);
		VEC3_Multiply(offset, -amount, offset);

	} else {
		VEC3_Set(offset, 0, -amount, 0);
	}

	inOutViewMat[12] += offset[0];
	inOutViewMat[13] += offset[1];
	inOutViewMat[14] += offset[2];
}

void FPSCAM_StrafeRight(float amount, float * inOutViewMat) {
	float offset[3] = { -amount, 0, 0 };

	inOutViewMat[12] += offset[0];
	inOutViewMat[13] += offset[1];
	inOutViewMat[14] += offset[2];
}

void FPSCAM_LookRight(float amount_deg, bool inFixedUp, float * inOutViewMat) {
    float rotMat[16];
	float rotVec[4];

	if( inFixedUp ) {
		MAT4_ExtractCol3(inOutViewMat, 2, rotVec);
	} else {
		VEC3_Set(rotVec, 0, 1, 0);
	}

	MAT4_MakeArbitraryRotation(rotVec, MTH_DegToRad(-amount_deg), rotMat);
	MAT4_Multiply(rotMat, inOutViewMat, inOutViewMat);
}

void FPSCAM_LookUp(float amount_deg, float * inOutViewMat) {

    float xRotation_mat4[16];
	MAT4_MakeXRotation(-amount_deg, xRotation_mat4);

	MAT4_Multiply(xRotation_mat4, inOutViewMat, inOutViewMat);
}

void GenerateSphere(float inRadius, vector<float> & vertices, vector<float> & normals) {
	const int Size = (int)((3.141592654 / .2 + 1)*(2*3.141592654 / .2 + 1)) * 6;
	vertices.resize(Size);
	normals.resize(Size);

	int i = 0;

	float radius = inRadius;
	for( float phi = 0.2f; phi <= 3.141592654f; phi += 0.2f ) {
		for( float theta = 0.0f; theta <= 2.0f*3.141592654f; theta += 0.2f ) {
			float px = radius * cosf(theta) * sinf(phi-0.2f);
			float py = radius * sinf(theta) * sinf(phi-0.2f);
			float pz = radius * cosf(phi-0.2f);
			float pmagnitude = sqrtf(px*px + py*py + pz*pz);
			
			normals[i] = px/pmagnitude;
			vertices[i++] = px;
			normals[i] = py/pmagnitude;
			vertices[i++] = py;
			normals[i] = pz/pmagnitude;
			vertices[i++] = pz;

            float x = radius * cosf(theta) * sinf(phi);
			float y = radius * sinf(theta) * sinf(phi);
			float z = radius * cosf(phi);
			float magnitude = sqrtf(x*x + y*y + z*z);

			normals[i] = x/magnitude;
			vertices[i++] = x;
			normals[i] = y/magnitude;
			vertices[i++] = y;
			normals[i] = z/magnitude;
			vertices[i++] = z;
		}
	}
	
}

void ExtractFrustumPoints(float * inVPInvertMat, float outFrustumPoints[3*8])
// frustum points are coordinates in world space in this order: NearBottomLeft, NTopL, NTRight, NBR, FarBL, FTL, FTR, FBR
{
	float ndcFrustumPts[] = {
		-1,-1,-1, 1,  -1, 1,-1, 1,  1, 1,-1, 1,  1,-1,-1,1, // near
		-1,-1, 1, 1,  -1, 1, 1, 1,  1, 1, 1, 1,  1,-1, 1, 1 // far
	};

	for( int i=0; i<8; i++ ) {
		float p[4];
		MAT4_VEC4_Multiply(inVPInvertMat, &ndcFrustumPts[i*4], p);
		VEC3_Divide(p, p[3], &outFrustumPoints[i*3]);
	}

}

bool RayPlaneIntersection(float * ptOnPlane, float * planeNormal, float * rayOrigin, float * rayDir, float * optOutT, float * outIntersection) {
	float * p = ptOnPlane;
	float * n = planeNormal;
	float * o = rayOrigin;
	float d[3];
	VEC3_Normalize(rayDir, d);

	
	float den = VEC3_Dot(n, d);
	if( den == 0 ) return false;

	float vec[3];
	VEC3_Subtract(o, p, vec);
	float num = -VEC3_Dot(n, vec);

	float t = num / den;
	if( optOutT != 0 ) {
		*optOutT = t;
	}
	
	VEC3_Multiply(d, t, vec);
	VEC3_Add(o, vec, outIntersection);

	return true;
}

void GL_DrawWireCube(WorldCube & cube) {
	
	WorldPt3 p0 = cube.Pointxyz();
	WorldPt3 p1 = cube.PointXYZ();


	glBegin(GL_LINE_LOOP);
	glVertex3d(p1.x, p1.y, p1.z);
	glVertex3d(p0.x, p1.y, p1.z);
	glVertex3d(p0.x, p0.y, p1.z);
	glVertex3d(p1.x, p0.y, p1.z);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3d(p1.x, p1.y, p0.z);
	glVertex3d(p0.x, p1.y, p0.z);
	glVertex3d(p0.x, p0.y, p0.z);
	glVertex3d(p1.x, p0.y, p0.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3d(p0.x, p0.y, p0.z);
	glVertex3d(p0.x, p0.y, p1.z);
	
	glVertex3d(p1.x, p0.y, p0.z);
	glVertex3d(p1.x, p0.y, p1.z);
	
	glVertex3d(p1.x, p1.y, p0.z);
	glVertex3d(p1.x, p1.y, p1.z);
	
	glVertex3d(p0.x, p1.y, p0.z);
	glVertex3d(p0.x, p1.y, p1.z);
	glEnd();

}

int DBG_glGetError() 
{
	int err = glGetError();
	char message[256];
	sprintf_s(message, "glError() == %d\n", err);
	DBG_Message(false, message);

	return err;
}

void DBG_Message(bool inAlert, char * inMessage, ...) 
// args are in printf format
{
	va_list args;
	va_start(args, inMessage);
	char dbgMessage[2048];
	vsprintf_s(dbgMessage, sizeof(dbgMessage), inMessage, args);
	
	OutputDebugStringA(dbgMessage);
	
	if( inAlert ) {
		MessageBoxA(GetActiveWindow(), dbgMessage, "Alert", MB_OK);
	}

	//DebugBreak();
	va_end(args);
}

void GL_SaveTextureAsRAW(GLuint inTex, GLenum inType, size_t inFilesize, const char * inFilename, ...) 
// saves an opengl texture as a .RAW image file 
// args are in printf format
{
	int tex;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &tex);

	glBindTexture(GL_TEXTURE_2D, inTex);

	static vector<unsigned char> bytes;
	bytes.resize(inFilesize);
	glGetTexImage(GL_TEXTURE_2D, 0, inType, GL_UNSIGNED_BYTE, &bytes[0]);
	char filename[128];
	va_list args;
	va_start(args, inFilename);
	vsprintf_s(filename, 128, inFilename, args);
	ofstream out(filename, ios_base::binary);
	out.write((const char *)&bytes[0], inFilesize);
	out.close();
	va_end(args);
	
	glBindTexture(GL_TEXTURE_2D, tex);
}

void GL_SaveBitsAsRAW(float * inBits, size_t inFilesize, const char * inFilename, ...) 
// saves an opengl texture as a .RAW image file 
// args are in printf format
{
	static vector<unsigned char> bytes;
	bytes.resize(inFilesize);
	for_i( inFilesize ) {
		unsigned char byte =unsigned char(max(min(inBits[i], 0.0f), 1.0f) * 255);
		bytes[i] = byte;
	}

	char filename[128];
	va_list args;
	va_start(args, inFilename);
	vsprintf_s(filename, 128, inFilename, args);
	ofstream out(filename, ios_base::binary);
	out.write((const char *)&bytes[0], inFilesize);
	out.close();
	va_end(args);
}

void GL_SaveBufferAsRAW(GLenum inFormat, const char * inFilename, ...) 
// saves an opengl buffer as a .RAW image file 
// args are in printf format
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	static vector<float> depths1;
	static vector<float> depths2;
	static vector<float> * oldDepths = &depths1;
	static vector<float> * newDepths = &depths2;

	newDepths->resize(viewport[2]*viewport[3]);
	glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], inFormat, GL_FLOAT, &(*newDepths)[0]);

	bool depthsDiffer = false;
	if( oldDepths->size() == newDepths->size() ) {
		for_i( newDepths->size() ) {
			float diff = (*newDepths)[i] - (*oldDepths)[i];
			if( diff != 0 ) {
				depthsDiffer = true;
				break;
			}
		}
	}
	swap(oldDepths, newDepths);

	

	static vector<unsigned char> bytes;
	bytes.resize(viewport[2]*viewport[3]);
	glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], inFormat, GL_UNSIGNED_BYTE, &bytes[0]);
	
	char filename[128];
	va_list args;
	va_start(args, inFilename);
	vsprintf_s(filename, 128, inFilename, args);
	ofstream out(filename, ios_base::binary);
	out.write((const char *)&bytes[0], bytes.size());
	out.close();
	va_end(args);
}

static bool FILE_Read(const char * inFilename, vector<char> & outFileContents) {

	ifstream file(inFilename, ios::binary);
	if( !file.is_open() ) {
		DBG_Message(true, "failed to open %s\n", inFilename);
		return false;
	} else {
		file.seekg(0, ios::end);
		std::streamoff size = file.tellg();
		file.seekg(0, ios::beg);
		outFileContents.resize(size_t(size+1));
		file.read(&outFileContents[0], size);
		outFileContents[size_t(size)] = 0;
		file.close();
		return true;
	}
}

static int GL_LoadGLSLProgramWithArgs(const char * inVS, const char * inFS, GLuint * outProgram, va_list & inArgs)
//  e.g.,
//  GFX_LoadGLSLProgram(vs, fs, &programID, 
//      eGLSLBindingAttribute, "inPosition", &position,
//      eGLSLBindingAttribute, "inNormal", &normal,
//      eGLSLBindingUniform, "kUniform1", &uniform1,
//      eGLSLBindingUniform, "kUniform2", &uniform2,
//      eGLSLBindingEnd);
{
    int errorCode = 0;
    
    char infoLog[1024];
    int len;
    GLint compileStatus;
    
    GLuint vShader = 0;
    GLuint fShader = 0;
    
    if( !errorCode ) {
        vShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vShader, 1, &inVS, 0);
        glCompileShader(vShader);
    
        glGetShaderInfoLog(vShader, 1024, &len, infoLog);
        glGetShaderiv(vShader, GL_COMPILE_STATUS, &compileStatus);
        if( compileStatus == GL_FALSE ) {
			DBG_Message(true, infoLog);
            errorCode = 1;
        }
    }
    
    if( !errorCode ) {
        fShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fShader, 1, &inFS, 0);
        glCompileShader(fShader);
    
        glGetShaderInfoLog(fShader, 1024, &len, infoLog);
        glGetShaderiv(fShader, GL_COMPILE_STATUS, &compileStatus);
        if( compileStatus == GL_FALSE ) {
			DBG_Message(true, infoLog);
            errorCode = 2;
        }
    }
    
    if( !errorCode ) {
        GLint linkStatus;
        
        *outProgram = glCreateProgram();
        glAttachShader(*outProgram, vShader);
        glAttachShader(*outProgram, fShader);
        
        glLinkProgram(*outProgram);
        
        glGetProgramInfoLog(*outProgram, 1024, &len, infoLog);
        glGetProgramiv(*outProgram, GL_LINK_STATUS, &linkStatus);
        if( linkStatus == GL_FALSE ) {
			DBG_Message(true, infoLog);
            errorCode = 3;  // failed to link
        }
    }
    if( !errorCode ) {
        typedef pair<char *, int *> GLSLIdentifier;
    
        vector<GLSLIdentifier> uniforms;
        vector<GLSLIdentifier> attribs;
    
        EGLSLBinding bindingType;
        while( (bindingType = va_arg(inArgs, EGLSLBinding)) != eGLSLBindingEnd && !errorCode ) 
        // extract the vertex attrib and uniform names so we can get their id later
        {
        
            vector<GLSLIdentifier> & array = (bindingType==eGLSLBindingAttribute) ? attribs : uniforms;
        
            char * name = va_arg(inArgs, char *);
            int * id = va_arg(inArgs, int *);
            
            if( !id ) {
                errorCode = 4;  // bad attrib parameter
            }
            array.push_back(GLSLIdentifier(name, id));
        }
            
        
        for( size_t i=0; i<attribs.size() && !errorCode; i++ ) {
            char * name = attribs[i].first;
            int * id = attribs[i].second;
            
            *id = glGetAttribLocation(*outProgram, name);
            //if( *id == -1 ) {
            //    errorCode = 5;  // attrib not found
            //}
            
        }

    
        for( size_t i=0; i<uniforms.size() && !errorCode; i++ ) {
            char * name = uniforms[i].first;
            int * id = uniforms[i].second;
        
            *id = glGetUniformLocation(*outProgram, name);
            //if( *id == -1 ) {
            //    errorCode = 6;  // uniform not found
            //}
        }
		        
    }
    
    if( errorCode > 0 ) {
        glDeleteProgram(*outProgram);
        *outProgram = 0;
        glDeleteShader(vShader);
        glDeleteShader(fShader);
    }
    
    return errorCode;
    
}

int GL_GenerateVA(GLuint * outVBO, GLuint * outVA, ...)
//  e.g.,
//  vector<float> vertices, normals, uvs;
//	vector<int> indices;
//	GenerateGeometry(vertices, normals, uvs, indices);
//
//	GLuint vbo, va;
//	int indicesOffset;
//  if( GL_GenerateVA(&vbo, &va, 
//			eVAAttrib, vertexPositionLoc, &vertices, 3,
//			eVAAttrib, vertexNormalLoc, &normals, 3,
//			eVAAttrib, vertexUVLoc, &uvs, 2,
//			eVAIndices, &indices, &indicesOffset,
//			eVAEnd) != 0 )
//		HandleError();
//	else
//		gDrawElements(GL_TRIANGLE_STRIP, indices.size(), GL_UNSIGNED_INT, (void *)indicesOffset);
{
    int errorCode = 0;
    
	va_list args;
	va_start(args, outVA);

	glGetError();

	if( glIsVertexArray(*outVA) ) {
		glDeleteVertexArrays(1, outVA);
	}

	glGenVertexArrays(1, outVA);
	glBindVertexArray(*outVA);

	if( glIsBuffer(*outVBO) ) {
		glDeleteBuffers(1, outVBO);
	}

	glGenBuffers(1, outVBO);
	
	if( glGetError() ) {
		DBG_Message(true, "GL_GenreateVA:  Error generating vertex array and/or vertex buffer objects");
		errorCode = 1;
	}

	struct ArrayInstruction {
		ArrayInstruction() {
			location = 0;
			array = 0;
			size = 0;
		}

		int location;
		vector<float> * array;
		int size;
	};

	struct IndexInstruction {
		IndexInstruction() {
			array = 0;
			bufferOffset = 0;
		}

		vector<int> * array;
		int size;
		unsigned int * bufferOffset;
	} indices;


	vector<ArrayInstruction> arrays;
	unsigned int vboSize = 0;

    if( !errorCode ) {
		EVA vaType;
        while( (vaType = va_arg(args, EVA)) != eVAEnd && !errorCode ) 
        // extract the data
        {
			if( vaType == eVAAttrib ) {
				ArrayInstruction instruction;
				instruction.location = va_arg(args, int);
				instruction.array = va_arg(args, vector<float> *);
				instruction.size = va_arg(args, int);
				arrays.push_back(instruction);
				vboSize += instruction.array->size() * sizeof(float);
			} else if( vaType == eVAIndices ) {
				if( indices.array != 0 ) {
					DBG_Message(true, "GL_GenerateVA: indices can only be specified once\n");
					errorCode = 2;
				}
				indices.array = va_arg(args, vector<int> *);
				indices.bufferOffset = va_arg(args, unsigned int *);
				vboSize += indices.array->size() * sizeof(int);
			}
        }
	}

	if( !errorCode ) {
		glBindBuffer(GL_ARRAY_BUFFER, *outVBO);
		glBufferData(GL_ARRAY_BUFFER, vboSize, 0, GL_STATIC_DRAW);

		if( glGetError() ) {
			DBG_Message(true, "could not allocate buffer object of size %d\n", vboSize);
			errorCode = 3;
		}
	}

	if( !errorCode ) {

		unsigned int offset = 0;
		
		for_i( arrays.size() ) {
			ArrayInstruction & ai = arrays[i];
			int size = ai.array->size()*sizeof(float);
			glBufferSubData(GL_ARRAY_BUFFER, offset, size, &(*ai.array)[0]);

			glEnableVertexAttribArray(ai.location);
			glVertexAttribPointer(ai.location, ai.size, GL_FLOAT, GL_FALSE, 0, (void *)offset);

			offset += size;
		}

		if( indices.array != 0 ) {
			int size = indices.array->size()*sizeof(int);
			glBufferSubData(GL_ARRAY_BUFFER, offset, size, (void *)&(*indices.array)[0]);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *outVBO);
			*indices.bufferOffset = offset;

			offset += size;
		}

		glBindVertexArray(0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

 	va_end(args);

    return errorCode;
}

int GL_LoadGLSLProgramFromFile(const char * inVSFilename, const char * inFSFilename, GLuint * outProgram, ...)
//  e.g.,
//  if( GFX_LoadGLSLProgram("shader.vs", "shader.fs", &programID, 
//			eGLSLBindingAttribute, "inPosition", &position,
//			eGLSLBindingAttribute, "inNormal", &normal,
//			eGLSLBindingUniform, "kUniform1", &uniform1,
//			eGLSLBindingUniform, "kUniform2", &uniform2,
//			eGLSLBindingEnd) ) != 0 )
//		HandleError();
{
	int errorCode = 0;

	vector<char> vs;
	vector<char> fs;

	if( !errorCode ) {
		if( !FILE_Read(inVSFilename, vs) ) {
			DBG_Message(true, "failed to open %s\n", inVSFilename);
			errorCode = 1;
		}
	}

	if( !errorCode ) {
		if( !FILE_Read(inFSFilename, fs) ) {
			DBG_Message(true, "failed to open %s\n", inFSFilename);
			errorCode = 2;
		}
	}

	if( !errorCode ) {
		va_list args;
		va_start(args, outProgram);
		errorCode = GL_LoadGLSLProgramWithArgs(&vs[0], &fs[0], outProgram, args);
		va_end(args);
	}

	return errorCode;
}

int GL_LoadGLSLProgram(const char * inVS, const char * inFS, GLuint * outProgram, ...)
//  e.g.,
//  GFX_LoadGLSLProgram(vs, fs, &programID, 
//      eGLSLBindingAttribute, "inPosition", &position,
//      eGLSLBindingAttribute, "inNormal", &normal,
//      eGLSLBindingUniform, "kUniform1", &uniform1,
//      eGLSLBindingUniform, "kUniform2", &uniform2,
//      eGLSLBindingEnd);
{
	va_list args;
    va_start(args, outProgram);
	int errorCode = GL_LoadGLSLProgramWithArgs(inVS, inFS, outProgram, args);
	va_end(args);

	return errorCode;
}

PerformanceCounter::PerformanceCounter() {
	if( gTickFrequency == 0.0 ) {
		LARGE_INTEGER tickFrequency;
		QueryPerformanceFrequency(&tickFrequency);
		gTickFrequency = (double)tickFrequency.QuadPart;
	}
	tick.QuadPart = 0;
}

double PingPerformanceCounter(PerformanceCounter & inOutPC) {
	LARGE_INTEGER tick;
	QueryPerformanceCounter(&tick);

	double dt = (tick.QuadPart-inOutPC.tick.QuadPart)/PerformanceCounter::gTickFrequency;	
	inOutPC.tick = tick;

	return dt;
}

void TIMER_Begin(Timer & inTimer) {
	QueryPerformanceCounter(&inTimer.tick);
}

float TIMER_End(Timer & inTimer) {
	LARGE_INTEGER tick;
	QueryPerformanceCounter(&tick);

	static float frequency = 0;
	if( frequency == 0 ) {
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		frequency = (float)freq.QuadPart;
	}

	float dSeconds = (tick.QuadPart - inTimer.tick.QuadPart) / frequency;
	return dSeconds;
}


std::stack<Timer> gProfileTimers;
std::map<string, std::pair<double, int> > gProfileData; // function name to <runtime, #calls>

void PROF_Begin() {
	Timer profileTimer;
	TIMER_Begin(profileTimer);
	gProfileTimers.push(profileTimer);
}

void PROF_End(const char * inFuncName) {
	// use __FUNCTION__ for inFuncName

	double runtime = TIMER_End(gProfileTimers.top());
	gProfileTimers.pop();

	auto it = gProfileData.find(inFuncName);
	if( it == gProfileData.end() ) {
		gProfileData[inFuncName] = pair<double, int>(0, 0);
		it = gProfileData.find(inFuncName);
	}

	it->second.first += runtime;
	it->second.second++;
}

int PROF_FormatStatistics(char * outStatistics) {
	int length = 0;
	for_it( gProfileData ) {
		int len = sprintf_s(outStatistics, 512, "%s %fs\n", it->first.c_str(), it->second.first/it->second.second);
		outStatistics += len;
		length += len;
	}
	return length;
}

int OS_ReadImageFile(const char * inFilename, ImageInfo * outImage) {

	//GDI_InitLazyGdiplus();

	size_t origsize = strlen(inFilename) + 1;
    size_t convertedChars = 0;
    wchar_t wide[256];
    mbstowcs_s(&convertedChars, wide, origsize, inFilename, _TRUNCATE);
    
	Gdiplus::Bitmap bm(wide);
	if( bm.GetLastStatus() != Gdiplus::Ok ) {
		return 1;	
	}

	Gdiplus::BitmapData bitmapData;
	bm.LockBits(0, Gdiplus::ImageLockModeRead, bm.GetPixelFormat(), &bitmapData);
	outImage->width = bitmapData.Width;
	outImage->height = bitmapData.Height;
	outImage->bitDepth = Gdiplus::GetPixelFormatSize(bitmapData.PixelFormat);
	outImage->rowBytes = bitmapData.Stride;
	outImage->image = new unsigned char[outImage->rowBytes*outImage->height];
	memcpy(outImage->image, bitmapData.Scan0, outImage->rowBytes*outImage->height);
	bm.UnlockBits(&bitmapData);

	return 0;
}

int GL_LoadTextureImage(const char * inFilename, GLuint & outTexID) {
	
	ImageInfo image;
	if( OS_ReadImageFile(inFilename, &image) != 0 ) {
		return 1;
	}


	int B = image.bitDepth/8;
	int w = image.width;
	int h = image.height;
	unsigned char * v = image.image;
	int rowByteMap[] = { 8, 1, 2, 1, 4, 1, 2, 1 };
	int rowBytesI = image.rowBytes%8;

	glPixelStorei(GL_UNPACK_ALIGNMENT, rowByteMap[rowBytesI]);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, w);

	glGenTextures(1, &outTexID);

	if( 0 == outTexID ) {
		return 2;
	}

	glBindTexture(GL_TEXTURE_2D, outTexID);
	
	GLuint format = GL_RGB;
	if( B == 1 ) {
		format = GL_LUMINANCE;
	} else if( B == 3 ) {
		format = GL_BGR_EXT;
	} else if( B == 4 ) {
		format = GL_BGRA_EXT;
	}
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);		
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);		
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	//glTexImage2D(GL_TEXTURE_2D, 0, B, w, h, 0, format, GL_UNSIGNED_BYTE, v);
	gluBuild2DMipmaps(GL_TEXTURE_2D, B, w, h, format, GL_UNSIGNED_BYTE, v);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // reset to default
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);

	return 0;
}