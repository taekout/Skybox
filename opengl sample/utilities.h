#ifndef utilities_h
#define utilities_h

#define GS_WIN 1
#define WIDGET_DEMO_APP 1


#include <Windows.h>
#include "VWMath.h"
#include "gl/glew.h"
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>

using namespace std;

#define for_i(size) for( size_t i=0; i<size; i++ )
#define for_j(size) for( size_t j=0; j<size; j++ )
#define for_k(size) for( size_t k=0; k<size; k++ )

#define for_it( container ) for( auto it = container.begin(); it != container.end(); it++ ) 
#define for_jt( container ) for( auto jt = container.begin(); jt != container.end(); jt++ ) 
#define for_kt( container ) for( auto kt = container.begin(); kt != container.end(); kt++ ) 

extern int gWidgetFrameworkErrorCode; // 0 means no error
bool GenerateErrorCode(int inErrorCode, bool inSerious = true);
float MTH_Round(float value);
double MTH_Round(double value);
double RND_DoubleInRange(double min, double max);
int RND_IntInRange(int min, int max);
void RND_GenerateNormalizedVector(double * vec3);
WorldPt3 RND_GenerateNormalizedVector();
TransformMatrix RND_MakeTransform(double minx, double miny, double minz, double maxx, double maxy, double maxz);
void VEC3_Equate(double * vec3, double * output_vec3);
void VEC3_Negate(double * vec3);
double VEC3_Dot(double * first_vec3, double * second_vec3);
double VEC3_Magnitude(double * vec3);
void VEC3_Normalize(double * vec3, double * normalized_vec3);
void VEC3_Normalize(double * vec3);
void VEC3_Cross(double * a_vec3, double * b_vec3, double * crossed_vec3);
void VEC3_Difference(double * a_vec3, double * b_vec3, double * difference_vec3);
void VEC3_Swap(double * a_vec3, double * b_vec3);
void VEC3_Add(double * first_vec3, double * second_vec3, double * output_vec3);
void VEC3_Subtract(double * first_vec3, double * second_vec3, double * output_vec3);
void VEC3_Multiply(double * first_vec3, double * second_vec3, double * output_vec3);
void VEC3_Multiply(double * vec3, double scalar, double * output_vec3);
double VEC3_AngleBetween(double * first_vec3, double * second_vec3);
void MAT4_ExtractCol3(double * mat4, int col, double * vec3);
void MAT4_ExtractRow3(double * mat4, int row, double * vec3);
void MAT4_SetCol3(int col, double * vec3, double * mat4);
void MAT4_LoadIdentity(double * mat4);
void MAT4_MatrixVectorMultiply(double * mat4, double * column_vec3, double * row_vec3);
void MAT4_VEC4_MatrixVectorMultiply(double * mat4, double * column_vec4, double * out_vec4);
void MAT4_VEC4_InverseMatrixVectorMultiply(double * mat4, double * column_vec4, double * out_vec4);
void MAT4_VectorMatrixMultiply(double * row_vec3, double * mat4, double * column_vec3);
void MAT4_Multiply(double * first_mat4, double * second_mat4, double * output_mat4);
void MAT4_TransposeUpper3(double * mat4, double * transpose_mat4);
void MAT4_SetEqual(double * mat4, double * output_mat4);
void MAT4_InvertRotation(double * rotation_mat4, double * invert_mat4);
double MTH_DegToRad(double deg);
double MTH_RadToDeg(double rad);
double MTH_Sign(double inValue);
TransformMatrix MTH_MakeArbitraryRotation(WorldPt3 & axis_vec3, double angle_rad);

void MAT4_MakeZRotation(double deg, double * zRotation_mat4);
void MAT4_MakeXRotation(double deg, double * xRotation_mat4);
void MAT4_MakeArbitraryRotation(double * axis_vec3, double angle_rad, double * rotation_mat4);
void MAT4_MakeArbitraryRotation(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * rotation_mat4);
void MAT4_PostRotate(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * inOutMat);
void MAT4_PreRotate(float inAxisX, float inAxisY, float inAxisZ, float inAngleRad, float * inOutMat);

void MAT4_MakeTranslation(float inX, float inY, float inZ, float * outMat);
void MAT4_PostTranslate(float inX, float inY, float inZ, float * inOutMat);
void MAT4_PreTranslate(float inX, float inY, float inZ, float * inOutMat);

void MAT4_MakeScale(float inScaleX, float inScaleY, float inScaleZ, float * outMat);
void MAT4_PreScale(float inScaleX, float inScaleY, float inScaleZ, float * inOutMat);
void MAT4_PostScale(float inScaleX, float inScaleY, float inScaleZ, float * inOutMat);

void MAT4_MakePerspective(double fovy_rad, double aspect, double zNear, double zFar, double * mat4);
void MAT4_MakeInversePerspective(double fovy_rad, double aspect, double zNear, double zFar, double * mat4);
void MAT4_MakeOrtho(double left, double right, double bottom, double top, double near, double far, double * mat4);
void MAT4_MakeFrustum(double left, double right, double bottom, double top, double zNear, double zFar, double * mat4);
void MAT4_LookAt(float eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ, double * viewMatrix);
void MAT4_Translate(double inDX, double inDY, double inDZ, double * inOutMat4);

void FPS_CAM_MoveForward(double amount, bool fixedZ, double * camera_mat4);
void FPS_CAM_StrafeRight(double amount, double * camera_mat4);
void FPS_CAM_LookRight(double amount_deg, double * camera_mat4);
void FPS_CAM_LookUp(double amount_deg, double * camera_mat4);
void RND_GenerateNormalizedVector(double * vec3);

 class TQuaternion {
 public:
	double x, y, z, w;
	
	TQuaternion() {
		x = 0;
		y = 0;
		z = 0;
		w = 1;
	}

	TQuaternion(double x, double y, double z, double w) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}


	TQuaternion(TransformMatrix & matrix) {
		double trace = matrix.v1.a00 + matrix.v1.a11 + matrix.v1.a22 + 1;
		double s = 0;

		// if trace is positive then perform instant calculation
		if( trace > 1e-8 ) {
			s = 0.5 / sqrt(trace);
			w = 0.25 / s;
			x = ( matrix.v1.a21 - matrix.v1.a12 ) * s;
			y = ( matrix.v1.a02 - matrix.v1.a20 ) * s;
			z = ( matrix.v1.a10 - matrix.v1.a01 ) * s;
		} else {
			int majorColumn = 0;
			if( matrix.v1.a11 > matrix.mat[majorColumn][majorColumn] ) majorColumn = 1;
			if( matrix.v1.a22 > matrix.mat[majorColumn][majorColumn] ) majorColumn = 2;

			switch( majorColumn ) {
				case 0:
					s = sqrt( 1.0 + matrix.v1.a00 - matrix.v1.a11 - matrix.v1.a22 ) * 2.0;
					x = 0.25 * s;
					y = (matrix.v1.a01 + matrix.v1.a10 ) / s;
					z = (matrix.v1.a02 + matrix.v1.a20 ) / s;
					w = (matrix.v1.a12 + matrix.v1.a21 ) / s;
					break;

				case 1:
					s = sqrt( 1.0 + matrix.v1.a11 - matrix.v1.a00 - matrix.v1.a22 ) * 2.0;
					x = (matrix.v1.a01 + matrix.v1.a10 ) / s;
					y = 0.25 * s;
					z = (matrix.v1.a12 + matrix.v1.a21 ) / s;
					w = (matrix.v1.a02 + matrix.v1.a20 ) / s;
					break;

				case 2:
					s = sqrt( 1.0 + matrix.v1.a22 - matrix.v1.a00 - matrix.v1.a11 ) * 2.0;
					x = (matrix.v1.a02 + matrix.v1.a20 ) / s;
					y = (matrix.v1.a12 + matrix.v1.a21 ) / s;
					z = 0.25 * s;
					w = (matrix.v1.a01 + matrix.v1.a10 ) / s;
					break;
			}
		}
	}

	TransformMatrix ToRotationMatrix() {
		TransformMatrix m;

		double xx = x * x;
		double xy = x * y;
		double xz = x * z;
		double xw = x * w;

		double yy = y * y;
		double yz = y * z;
		double yw = y * w;

		double zz = z * z;
		double zw = z * w;

		m.v1.a00 = 1 - 2 * ( yy + zz );
		m.v1.a01 = 2 * ( xy - zw );
		m.v1.a02 = 2 * ( xz + yw );
		
		m.v1.a10 = 2 * ( xy + zw );
		m.v1.a11 = 1 - 2 * ( xx + zz );
		m.v1.a12 = 2 * ( yz - xw );
		
		m.v1.a20 = 2 * ( xz - yw );
		m.v1.a21 = 2 * ( yz + xw );
		m.v1.a22 = 1 - 2 * ( xx + yy );
		
		m.P() = WorldPt3(0, 0, 0);

		return m;
	}

	void AxisAngle(WorldPt3 & axis, double & angle /*rad*/) {
		// compute the axis and angle representation of 'this' rotation

		angle = 2.0*acos(w);
		double sinHalfAngle = sin(angle/2.0);
		axis = WorldPt3(x / sinHalfAngle, y / sinHalfAngle, z / sinHalfAngle);
	}
	
	TQuaternion Slerp(TQuaternion end, double t) {
		// Spherical Linear intERPolation between 'this' and 'end' at time 't' which 
		// ranges from [0, 1]
		// end - The end orientation
		// t - tnterpolation time in [0, 1]
		double to1[4];
        double omega, cosom, sinom, scale0, scale1;

        // calc cosine
        cosom = x*end.x + y*end.y + z*end.z + w*end.w;

        // adjust signs (if necessary)
        if ( cosom < 0.0 ) {
			cosom = -cosom;
			to1[0] = - end.x;
			to1[1] = - end.y;
			to1[2] = - end.z;
			to1[3] = - end.w;
        } else {
			to1[0] = end.x;
			to1[1] = end.y;
			to1[2] = end.z;
			to1[3] = end.w;
        }


        // calculate coefficients

       if ( (1.0 - cosom) > 1e-6 ) {
			// standard case (slerp)
			omega = acos(cosom);
			sinom = sin(omega);
			scale0 = sin((1.0 - t) * omega) / sinom;
			scale1 = sin(t * omega) / sinom;

	   } else {
			// "from" and "to" quaternions are very close 
			//  ... so we can do a linear interpolation
               scale0 = 1.0 - t;
               scale1 = t;
       }

		// calculate final values
		return TQuaternion(
			scale0 * x + scale1 * to1[0],
			scale0 * y + scale1 * to1[1],
			scale0 * z + scale1 * to1[2],
			scale0 * w + scale1 * to1[3]);
	}
};

void MAT4_VWViewMatToGLViewMat(TransformMatrix & inMat, double inCrsPersDistance, float * outMat);
void MAT4_VWToGLMatrix(TransformMatrix & inMat, float * outMat);
void MAT4_VWToGLMatrix(TransformMatrix & mat, double * glMatrix);
void MAT4_GLToVWMatrix(double * glMatrix, TransformMatrix & mat);

void GL_Frustum(double left, double right, double bottom, double top, double zNear, double zFar, double * projectionMatrix);
void MAT4_MakeOrtho(double left, double right, double bottom, double top, double zNear, double zFar, double * mat4);
void GLU_LookAt(double eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ, double * viewMatrix);

void MAT4_VEC4_Multiply(double * mat4, double * column_vec4, double * row_vec4);

bool RaySphereIntersection(WorldPt3 & sphereCenter, double radius, WorldPt3 & rayOrigin, WorldPt3 & rayDir, /*out*/ double * t, /*out*/ WorldPt3 * intersection);
bool RayCylinderIntersection(double radius, double * minZ, double * maxZ, WorldPt3 & rayOrigin, WorldPt3 & rayDir, double * t, WorldPt3 * intersection, bool findNegativeIntersections);
bool GEO_LineLineIntersection(WorldPt & inP1, WorldPt & inV1, WorldPt & inP2, WorldPt & inV2, bool inSegmentIntersection, bool & outParallel, double & outT1, double & outT2, WorldPt & outIntersection);
void GEO_GenerateSphere(int inPrimitiveType, float * inOrigin, float inRadius, int inSamples, vector<float> & outVertices, vector<float> & outNormals, vector<float> & outTexCoords, vector<int> & outIndices, int * outNumTriangles);

vector<double> MTH_QuadraticEquation(double a, double b, double c);
void MAT4_RotationFromSphericalCoords(double theta, double phi, TransformMatrix & cameraMat);
bool MAT4_Invert(const double * mat4, double * out4);
bool MAT4_Transpose(double * mat4, double * out4);
bool GL_DrawString(const char * text, unsigned int textLen, int x, int y, int columnWidth, double * textWidth, double * textHeight);
void GL_DrawDiamond(const WorldPt & center, double width, double height);
void GL_FillDiamond(const WorldPt & center, double width, double height);
void GL_DrawCircle(double radius, const WorldPt3 & normal, int segments);
void GL_DrawCircle(const WorldPt & center, double radius, int segments);
void GL_Circle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments);
void GL_DrawRect(const WorldPt3 & inOrigin, const WorldPt3 & inNormal, const WorldRect & inRect);
void GL_FillDisc(const WorldPt & center, double r1, double r2, int segments);
void GL_DrawCircle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments);
void GL_FillCircle(const WorldPt3 & center, double radius, const WorldPt3 & normal, int segments);
void GL_FillCircle(const WorldPt & center, double radius, int segments);
void GL_DrawDepthCuedLine(const WorldPt3 & origin, const WorldPt3 & dir, const WorldPt3 & localViewer, double maxAlpha, double minAlpha, double maxDist, double maxWidth, int segments);
void GL_DrawWireCube(WorldCube & cube);
WorldCube TransformAxesAlignedCube(WorldCube & worldCube, TransformMatrix & transform);
WorldCube InverseTransformAxesAlignedCube(WorldCube & worldCube, TransformMatrix & transform);
void CheckGLError();



float VEC3_DistanceBetween(float * inV1, float * inV2);
float VEC3_DistanceBetween_Sq(float * inV1, float * inV2);

void VEC3_Set(float * outVec, float x, float y, float z);
void VEC4_Set(float * outVec, float x, float y, float z, float w);

void VEC3_Set(float x, float y, float z, float * outVec);
void VEC4_Set(float x, float y, float z, float w, float * outVec);

void VEC4_Equate(float * inVec, float * outVec);
void VEC3_Equate(float * vec3, float * output_vec3);
void VEC3_Negate(float * vec3);
float VEC3_Dot(float * first_vec3, float * second_vec3);
float VEC3_Magnitude_Sq(float * vec3);
float VEC3_Magnitude(float * vec3);
void VEC3_Normalize(float * vec3, float * normalized_vec3);
void VEC3_Normalize(float * vec3);
void VEC3_Cross(float * a_vec3, float * b_vec3, float * crossed_vec3);
void VEC3_Difference(float * a_vec3, float * b_vec3, float * difference_vec3);
void VEC3_Swap(float * a_vec3, float * b_vec3);
void VEC3_Add(float * first_vec3, float * second_vec3, float * output_vec3);
void VEC3_Subtract(float * first_vec3, float * second_vec3, float * output_vec3);
void VEC3_Multiply(float * first_vec3, float * second_vec3, float * output_vec3);
void VEC3_Multiply(float * vec3, float scalar, float * output_vec3);
void VEC4_Multiply(float * inVec4, float scalar, float * outVec);
void VEC3_Divide(float * inVec3, float scalar, float * outVec);
void VEC4_Divide(float * inVec4, float scalar, float * outVec);
float VEC3_AngleBetween(float * first_vec3, float * second_vec3);
void VEC3_MAT4_Multiply(float * row_vec3, float * mat4, float * column_vec3);
void MAT4_Set(float in0, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float in8, float in9, float in10, float in11, float in12, float in13, float in14, float in15, float * outMat);
void MAT4_ExtractCol3(float * mat4, int col, float * vec3);
void MAT4_ExtractRow3(float * mat4, int row, float * vec3);
void MAT4_Equate(float * inMat, float * outMat);
void MAT4_SetCol3(int col, float * vec3, float * mat4);
void MAT4_LoadIdentity(float * mat4);
void MAT3_VEC3_Multiply(float * mat4, float * column_vec3, float * row_vec3);
void MAT4_VEC3_Multiply(float * mat4, float * column_vec3, float * row_vec3);
void MAT4_VEC4_MatrixVectorMultiply(float * mat4, float * column_vec4, float * out_vec4);
void MAT4_Multiply(float * first_mat4, float * second_mat4, float * output_mat4);
void MAT4_TransposeUpper3(float * mat4, float * transpose_mat4);
void MAT4_Transpose(float * mat4);
void MAT4_SetEqual(float * mat4, float * output_mat4);
void MAT4_InvertQuick(float * inMat, float * outMat);
float MTH_DegToRad(float deg);
float MTH_RadToDeg(float rad);
void MAT4_MakeZRotation(float deg, float * zRotation_mat4);
void MAT4_MakeXRotation(float deg, float * xRotation_mat4);
void MAT4_MakeYRotation(float deg, float * xRotation_mat4);
void MAT4_MakeArbitraryRotation(float * axis_vec3, float angle_rad, float * rotation_mat4);
void MAT4_MakeInversePerspective(float fovy_rad, float aspect, float zNear, float zFar, float * mat4);
void MAT4_MakePerspective(float fovy_rad, float aspect, float zNear, float zFar, float * mat4);
void MAT4_MakeOrtho(float left, float right, float bottom, float top, float zNear, float zFar, float * mat4);
void MAT4_MakeFrustum(float left, float right, float bottom, float top, float zNear, float zFar, float * mat4);
void MAT4_LookAt(float eyeX, float eyeY, float eyeZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ, float * viewMatrix);
void MAT4_MakeArbitraryRotation(float * axis_vec3, float angle_rad, float * rotation_mat4);
void MAT4_Translate(float inX, float inY, float inZ, float * inOutMat4);
void MAT4_SetUpper3x3Rowwise(float * outMat, float * inR0, float * inR1, float * inR2);
void MAT4_SetUpper3x3(float * outMat, float * inR0, float * inR1, float * inR2);
void MAT4_AxisAngleToMatrix(float * inAxis, float inAngleDeg, float * outMat);
void GLU_LookAt(float eyeX, float eyeY, float eyeZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ, float * viewMatrix);
void MAT4_VEC4_Multiply(float * mat4, float * column_vec4, float * row_vec4);
bool MAT4_Invert(const float * mat4, float * outMat);
void MAT4_Transpose(float * inMat4, float * outMat4);
void MAT4_ExtractUpper3x3(float * inMat4, float * outMat3);
void MAT4_ExtractFrustumPoints(float * inVPInvertMat, float outFrustumPoints[3*8]);

float MTH_GetFraction(float inF);

void QT_Equate(float * inQ, float * outQ);
void QT_Set(float inX, float inY, float inZ, float inW, float * outQ);
void QT_MAT4_MakeQuaternion(float * inMat, float * outQ);
void QT_MAT4_MakeMatrix(float * inQ, float * outMat);
void QT_GetAxisAngle(float * inQ, float * outAxis, float outAngle /*rad*/);
void QT_SLERP(float * inQ0, float * inQ1, float inT, float * outQ);


int ExtractFrustumCorners(float * inProjectionMat, float * inModelviewMat, float * outCorners);
void FPSCAM_MoveForward(float amount, bool fixedZ, float * inOutViewMat);
void FPSCAM_MoveUp(float amount, bool fixedZ, float * inOutViewMat);
void FPSCAM_StrafeRight(float amount, float * inOutViewMat);
void FPSCAM_LookRight(float amount_deg, bool inFixedUp, float * inOutViewMat);
void FPSCAM_LookUp(float amount_deg, float * inOutViewMat);
void GenerateSphere(float inRadius, vector<float> & vertices, vector<float> & normals);
void ExtractFrustumPoints(float * inVPInvertMat, float outFrustumPoints[3*8]);
bool RayPlaneIntersection(float * ptOnPlane, float * planeNormal, float * rayOrigin, float * rayDir, float * optOutT, float * outIntersection);



template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1) {
	inOutList.push_back(in1);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5, const T2 & in6) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
	inOutList.push_back(in6);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5, const T2 & in6, const T2 & in7) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
	inOutList.push_back(in6);
	inOutList.push_back(in7);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5, const T2 & in6, const T2 & in7, const T2 & in8) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
	inOutList.push_back(in6);
	inOutList.push_back(in7);
	inOutList.push_back(in8);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5, const T2 & in6, const T2 & in7, const T2 & in8, const T2 & in9) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
	inOutList.push_back(in6);
	inOutList.push_back(in7);
	inOutList.push_back(in8);
	inOutList.push_back(in9);
}

template<typename T1, typename T2>
void push_back(T1 & inOutList, const T2 & in1, const T2 & in2, const T2 & in3, const T2 & in4, const T2 & in5, const T2 & in6, const T2 & in7, const T2 & in8, const T2 & in9, const T2 & in10) {
	inOutList.push_back(in1);
	inOutList.push_back(in2);
	inOutList.push_back(in3);
	inOutList.push_back(in4);
	inOutList.push_back(in5);
	inOutList.push_back(in6);
	inOutList.push_back(in7);
	inOutList.push_back(in8);
	inOutList.push_back(in9);
	inOutList.push_back(in10);
}

template<typename T>
class SimpleTable {
	vector<T *> data;
	list<int> freeIDs;

public:

	SimpleTable() {
		data.resize(1);
	}

	~SimpleTable() {
		clear();
	}

	
	class iterator {
	public:

		struct IDValue {
			IDValue(int inID, T & inValue) : first(inID), second(inValue) {}
			void operator =(IDValue &) {}

			int first;
			T & second;
		};
		
		IDValue operator *() {
			return &idValue;
		}

		IDValue * operator ->() {
			idValue.first = index;
			idValue.second = table->data[index];
			return (IDValue *)&idValue; // cast magic to set reference
		}

		void operator ++(int) {
			do {
				index++;
			} while( index<(int)table->data.size() && table->data[index]==0 );
		}

		void operator --() {
			do {
				index--;
			} while( index >= 0 && table->data[index]==0 );			
		}
		
		bool operator !=(const iterator & it) {
			return index != it.index;
		}

	
	private:
		friend class SimpleTable;
		iterator(SimpleTable<T> * inTable, int inIndex) {
			table = inTable;
			index = inIndex;
		}

		int index;
		SimpleTable<T> * table;
		
		struct {
			int first;
			T * second;
		} idValue;
	
	};

	iterator begin() {
		if( empty() ) {
			return end();
		}

		int index = 0;
		while( !data[index] ) {
			index++;
		}
		iterator it(this, index);

		return it;
	}

	iterator end() {
		int index = (int)data.size();
		iterator it(this, index);

		return it;
	}

	bool insert(int inID, const T & inValue = T()) {

		if( inID <= 0 ) {
			return false;
		}

		if( (int)data.size() < inID+1 ) {
			for( int i=data.size(); i<inID; i++ ) {
				freeIDs.push_back(i);
			}
			data.resize(inID+1);

			data[inID] = new T();
			
		} else if( data[inID] == 0 ) {
			list<int>::iterator it = ::find(freeIDs.begin(), freeIDs.end(), inID);
			if( it != freeIDs.end() ) {
				freeIDs.erase(it);
			}

			data[inID] = new T();
		}

		*data[inID] = inValue;

		return true;
	}

	int insert(const T & inValue = T()) {

		int newID = 0;
		if( !freeIDs.empty() ) {
			newID = (int)freeIDs.front();
			freeIDs.pop_front();
			data[newID] = new T();

		} else {
			newID = (int)data.size();
			data.push_back(new T());
		}

		*data[newID] = inValue;

		return newID;
	}

	int erase(int inID) {
		int errorCode = 0;

		delete data[inID];
		data[inID] = 0;
		freeIDs.push_back(inID);

		return errorCode;
	}

	int clear() {
		int errorCode = 0;

		for( size_t i=0; i<data.size(); i++ ) {
			if( data[i] ) {
				delete data[i];
			}
		}

		data.resize(1);
		freeIDs.clear();

		return errorCode;
	}

	T & operator [](int inID) {
		if( !DoesExist(inID) ) {
			insert(inID, T());
		}
		return *data[inID];
	}

	int size() {
		return data.size() - freeIDs.size() - 1;
	}
	bool empty() {
		return size() == 0;
	}

	iterator find(int inID) {
		if( DoesExist(inID) ) {
			iterator it(this, inID);
			return it;
		}

		return end();
	}

	bool DoesExist(int inID) {
		bool found = false;

		if( inID>0 && inID<(int)data.size() && data[inID]!=0 ) {
			found = true;
		}

		return found;
	}
};

template<typename T>
class SimpleTable2 { // IDs dealt out are never redealt (unless DeleteAllEntries is called)
	vector<T *> data;
	size_t size;

public:

	SimpleTable2() {
		size = 0;
		data.push_back(0); // unused id 0
	}

	~SimpleTable2() {
		DeleteAllEntries();
	}

	int NewEntry(int * outID) {
		int errorCode = 0;

		*outID = (int)data.size();
		data.push_back(new T());
		size++;

		return errorCode;
	}

	int DeleteEntry(int & inID) {
		int errorCode = 0;

		delete data[inID];
		data[inID] = 0;
		size--;

		return errorCode;
	}

	int DeleteAllEntries() {
		int errorCode = 0;

		for( size_t i=0; i<data.size(); i++ ) {
			if( data[i] ) {
				delete data[i];
			}
		}

		data.clear();
		data.push_back(0);
		size = 0;

		return errorCode;
	}

	int GetEntry(int & inID, T ** outValue) {
		int errorCode = 0;

		*outValue = data[inID];

		return errorCode;
	}

	int SetEntry(int & inID, const T & inValue) {
		int errorCode = 0;

		*data[inID] = inValue;

		return errorCode;
	}
};

int DBG_glGetError();
void DBG_Message(bool inAlert, char * inMessage, ...);

void GL_SaveTextureAsRAW(GLuint inTex, GLenum inType, size_t inFilesize, const char * inFilename, ...);
void GL_SaveBitsAsRAW(float * inBits, size_t inFilesize, const char * inFilename, ...);
void GL_SaveBufferAsRAW(GLenum inFormat, const char * inFilename, ...);

enum EGLSLBinding {  // for use with GL_LoadGLSLProgram{FromFile}
    eGLSLBindingUniform,
    eGLSLBindingAttribute,
    eGLSLBindingEnd
};

int GL_LoadGLSLProgram(const char * inVS, const char * inFS, GLuint * outProgram, ...);
int GL_LoadGLSLProgramFromFile(const char * inVSFilename, const char * inFSFilename, GLuint * outProgram, ...);

enum EVA {  // for use with GL_GenerateVA()
	eVAAttrib,
	eVAIndices,
	eVAEnd
};

int GL_GenerateVA(GLuint * outVBO, GLuint * outVA, ...);

struct PerformanceCounter {
	PerformanceCounter();

	LARGE_INTEGER tick;
	static double gTickFrequency;
};

double PingPerformanceCounter(PerformanceCounter & inOutPC);

struct Timer {
	LARGE_INTEGER tick;
};

void TIMER_Begin(Timer & inTimer);
float TIMER_End(Timer & inTimer);


void PROF_Begin();
void PROF_End(const char * inFuncName);
int PROF_FormatStatistics(char * outStatistics);

struct ImageInfo {
	ImageInfo() {
		image = 0;
		width = 0;
		height = 0;
		bitDepth = 0;
		rowBytes = 0;
	}
	ImageInfo(int inWidth, int inHeight, int inBitDepth, int inRowBytes, unsigned char * inImage) {
		image = inImage;
		width = inWidth;
		height = inHeight;
		bitDepth = inBitDepth;
		rowBytes = inRowBytes;
	}
	int width, height, bitDepth, rowBytes;
	unsigned char * image;
};

int GL_LoadTextureImage(const char * inFilename, GLuint & outTexID);


#endif
