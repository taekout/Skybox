#include "VWMath.h"
#include <math.h>


double xpwr2(double x) { return x*x; }

Boolean DoublesAreNearlyEqual(double n1, double n2, double epsilon)
{
	double abs_n1 = fabs(n1);
	double abs_n2 = fabs(n2);
	double abs_n1n2 = fabs(n1 - n2);

	if ((abs_n1 <= epsilon) || (abs_n2 <= epsilon))   // if either is very nearly zero, don't take the magnitude into account
		return (abs_n1n2 <= epsilon);  
	else

	return (abs_n1n2 <= ( abs_n1 > abs_n2 ? abs_n1 : abs_n2 ) * epsilon);
}


Boolean WorldCoordsAreNearlyEqual(double n1, double n2)
{
	return DoublesAreNearlyEqual(n1, n2, kNearlyEqualEpsilonForWorldCoords);
}

Boolean DoubleIsNearlyZero(double n)
{
	return DoublesAreNearlyEqual(n, 0);
}


Boolean WorldCoordIsNearlyZero(double n)
{
	return WorldCoordsAreNearlyEqual(n, 0);
}

Boolean NormalizedValuesAreNearlyEqual(double n1, double n2)
{ 
	// Factor out the division for speed and div0. Also, always use the 
	// larger value * epsilon so that DoublesAreNearlyEqual(x,y) == DoublesAreNearlyEqual(y,x)
	if ((n1 == 0) || (n2 == 0))
		return (fabs(n1 - n2) <= kNearlyEqualEpsilonForNormalizedValues);
	else
		return (fabs(n1 - n2) <= ( fabs(n1) > fabs(n2) ? fabs(n1) : fabs(n2) ) * kNearlyEqualEpsilonForNormalizedValues);
}


Boolean NormalizedValueIsNearlyOne(double n)
{
	return NormalizedValuesAreNearlyEqual(n, 1);
}



WorldPt3 WorldPt3::Normal() const				{ WorldPt3 r = *this; r.Normalize(); return r; }
//--------------------------------------------------------------------
WorldPt3 CrossProduct(const WorldPt3& p1, const WorldPt3& p2)
	{ return WorldPt3(p2.z*p1.y - p1.z*p2.y, p2.x*p1.z - p1.x*p2.z, p2.y*p1.x - p1.y*p2.x); }

//--------------------------------------------------------------------

bool TransformMatrix::IsOrthogonal() const
//
// This function checks to see if the upper left 3x3 submatrix is orthogonal. 
// An orthogonal matrix does rigid psBody transformations (no scaling or skewing)
// and has nice properties like M inverse == M transpose. VectorWorks objects 
// should all have orthogonal matrices. 
//
/////// John Williams, 2 October 1999 //////////////////////////////////////////
{
	WorldPt3	row1(v2.i);
	WorldPt3	row2(v2.j);
	WorldPt3	row3(v2.k);

	if (row1.IsNormalized() && 
		row2.IsNormalized() && 
		row3.IsNormalized() &&
		NormalizedValueIsNearlyOne(abs(DotProduct(CrossProduct(row1, row2),  row3)))
	)
		return true;
	else
		return false;
}





////////////////////////////////////////////////////////////////////////////////
void TransformMatrix::SetToIdentity()
// Sets the matrix to the identity matrix.
/////// John Williams, 4 October 1999 //////////////////////////////////////////
// Optimized out for loops.  -DLD 8/15/2001
{
	static bool gIdentityTransformDefined = false;
	static TransformMatrix kIdentityTransform;
	
	if (!gIdentityTransformDefined) {

		int j;

		for (int i = 0; i <= 3; i++) {
			for (j = 0; j <= 2; j++)
				kIdentityTransform.mat[i][j] = (j == i) ? 1.0 : 0.0;
		}
		
		gIdentityTransformDefined = true;
	}

	// On the Mac anyway, changing from = to memcpy is way faster.  -DLD 8/16/2001
	//*this = kIdentityTransform;

	memcpy(this, &kIdentityTransform, sizeof(TransformMatrix));
}


void TransformXMatrix::SetToIdentity()
// Sets the matrix to the identity matrix.
/////// John Williams, 5 October 1999 //////////////////////////////////////////
// Optimized out for loops.  -DLD 8/15/2001
{
	static bool gIdentityTransformXDefined = false;
	static TransformXMatrix kIdentityTransformX;
	
	if (!gIdentityTransformXDefined) {

		int j;

		for (int i = 0; i <= 3; i++) {
			for (j = 0; j <= 3; j++)
				kIdentityTransformX.mat[i][j] = (j == i) ? 1.0 : 0.0;
		}
		
		gIdentityTransformXDefined = true;
	}

	// On the Mac anyway, changing from = to memcpy is way faster.  -DLD 8/16/2001
	//*this = kIdentityTransformX;

	memcpy(this, &kIdentityTransformX, sizeof(TransformXMatrix));
}




////////////////////////////////////////////////////////////////////////////////
bool TransformMatrix::IsIdentity() const
// Checks to see if the matrix is an identity matrix.
/////// John Williams, 4 October 1999 //////////////////////////////////////////
{
	bool isIdentity = true;
	for (int i = 0; i <= 3 && isIdentity; i++)
		for (int j = 0; j <= 2 && isIdentity; j++)
			if (mat[i][j] != ((j == i) ? 1.0 : 0.0))
				isIdentity = false;
	return isIdentity;
}


/////////////////////////////
double DotProduct(const WorldPt3& p1, const WorldPt3& p2)
	{ return (p2.x * p1.x) + (p2.y * p1.y) + (p2.z * p1.z); }

double DotProductX(const WorldPt3& p)
//	{ return (p.x * 1) + (p.y * 0) + (p.z * 0); }
	{ return p.x; }

double DotProductY(const WorldPt3& p)
//	{ return (p.x * 0) + (p.y * 1) + (p.z * 0); }
	{ return p.y; }

double DotProductZ(const WorldPt3& p)
//	{ return (p.x * 0) + (p.y * 0) + (p.z * 1); }
	{ return p.z; }





void MatrixToXMatrix(const TransformMatrix &source, TransformXMatrix &dest)
{
	short						a;
	short						b;

	for (a = 0; a < 3; a++) {
		for (b = 0; b < 3; b++) {
			dest.mat[a][b] = source.mat[a][b];
		}
	}
	for (b = 0; b < 3; b++) {
		dest.mat[3][b] = source.mat[3][b];
	}
	for (b = 0; b < 3; b++) {
		dest.mat[b][3] = 0.0;
	}
	dest.mat[3][3] = 1.0;
}

void XMatrixToMatrix(const TransformXMatrix &source, TransformMatrix &dest)
{
	short						a;
	short						b;

	for (a = 0; a < 3; a++) {
		for (b = 0; b < 3; b++) {
			dest.mat[a][b] = source.mat[a][b];
		}
	}
	for (b = 0; b < 3; b++) {
		dest.mat[3][b] = source.mat[3][b];
	}
}

static void Invert4by4_ByeBye(void) {
}


void InvertXMatrix(const TransformXMatrix &source, TransformXMatrix &dest)
{
	TransformXMatrix			a;
	double					det;
	double					temp;
	double					pos;
	double					neg;
	double					mag;

	a = source;
	pos = 0.0;
	neg = 0.0;
	
	#if 0
	temp = a.mat[0][0]*a.mat[1][1]*a.mat[2][2];

	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;
	
	temp = a.mat[0][1]*a.mat[1][2]*a.mat[2][0];

	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;
	
	temp = a.mat[0][2]*a.mat[1][0]*a.mat[2][1];
	
	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;

	temp = -(a.mat[2][0]*a.mat[1][1]*a.mat[0][2]);
	
	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;

	temp = -(a.mat[2][1]*a.mat[1][2]*a.mat[0][0]);
	
	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;

	temp = -(a.mat[2][2]*a.mat[1][0]*a.mat[0][1]);
	
	if (temp >= 0) 
		pos = pos+temp;
	else 
		neg = neg+temp;
	
	#else
	
	temp = a.mat[0][0]*a.mat[1][1]*a.mat[2][2];		(temp >= 0 ? pos : neg) += temp;
	temp = a.mat[0][1]*a.mat[1][2]*a.mat[2][0];		(temp >= 0 ? pos : neg) += temp;
	temp = a.mat[0][2]*a.mat[1][0]*a.mat[2][1];		(temp >= 0 ? pos : neg) += temp;
	temp = -(a.mat[2][0]*a.mat[1][1]*a.mat[0][2]);	(temp >= 0 ? pos : neg) += temp;
	temp = -(a.mat[2][1]*a.mat[1][2]*a.mat[0][0]);	(temp >= 0 ? pos : neg) += temp;
	temp = -(a.mat[2][2]*a.mat[1][0]*a.mat[0][1]);	(temp >= 0 ? pos : neg) += temp;

	#endif
	
	det = pos+neg;
	mag = pos-neg;
	
	if (mag == 0.0) {
		Invert4by4_ByeBye();
		return;
	}
	else {
		temp = det/mag;
		if (DoubleIsNearlyZero(temp)) { 
			Invert4by4_ByeBye();
			return;
		}
	}
	
	dest.mat[0][0] = (a.mat[1][1]*a.mat[2][2]-a.mat[2][1]*a.mat[1][2])/det;
	dest.mat[0][1] = -((a.mat[0][1]*a.mat[2][2]-a.mat[2][1]*a.mat[0][2])/det);
	dest.mat[0][2] = (a.mat[0][1]*a.mat[1][2]-a.mat[1][1]*a.mat[0][2])/det;
	dest.mat[1][0] = -((a.mat[1][0]*a.mat[2][2]-a.mat[2][0]*a.mat[1][2])/det);
	dest.mat[1][1] = (a.mat[0][0]*a.mat[2][2]-a.mat[2][0]*a.mat[0][2])/det;
	dest.mat[1][2] = -((a.mat[0][0]*a.mat[1][2]-a.mat[1][0]*a.mat[0][2])/det);
	dest.mat[2][0] = (a.mat[1][0]*a.mat[2][1]-a.mat[2][0]*a.mat[1][1])/det;
	dest.mat[2][1] = -((a.mat[0][0]*a.mat[2][1]-a.mat[2][0]*a.mat[0][1])/det);
	dest.mat[2][2] = (a.mat[0][0]*a.mat[1][1]-a.mat[1][0]*a.mat[0][1])/det;
	dest.mat[3][0] = -(a.mat[3][0]*dest.mat[0][0]+a.mat[3][1]*dest.mat[1][0]+a.mat[3][2]*dest.mat[2][0]);
	dest.mat[3][1] = -(a.mat[3][0]*dest.mat[0][1]+a.mat[3][1]*dest.mat[1][1]+a.mat[3][2]*dest.mat[2][1]);
	dest.mat[3][2] = -(a.mat[3][0]*dest.mat[0][2]+a.mat[3][1]*dest.mat[1][2]+a.mat[3][2]*dest.mat[2][2]);
	dest.mat[0][3] = 0.0;
	dest.mat[1][3] = 0.0;
	dest.mat[2][3] = 0.0;
	dest.mat[3][3] = 1.0;
	
}


void InvertMatrix(const TransformMatrix &source, TransformMatrix &dest)
{
	TransformXMatrix					tempMat;

	MatrixToXMatrix(source, tempMat);
	InvertXMatrix(tempMat, tempMat);
	XMatrixToMatrix(tempMat, dest);
}


void BasePointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b)
//
// The base transformation code, so we can have
// one version that checks the matrix is a linear
// transformation, and one that doesn't.
//
{
	// presumably this copy is done so a and b can be the same object; we could simply pass
	// a by value though.  In fact, because of aliasing (the possibility that mat and b are
	// overlapping), we really should probably compute into d and then assign to b at the
	// end.  I don't want to have to deal with testing this at the moment, though. -- ACB
	WorldPt3 d = a;

	b.x =	d.x * mat.v1.a00 + d.y * mat.v1.a10 + d.z * mat.v1.a20 + mat.v1.xOff;
	b.y =	d.x * mat.v1.a01 + d.y * mat.v1.a11 + d.z * mat.v1.a21 + mat.v1.yOff;
	b.z =	d.x * mat.v1.a02 + d.y * mat.v1.a12 + d.z * mat.v1.a22 + mat.v1.zOff;
}


void PointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b)
//
// This function performs what is assumed to be a
// linear transformation on a and puts the result
// in b.
//
{
	BasePointTransformN(a, mat, b);
}

void NonLinearPointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b)
//
// This function performs what is assumed to be a
// linear transformation on a and puts the result
// in b.
//
{
	BasePointTransformN(a, mat, b);
}

void InversePointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b)
{
	WorldPt3				dOffset;

	dOffset.x = a.x - mat.v1.xOff;
	dOffset.y = a.y - mat.v1.yOff;
	dOffset.z = a.z - mat.v1.zOff;
	b.x =	dOffset.x * mat.v1.a00 + dOffset.y * mat.v1.a01 + dOffset.z * mat.v1.a02;
	b.y =	dOffset.x * mat.v1.a10 + dOffset.y * mat.v1.a11 + dOffset.z * mat.v1.a12;
	b.z =	dOffset.x * mat.v1.a20 + dOffset.y * mat.v1.a21 + dOffset.z * mat.v1.a22;
}

void VectorTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b) 
{
	WorldPt3 d = a;  // need this for VectorTransformN(v, t, v);
	
	b.x =	d.x * mat.v1.a00 + d.y * mat.v1.a10 + d.z * mat.v1.a20;
	b.y =	d.x * mat.v1.a01 + d.y * mat.v1.a11 + d.z * mat.v1.a21;
	b.z =	d.x * mat.v1.a02 + d.y * mat.v1.a12 + d.z * mat.v1.a22;
}

void InverseVectorTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b)
{
	WorldPt3 d = a;  // need this for InverseVectorTransformN(v, t, v);
	
	b.x =	d.x * mat.v1.a00  +  d.y * mat.v1.a01  +  d.z * mat.v1.a02;
	b.y =	d.x * mat.v1.a10  +  d.y * mat.v1.a11  +  d.z * mat.v1.a12;
	b.z =	d.x * mat.v1.a20  +  d.y * mat.v1.a21  +  d.z * mat.v1.a22;
}



WorldPt3 PointTransformN(const WorldPt3 &a, const TransformMatrix &mat)
{
	WorldPt3 result;
	PointTransformN(a, mat, result);
	return result;
}

WorldPt3 NonLinearPointTransformN(const WorldPt3 &a, const TransformMatrix &mat)
{
	WorldPt3 result;
	PointTransformN(a, mat, result);
	return result;
}

WorldPt3 InversePointTransformN(const WorldPt3 &a, const TransformMatrix &mat)
{
	WorldPt3 result;
	InversePointTransformN(a, mat, result);
	return result;
}

WorldPt3 VectorTransformN(const WorldPt3 &a, const TransformMatrix &mat)
{
	WorldPt3 result;
	VectorTransformN(a, mat, result);
	return result;
}

WorldPt3 InverseVectorTransformN(const WorldPt3 &a, const TransformMatrix &mat)
{
	WorldPt3 result;
	InverseVectorTransformN(a, mat, result);
	return result;
}

TransformMatrix MatrixMultiply(const TransformMatrix &mat1, const TransformMatrix &mat2)
{
	TransformMatrix				matTemp;

	matTemp.v1.a00 =	mat2.v1.a00 * mat1.v1.a00 +
						mat2.v1.a10 * mat1.v1.a01 +
						mat2.v1.a20 * mat1.v1.a02;
	matTemp.v1.a01 =	mat2.v1.a01 * mat1.v1.a00 +
						mat2.v1.a11 * mat1.v1.a01 +
						mat2.v1.a21 * mat1.v1.a02;
	matTemp.v1.a02 =	mat2.v1.a02 * mat1.v1.a00 +
						mat2.v1.a12 * mat1.v1.a01 +
						mat2.v1.a22 * mat1.v1.a02;
	matTemp.v1.a10 =	mat2.v1.a00 * mat1.v1.a10 +
						mat2.v1.a10 * mat1.v1.a11 +
						mat2.v1.a20 * mat1.v1.a12;
	matTemp.v1.a11 =	mat2.v1.a01 * mat1.v1.a10 +
						mat2.v1.a11 * mat1.v1.a11 +
						mat2.v1.a21 * mat1.v1.a12;
	matTemp.v1.a12 =	mat2.v1.a02 * mat1.v1.a10 +
						mat2.v1.a12 * mat1.v1.a11 +
						mat2.v1.a22 * mat1.v1.a12;
	matTemp.v1.a20 =	mat2.v1.a00 * mat1.v1.a20 +
						mat2.v1.a10 * mat1.v1.a21 +
						mat2.v1.a20 * mat1.v1.a22;
	matTemp.v1.a21 =	mat2.v1.a01 * mat1.v1.a20 +
						mat2.v1.a11 * mat1.v1.a21 +
						mat2.v1.a21 * mat1.v1.a22;
	matTemp.v1.a22 =	mat2.v1.a02 * mat1.v1.a20 +
						mat2.v1.a12 * mat1.v1.a21 +
						mat2.v1.a22 * mat1.v1.a22;
	matTemp.v1.xOff =	mat1.v1.xOff * mat2.v1.a00 +
						mat1.v1.yOff * mat2.v1.a10 +
						mat1.v1.zOff * mat2.v1.a20 + mat2.v1.xOff;
	matTemp.v1.yOff =	mat1.v1.xOff * mat2.v1.a01 +
						mat1.v1.yOff * mat2.v1.a11 +
						mat1.v1.zOff * mat2.v1.a21 + mat2.v1.yOff;
	matTemp.v1.zOff =	mat1.v1.xOff * mat2.v1.a02 +
						mat1.v1.yOff * mat2.v1.a12 +
						mat1.v1.zOff * mat2.v1.a22 + mat2.v1.zOff;
	return matTemp;
}

WorldPt3 operator*(const double scalar, const WorldPt3 thePt)
	{ return WorldPt3(thePt.x * scalar, thePt.y * scalar, thePt.z * scalar); }

WorldPt3 operator*(const WorldPt3 inPt, const double inScalar)
	{ return WorldPt3(inPt.x * inScalar, inPt.y * inScalar, inPt.z * inScalar); }

void TranslateMatrix(TransformMatrix &f, double x, double y, double z)
{
	f.v1.xOff += x;
	f.v1.yOff += y;
	f.v1.zOff += z;
}

void SetAxisRotationXMatrix(MajorAxisSpec axis, double degrees, TransformXMatrix &mat)

{
	double					myCos;
	double					mySin;
	double					radians;

	mat.SetToIdentity();
	radians = degrees * kRadiansPerDegree;

	////VW901 fix, bug#19245///////////////////////////////////////
	mySin = WorldCoordsAreNearlyEqual(abs(sin(radians)), 0)? 0.0 : sin(radians);//sin(radians);
	myCos = WorldCoordsAreNearlyEqual(abs(cos(radians)), 0)? 0.0 : cos(radians);//cos(radians);
	/////////////////////////////////////////////////////////////

	switch (axis) {
		case kXAxis: {
			mat.mat[0][0] = 1.0;
			mat.mat[1][1] = myCos;
			mat.mat[1][2] = mySin;
			mat.mat[2][1] = -mySin;
			mat.mat[2][2] = myCos;
			break;
		}
		case kYAxis: {
			mat.mat[2][0] = mySin;
			mat.mat[2][2] = myCos;
			mat.mat[0][0] = myCos;
			mat.mat[0][2] = -mySin;
			mat.mat[1][1] = 1.0;
			break;
		}
		case kZAxis: {
			mat.mat[0][0] = myCos;
			mat.mat[0][1] = mySin;
			mat.mat[1][0] = -mySin;
			mat.mat[1][1] = myCos;
			mat.mat[2][2] = 1.0;
			break;
		}
	}
}


void SetAxisRotationMatrix(MajorAxisSpec axis, double degrees, TransformMatrix &mat)
{
	TransformXMatrix					tmpMat;

	SetAxisRotationXMatrix(axis, degrees, tmpMat);
	XMatrixToMatrix(tmpMat, mat);
}



void AngleBetween(const WorldPt3 v1, const WorldPt3 v2, double & rdAngleRadians) {
	double vec1[] = { v1.x, v1.y, v1.z };
	double vec2[] = { v2.x, v2.y, v2.z };
	double VEC3_AngleBetween(double * first_vec3, double * second_vec3);
	rdAngleRadians = VEC3_AngleBetween(vec1, vec2);
}



double CrossProductMagnitude(const WorldPt& a, const WorldPt& b)  { return a.x * b.y - a.y * b.x; }


double DotProduct(const WorldPt& a, const WorldPt& b)   		     { return a.x * b.x + a.y * b.y; }

Boolean ApproxEqualVectors(const WorldPt3 &a, const WorldPt3 &b)
{
	return NormalizedValueIsNearlyOne(DotProduct(a, b));
}

bool ApproxEqualMatrix3by3(const TransformMatrix &a, const TransformMatrix &b)
// Checks to see if the upper left 3x3 submatrices are the same within a
// floating point rounding tolerance.
{
	return (ApproxEqualVectors(a.v2.i, b.v2.i) && 
			ApproxEqualVectors(a.v2.j, b.v2.j) && 
			ApproxEqualVectors(a.v2.k, b.v2.k)	);
}

/* AxisAngleToMatrix( ) - Generates a 4x4 rotation matrix representing rotation about 'axis' 
		'axisDeg' degrees.
	axis - axis of rotation.  Does not need to be normalized.  Should have magnitude > 0.
	angleDeg - angle of rotation about 'axis' in degrees.  Any degree is valid.
	return - a 4x4 rotation matrix representing rotation about 'axis' 'axisDeg' degrees.
*/		
TransformMatrix AxisAngleToMatrix(WorldPt3 axis, double angleDeg) {
	TransformMatrix result;
	result.SetToIdentity();
	
	if( axis.IsZero() ) {
		return result;
	}

	axis.Normalize();

	double angleRad = angleDeg * PI/180.0;
	double c = cos(angleRad);
	double s = sin(angleRad);
	double t = 1-c;

	// variable rename
	double & x = axis.x;
	double & y = axis.y;
	double & z = axis.z;

	// equation from http://www.gamedev.net/reference/articles/article1199.asp
	result.v1.a00 = t*x*x + c;
	result.v1.a10 = t*x*y + s*z;
	result.v1.a20 = t*x*z - s*y;
	result.v1.a01 = t*x*y - s*z;
	result.v1.a11 = t*y*y + c;
	result.v1.a21 = t*y*z + s*x;
	result.v1.a02 = t*x*z + s*y;
	result.v1.a12 = t*y*z - s*x;
	result.v1.a22 = t*z*z + c;

	return result;
}

bool RayPlaneIntersection(const WorldPt3 & ptOnPlane, const WorldPt3 & planeNormal, const WorldPt3 & rayOrigin, const WorldPt3 & rayDir, double * t, WorldPt3 * intersection) {
	WorldPt3 p = ptOnPlane;
	WorldPt3 n = planeNormal;
	WorldPt3 o = rayOrigin;
	WorldPt3 d = rayDir.Normal();

	
	double den = n.Dot(d);
	if( den == 0 ) return false;

	double num = -n.Dot(o - p);

	double tt = num / den;
	if( t != 0 ) {
		*t = tt;
	}
	if( intersection != 0 ) {
		*intersection = o + d * tt;
	}
	return true;
}




int NearestPowerOfTwo(int positiveNum, int minNum, int maxNum) {
	// returns the closest power of two number to 'positiveNum' in [minNum, maxNum].  
	// e.g., NearestPowerOfTwo(20, 1, 512) == 16, NearestPowerOfTwo(30, 1, 512) == 32

	int powerOfTwo = 0;

	double nearest = log((double)positiveNum) / log(2.0);
	if( nearest-(int)nearest >= 0.5 ) {				// if nearest power of two is larger than positiveNum
		powerOfTwo = (int)pow(2.0, ((int)nearest) + 1);
	} else {										// if nearest power of two is smaller than positiveNum
		powerOfTwo = (int)pow(2.0, (int)nearest);
	}

	return min(max(powerOfTwo, minNum), maxNum);
}

TransformMatrix MakeScaleMatrix(double scale) {
	/* Simply makes and returns a scale matrix.
		scale - scale of the matrix; this value runs down the diagonal of the rotation part of the matrix
		return - a scale matrix
	*/

	TransformMatrix scaleMatrix;
	scaleMatrix.SetToIdentity();

	scaleMatrix.v1.a00 = scale;
	scaleMatrix.v1.a11 = scale;
	scaleMatrix.v1.a22 = scale;

	return scaleMatrix;
}

bool WorldPt3NearlyEqual(const WorldPt3 & p1, const WorldPt3 & p2) {
	return WorldCoordsAreNearlyEqual(p1.x, p2.x)
		&& WorldCoordsAreNearlyEqual(p1.y, p2.y)
		&& WorldCoordsAreNearlyEqual(p1.z, p2.z);
}
bool IsPointNearLine(const WorldPt & inPoint, const WorldPt & inP0, const WorldPt & inP1, double inNearDist, int inLineType)
// returns whether point is near a line.  line may be an infinite line, a ray or a segment.
//		inLineType - 1 == infinite line, 2 == ray from inP0 toward inP1, 3 == segment from inP0 to inP1
{
	WorldPt vec1 = inP1 - inP0;
	WorldPt vec1Norm = vec1.Normal();
	WorldPt vec2 = inPoint - inP0;
	double dot = vec2.Dot(vec1Norm);
	
	if( inLineType == 2 ) {												// if line is a ray
		if( dot < -inNearDist ) {										// if pt on wrong side of ray
			return false;
		}
	} else if( inLineType == 3 ) {										// if line is a segment
		if( dot < -inNearDist || dot > vec1.Magnitude()+inNearDist ) {	// if pt outside of segment
			return false;
		}
	}

	double perpDot = vec2.CrossMagnitude(vec1Norm);
	return abs(perpDot) <= inNearDist;
}
