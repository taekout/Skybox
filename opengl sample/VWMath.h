#ifndef VWMath_h
#define VWMath_h

#include <algorithm>
#include <math.h>
#include <float.h>
#include <string.h>

using namespace std;

typedef int Boolean;

const double kMaxLegalWorldCoord =  1.0E100;
const double kMaxOutOfRangeWorldCoord =  DBL_MAX;
const double kNearlyEqualEpsilonForDoubles          = 1e-13;
const double kNearlyEqualEpsilonForNormalizedValues = 1e-06;
const double kNearlyEqualEpsilonForWorldCoords = kNearlyEqualEpsilonForDoubles;
const double PI = 3.141592654;
const double kDegreesPerRadian	= 180.0 / PI;
const double kRadiansPerDegree	= PI / 180.0;
enum MajorAxisSpec {
	kNoAxis,
	kXAxis,
	kYAxis,
	kZAxis
};

struct WorldPt3;
class TransformMatrix;

WorldPt3 CrossProduct(const WorldPt3& p1, const WorldPt3& p2);
double DotProduct(const WorldPt3& p1, const WorldPt3& p2);

double xpwr2(double x);
Boolean DoublesAreNearlyEqual(double n1, double n2, double epsilon = kNearlyEqualEpsilonForDoubles);
Boolean WorldCoordsAreNearlyEqual(double n1, double n2);
Boolean DoubleIsNearlyZero(double n);
Boolean WorldCoordIsNearlyZero(double n);
Boolean NormalizedValuesAreNearlyEqual(double n1, double n2);
Boolean NormalizedValueIsNearlyOne(double n);
inline Boolean WorldCoord1_GE_WorldCoord2(double n1, double n2)
{
	return ((n1 > n2) || WorldCoordsAreNearlyEqual(n1, n2));
}

inline Boolean WorldCoord1_GT_WorldCoord2(double n1, double n2)
{
	return ((n1 > n2) && (! WorldCoordsAreNearlyEqual(n1, n2)));
}

inline Boolean WorldCoord1_LE_WorldCoord2(double n1, double n2)
{
	return ((n1 < n2) || WorldCoordsAreNearlyEqual(n1, n2));
}

inline Boolean WorldCoord1_LT_WorldCoord2(double n1, double n2)
{
	return ((n1 < n2) && (! WorldCoordsAreNearlyEqual(n1, n2)));
}

class WorldPt;

double CrossProductMagnitude(const WorldPt& a, const WorldPt& b);
double DotProduct(const WorldPt& a, const WorldPt& b);

#pragma pack(push, 2) // BEGIN 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)


class _WorldPt 
{
public:
	double y, x;

	double X() const { return x; }
	double Y() const { return y; }

//	operator const WorldPt&() const			{ return *((const WorldPt *) this); }
//	operator WorldPt&()						{ return *((WorldPt *) this); }

	// As Mike H. points out, these are dangerous if other classes try to derive from _WorldPt (which of course they do).
	// We should move these to WorldPt. [MAF 10/5/01]
	WorldPt*       operator&(void)			{ return (WorldPt *) this; }
	const WorldPt* operator&(void) const 	{ return (const WorldPt *) this; }
	
	Boolean operator==(const _WorldPt& p) const { return WorldCoordsAreNearlyEqual(y, p.y) 
		                                              && WorldCoordsAreNearlyEqual(x, p.x); }
	Boolean operator!=(const _WorldPt& p) const { return !(*this == p);  }

	void ByteSwap();
};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------


//--------------------------------------------------------------------
#pragma pack(push, 2) // BEGIN 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)

class WorldPt : public _WorldPt {
public:

	WorldPt(void)				{ this->x = kMaxLegalWorldCoord; this->y = kMaxLegalWorldCoord;}
	WorldPt(double x, double y)	{ this->x = x; this->y = y; }

//	WorldPt(const WorldPt& p)	{ this->y = p.y; this->x = p.x; }
	WorldPt(const WorldPt& p)	{ memcpy(this, &p, sizeof(WorldPt)); }   // faster for misaligned doubles [MAF FP Speedup, 2/21/01]
	
	
	WorldPt(const _WorldPt& p)	{ this->y = p.y; this->x = p.x; }
	
	void Set(double x, double y) { this->x = x; this->y = y; }
	void SetAngle(double radians) { x = cos(radians); y = sin(radians); }

	Boolean operator==(const WorldPt& p) const { return WorldCoordsAreNearlyEqual(y, p.y) 
		                                             && WorldCoordsAreNearlyEqual(x, p.x); }
	Boolean operator!=(const WorldPt& p) const { return !(*this == p); }

	WorldPt operator+(const WorldPt& p) const	{ return WorldPt(x + p.x, y + p.y); }
	WorldPt operator-(const WorldPt& p) const	{ return WorldPt(x - p.x, y - p.y); }
	WorldPt	operator*(double scalar) const	{ return WorldPt((x * scalar), (y * scalar)); }
	WorldPt	operator/(double scalar) const	{ return WorldPt((x / scalar), (y / scalar)); }
	
	WorldPt operator-(void) const				{ return WorldPt(-x, -y); }
	
	WorldPt& operator=(const WorldPt& p) { memcpy(this, &p, sizeof(WorldPt)); return *this; }
	
	WorldPt& operator+=(const WorldPt& p) { x += p.x;    y += p.y;    return *this; }
	WorldPt& operator-=(const WorldPt& p) { x -= p.x;    y -= p.y;    return *this; }
	WorldPt& operator*=(double scalar) { x *= scalar; y *= scalar; return *this; }
	WorldPt& operator/=(double scalar) { x /= scalar; y /= scalar; return *this; }

	Boolean IsZero(void) const	{ return WorldCoordIsNearlyZero(x) && WorldCoordIsNearlyZero(y); } 

	double Magnitude(void)        const { return sqrt(x*x + y*y); }
	double MagnitudeSquared(void) const { return x*x + y*y; }
	void       Normalize(void) { if (!this->IsZero()) *this /= this->Magnitude(); }

	// This function assumes that the WorldPt is a vector. 
	// This poorly named function returns a unit vector, or a normalized vector, NOT a perpindicular vector.
	WorldPt    Normal(void) const { WorldPt r = *this; r.Normalize(); return r; }

	// returns a perpendicular 2d vector turning clockwise or counter-clockwise depending on the parameter
	WorldPt Perpendicular(bool clockwise = false) const {
		return clockwise ? WorldPt(y, -x) : WorldPt(-y, x);
	}

	bool Parallel(const WorldPt & p) const {
		const double Epsilon = 1e-9;

		WorldPt p1 = Normal();
		WorldPt p2 = p.Normal();

		return fabs(p1.x-p2.x)<Epsilon && fabs(p1.y-p2.y)<Epsilon;
	}
	// crossing 2d vectors is useful in that it revels a winding order.
	// This function returns the z value of the cross product.
	double CrossMagnitude(const WorldPt & p) {
		return CrossProductMagnitude(*this, p);
	}
	// Perform dot product
	double Dot(const WorldPt & p) {
		return DotProduct(*this, p);
	}
	double DistanceFrom(WorldPt otherPoint)        const { return (*this - otherPoint).Magnitude(); }
	double SquaredDistanceFrom(WorldPt otherPoint) const { return (*this - otherPoint).MagnitudeSquared(); }

	double AngleInDegrees(void) const { if (this->IsZero()) return 0.0; else return kDegreesPerRadian * atan2(y, x); }
	double DegreeAngleFrom(const WorldPt inOrigin) const	 { WorldPt p = *this - inOrigin; return p.AngleInDegrees(); }

	WorldPt Abs(void) const		{ return WorldPt(abs(x), abs(y)); }
	WorldPt Half(void) const	{ return WorldPt(x / 2, y / 2); }

};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------



struct WorldCoordBase { 
	double x, y, z;
	double&		operator[](int i)       { return ((double *) this)[i]; }
	double		operator[](int i) const { return ((const double *) this)[i]; }

	Boolean			operator==(const WorldCoordBase& p) const { return WorldCoordsAreNearlyEqual(x, p.x) 
																	&& WorldCoordsAreNearlyEqual(y, p.y)
																	&& WorldCoordsAreNearlyEqual(z, p.z);  }
	Boolean			operator!=(const WorldCoordBase& p) const { return !(*this == p); }
};


#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>
struct _WorldPt3 : public WorldCoordBase {
//	operator WorldPt3&()					{ return *((WorldPt3 *) this); }
//	operator const WorldPt3&() const		{ return *((const WorldPt3 *) this); }
	typedef double WorldCoordType;
	
	double X() const { return x; }	// XXX_JDW_MISC const ?
	double Y() const { return y; }
	double Z() const { return z; }

	WorldPt3* operator&(void)				{ return (WorldPt3 *) this; }
	const WorldPt3* operator&(void) const	{ return (const WorldPt3 *) this; }
};
#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------

//--------------------------------------------------------------------
#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>
struct WorldPt3 : public _WorldPt3 {

	WorldPt3(void)						{ this->x = 0; this->y = 0; this->z = 0;}
	
	WorldPt3(const WorldPt3& p)	  { memcpy(this, &p, sizeof(WorldPt3)); }   // faster for misaligned doubles [MAF FP Speedup, 2/21/01]
	WorldPt3(const _WorldPt3& p)		{ this->y = p.y; this->x = p.x; this->z = p.z; }

	WorldPt3(WorldPt p, double z)		{ this->x = p.x; this->y = p.y; this->z = z; }
	WorldPt3(double x, double y, double z)	{ this->x = x; this->y = y; this->z = z; }

	void Set(double x, double y, double z) { this->x = x; this->y = y; this->z = z; }
	
	Boolean operator==(const WorldPt3& p) const { return WorldCoordsAreNearlyEqual(x, p.x) 
	                                                  && WorldCoordsAreNearlyEqual(y, p.y)
	                                                  && WorldCoordsAreNearlyEqual(z, p.z); }
	Boolean operator!=(const WorldPt3& p) const { return !(*this == p); }

	WorldPt3	operator+(const WorldPt3& p) const { return WorldPt3(x + p.x, y + p.y, z + p.z); }
	WorldPt3	operator-(const WorldPt3& p) const { return WorldPt3(x - p.x, y - p.y, z - p.z); }
	WorldPt3	operator-(void) const { return WorldPt3(-x, -y, -z); }

	WorldPt3&	operator=(const WorldPt3& p)  { memcpy(this, &p, sizeof(WorldPt3)); return *this; }
	WorldPt3&	operator+=(const WorldPt3& p) { x += p.x; y += p.y; z += p.z; return *this; }
	WorldPt3&	operator-=(const WorldPt3& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }

	WorldPt3& operator*=(double s)			{ x *= s; y *= s; z *= s; return *this; }
	WorldPt3& operator/=(double s)			{ x /= s; y /= s; z /= s; return *this; }
	
	//WorldPt3 operator*(double scalar) const;
	WorldPt3 operator/(double scalar) const;

	Boolean IsZero(void) const	{ return WorldCoordIsNearlyZero(x) && WorldCoordIsNearlyZero(y) && WorldCoordIsNearlyZero(z); } 
	
	WorldPt3 abs(void) const;
	WorldPt3 Half(void) const { return WorldPt3(x/2, y/2, z/2); }		// XXX_JDW_MISC - maybe remove

	bool Parallel(const WorldPt3 & p) const {
		const double Epsilon = 1e-9;

		WorldPt3 p1 = Normal();
		WorldPt3 p2 = p.Normal();

		return fabs(p1.x-p2.x)<Epsilon && fabs(p1.y-p2.y)<Epsilon && fabs(p1.z-p2.z)<Epsilon;
	}

	// compute the distance from 'this' WorldPt3 to parameter 'p'
	double DistanceFrom(const WorldPt3 & p) const {
        return (*this-p).Magnitude();
	}

	// It's faster to compute squared distances
	double SquaredDistanceFrom(const WorldPt3 & p) const {
        return (*this-p).MagnitudeSquared();
	}

	double Magnitude(void)        const { return (sqrt(xpwr2(x) + xpwr2(y) + xpwr2(z))); }
	double MagnitudeSquared(void) const { return xpwr2(x) + xpwr2(y) + xpwr2(z); }

	// Perform cross product
	WorldPt3 Cross(const WorldPt3 & p) {
		return CrossProduct(*this, p);
	}
	// Perform dot product
	double Dot(const WorldPt3 & p) {
		return DotProduct(*this, p);
	}
	void            Normalize(void) { if (!this->IsZero()) *this /= this->Magnitude(); }
	WorldPt3 Normal(void) const;
	Boolean         IsNormalized(void) const { return NormalizedValueIsNearlyOne(this->MagnitudeSquared()); }


	
};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------

WorldPt3 CrossProduct(const WorldPt3& p1, const WorldPt3& p2);

//--------------------------------------------------------------------
#pragma pack(push, 2) // BEGIN 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)

class _TransformMatrix 	
{
public:	
	union {
		
		double	mat[4][3];
		struct {
			double a00;    double a01;    double a02;
			double a10;    double a11;    double a12;
			double a20;    double a21;    double a22;
			double xOff;   double yOff;   double zOff;
		} v1;
		struct {
			_WorldPt3		i, j, k;
			_WorldPt3	offset;
		} v2;
	};

	TransformMatrix* operator&(void)				{ return (TransformMatrix *) this; }
	const TransformMatrix* operator&(void) const	{ return (const TransformMatrix *) this; }
};

#pragma pack(pop) // END 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)
//--------------------------------------------------------------------


//--------------------------------------------------------------------
#pragma pack(push, 2) // BEGIN 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)

class TransformMatrix : public _TransformMatrix
// This is a standard 4x4 matrix with the right column vector (0 0 0 1) removed
// (perspective transformation & w).  It should always be orthogonal (all basis 
// vectors are orthonormal).
{
public:
	TransformMatrix()  { }

	TransformMatrix& operator=(const TransformMatrix& tm)	{ memcpy(this, &tm, sizeof(TransformMatrix)); return *this; }
	TransformMatrix(const _TransformMatrix& t)	{ this->v2.i = t.v2.i; this->v2.j = t.v2.j; this->v2.k = t.v2.k; this->v2.offset = t.v2.offset; }

	// XXX_JDW_MISC - maybe R0() R1() R2() P() C0() C1() C2() ?  XDir() YDir() ZDir() Pt() IDir() JDir() KDir() ? Review.

	//const WorldPt3& R0() const	{ return reinterpret_cast<const WorldPt3&>(v2.i); }	// mat[0][*]
	//const WorldPt3& R1() const	{ return reinterpret_cast<const WorldPt3&>(v2.j); }	// mat[1][*]
	//const WorldPt3& R2() const	{ return reinterpret_cast<const WorldPt3&>(v2.k); }	// mat[2][*]
	//const WorldPt3& P() const		{ return reinterpret_cast<const WorldPt3&>(v2.offset); }// mat[3][*]  rename? XXX_JDW_MISC
	
	WorldPt3 C0() const	{ return WorldPt3(v1.a00, v1.a10, v1.a20); }	// mat[*][0]
	WorldPt3 C1() const	{ return WorldPt3(v1.a01, v1.a11, v1.a21); }	// mat[*][1]
	WorldPt3 C2() const	{ return WorldPt3(v1.a02, v1.a12, v1.a22); }	// mat[*][2]

	// I'm opposed to non-const ref returning functions like this -- it's basically the same as having
	// a public member variable -- but we have direct member access now anyway
	// so I'll live with it. -- ACB 1/28/00
	WorldPt3& R0()				{ return reinterpret_cast<WorldPt3&>(v2.i); }		// mat[0][*]
	WorldPt3& R1()				{ return reinterpret_cast<WorldPt3&>(v2.j); }		// mat[1][*]
	WorldPt3& R2()				{ return reinterpret_cast<WorldPt3&>(v2.k); }		// mat[2][*]
	WorldPt3& P()					{ return reinterpret_cast<WorldPt3&>(v2.offset); }		// mat[3][*]  rename? XXX_JDW_MISC

	// We might want to assert that the vectors are normalized, but the plan is to make
	// WorldPt3 do some of that and not just be a typedef.
	void SetRow0(const WorldPt3& vec)
		{ mat[0][0] = vec.x; mat[0][1] = vec.y; mat[0][2] = vec.z; }
	void SetRow1(const WorldPt3& vec)
		{ mat[1][0] = vec.x; mat[1][1] = vec.y; mat[1][2] = vec.z; }
	void SetRow2(const WorldPt3& vec)
		{ mat[2][0] = vec.x; mat[2][1] = vec.y; mat[2][2] = vec.z; }
	
	bool IsOrthogonal() const;
	bool 		Is2DOrHybrid() const	{ return (abs(v2.k.z) == 1.0); };
	bool IsIdentity() const;
	
	void SetToIdentity();
	
	void ByteSwap();
};

#pragma pack(pop) // END 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)
//--------------------------------------------------------------------

struct TransformXMatrix
// 4x4 matrix of doubles. Unless you need the fourth column for something like perspective
// calculations, you should use the 3x4 TransformMatrix instead. 
{	
	double	mat[4][4];
	void SetToIdentity();
};



/////////////////////////////
double DotProduct(const WorldPt3& p1, const WorldPt3& p2);
double DotProductX(const WorldPt3& p);
double DotProductY(const WorldPt3& p);
double DotProductZ(const WorldPt3& p);
void MatrixToXMatrix(const TransformMatrix &source, TransformXMatrix &dest);
void XMatrixToMatrix(const TransformXMatrix &source, TransformMatrix &dest);
static void Invert4by4_ByeBye(void);
void InvertXMatrix(const TransformXMatrix &source, TransformXMatrix &dest);
void InvertMatrix(const TransformMatrix &source, TransformMatrix &dest);
void BasePointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
void PointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
void NonLinearPointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
void InversePointTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
void VectorTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
void InverseVectorTransformN(const WorldPt3 &a, const TransformMatrix &mat, WorldPt3 &b);
WorldPt3 PointTransformN(const WorldPt3 &a, const TransformMatrix &mat);
WorldPt3 NonLinearPointTransformN(const WorldPt3 &a, const TransformMatrix &mat);
WorldPt3 InversePointTransformN(const WorldPt3 &a, const TransformMatrix &mat);
WorldPt3 VectorTransformN(const WorldPt3 &a, const TransformMatrix &mat);
WorldPt3 InverseVectorTransformN(const WorldPt3 &a, const TransformMatrix &mat);
TransformMatrix MatrixMultiply(const TransformMatrix &mat1, const TransformMatrix &mat2);
WorldPt3 operator*(const double scalar, const WorldPt3 thePt);
WorldPt3 operator*(const WorldPt3 inPt, const double inScalar);
void TranslateMatrix(TransformMatrix &f, double x, double y, double z);
void SetAxisRotationMatrix(MajorAxisSpec axis, double degrees, TransformMatrix &mat);

void MAT4_VWToGLMatrix(TransformMatrix & mat, double * glMatrix);
void MAT4_GLToVWMatrix(double * glMatrix, TransformMatrix & mat);

void AngleBetween(const WorldPt3 v1, const WorldPt3 v2, double & rdAngleRadians);


Boolean ApproxEqualVectors(const WorldPt3 &a, const WorldPt3 &b);

bool ApproxEqualMatrix3by3(const TransformMatrix &a, const TransformMatrix &b);
TransformMatrix AxisAngleToMatrix(WorldPt3 axis, double angleDeg);
bool RayPlaneIntersection(const WorldPt3 & ptOnPlane, const WorldPt3 & planeNormal, const WorldPt3 & rayOrigin, const WorldPt3 & rayDir, double * t, WorldPt3 * intersection);

int NearestPowerOfTwo(int positiveNum, int minNum, int maxNum);

bool MTH_SnapToMultiple(double & value, double multiple, double snapTolerance);

#if TEST
	#define WORLDRECTISEMPTY(rect, assert)	(rect.IsEmpty(assert))
	#define WORLDCUBEISEMPTY(cube, assert)	(cube.IsEmpty(assert))
#else
	#define WORLDRECTISEMPTY(rect, assert)	((rect.left > rect.right) || (rect.bottom > rect.top))
	#define WORLDCUBEISEMPTY(cube, assert)	((cube.left > cube.right) || (cube.bottom > cube.top) || (cube.back > cube.front))
#endif

class WorldRect;
class WorldCube;

typedef enum { eAssert, eDontAssert } EAssertStatus; // controlled assert behavior when manipulating WorldRects - TDW
typedef enum { eStandardValue, eSpecialValue } EEmptyKind;

//--------------------------------------------------------------------
#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>

class _WorldRect {
public:
//protected:
	double top, left, bottom, right;   // this should be protected, not public [MAF 10/5/99]
	
//	operator const WorldRect&() const		{ return *((const WorldRect *) this); }
//	operator WorldRect& ()					{ return *((WorldRect *) this); }

	WorldRect* operator&(void)				{ return (WorldRect *) this; }
	const WorldRect* operator&(void) const	{ return (const WorldRect *) this; }

	friend class WorldRect;	
	friend class WorldCube;
};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------

//--------------------------------------------------------------------
#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>

class WorldRect : public _WorldRect
//
// See comments above
//
{
public:

	WorldRect(void)                               { MakeEmpty(); }
	WorldRect(double l, double t, 
		      double r, double b)         { this->SetLTRB(l, t, r, b); }
	WorldRect(const WorldPt& pt)                  { this->SetLTRB(pt.x, pt.y, pt.x, pt.y); }

	// these constructors must not aggressively validate their data because
	// they are used in very general purpose (and implicit) situations.
	WorldRect(const WorldPt& a, const WorldPt& b) { this->SetLTRB(min(a.x, b.x), max(a.y, b.y), max(a.x, b.x), min(a.y, b.y)); }

//	WorldRect(const WorldRect& r)	   { this->SetLTRB(r.left, r.top, r.right, r.bottom); }
	WorldRect(const WorldRect& r)	   { memcpy(this, &r, sizeof(WorldRect)); }   // faster for misaligned doubles [MAF FP Speedup, 2/21/01]
		
	WorldRect(const _WorldRect& r)						{ this->SetLTRB(r.left, r.top, r.right, r.bottom); }
	WorldRect(const WorldPt& cen, double radius)	{ this->Set(cen, radius); }

	void Set(double l, double t,             // set all rect dimensions
		     double r, double b)          { SetLTRB(l, t, r, b); }
	void Set(const WorldPt& a, const WorldPt& b);    // set based on point locations
	void SetCornersTLBR(WorldPt tl, WorldPt br);     // set the topLeft and bottomRight of the rect without adjusting to make rect correct
	void Set(const WorldPt& cen, double radius); // set the center and outset distance
	void SetLTRB(double l, double t, 
		         double r, double b);        // set all rect dimensions
	void AdjustLTRB(double l, double t, 
		         double r, double b);        // adjust all rect dimensions
	void SetLTRBWithoutValidating(double l, double t, 
				  				  double r, double b); 
	
	// Gets for the individual fields. There are no individual Sets because setting left and right
	// individually could trip the IsEmpty checks overaggressively. The multi-field Sets above should
	// handle all cases (see me if they don't). [MAF 11/17/99]

	const double Top()    const { return WORLDRECTISEMPTY((*this), eAssert) ? 0 : top;   }   
	const double Left()   const { return WORLDRECTISEMPTY((*this), eAssert) ? 0 : left;  }   
	const double Bottom() const { return WORLDRECTISEMPTY((*this), eAssert) ? 0 : bottom;}   
	const double Right()  const { return WORLDRECTISEMPTY((*this), eAssert) ? 0 : right; }   

#if GS_WIN
	friend class  WorldRect;  // huh? [MAF VW10 Revisit, 11/23/99]
#endif
	friend class  WorldCube;

	const WorldPt TopLeft(void)  const { return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(left,  top); }
	const WorldPt BotRight(void) const { return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(right, bottom); }
	const WorldPt BotLeft(void)  const { return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(left,  bottom); }
	const WorldPt TopRight(void) const { return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(right, top); }

	WorldRect& operator=(const WorldRect& r) { memcpy(this, &r, sizeof(WorldRect)); return *this; }

	Boolean   operator==(const WorldRect& r) const;
	Boolean   operator!=(const WorldRect& r) const { return !(*this == r); }
	
	void      Inset(double dx,  double dy); 
	void      Outset(double dx, double dy); 

	void      Offset(const WorldPt& p); 

	void      Unite(const WorldPt&   p);	
	void	  UniteNoCheck(const WorldPt& p);
	void      Unite(const WorldRect& r);

	WorldRect Union(const WorldPt&   p) const;	
	WorldRect Union(const WorldRect& r) const;	

	void      Intersect(const WorldRect& r);        	
	WorldRect Intersection(const WorldRect& r) const;	
	
	WorldPt		Center(void) const	{ return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt((left+right)/2,	(top+bottom)/2); }
	WorldPt		CenterLeft() const	{ return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(left,				(top+bottom)/2); }
	WorldPt		CenterTop() const	{ return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt((left+right)/2,	top); }
	WorldPt		CenterRight() const	{ return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt(right,			(top+bottom)/2); }
	WorldPt		CenterBottom() const{ return WORLDRECTISEMPTY((*this), eAssert) ? WorldPt(0, 0) : WorldPt((left+right)/2,	bottom); }
	double	Height(void) const	{ return WORLDRECTISEMPTY((*this), eAssert) ? 0 : top   - bottom; }
	double	Width(void) const	{ return WORLDRECTISEMPTY((*this), eAssert) ? 0 : right - left; }

	double Perimeter()  const { return ((this->Width() * 2) + (this->Height() * 2));}


	Boolean IsEmpty(EAssertStatus as = eDontAssert) const; 
	Boolean IsAProperlyDefinedEmptyRect() const;
	void    MakeEmpty(EEmptyKind kind = eStandardValue);
	Boolean IsAPoint() const { return ((bottom == top) && (left == right)); }
	Boolean IsNotEmpty(EAssertStatus /*as = eDontAssert*/) const { return !WORLDRECTISEMPTY((*this), eDontAssert); }
	Boolean IsNotAPoint() const { return !this->IsAPoint(); }


	Boolean Intersects(const WorldRect& r) const;    // True even if the two rects just share a point or a line.
	Boolean Overlaps(const WorldRect& r) const;      // True if the two rects share more than a point or a line.
	Boolean ContainsRect(const WorldRect& r) const;
	Boolean ContainsPoint(const WorldPt& p) const;

	
	WorldRect& operator*=(double scale) 
	{
		if (!WORLDRECTISEMPTY((*this), eAssert)) {
			left   *= scale; 
			top    *= scale; 
			right  *= scale; 
			bottom *= scale;
		}
			
		return *this;
	}

	WorldRect& operator/=(double scale) 
	{
		if (!WORLDRECTISEMPTY((*this), eAssert))
			*this *= 1 / scale;
			
		return *this;
	}

	WorldRect operator*(double s) const 
	{
		WorldRect ret = *this; ret *= s; return ret;
	}

	WorldRect operator/(double s) const 
	{
		WorldRect ret = *this; ret /= s; return ret;
	}

	void ByteSwap();
};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------

inline void WorldRect::MakeEmpty(EEmptyKind kind)
// by default, this will use the standard empty value. If a special
// empty value is requested, a value will be assigned in the TEST and BUG
// version which will guarantee that the WorldRect will fail a test
// by IsAProperlyDefinedEmptyRect();

{
	kind;
	#if TEST

	if (kind == eSpecialValue) {
		bottom = left = kMaxOutOfRangeWorldCoord/10.0;
		top = right = -kMaxOutOfRangeWorldCoord/10.0;
	}
	else
	#endif
	{
		bottom = left = kMaxOutOfRangeWorldCoord;
		top = right = -kMaxOutOfRangeWorldCoord;
	}
}

inline void WorldRect::Set(const WorldPt& a, const WorldPt& b)
// DOES swap values as necessary to make rect non-empty
{ 
	SetLTRB(min(a.x, b.x), max(a.y, b.y), max(a.x, b.x), min(a.y, b.y));
}

inline void WorldRect::SetCornersTLBR(WorldPt tl, WorldPt br)
// DOES NOT swap values as necessary to make rect non-empty
{
	SetLTRB(tl.x, tl.y, br.x, br.y);
}

inline void WorldRect::Set(const WorldPt& cen, double radius)
{
	SetLTRB(cen.x - abs(radius), cen.y + abs(radius), cen.x + abs(radius), cen.y - abs(radius)); 
}

inline void WorldRect::SetLTRBWithoutValidating(double l, double t, double r, double b)
{
	left   = l; 
	top    = t; 
	right  = r; 
	bottom = b; 
}

inline void WorldRect::SetLTRB(double l, double t, double r, double b)
{ 
	left   = l; 
	top    = t; 
	right  = r; 
	bottom = b; 
	
#if TEST
	static int assertCount = 3;	// if you want more, change value in debugger to large number

	if ((l <= r) && (b <= t)) {   // valid non-empty rect

		// Is it legal?
		if (assertCount > 0) {
			if (!VERIFYN(kMark, ((left   >= -kMaxLegalWorldCoord)
							  && (bottom >= -kMaxLegalWorldCoord)
							  && (top    <=  kMaxLegalWorldCoord)
							  && (right  <=  kMaxLegalWorldCoord)))
			)
				assertCount--;
		}
	}
	else

		// The rect is empty

		if (!this->IsAProperlyDefinedEmptyRect()) {
			if (assertCount > 0) {
				DSTOP((kMark, "Invalid partially empty WorldRect specified. Check for uninitialized values."));
				assertCount--;
			}
		}

#endif
}

inline void WorldRect::AdjustLTRB(double l, double t, double r, double b)
{ 
	left   += l; 
	top    += t; 
	right  += r; 
	bottom += b; 
	
	#if TEST
	static int assertCount = 3;	// if you want more, change value in debugger to large number
	#endif

	if ((left <= right) && (bottom <= top)) {   // valid non-empty rect
		;

		#if TEST
		// Is it legal?
		if (assertCount > 0) {
			if (!VERIFYN(kMark, ((left   >= -kMaxLegalWorldCoord)
							  && (bottom >= -kMaxLegalWorldCoord)
							  && (top    <=  kMaxLegalWorldCoord)
							  && (right  <=  kMaxLegalWorldCoord)))
			)
				assertCount--;
		}
		#endif
	}
	else
		if (this->IsAProperlyDefinedEmptyRect())
			; // valid empty rect
		else {
			#if TEST
			if (assertCount > 0) {
				DSTOP((kMark, "Invalid partially empty WorldRect specified. Check for uninitialized values."));
				assertCount--;
			}
			#endif
			
			this->MakeEmpty();
		}
}

inline Boolean WorldRect::operator==(const WorldRect& r) const 
{ 
	if (WORLDRECTISEMPTY((*this), eDontAssert) && WORLDRECTISEMPTY(r, eDontAssert))
		return true;
	else
		if ((WORLDRECTISEMPTY((*this), eDontAssert)) ^ (WORLDRECTISEMPTY(r, eDontAssert)))
			return false;
		else
			return ((WorldCoordsAreNearlyEqual(top,    r.top)) 
			     && (WorldCoordsAreNearlyEqual(right,  r.right))
				 && (WorldCoordsAreNearlyEqual(left,   r.left))
				 && (WorldCoordsAreNearlyEqual(bottom, r.bottom))); 
}


inline void WorldRect::Inset(double dx, double dy) 
{ 
	if (!WORLDRECTISEMPTY((*this), eAssert)) {
		left   += dx; 
		right  -= dx; 
		bottom += dy; 
		top    -= dy; 
	}
}


inline void WorldRect::Outset(double dx, double dy) 
{ 
	if (!WORLDRECTISEMPTY((*this), eAssert)) {
		left   -= dx; 
		right  += dx; 
		bottom -= dy; 
		top    += dy; 
	}
}


inline void WorldRect::Offset(const WorldPt& p) 
{ 
	if (!WORLDRECTISEMPTY((*this), eAssert)) {
		left   += p.x; 
		right  += p.x; 
		bottom += p.y; 
		top    += p.y; 
	}
}



inline void WorldRect::Unite(const WorldPt& p)
{
	if (WORLDRECTISEMPTY((*this), eDontAssert)) 
		this->Set(p, p);
	else {
		left   = min(left,   p.x);
		right  = max(right,  p.x);
		bottom = min(bottom, p.y);
		top    = max(top,    p.y);
	}
}

inline void WorldRect::UniteNoCheck(const WorldPt& p)
{
	left = min(left, p.x);
	right = max(right, p.x);
	bottom = min(bottom, p.y);
	top = max(top, p.y);
}

inline WorldRect WorldRect::Union(const WorldPt& p) const
{
	WorldRect r = *this; r.Unite(p); return r;
}

inline void WorldRect::Unite(const WorldRect& r)
{
	if (!WORLDRECTISEMPTY(r, eDontAssert)) {
		if (WORLDRECTISEMPTY((*this), eDontAssert))
			*this = r;
		else {
			left 	= min(left,   r.left);
			right 	= max(right,  r.right);
			bottom 	= min(bottom, r.bottom);
			top 	= max(top,    r.top);
		}
	}
}

inline WorldRect WorldRect::Union(const WorldRect& r) const
{
	WorldRect ret = *this; ret.Unite(r); return ret;
}

inline void WorldRect::Intersect(const WorldRect& r)
{
	if (this->Intersects(r)) {
		top		= min(top,    r.top);	
		left	= max(left,   r.left);
		bottom	= max(bottom, r.bottom);
		right	= min(right,  r.right);

		// possible because WorldCoord1_GE_WorldCoord2 will return true if 
		// WorldCoord2 is infinitesimally greater than WorldCoord1, as
		// we consider them equal in that case
		if (bottom > top) 
			top = bottom;
		if (left > right)
			right = left;
	}
	else
		this->MakeEmpty();
}

inline WorldRect WorldRect::Intersection(const WorldRect& r) const
{
	WorldRect ret = *this; ret.Intersect(r); return ret;
}

inline Boolean WorldRect::Intersects(const WorldRect& r) const
// True even if the two rects just share a point.
{ 
	if (WORLDRECTISEMPTY((*this), eDontAssert) || WORLDRECTISEMPTY(r, eDontAssert))
		return false;
	else
		return (WorldCoord1_GE_WorldCoord2(top, r.bottom) 
		     && WorldCoord1_GE_WorldCoord2(right, r.left) 
		     && WorldCoord1_LE_WorldCoord2(left, r.right) 
		     && WorldCoord1_LE_WorldCoord2(bottom, r.top)); 
}

inline Boolean WorldRect::Overlaps(const WorldRect& r) const
// True if the two rects share more than a point or a line.
{ 
	if (WORLDRECTISEMPTY((*this), eDontAssert) || WORLDRECTISEMPTY(r, eDontAssert))
		return false;
	else
		return (WorldCoord1_GT_WorldCoord2(top, r.bottom)
		     && WorldCoord1_GT_WorldCoord2(right, r.left)
			 && WorldCoord1_LT_WorldCoord2(left, r.right) 
			 && WorldCoord1_LT_WorldCoord2(bottom, r.top)); 
}

inline Boolean WorldRect::ContainsPoint(const WorldPt& p) const
{ 
	if (WORLDRECTISEMPTY((*this), eDontAssert))
		return false;
	else
		return (WorldCoord1_LE_WorldCoord2(left, p.x)
		     && WorldCoord1_GE_WorldCoord2(right, p.x)
			 && WorldCoord1_LE_WorldCoord2(bottom, p.y)
			 && WorldCoord1_GE_WorldCoord2(top, p.y));
}

inline Boolean WorldRect::ContainsRect(const WorldRect& r) const
{ 
	if (WORLDRECTISEMPTY((*this), eDontAssert) || WORLDRECTISEMPTY(r, eDontAssert))
		return false;
	else
		return (WorldCoord1_LE_WorldCoord2(left, r.left)
		    &&  WorldCoord1_GE_WorldCoord2(right, r.right)
			&&  WorldCoord1_LE_WorldCoord2(bottom, r.bottom)
			&&  WorldCoord1_GE_WorldCoord2(top, r.top)); 
}

inline Boolean WorldRect::IsEmpty(EAssertStatus as) const 
{
	as;
	#if TEST
	static int assertCount = 4;	// if you want more, change value in debugger to large number
	if (as == eAssert && assertCount > 0) {
		if ((top < bottom) || (right < left)) {
			--assertCount;
			DSTOP((kMark, "Trying to operate on or access individual fields of an empty rect."));
		}
		if (((top < bottom) ^ (right < left)) && assertCount > 0) {
			--assertCount;
			DSTOP((kMark, "If Rect is empty in one dimension it should be empty in the other dimension too"));
		}
	}
	#endif
	
	return ((top < bottom) || (right < left)); 
}

inline Boolean WorldRect::IsAProperlyDefinedEmptyRect() const
{
	return ((top    == -kMaxOutOfRangeWorldCoord)
	&&      (right  == -kMaxOutOfRangeWorldCoord)
	&&      (bottom ==  kMaxOutOfRangeWorldCoord)
	&&      (left   ==  kMaxOutOfRangeWorldCoord));
}

const extern WorldRect kEmptyWorldRect;

inline WorldRect EmptyWorldRect()
{
//	return WorldRect();
	return kEmptyWorldRect;
}

//--------------------------------------------------------------------
#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>

struct _WorldCube
//
// This structure defines a bounding Grid-Aligned
// Rectilinear Parallelepiped in 3-D space (I
// suppose one could call that a GARP, but somehow
// I doubt that name would catch on.)  This base
// version exists so we can put them in unions as
// it has no constructor.
//
{
	friend class WorldCube; // updated for compatibilty to CW 6 project
	//protected:
	double				right, top, front;
	double				left, bottom, back;
	public:
	
//	operator const WorldCube&() const		{ return *((const WorldCube *) this); }
//	operator WorldCube& ()					{ return *((WorldCube *) this); }
	double				MinX() const { return left; }
	double				MinY() const { return bottom; }
	double				MinZ() const { return back; }
	double				MaxX() const { return right; }
	double				MaxY() const { return top; }
	double				MaxZ() const { return front; }

	//WorldCube* operator&(void)				{ return (WorldCube *) this; }
	//const WorldCube* operator&(void) const	{ return (const WorldCube *) this; }

	//void ByteSwap();
};

#pragma pack(pop)	// END 2 BYTE STUCTURE AIGNMENT <.vwx file format>
//--------------------------------------------------------------------


//--------------------------------------------------------------------
#pragma pack(push, 2)	// BEGIN 2 BYTE STUCTURE AIGNMENT <.vwx file format>

class WorldCube : public _WorldCube
//
// This derivation of the base _WorldCube adds
// constructors, accessors, etc.  The idea is that
// this class should always be maintained in a
// consistent state (max >= min for X,Y, and Z),
// unless it is empty.  Currently, emptiness is
// handled by setting all values to 0, which is a
// poor way to represent emptiness -- it's the
// same as a point-sized cube at the origin, and
// shouldn't be.
//
{
public:
	WorldCube(void) { MakeEmpty(); }
	WorldCube(const WorldCube& c);
	WorldCube(double top, double bottom, double right, double left, double front, double back)
		{ this->SetYyXxZz(top, bottom, right, left, front, back); }
	WorldCube(const WorldPt3& a, const WorldPt3& b);
	WorldCube(const _WorldCube& c);
		
	WorldPt3 High(void) const;
	WorldPt3 Low(void) const;
	double				MinX() const { return this->_WorldCube::left; }
	double				MinY() const { return this->_WorldCube::bottom; }
	double				MinZ() const { return this->_WorldCube::back; }
	double				MaxX() const { return this->_WorldCube::right; }
	double				MaxY() const { return this->_WorldCube::top; }
	double				MaxZ() const { return this->_WorldCube::front; }
	double				Left() const { return this->_WorldCube::left; }
	double				Bottom() const { return this->_WorldCube::bottom; }
	double				Back() const { return this->_WorldCube::back; }
	double				Right() const { return this->_WorldCube::right; }
	double				Top() const { return this->_WorldCube::top; }
	double				Front() const { return this->_WorldCube::front; }

	//void ByteSwap() { this->_WorldCube::ByteSwap(); }

	void Clear() { MakeEmpty(); } 
	Boolean IsAProperlyDefinedEmptyCube() const;
	void MakeEmpty(EEmptyKind kind = eStandardValue);
	void Set(double top, double bottom, double right, double left, double front, double back)
		{ this->SetYyXxZz(top, bottom, right, left, front, back); }
	void SetYyXxZz(double top, double bottom, double right, double left, double front, double back);
	void SetYyXxZzWithoutValidating(double top, double bottom, double right, double left, double front, double back);
	void Set(const WorldPt3& a, const WorldPt3& b);
	// This is weird, a "unit" cube with edge lengths of two.  This is only used one place in
	// the code, and thus really shouldn't be a member function.
	void SetToUnitCube()
		//{ top = 1.0; bottom = -1.0; right = 1.0; left = -1.0; front = 1.0; back = -1.0; }
		{ this->SetYyXxZz(1.0, -1.0, 1.0, -1.0, 1.0, -1.0); }
	void SetToPoint(const WorldPt3& inP)
		{ top = bottom = inP.y; right = left = inP.x; front = back = inP.z; }
	void SetToPoint(const double& x, const double& y, const double& z)
		{ top = bottom = y; right = left = x; front = back = z; }
	void SetHighLow(WorldPt3 high, WorldPt3 low)
		{ SetYyXxZz(high.y, low.y, high.x, low.x, high.z, low.z); }
		
	WorldCube& operator=(const WorldCube& c) { memcpy(this, &c, sizeof(WorldCube)); return *this; }

	Boolean operator==(const WorldCube& c) const;
	Boolean operator!=(const WorldCube& c) const;
	WorldCube&	operator*=(double s);
	WorldCube&	operator/=(double s);  
	WorldCube	operator*(double s) const; 
	WorldCube	operator/(double s) const;

	void Offset(const WorldPt3 &p);
	WorldCube GetOffset(const WorldPt3& p) const;
	void UniteOffset(const WorldPt3& p);
	void Inset(double dx, double dy, double dz);

	void Unite(const WorldPt3 &p);
	void Unite(const WorldCube &c);
	void Unite(const WorldRect &r, const double z);
	void Unite(const double& x, const double& y, const double& z);
	
	Boolean Intersects(const WorldCube& c) const;
	void IntersectWith(const WorldCube& c);	
	WorldCube IntersectionWith(const WorldCube& c) const;
	WorldCube Union(const WorldCube &c) const;
	
	Boolean IsEmpty(EAssertStatus as = eDontAssert) const
		{ 	
			if (as == eAssert) {
				#if TEST
				static int assertCount = 4;	// if you want more, change value in debugger to large number
				if (as == eAssert && assertCount > 0) {
					if ((left > right || bottom > top || back > front)) 
					{
						--assertCount;
						DSTOP((kMCCoordTypes, "Trying to operate on or access individual fields of an empty cube."));
					}
					if (((top < bottom) ^ (right < left)) && assertCount > 0) 
					{
						--assertCount;
						DSTOP((kMCCoordTypes, "If cube is empty in one dimension it should be empty in the other dimensions too"));
					}
				}
				#endif
			}
			return (left > right || bottom > top || back > front);
		}
	Boolean IsZero(void) const
		{ return WORLDCUBEISEMPTY((*this), eDontAssert); }

	WorldPt3 Center(void) const
		{ return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3((left + right) / 2.0, (bottom + top) / 2.0, (back + front) / 2.0); }
	double CenterX(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : (left + right) / 2.0; }
	double CenterY(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : (bottom + top) / 2.0; }
	double CenterZ(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : (back + front) / 2.0; }

	Boolean ContainsPoint(const WorldPt3& p) const
		{ return WorldCoord1_GE_WorldCoord2(p.x, left) && 
				WorldCoord1_LE_WorldCoord2(p.x, right) && 
				WorldCoord1_GE_WorldCoord2(p.y, bottom) && 
				WorldCoord1_LE_WorldCoord2(p.y, top) && 
				WorldCoord1_GE_WorldCoord2(p.z, back) && 
				WorldCoord1_LE_WorldCoord2(p.z, front); }
	Boolean ContainsCube(const WorldCube& c) const
		{ WorldCube temp(*this); temp.Unite(c); return (temp == *this); }

	double Height(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : top-bottom; }
	double Width(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : right-left; }
	double Depth(void) const	{ return WORLDCUBEISEMPTY((*this), eAssert) ? 0 : front-back; }

	// lower case means low value, upper case high value
	WorldPt3 Pointxyz() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MinX(), MinY(), MinZ()); } // left, bottom, back
	WorldPt3 PointxyZ() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MinX(), MinY(), MaxZ()); } // left, bottom, front
	WorldPt3 PointxYz() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MinX(), MaxY(), MinZ()); } // left, top, back
	WorldPt3 PointxYZ() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MinX(), MaxY(), MaxZ()); } // left, top, front
	WorldPt3 PointXyz() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MaxX(), MinY(), MinZ()); } // right, bottom, back
	WorldPt3 PointXyZ() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MaxX(), MinY(), MaxZ()); } // right, bottom, front
	WorldPt3 PointXYz() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MaxX(), MaxY(), MinZ()); } // right, top, back
	WorldPt3 PointXYZ() const { return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(MaxX(), MaxY(), MaxZ()); } // right, top, front

	WorldRect FrontWorldRect(void) const
		{ WorldRect r; r.top = top; r.left = left; r.bottom = bottom; r.right = right; return r; }
};

#pragma pack(pop) // END 2 BYTE STRUCTURE ALIGNMENT (.vwx file format)
//--------------------------------------------------------------------

inline WorldCube::WorldCube(const WorldCube& c)
{
//	right = c.right; top = c.top; front = c.front;
//	left = c.left; bottom = c.bottom; back = c.back; 

	memcpy(this, &c, sizeof(WorldCube));    // faster for misaligned doubles [MAF FP Speedup, 2/21/01]
}

inline WorldCube::WorldCube(const _WorldCube& c)
{
	right = c.right; top = c.top; front = c.front;
	left = c.left; bottom = c.bottom; back = c.back; 
}
		
inline WorldPt3 WorldCube::High(void) const	
{ 
	return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(right, top, front);
}

inline WorldPt3 WorldCube::Low(void) const		
{ 
	return WORLDCUBEISEMPTY((*this), eAssert) ? WorldPt3(0,0,0) : WorldPt3(left, bottom, back); 
}

inline Boolean WorldCube::IsAProperlyDefinedEmptyCube() const
{
	return ((top    == -kMaxOutOfRangeWorldCoord)
	&&      (right  == -kMaxOutOfRangeWorldCoord)
	&&      (front  == -kMaxOutOfRangeWorldCoord)
	&&      (bottom ==  kMaxOutOfRangeWorldCoord)
	&&      (left   ==  kMaxOutOfRangeWorldCoord)
	&&      (back   ==  kMaxOutOfRangeWorldCoord));
}

inline void WorldCube::MakeEmpty(EEmptyKind kind) 
{
	kind;
	#if TEST

	if (kind == eSpecialValue) {
		left = bottom = back = kMaxOutOfRangeWorldCoord/10.0;
		right = top = front = -kMaxOutOfRangeWorldCoord/10.0;
	}
	else
	#endif
	{
		left = bottom = back = kMaxOutOfRangeWorldCoord;
		right = top = front = -kMaxOutOfRangeWorldCoord;
	}
}

inline void WorldCube::SetYyXxZz(double top, double bottom, double right, double left, double front, double back)
{ 
	if( (left <= right && bottom <= top && back <= front)) {
		this->top = top; this->bottom = bottom; this->right = right;
		this->left = left; this->front = front; this->back = back; 
	}
}

inline void WorldCube::SetYyXxZzWithoutValidating(double top, double bottom, double right, double left, double front, double back)
{ 
	this->top = top; this->bottom = bottom; this->right = right;
	this->left = left; this->front = front; this->back = back; 
}

inline Boolean WorldCube::operator==(const WorldCube& c) const
{ 
	return WorldCoordsAreNearlyEqual(right,  c.right) 
		&& WorldCoordsAreNearlyEqual(top,    c.top) 
	    && WorldCoordsAreNearlyEqual(front,  c.front) 
		&& WorldCoordsAreNearlyEqual(left,   c.left) 
		&& WorldCoordsAreNearlyEqual(bottom, c.bottom) 
		&& WorldCoordsAreNearlyEqual(back,   c.back); 
}

inline Boolean WorldCube::operator!=(const WorldCube& c) const
{ 
	return !(*this == c); 
}

inline WorldCube&	WorldCube::operator*=(double s) 
{
	if ( !WORLDCUBEISEMPTY((*this), eAssert) )
		left *= s; right *= s; top *= s; bottom *= s; front *= s; back *= s;
	return *this;
}

inline WorldCube&	WorldCube::operator/=(double s)  
{
	if ( !WORLDCUBEISEMPTY((*this), eAssert) )
		left /= s; right /= s; top /= s; bottom /= s; front /= s; back /= s;
	return *this;
}

inline WorldCube	WorldCube::operator*(double s) const 
{
	if ( WORLDCUBEISEMPTY((*this), eAssert))
		return *this;
	return WorldCube(top * s, bottom * s, right * s, left * s, front * s, back * s);
}

inline WorldCube	WorldCube::operator/(double s) const
{
	if ( WORLDCUBEISEMPTY((*this), eAssert))
		return *this;
	return WorldCube(top / s, bottom / s, right / s, left / s, front / s, back / s);
}

inline void WorldCube::Offset(const WorldPt3 &p)
{
	if ( !WORLDCUBEISEMPTY((*this), eAssert))
		 left += p.x; right += p.x; bottom += p.y; top += p.y; back += p.z; front += p.z;
}

inline WorldCube WorldCube::GetOffset(const WorldPt3& p) const
{
	if ( WORLDCUBEISEMPTY((*this), eAssert))
		return *this;
	 return WorldCube(top + p.y, bottom + p.y, right + p.x, left + p.x, front + p.z, back + p.z);
}

inline void WorldCube::UniteOffset(const WorldPt3& p)
{
	if ( !WORLDCUBEISEMPTY((*this), eAssert))
		this->Unite(this->GetOffset(p));
}

inline void WorldCube::Inset(double dx, double dy, double dz)
{
	if ( !WORLDCUBEISEMPTY((*this), eAssert)) {
		left += dx; right -= dx; bottom += dy; top -= dy; back += dz; front -= dz;
	}
}

inline void WorldCube::Unite(const WorldPt3 &p)
{ 
	this->Unite(p.x, p.y, p.z); 
}	

inline void WorldCube::Unite(const WorldCube &c)
{ 
	if (WORLDCUBEISEMPTY((*this), eDontAssert)) 
		*this = c;
	else if (!WORLDCUBEISEMPTY(c, eDontAssert)) {
		right = max(right, c.right); top = max(top, c.top); front = max(front, c.front);
		left = min(left, c.left); bottom = min(bottom, c.bottom); back = min(back, c.back); 
	}
}

inline void WorldCube::Unite(const WorldRect &r, const double z)
{
	if (!WORLDRECTISEMPTY(r, eDontAssert)) {
		if (WORLDCUBEISEMPTY((*this), eDontAssert)) 
			this->SetYyXxZz(r.Top(), r.Bottom(), r.Right(), r.Left(), z, z);
		else {
			right = max(right, r.Right()); top = max(top, r.Top()); front = max(front, z);
			left = min(left, r.Left()); bottom = min(bottom, r.Bottom()); back = min(back, z); 
		}
	}
}

inline void WorldCube::Unite(const double& x, const double& y, const double& z)
{ 
	if (WORLDCUBEISEMPTY((*this), eDontAssert)) 
		this->SetToPoint(x, y, z);
	else {
		right = max(right, x); top = max(top, y); front = max(front, z);
		left  = min(left, x); bottom = min(bottom, y); back = min(back, z); 
	}
}


inline void WorldCube::IntersectWith(const WorldCube& c)
{ 
	if (this->Intersects(c)) {
		right = min(right, c.right); top = min(top, c.top); front = min(front, c.front);
		left = max(left, c.left); bottom = max(bottom, c.bottom); back = max(back, c.back); 
	}
	else
		this->MakeEmpty();
}

inline WorldCube WorldCube::IntersectionWith(const WorldCube& c) const
{
	WorldCube rc = *this;
	rc.IntersectWith(c);
	return rc;
}

inline Boolean WorldCube::Intersects(const WorldCube& c) const
// True even if the two cubes just share a point.
{ 
	if (WORLDCUBEISEMPTY((*this), eDontAssert) || WORLDCUBEISEMPTY(c, eDontAssert))
		return false;
	else
		return (WorldCoord1_GE_WorldCoord2(top, c.bottom) 
		     && WorldCoord1_LE_WorldCoord2(bottom, c.top)
		     && WorldCoord1_GE_WorldCoord2(right, c.left) 
		     && WorldCoord1_LE_WorldCoord2(left, c.right) 
		     && WorldCoord1_GE_WorldCoord2(front, c.back) 
		     && WorldCoord1_LE_WorldCoord2(back, c.front)); 
}

TransformMatrix MakeScaleMatrix(double scale);
bool WorldPt3NearlyEqual(const WorldPt3 & p1, const WorldPt3 & p2);
bool IsPointNearLine(const WorldPt & inPoint, const WorldPt & inP0, const WorldPt & inP1, double inNearDist, int inLineType);

#endif
