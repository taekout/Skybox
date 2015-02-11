#ifndef GLDebug_h
#define GLDebug_h

#define _WINDOWS 1
#define GL_DEBUG 1
#define EXPORT_GEO_FOR_CULLING_DEMO 0 // uncomment this line to use - hold CTRL when turning on opengl to write geometry to c:\a\model.txt

#if _WINDOWS 
	#include <GL/gl.h>
	#include <GL/glu.h>
	#include <gl/glext.h>
	//#include "OPGLRenderer.Win.h"
#else
	#include <AGL/glu.h>
	#include <AGL/glm.h>
	#include <AGL/aglContext.h>
	#include "glext.h"
	#include "OPGLRenderer.Mac.h"
#endif
#include <vector>
#include <fstream>

using namespace std;

#if GL_DEBUG

void GL_DebugPrintState2();
void GL_DebugPrintState(ostream * outStream = 0);
void GL_BeginLoggingCalls();
void GL_EndLoggingCalls();
void GL_SetFunctionLogability(const char * inFunctionName, bool inLog);

void _glPixelStorei(GLenum inEnum, int inValue);
void _glGetFloatv(GLenum inEnum, float * outFloat);
void _glDrawBuffer(GLenum inEnum);
void _glClearColor(float inR, float inG, float inB, float inA);
void _glClear(GLenum inEnum);
void _glReadBuffer(GLenum inEnum);
void _glColorMask(GLboolean inR, GLboolean inG, GLboolean inB, GLboolean inA);
void _glShadeModel(GLenum inEnum);
void _glPolygonOffset(float inFactor, float inBias);
void _glMatrixMode(GLenum inEnum);
void _glLoadMatrixf(GLfloat * inMat);
void _glEnable(GLenum inEnum);
void _glDisable(GLenum inEnum);
void _glBindBuffer(GLenum target, GLuint buffer);
void _glUnmapBuffer(GLenum target);
void _glBindFramebuffer(GLenum inEnum, GLuint inFBO);
void _glFramebufferRenderbuffer(GLenum inEnum1, GLenum inEnum2, GLenum inEnum3, GLuint inRB);
void _glBlitFramebuffer(GLint inX0, GLint inY0, GLint inW0, GLint inH0, GLint inX1, GLint inY1, GLint inW1, GLint inH1, GLbitfield inBitfield, GLenum inEnum);
void _glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels);
void _glFlush();
void _glFinish();
void _glLoadIdentity();
void _glViewport(GLint x, GLint y, GLsizei width, GLsizei height);
void _glFramebufferTexture2D(GLenum inEnum1, GLenum inEnum2, GLenum inEnum3, GLuint inUint, GLint inInt);
void _glBindTexture(GLenum target, GLuint texture);
void _glGetTexImage(GLenum target, GLint level, GLenum format, GLenum type, GLvoid *pixels);
void _glActiveTexture(GLenum texture);
void _glTexSubImage2D(GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum format, GLenum type, const GLvoid *pixels);
void _glDrawBuffers(GLsizei n, const GLenum *bufs);
void _glUseProgram(GLuint program);
void _glUniform3fv(GLint location, GLsizei count, const GLfloat *value);
void _glUniform1i(GLint location, GLint v0);
void _glUniform1fv(GLint location, GLsizei count, const GLfloat *value);
void _glUniform2fv(GLint location, GLsizei count, const GLfloat *value);
void _glUniform1iv(GLint location, GLsizei count, const GLint *value);
void _glBlendFunc(GLenum sfactor, GLenum dfactor);
void _glDepthMask(GLboolean flag);
GLenum _glCheckFramebufferStatus(GLenum inEnum);
void _glGetFramebufferAttachmentParameterivEXT(GLenum target, GLenum attachment, GLenum pname, GLint *params);
void _glGetBooleanv(GLenum pname, GLboolean * params);
void _glGetDoublev(GLenum pname, GLdouble * params);
void _glGetFloatv(GLenum pname, GLfloat * params);
void _glGetIntegerv(GLenum pname, GLint * params);
void _glCallList(GLuint inList);
void _glBegin(GLenum type);
void _glEnd();
void _glNormal3d(GLdouble x, GLdouble y, GLdouble z);
void _glNormal3dv(GLdouble * verts);
void _glVertex3d(GLdouble x, GLdouble y, GLdouble z);
void _glVertex3dv(GLdouble * verts);
void _glEdgeFlag(GLboolean flag);
void _glEdgeFlagv(GLboolean * flag);

#define glEdgeFlagv _glEdgeFlag
#define glEdgeFlag _glEdgeFlag
#define glVertex3dv _glVertex3dv
#define glVertex3d _glVertex3d
#define glNormal3d _glNormal3d
#define glNormal3d _glNormal3d
#define glEnd _glEnd
#define glBegin _glBegin
#define glCallList _glCallList
#define glGetIntegerv _glGetIntegerv
#define glGetFloatv _glGetFloatv
#define glGetDoublev _glGetDoublev
#define glGetBooleanv _glGetBooleanv
#define glGetFramebufferAttachmentParameterivEXT _glGetFramebufferAttachmentParameterivEXT
#define glCheckFramebufferStatus _glCheckFramebufferStatus
#define glDepthMask _glDepthMask
#define glBlendFunc _glBlendFunc
#define glUniform1iv _glUniform1iv
#define glUniform2fv _glUniform2fv
#define glUniform1fv _glUniform1fv
#define glUniform1i _glUniform1i
#define glUniform3fv _glUniform3fv
#define glUseProgram _glUseProgram
#define glDrawBuffers _glDrawBuffers
#define glTexSubImage2D _glTexSubImage2D
#define glActiveTexture _glActiveTexture
#define glBindTexture _glBindTexture
#define glGetTexImage _glGetTexImage
#define glFramebufferTexture2D _glFramebufferTexture2D
#define glViewport _glViewport
#define glLoadIdentity _glLoadIdentity
#define glFinish _glFinish
#define glFlush _glFlush
#define glReadPixels _glReadPixels
#define glBlitFramebuffer _glBlitFramebuffer
#define glFramebufferRenderbuffer _glFramebufferRenderbuffer
#define glBindFramebuffer _glBindFramebuffer
#define glUnmapBuffer _glUnmapBuffer
#define glBindBuffer _glBindBuffer
#define glDisable _glDisable
#define glEnable _glEnable
#define glLoadMatrixf _glLoadMatrixf
#define glMatrixMode _glMatrixMode
#define glPolygonOffset _glPolygonOffset
#define glShadeModel _glShadeModel
#define glColorMask _glColorMask	
#define glReadBuffer _glReadBuffer	
#define glClear _glClear
#define glClearColor _glClearColor
#define glDrawBuffer _glDrawBuffer
#define glGetFloatv _glGetFloatv
#define glPixelStorei _glPixelStorei

#if EXPORT_GEO_FOR_CULLING_DEMO

struct PipelineEnv;

struct SymGeo {
	vector<double> normals;
	vector<double> vertices;
	vector<int> newObject; // index into vertices of the start of a new object
};
extern SymGeo * gSymGeo;
extern map<BlankHandle, SymGeo *> gSymGeoMap;
extern bool gWritingMeshData;
extern bool gWritingVertexData;
extern ofstream geometryStream;
extern PipelineEnv * gEnv;

void CULL_glBegin(GLenum type);
void CULL_glEnd();
void CULL_glNormal3d(GLdouble x, GLdouble y, GLdouble z);
void CULL_glNormal3dv(GLdouble * verts);
void CULL_glVertex3d(GLdouble x, GLdouble y, GLdouble z);
void CULL_glVertex3dv(GLdouble * verts);

#endif


#endif

#endif
