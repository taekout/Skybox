#include "StdHeaders.h"	

#pragma warning(disable:4996)

#include "GLDebug.h"
#include <map>
#include <string>
#include <time.h>

using namespace std;

#if GL_DEBUG

#undef glEdgeFlag 
#undef glEdgeFlagv
#undef glVertex3dv
#undef glVertex3d 
#undef glNormal3d 
#undef glNormal3d 
#undef glEnd
#undef glBegin
#undef glCallList
#undef glGetIntegerv
#undef glGetFloatv
#undef glGetDoublev
#undef glGetBooleanv
#undef glGetFramebufferAttachmentParameterivEXT
#undef glCheckFramebufferStatus
#undef glDepthMask
#undef glBlendFunc
#undef glUniform1iv
#undef glUniform2fv
#undef glUniform1fv
#undef glUniform1i
#undef glUniform3fv
#undef glUseProgram
#undef glDrawBuffers
#undef glTexSubImage2D
#undef glActiveTexture
#undef glBindTexture
#undef glGetTexImage
#undef glFramebufferTexture2D
#undef glViewport
#undef glLoadIdentity
#undef glFinish
#undef glFlush
#undef glReadPixels
#undef glBlitFramebuffer
#undef glFramebufferRenderbuffer
#undef glBindFramebuffer
#undef glUnmapBuffer
#undef glBindBuffer
#undef glDisable
#undef glEnable
#undef glLoadMatrixf
#undef glMatrixMode
#undef glPolygonOffset
#undef glShadeModel
#undef glColorMask
#undef glReadBuffer
#undef glClear
#undef glClearColor
#undef glDrawBuffer
#undef glGetFloatv
#undef glPixelStorei

// OpenGL 2.1 functions (dynamically located via wglGetProcAdderss)
PFNGLCREATEPROGRAMPROC				glCreateProgram						= 0;
PFNGLDELETEPROGRAMPROC				glDeleteProgram						= 0;
PFNGLISPROGRAMPROC					glIsProgram							= 0;
PFNGLATTACHSHADERPROC				glAttachShader						= 0;
PFNGLDETACHSHADERPROC				glDetachShader						= 0;
PFNGLISSHADERPROC					glIsShader							= 0;
PFNGLLINKPROGRAMPROC				glLinkProgram						= 0;
PFNGLUSEPROGRAMPROC					glUseProgram						= 0;
PFNGLGETATTACHEDSHADERSPROC			glGetAttachedShaders				= 0;
PFNGLGETPROGRAMIVPROC				glGetProgramiv						= 0;
PFNGLGETPROGRAMINFOLOGPROC			glGetProgramInfoLog					= 0;
PFNGLVALIDATEPROGRAMPROC			glValidateProgram					= 0;
PFNGLGETSHADERINFOLOGPROC			glGetShaderInfoLog					= 0;
PFNGLGETUNIFORMLOCATIONPROC			glGetUniformLocation				= 0;
PFNGLGETUNIFORMIVPROC				glGetUniformiv						= 0;
PFNGLUNIFORM1IPROC					glUniform1i							= 0;
PFNGLUNIFORM1IVPROC					glUniform1iv						= 0;
PFNGLUNIFORM2IVPROC					glUniform2iv						= 0;
PFNGLUNIFORM3IVPROC					glUniform3iv						= 0;			
PFNGLUNIFORM4IVPROC					glUniform4iv						= 0;
PFNGLUNIFORM1FPROC					glUniform1f							= 0;
PFNGLUNIFORM2FPROC					glUniform2f							= 0;
PFNGLUNIFORM1FVPROC					glUniform1fv						= 0;
PFNGLUNIFORM2FVPROC					glUniform2fv						= 0;
PFNGLUNIFORM3FVPROC					glUniform3fv						= 0;
PFNGLUNIFORM4FVPROC					glUniform4fv						= 0;
PFNGLUNIFORMMATRIX4FVPROC			glUniformMatrix3fv					= 0;
PFNGLUNIFORMMATRIX4FVPROC			glUniformMatrix4fv					= 0;
PFNGLCREATESHADERPROC				glCreateShader						= 0;
PFNGLDELETESHADERPROC				glDeleteShader						= 0;
PFNGLSHADERSOURCEPROC				glShaderSource						= 0;
PFNGLCOMPILESHADERPROC				glCompileShader						= 0;
PFNGLGETSHADERIVPROC				glGetShaderiv						= 0;
PFNGLACTIVETEXTUREPROC				glActiveTexture						= 0;
PFNGLBLENDFUNCSEPARATEPROC			glBlendFuncSeparate					= 0;

// from the extionsion ARB_framebuffer_object
PFNGLGENFRAMEBUFFERSPROC			glGenFramebuffers					= 0;
PFNGLDELETEFRAMEBUFFERSPROC			glDeleteFramebuffers				= 0;
PFNGLDRAWBUFFERSPROC				glDrawBuffers						= 0;
PFNGLBINDFRAMEBUFFERPROC			glBindFramebuffer					= 0; 
PFNGLCHECKFRAMEBUFFERSTATUSPROC		glCheckFramebufferStatus			= 0;
PFNGLFRAMEBUFFERTEXTURE2DPROC		glFramebufferTexture2D				= 0;
PFNGLISFRAMEBUFFERPROC				glIsFramebuffer						= 0;
PFNGLBLENDEQUATIONPROC				glBlendEquation						= 0;
PFNGLISRENDERBUFFERPROC				glIsRenderbuffer					= 0;
PFNGLBINDRENDERBUFFERPROC			glBindRenderbuffer					= 0;
PFNGLDELETERENDERBUFFERSPROC		glDeleteRenderbuffers				= 0;
PFNGLGENRENDERBUFFERSPROC			glGenRenderbuffers					= 0;
PFNGLRENDERBUFFERSTORAGEPROC		glRenderbufferStorage				= 0;
PFNGLGETRENDERBUFFERPARAMETERIVPROC glGetRenderbufferParameteriv		= 0;
PFNGLFRAMEBUFFERRENDERBUFFERPROC	glFramebufferRenderbuffer			= 0;
PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVEXTPROC glGetFramebufferAttachmentParameterivEXT = 0;
PFNGLBLITFRAMEBUFFERPROC			glBlitFramebuffer					= 0;
PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC glRenderbufferStorageMultisample = 0;
PFNGLSAMPLECOVERAGEPROC				glSampleCoverage					= 0;

bool gGLExtensionSupported_GL_EXT_draw_buffers2 = false;
PFNGLENABLEINDEXEDEXTPROC			glEnableIndexedEXT					= 0;
PFNGLDISABLEINDEXEDEXTPROC			glDisableIndexedEXT					= 0;
	
bool gGLExtensionSupported_GL_ARB_draw_buffers_blend = false;
PFNGLBLENDEQUATIONIARBPROC			glBlendEquationiARB					= 0;
PFNGLBLENDEQUATIONSEPARATEIARBPROC	glBlendEquationSeparateiARB			= 0;
PFNGLBLENDFUNCIARBPROC				glBlendFunciARB						= 0;
PFNGLBLENDFUNCSEPARATEIARBPROC		glBlendFuncSeparateiARB				= 0;

PFNGLBINDBUFFERPROC					glBindBuffer						= 0;
PFNGLDELETEBUFFERSPROC				glDeleteBuffers						= 0;
PFNGLGENBUFFERSPROC					glGenBuffers						= 0;
PFNGLISBUFFERPROC					glIsBuffer							= 0;
PFNGLBUFFERDATAPROC					glBufferData						= 0;
PFNGLBUFFERSUBDATAPROC				glBufferSubData						= 0;
PFNGLGETBUFFERSUBDATAPROC			glGetBufferSubData					= 0;
PFNGLMAPBUFFERPROC					glMapBuffer							= 0;
PFNGLUNMAPBUFFERPROC				glUnmapBuffer						= 0;
PFNGLGETBUFFERPARAMETERIVPROC		glGetBufferParameteriv				= 0;
PFNGLGETBUFFERPOINTERVPROC			glGetBufferPointerv					= 0;


char gDebugMessage[1024];
#define kEveryone 1
#define kJHutchison 2
void DMSG2(int, char * message, ...) {
	va_list args;
	va_start(args, message);

	sprintf_s(gDebugMessage, 1024, message, args);
	OutputDebugString(gDebugMessage);

	va_end(args);
}

#define DMSG(params) DMSG2 params


// logging gl function calls file
static ofstream gGLDebugOut;
static int		gTimeOfLastLog = 0;
static int		gTimeOfLastStateDump = 0;
static int		gCallsSinceLastStateDump = 0;
static map<string, bool> gFunctionLogability; // function names not in here are logged; otherwise the value says whether the function is logable

void GL_SetFunctionLogability(const char * inFunctionName, bool inLog) {
	gFunctionLogability[inFunctionName] = inLog;
}

void OutputGLStateDump() {
	int tick = GetTickCount();
	float dt = (tick - gTimeOfLastStateDump) / 1000.0f;

	if( !gGLDebugOut.is_open() || dt<30 || gCallsSinceLastStateDump==0 ) {
		return;
	}

	gTimeOfLastStateDump = tick;
	gCallsSinceLastStateDump = 0;

	gGLDebugOut << endl << endl << "**********GL State Dump**************\n";

	struct tm *local;
	time_t t;

	t = time(NULL);
	local = localtime(&t);
	
	gGLDebugOut << asctime(local) << endl;
	GL_DebugPrintState(&gGLDebugOut);
	gGLDebugOut << "*************************************\n\n";
	gGLDebugOut.flush();

}

void GL_BeginLoggingCalls() 
// To use the OpenGL function logger call GL_BeginLoggingCalls from the debugger, let opengl
// run, then call GL_EndLoggingCalls from the debugger to stop logging and view the results at c:\a\glDebug.txt
{
	gTimeOfLastLog = GetTickCount();
	gTimeOfLastStateDump = 0;
	gCallsSinceLastStateDump = 0;

	if( !gGLDebugOut.is_open() ) {
		gGLDebugOut.open("c:\\a\\glDebug.txt", ios_base::app);
	}

	OutputGLStateDump();
}

void GL_EndLoggingCalls() 
// Stop logging calls
{
	if( gGLDebugOut.is_open() ) {
		gGLDebugOut.close();
	}
}

void GL_DebugPrintState2() {
	GL_DebugPrintState(0);
}

void GL_DebugPrintState(ostream * outStream) {
	int i;
#if GS_WIN
	i = (int)wglGetCurrentContext();
	DMSG((kEveryone, "current gl context = %d\n", i));
#endif
	
	char text[2048];
	int len = 0;

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	len += sprintf(text+len, "viewport = %d, %d, %d, %d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
	
	float mat[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mat);
	len += sprintf(text+len, "modelview mat = \n\t%f, %f, %f, %f,\n\t%f, %f, %f, %f, \n\t%f, %f, %f, %f, \n\t%f, %f, %f, %f\n", mat[0], mat[4], mat[8], mat[12], mat[1], mat[5], mat[9], mat[13], mat[2], mat[6], mat[10], mat[14], mat[3], mat[7], mat[11], mat[15]);

	glGetFloatv(GL_PROJECTION_MATRIX, mat);
	len += sprintf(text+len, "projection mat = \n\t%f, %f, %f, %f,\n\t%f, %f, %f, %f, \n\t%f, %f, %f, %f, \n\t%f, %f, %f, %f\n", mat[0], mat[4], mat[8], mat[12], mat[1], mat[5], mat[9], mat[13], mat[2], mat[6], mat[10], mat[14], mat[3], mat[7], mat[11], mat[15]);

	GLboolean b;
	b = glIsEnabled(GL_LIGHTING);
	len += sprintf(text+len, "lighting enabled = %d\n", b);

	b = glIsEnabled(GL_BLEND);
	len += sprintf(text+len, "blending enabled = %d\n", b);

	b = glIsEnabled(GL_TEXTURE_2D);
	len += sprintf(text+len, "texturing enabled = %d\n", b);

	b = glIsEnabled(GL_DEPTH_TEST);
	len += sprintf(text+len, "depth test enabled = %d\n", b);

	b = glIsEnabled(GL_STENCIL_TEST);
	len += sprintf(text+len, "stencil test enabled = %d\n", b);

	b = glIsEnabled(GL_ALPHA_TEST);
	len += sprintf(text+len, "alpha test enabled = %d\n", b);

	b = glIsEnabled(GL_SCISSOR_TEST);
	len += sprintf(text+len, "scissor test enabled = %d\n", b);

	b = glIsEnabled(GL_POLYGON_OFFSET_FILL);
	len += sprintf(text+len, "polygon offset fill enabled = %d\n", b);

	b = glIsEnabled(GL_POLYGON_OFFSET_LINE);
	len += sprintf(text+len, "polygon offset line enabled = %d\n", b);

	glGetIntegerv(GL_PACK_ALIGNMENT, &i);
	len += sprintf(text+len, "pack alignment = %d\n", i);

	glGetIntegerv(GL_PACK_IMAGE_HEIGHT, &i);
	len += sprintf(text+len, "pack image height = %d\n", i);

	glGetIntegerv(GL_PACK_ROW_LENGTH, &i);
	len += sprintf(text+len, "pack row length = %d\n", i);

	glGetIntegerv(GL_PACK_SKIP_PIXELS, &i);
	len += sprintf(text+len, "pack skip pixels = %d\n", i);

	glGetIntegerv(GL_PACK_SKIP_ROWS, &i);
	len += sprintf(text+len, "pack skip rows = %d\n", i);

	glGetIntegerv(GL_UNPACK_ALIGNMENT, &i);
	len += sprintf(text+len, "unpack alignment = %d\n", i);

	glGetIntegerv(GL_UNPACK_IMAGE_HEIGHT, &i);
	len += sprintf(text+len, "unpack image height = %d\n", i);

	glGetIntegerv(GL_UNPACK_ROW_LENGTH, &i);
	len += sprintf(text+len, "unpack row length = %d\n", i);

	glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &i);
	len += sprintf(text+len, "unpack skip pixels = %d\n", i);

	glGetIntegerv(GL_UNPACK_SKIP_ROWS, &i);
	len += sprintf(text+len, "unpack skip rows = %d\n", i);

	float f;
	glGetFloatv(GL_LINE_WIDTH, &f);
	len += sprintf(text+len, "line width = %f\n", f);

	if( outStream ) {
		(*outStream) << text;
	} else {
		DMSG((kEveryone, text));
	}

}

static void GL_ReportError(const char * inCallingFunc) {
	GLenum err = glGetError();
	
	static map<GLenum, string> errorMap;
	if( errorMap.empty() ) {
		errorMap[GL_INVALID_ENUM]		= "GL_INVALID_ENUM";
		errorMap[GL_INVALID_VALUE]		= "GL_INVALID_VALUE";
		errorMap[GL_INVALID_OPERATION]	= "GL_INVALID_OPERATION";
		errorMap[GL_STACK_OVERFLOW]		= "GL_STACK_OVERFLOW";
		errorMap[GL_STACK_UNDERFLOW]	= "GL_STACK_UNDERFLOW";
		errorMap[GL_OUT_OF_MEMORY]		= "GL_OUT_OF_MEMORY";
		errorMap[GL_TABLE_TOO_LARGE]	= "GL_TABLE_TOO_LARGE";
		errorMap[GL_INVALID_FRAMEBUFFER_OPERATION] = "GL_INVALID_FRAMEBUFFER_OPERATION";
	}
	
	char text[512] = { 0 };

	map<GLenum, string>::iterator it = errorMap.find(err);
	if( it != errorMap.end() ) {
		sprintf(text, "OpenGL Error: %s caused %s\n", inCallingFunc, it->second.c_str());
		DMSG((kEveryone, text));
		
	} else if( err != GL_NO_ERROR ) {
		sprintf(text, "OpenGL Error: %s caused %d\n", inCallingFunc, err);
		DMSG((kJHutchison, text));		
	}

	if( gGLDebugOut.is_open() && text[0] != 0 ) {
		map<string, bool>::iterator it = gFunctionLogability.find(inCallingFunc);
		if( it==gFunctionLogability.end() || it->second ) {
			gGLDebugOut << text << endl;
			gGLDebugOut.flush();
		}
	}

}

#define LogGLFunctionBegin(functionName) \
	if( gGLDebugOut.is_open() ) {\
		map<string, bool>::iterator it = gFunctionLogability.find(functionName);\
		if( it==gFunctionLogability.end() || it->second ) {\
			gCallsSinceLastStateDump++;\
			OutputGLStateDump();

#define LogGLFunctionEnd \
			gGLDebugOut.flush();\
		}\
	}
	

#define LogGLFunction0(functionName) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "()\n";\
	LogGLFunctionEnd

#define LogGLFunction1(functionName, p1) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction2(functionName, p1, p2) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ")\n";\
	LogGLFunctionEnd
	
#define LogGLFunction3(functionName, p1, p2, p3) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction4(functionName, p1, p2, p3, p4) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction5(functionName, p1, p2, p3, p4, p5) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction6(functionName, p1, p2, p3, p4, p5, p6) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ", " << p6 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction7(functionName, p1, p2, p3, p4, p5, p6, p7) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ", " << p6 << ", " << p7 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction8(functionName, p1, p2, p3, p4, p5, p6, p7, p8) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ", " << p6 << ", " << p7 << ", " << p8 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction9(functionName, p1, p2, p3, p4, p5, p6, p7, p8, p9) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ", " << p6 << ", " << p7 << ", " << p8 << ", " << p9 << ")\n";\
	LogGLFunctionEnd

#define LogGLFunction10(functionName, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) \
	LogGLFunctionBegin(functionName)\
	gGLDebugOut << functionName << "(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ", " << p5 << ", " << p6 << ", " << p7 << ", " << p8 << ", " << p9 << ", " << p10 << ")\n";\
	LogGLFunctionEnd

char * GL_BooleanToString(GLboolean inB) {
	static char trueText[] = "true";
	static char falseText[] = "false";
	
	if( inB ) {
		return trueText;
	} else {
		return falseText;
	}

}

char * GL_EnumToString(int inEnum) {
	static map<int, string> enumToStringMap;
	if( enumToStringMap.empty() ) {
		enumToStringMap[GL_PACK_SWAP_BYTES] = "GL_PACK_SWAP_BYTES";
		enumToStringMap[GL_PACK_LSB_FIRST] = "GL_PACK_LSB_FIRST";
		enumToStringMap[GL_PACK_ROW_LENGTH] = "GL_PACK_ROW_LENGTH";
		enumToStringMap[GL_PACK_IMAGE_HEIGHT] = "GL_PACK_IMAGE_HEIGHT";
		enumToStringMap[GL_PACK_SKIP_PIXELS] = "GL_PACK_SKIP_PIXELS";
		enumToStringMap[GL_PACK_SKIP_ROWS] = "GL_PACK_SKIP_ROWS";
		enumToStringMap[GL_PACK_SKIP_IMAGES] = "GL_PACK_SKIP_IMAGES";
		enumToStringMap[GL_PACK_ALIGNMENT] = "GL_PACK_ALIGNMENT";
		enumToStringMap[GL_UNPACK_SWAP_BYTES] = "GL_UNPACK_SWAP_BYTES";
		enumToStringMap[GL_UNPACK_LSB_FIRST] = "GL_UNPACK_LSB_FIRST";
		enumToStringMap[GL_UNPACK_ROW_LENGTH] = "GL_UNPACK_ROW_LENGTH";
		enumToStringMap[GL_UNPACK_IMAGE_HEIGHT] = "GL_UNPACK_IMAGE_HEIGHT";
		enumToStringMap[GL_UNPACK_SKIP_PIXELS] = "GL_UNPACK_SKIP_PIXELS";
		enumToStringMap[GL_UNPACK_SKIP_ROWS] = "GL_UNPACK_SKIP_ROWS";
		enumToStringMap[GL_UNPACK_SKIP_IMAGES] = "GL_UNPACK_SKIP_IMAGES";
		enumToStringMap[GL_UNPACK_ALIGNMENT] = "GL_UNPACK_ALIGNMENT";

		enumToStringMap[GL_NONE] = "GL_NONE";
		enumToStringMap[GL_FRONT_LEFT] = "GL_FRONT_LEFT";
		enumToStringMap[GL_FRONT_RIGHT] = "GL_FRONT_RIGHT";
		enumToStringMap[GL_BACK_LEFT] = "GL_BACK_LEFT";
		enumToStringMap[GL_BACK_RIGHT] = "GL_BACK_RIGHT";
		enumToStringMap[GL_FRONT] = "GL_FRONT";
		enumToStringMap[GL_BACK] = "GL_BACK";
		enumToStringMap[GL_LEFT] = "GL_LEFT";
		enumToStringMap[GL_RIGHT] = "GL_RIGHT";
		enumToStringMap[GL_FRONT_AND_BACK] = "GL_FRONT_AND_BACK";
		enumToStringMap[GL_AUX0] = "GL_AUX0";
		enumToStringMap[GL_AUX1] = "GL_AUX1";
		enumToStringMap[GL_AUX2] = "GL_AUX2";
		enumToStringMap[GL_AUX3] = "GL_AUX3";

		enumToStringMap[GL_COLOR_BUFFER_BIT] = "GL_COLOR_BUFFER_BIT";
		enumToStringMap[GL_DEPTH_BUFFER_BIT] = "GL_DEPTH_BUFFER_BIT";
		enumToStringMap[GL_ACCUM_BUFFER_BIT] = "GL_ACCUM_BUFFER_BIT";
		enumToStringMap[GL_STENCIL_BUFFER_BIT] = "GL_STENCIL_BUFFER_BIT";

		enumToStringMap[GL_FLAT] = "GL_FLAT";
		enumToStringMap[GL_SMOOTH] = "GL_SMOOTH";

		enumToStringMap[GL_MODELVIEW] = "GL_MODELVIEW";
		enumToStringMap[GL_PROJECTION] = "GL_PROJECTION";
		enumToStringMap[GL_TEXTURE] = "GL_TEXTURE";
		enumToStringMap[GL_COLOR] = "GL_COLOR";

		enumToStringMap[GL_ALPHA_TEST] = "GL_ALPHA_TEST";
		enumToStringMap[GL_AUTO_NORMAL] = "GL_AUTO_NORMAL";
		enumToStringMap[GL_BLEND] = "GL_BLEND";
		enumToStringMap[GL_CLIP_PLANE0] = "GL_CLIP_PLANE0";
		enumToStringMap[GL_CLIP_PLANE1] = "GL_CLIP_PLANE1";
		enumToStringMap[GL_CLIP_PLANE2] = "GL_CLIP_PLANE2";
		enumToStringMap[GL_CLIP_PLANE3] = "GL_CLIP_PLANE3";
		enumToStringMap[GL_CLIP_PLANE4] = "GL_CLIP_PLANE4";
		enumToStringMap[GL_CLIP_PLANE5] = "GL_CLIP_PLANE5";
		enumToStringMap[GL_COLOR_LOGIC_OP] = "GL_COLOR_LOGIC_OP";
		enumToStringMap[GL_COLOR_MATERIAL] = "GL_COLOR_MATERIAL";
		enumToStringMap[GL_COLOR_SUM] = "GL_COLOR_SUM";
		enumToStringMap[GL_COLOR_TABLE] = "GL_COLOR_TABLE";
		enumToStringMap[GL_CONVOLUTION_1D] = "GL_CONVOLUTION_1D";
		enumToStringMap[GL_CONVOLUTION_2D] = "GL_CONVOLUTION_2D";
		enumToStringMap[GL_CULL_FACE] = "GL_CULL_FACE";
		enumToStringMap[GL_DEPTH_TEST] = "GL_DEPTH_TEST";
		enumToStringMap[GL_DITHER] = "GL_DITHER";
		enumToStringMap[GL_FOG] = "GL_FOG";
		enumToStringMap[GL_HISTOGRAM] = "GL_HISTOGRAM";
		enumToStringMap[GL_INDEX_LOGIC_OP] = "GL_INDEX_LOGIC_OP";
		enumToStringMap[GL_LIGHT0] = "GL_LIGHT0";
		enumToStringMap[GL_LIGHT1] = "GL_LIGHT1";
		enumToStringMap[GL_LIGHT2] = "GL_LIGHT2";
		enumToStringMap[GL_LIGHT3] = "GL_LIGHT3";
		enumToStringMap[GL_LIGHT4] = "GL_LIGHT4";
		enumToStringMap[GL_LIGHT5] = "GL_LIGHT5";
		enumToStringMap[GL_LIGHT6] = "GL_LIGHT6";
		enumToStringMap[GL_LIGHT7] = "GL_LIGHT7";
		enumToStringMap[GL_LIGHTING] = "GL_LIGHTING";
		enumToStringMap[GL_LINE_SMOOTH] = "GL_LINE_SMOOTH";
		enumToStringMap[GL_LINE_STIPPLE] = "GL_LINE_STIPPLE";
		enumToStringMap[GL_MAP1_COLOR_4] = "GL_MAP1_COLOR_4";
		enumToStringMap[GL_MAP1_INDEX] = "GL_MAP1_INDEX";
		enumToStringMap[GL_MAP1_NORMAL] = "GL_MAP1_NORMAL";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_1] = "GL_MAP1_TEXTURE_COORD_1";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_2] = "GL_MAP1_TEXTURE_COORD_2";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_3] = "GL_MAP1_TEXTURE_COORD_3";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_4] = "GL_MAP1_TEXTURE_COORD_4";
		enumToStringMap[GL_MAP1_VERTEX_3] = "GL_MAP1_VERTEX_3";
		enumToStringMap[GL_MAP1_VERTEX_4] = "GL_MAP1_VERTEX_4";
		enumToStringMap[GL_MAP2_COLOR_4] = "GL_MAP2_COLOR_4";
		enumToStringMap[GL_MAP2_INDEX] = "GL_MAP2_INDEX";
		enumToStringMap[GL_MAP2_NORMAL] = "GL_MAP2_NORMAL";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_1] = "GL_MAP2_TEXTURE_COORD_1";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_2] = "GL_MAP2_TEXTURE_COORD_2";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_3] = "GL_MAP2_TEXTURE_COORD_3";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_4] = "GL_MAP2_TEXTURE_COORD_4";
		enumToStringMap[GL_MAP2_VERTEX_3] = "GL_MAP2_VERTEX_3";
		enumToStringMap[GL_MAP2_VERTEX_4] = "GL_MAP2_VERTEX_4";
		enumToStringMap[GL_MINMAX] = "GL_MINMAX";
		enumToStringMap[GL_MULTISAMPLE] = "GL_MULTISAMPLE";
		enumToStringMap[GL_NORMALIZE] = "GL_NORMALIZE";
		enumToStringMap[GL_POINT_SMOOTH] = "GL_POINT_SMOOTH";
		enumToStringMap[GL_POINT_SPRITE] = "GL_POINT_SPRITE";
		enumToStringMap[GL_POLYGON_OFFSET_FILL] = "GL_POLYGON_OFFSET_FILL";
		enumToStringMap[GL_POLYGON_OFFSET_LINE] = "GL_POLYGON_OFFSET_LINE";
		enumToStringMap[GL_POLYGON_OFFSET_POINT] = "GL_POLYGON_OFFSET_POINT";
		enumToStringMap[GL_POLYGON_SMOOTH] = "GL_POLYGON_SMOOTH";
		enumToStringMap[GL_POLYGON_STIPPLE] = "GL_POLYGON_STIPPLE";
		enumToStringMap[GL_POST_COLOR_MATRIX_COLOR_TABLE] = "GL_POST_COLOR_MATRIX_COLOR_TABLE";
		enumToStringMap[GL_POST_CONVOLUTION_COLOR_TABLE] = "GL_POST_CONVOLUTION_COLOR_TABLE";
		enumToStringMap[GL_RESCALE_NORMAL] = "GL_RESCALE_NORMAL";
		enumToStringMap[GL_SAMPLE_ALPHA_TO_COVERAGE] = "GL_SAMPLE_ALPHA_TO_COVERAGE";
		enumToStringMap[GL_SAMPLE_ALPHA_TO_ONE] = "GL_SAMPLE_ALPHA_TO_ONE";
		enumToStringMap[GL_SAMPLE_COVERAGE] = "GL_SAMPLE_COVERAGE";
		enumToStringMap[GL_SEPARABLE_2D] = "GL_SEPARABLE_2D";
		enumToStringMap[GL_SCISSOR_TEST] = "GL_SCISSOR_TEST";
		enumToStringMap[GL_STENCIL_TEST] = "GL_STENCIL_TEST";
		enumToStringMap[GL_TEXTURE_1D] = "GL_TEXTURE_1D";
		enumToStringMap[GL_TEXTURE_2D] = "GL_TEXTURE_2D";
		enumToStringMap[GL_TEXTURE_3D] = "GL_TEXTURE_3D";
		enumToStringMap[GL_TEXTURE_CUBE_MAP] = "GL_TEXTURE_CUBE_MAP";
		enumToStringMap[GL_TEXTURE_GEN_Q] = "GL_TEXTURE_GEN_Q";
		enumToStringMap[GL_TEXTURE_GEN_R] = "GL_TEXTURE_GEN_R";
		enumToStringMap[GL_TEXTURE_GEN_S] = "GL_TEXTURE_GEN_S";
		enumToStringMap[GL_TEXTURE_GEN_T] = "GL_TEXTURE_GEN_T";
		enumToStringMap[GL_VERTEX_PROGRAM_POINT_SIZE] = "GL_VERTEX_PROGRAM_POINT_SIZE";
		enumToStringMap[GL_VERTEX_PROGRAM_TWO_SIDE] = "GL_VERTEX_PROGRAM_TWO_SIDE";

		enumToStringMap[GL_ARRAY_BUFFER] = "GL_ARRAY_BUFFER";
		enumToStringMap[GL_ELEMENT_ARRAY_BUFFER] = "GL_ELEMENT_ARRAY_BUFFER";
		enumToStringMap[GL_PIXEL_PACK_BUFFER] = "GL_PIXEL_PACK_BUFFER";
		enumToStringMap[GL_PIXEL_UNPACK_BUFFER] = "GL_PIXEL_UNPACK_BUFFER";

		enumToStringMap[GL_DRAW_FRAMEBUFFER] = "GL_DRAW_FRAMEBUFFER";
		enumToStringMap[GL_READ_FRAMEBUFFER] = "GL_READ_FRAMEBUFFER";
		enumToStringMap[GL_FRAMEBUFFER] = "GL_FRAMEBUFFER";

		enumToStringMap[GL_RENDERBUFFER] = "GL_RENDERBUFFER";

		enumToStringMap[GL_COLOR_INDEX] = "GL_COLOR_INDEX";
		enumToStringMap[GL_STENCIL_INDEX] = "GL_STENCIL_INDEX";
		enumToStringMap[GL_DEPTH_COMPONENT] = "GL_DEPTH_COMPONENT";
		enumToStringMap[GL_RED] = "GL_RED";
		enumToStringMap[GL_GREEN] = "GL_GREEN";
		enumToStringMap[GL_BLUE] = "GL_BLUE";
		enumToStringMap[GL_ALPHA] = "GL_ALPHA";
		enumToStringMap[GL_RGB] = "GL_RGB";
		enumToStringMap[GL_BGR] = "GL_BGR";
		enumToStringMap[GL_RGBA] = "GL_RGBA";
		enumToStringMap[GL_BGRA] = "GL_BGRA";
		enumToStringMap[GL_LUMINANCE_ALPHA] = "GL_LUMINANCE_ALPHA";

		enumToStringMap[GL_UNSIGNED_BYTE] = "GL_UNSIGNED_BYTE";
		enumToStringMap[GL_BYTE] = "GL_BYTE";
		enumToStringMap[GL_BITMAP] = "GL_BITMAP";
		enumToStringMap[GL_UNSIGNED_SHORT] = "GL_UNSIGNED_SHORT";
		enumToStringMap[GL_SHORT] = "GL_SHORT";
		enumToStringMap[GL_UNSIGNED_INT] = "GL_UNSIGNED_INT";
		enumToStringMap[GL_INT] = "GL_INT";
		enumToStringMap[GL_FLOAT] = "GL_FLOAT";
		enumToStringMap[GL_UNSIGNED_BYTE_3_3_2] = "GL_UNSIGNED_BYTE_3_3_2";
		enumToStringMap[GL_UNSIGNED_BYTE_2_3_3_REV] = "GL_UNSIGNED_BYTE_2_3_3_REV";
		enumToStringMap[GL_UNSIGNED_SHORT_5_6_5] = "GL_UNSIGNED_SHORT_5_6_5";
		enumToStringMap[GL_UNSIGNED_SHORT_5_6_5_REV] = "GL_UNSIGNED_SHORT_5_6_5_REV";
		enumToStringMap[GL_UNSIGNED_SHORT_4_4_4_4] = "GL_UNSIGNED_SHORT_4_4_4_4";
		enumToStringMap[GL_UNSIGNED_SHORT_4_4_4_4_REV] = "GL_UNSIGNED_SHORT_4_4_4_4_REV";
		enumToStringMap[GL_UNSIGNED_SHORT_5_5_5_1] = "GL_UNSIGNED_SHORT_5_5_5_1";
		enumToStringMap[GL_UNSIGNED_SHORT_1_5_5_5_REV] = "GL_UNSIGNED_SHORT_1_5_5_5_REV";
		enumToStringMap[GL_UNSIGNED_INT_8_8_8_8] = "GL_UNSIGNED_INT_8_8_8_8";
		enumToStringMap[GL_UNSIGNED_INT_8_8_8_8_REV] = "GL_UNSIGNED_INT_8_8_8_8_REV";
		enumToStringMap[GL_UNSIGNED_INT_10_10_10_2] = "GL_UNSIGNED_INT_10_10_10_2";
		enumToStringMap[GL_UNSIGNED_INT_2_10_10_10_REV] = "GL_UNSIGNED_INT_2_10_10_10_REV";

		enumToStringMap[GL_COLOR_ATTACHMENT0] = "GL_COLOR_ATTACHMENT0";
		enumToStringMap[GL_COLOR_ATTACHMENT1] = "GL_COLOR_ATTACHMENT1";
		enumToStringMap[GL_COLOR_ATTACHMENT2] = "GL_COLOR_ATTACHMENT2";
		enumToStringMap[GL_COLOR_ATTACHMENT3] = "GL_COLOR_ATTACHMENT3";
		enumToStringMap[GL_COLOR_ATTACHMENT4] = "GL_COLOR_ATTACHMENT4";
		enumToStringMap[GL_COLOR_ATTACHMENT5] = "GL_COLOR_ATTACHMENT5";
		enumToStringMap[GL_COLOR_ATTACHMENT6] = "GL_COLOR_ATTACHMENT6";
		enumToStringMap[GL_COLOR_ATTACHMENT7] = "GL_COLOR_ATTACHMENT7";
		enumToStringMap[GL_COLOR_ATTACHMENT8] = "GL_COLOR_ATTACHMENT8";
		enumToStringMap[GL_COLOR_ATTACHMENT9] = "GL_COLOR_ATTACHMENT9";
		enumToStringMap[GL_COLOR_ATTACHMENT10] = "GL_COLOR_ATTACHMENT10";
		enumToStringMap[GL_COLOR_ATTACHMENT11] = "GL_COLOR_ATTACHMENT11";
		enumToStringMap[GL_COLOR_ATTACHMENT12] = "GL_COLOR_ATTACHMENT12";
		enumToStringMap[GL_COLOR_ATTACHMENT13] = "GL_COLOR_ATTACHMENT13";
		enumToStringMap[GL_COLOR_ATTACHMENT14] = "GL_COLOR_ATTACHMENT14";
		enumToStringMap[GL_COLOR_ATTACHMENT15] = "GL_COLOR_ATTACHMENT15";
		enumToStringMap[GL_DEPTH_ATTACHMENT] = "GL_DEPTH_ATTACHMENT";
		enumToStringMap[GL_STENCIL_ATTACHMENT] = "GL_STENCIL_ATTACHMENT";
		enumToStringMap[GL_DEPTH_STENCIL_ATTACHMENT] = "GL_DEPTH_STENCIL_ATTACHMENT";

		enumToStringMap[GL_TEXTURE_1D] = "GL_TEXTURE_1D";
		enumToStringMap[GL_TEXTURE_2D] = "GL_TEXTURE_2D";
		enumToStringMap[GL_TEXTURE_3D] = "GL_TEXTURE_3D";
		enumToStringMap[GL_TEXTURE_CUBE_MAP] = "GL_TEXTURE_CUBE_MAP";

		enumToStringMap[GL_TEXTURE_1D] = "GL_TEXTURE_1D";
		enumToStringMap[GL_TEXTURE_2D] = "GL_TEXTURE_2D";
		enumToStringMap[GL_TEXTURE_3D] = "GL_TEXTURE_3D";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_POSITIVE_X] = "GL_TEXTURE_CUBE_MAP_POSITIVE_X";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_NEGATIVE_X] = "GL_TEXTURE_CUBE_MAP_NEGATIVE_X";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_POSITIVE_Y] = "GL_TEXTURE_CUBE_MAP_POSITIVE_Y";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_NEGATIVE_Y] = "GL_TEXTURE_CUBE_MAP_NEGATIVE_Y";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_POSITIVE_Z] = "GL_TEXTURE_CUBE_MAP_POSITIVE_Z";
		enumToStringMap[GL_TEXTURE_CUBE_MAP_NEGATIVE_Z] = "GL_TEXTURE_CUBE_MAP_NEGATIVE_Z";

		enumToStringMap[GL_TEXTURE0] = "GL_TEXTURE0";
		enumToStringMap[GL_TEXTURE1] = "GL_TEXTURE1";
		enumToStringMap[GL_TEXTURE2] = "GL_TEXTURE2";
		enumToStringMap[GL_TEXTURE3] = "GL_TEXTURE3";
		enumToStringMap[GL_TEXTURE4] = "GL_TEXTURE4";
		enumToStringMap[GL_TEXTURE5] = "GL_TEXTURE5";
		enumToStringMap[GL_TEXTURE6] = "GL_TEXTURE6";
		enumToStringMap[GL_TEXTURE7] = "GL_TEXTURE7";
		enumToStringMap[GL_TEXTURE8] = "GL_TEXTURE8";
		enumToStringMap[GL_TEXTURE9] = "GL_TEXTURE9";
		enumToStringMap[GL_TEXTURE10] = "GL_TEXTURE10";
		enumToStringMap[GL_TEXTURE11] = "GL_TEXTURE11";
		enumToStringMap[GL_TEXTURE12] = "GL_TEXTURE12";
		enumToStringMap[GL_TEXTURE13] = "GL_TEXTURE13";
		enumToStringMap[GL_TEXTURE14] = "GL_TEXTURE14";
		enumToStringMap[GL_TEXTURE15] = "GL_TEXTURE15";
		enumToStringMap[GL_TEXTURE16] = "GL_TEXTURE16";
		enumToStringMap[GL_TEXTURE17] = "GL_TEXTURE17";
		enumToStringMap[GL_TEXTURE18] = "GL_TEXTURE18";
		enumToStringMap[GL_TEXTURE19] = "GL_TEXTURE19";
		enumToStringMap[GL_TEXTURE20] = "GL_TEXTURE20";
		enumToStringMap[GL_TEXTURE21] = "GL_TEXTURE21";
		enumToStringMap[GL_TEXTURE22] = "GL_TEXTURE22";
		enumToStringMap[GL_TEXTURE23] = "GL_TEXTURE23";
		enumToStringMap[GL_TEXTURE24] = "GL_TEXTURE24";
		enumToStringMap[GL_TEXTURE25] = "GL_TEXTURE25";
		enumToStringMap[GL_TEXTURE26] = "GL_TEXTURE26";
		enumToStringMap[GL_TEXTURE27] = "GL_TEXTURE27";
		enumToStringMap[GL_TEXTURE28] = "GL_TEXTURE28";
		enumToStringMap[GL_TEXTURE29] = "GL_TEXTURE29";
		enumToStringMap[GL_TEXTURE30] = "GL_TEXTURE30";
		enumToStringMap[GL_TEXTURE31] = "GL_TEXTURE31";

		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE] = "GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_RED_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_RED_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_GREEN_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_GREEN_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_BLUE_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_BLUE_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_ALPHA_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_ALPHA_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE] = "GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_COMPONENT_TYPE] = "GL_FRAMEBUFFER_ATTACHMENT_COMPONENT_TYPE";
		enumToStringMap[GL_FLOAT] = "GL_FLOAT";
		enumToStringMap[GL_INT] = "GL_INT";
		enumToStringMap[GL_UNSIGNED_INT] = "GL_UNSIGNED_INT";
		enumToStringMap[GL_SIGNED_NORMALIZED] = "GL_SIGNED_NORMALIZED";
		enumToStringMap[GL_UNSIGNED_NORMALIZED] = "GL_UNSIGNED_NORMALIZED";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING] = "GL_FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING";
		enumToStringMap[GL_LINEAR] = "GL_LINEAR";
		enumToStringMap[GL_SRGB] = "GL_SRGB";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE] = "GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE";
		enumToStringMap[GL_RENDERBUFFER] = "GL_RENDERBUFFER";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME] = "GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME] = "GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL] = "GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE] = "GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LAYER] = "GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LAYER";
		enumToStringMap[GL_FRAMEBUFFER_ATTACHMENT_LAYERED] = "GL_FRAMEBUFFER_ATTACHMENT_LAYERED";

		enumToStringMap[GL_ACCUM_ALPHA_BITS] = "GL_ACCUM_ALPHA_BITS";
		enumToStringMap[GL_ACCUM_BLUE_BITS] = "GL_ACCUM_BLUE_BITS";
		enumToStringMap[GL_ACCUM_CLEAR_VALUE] = "GL_ACCUM_CLEAR_VALUE";
		enumToStringMap[GL_ACCUM_GREEN_BITS] = "GL_ACCUM_GREEN_BITS";
		enumToStringMap[GL_ACCUM_RED_BITS] = "GL_ACCUM_RED_BITS";
		enumToStringMap[GL_ACTIVE_TEXTURE] = "GL_ACTIVE_TEXTURE";
		enumToStringMap[GL_ALIASED_POINT_SIZE_RANGE] = "GL_ALIASED_POINT_SIZE_RANGE";
		enumToStringMap[GL_ALIASED_LINE_WIDTH_RANGE] = "GL_ALIASED_LINE_WIDTH_RANGE";
		enumToStringMap[GL_ALPHA_BIAS] = "GL_ALPHA_BIAS";
		enumToStringMap[GL_ALPHA_BITS] = "GL_ALPHA_BITS";
		enumToStringMap[GL_ALPHA_SCALE] = "GL_ALPHA_SCALE";
		enumToStringMap[GL_ALPHA_TEST] = "GL_ALPHA_TEST";
		enumToStringMap[GL_ALPHA_TEST_FUNC] = "GL_ALPHA_TEST_FUNC";
		enumToStringMap[GL_ALPHA_TEST_REF] = "GL_ALPHA_TEST_REF";
		enumToStringMap[GL_ARRAY_BUFFER_BINDING] = "GL_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_ATTRIB_STACK_DEPTH] = "GL_ATTRIB_STACK_DEPTH";
		enumToStringMap[GL_AUTO_NORMAL] = "GL_AUTO_NORMAL";
		enumToStringMap[GL_AUX_BUFFERS] = "GL_AUX_BUFFERS";
		enumToStringMap[GL_BLEND] = "GL_BLEND";
		enumToStringMap[GL_BLEND_COLOR] = "GL_BLEND_COLOR";
		enumToStringMap[GL_BLEND_DST_ALPHA] = "GL_BLEND_DST_ALPHA";
		enumToStringMap[GL_BLEND_DST_RGB] = "GL_BLEND_DST_RGB";
		enumToStringMap[GL_BLEND_EQUATION_RGB] = "GL_BLEND_EQUATION_RGB";
		enumToStringMap[GL_BLEND_EQUATION_ALPHA] = "GL_BLEND_EQUATION_ALPHA";
		enumToStringMap[GL_BLEND_SRC_ALPHA] = "GL_BLEND_SRC_ALPHA";
		enumToStringMap[GL_BLEND_SRC_RGB] = "GL_BLEND_SRC_RGB";
		enumToStringMap[GL_BLUE_BIAS] = "GL_BLUE_BIAS";
		enumToStringMap[GL_BLUE_BITS] = "GL_BLUE_BITS";
		enumToStringMap[GL_BLUE_SCALE] = "GL_BLUE_SCALE";
		enumToStringMap[GL_CLIENT_ACTIVE_TEXTURE] = "GL_CLIENT_ACTIVE_TEXTURE";
		enumToStringMap[GL_CLIENT_ATTRIB_STACK_DEPTH] = "GL_CLIENT_ATTRIB_STACK_DEPTH";
		enumToStringMap[GL_COLOR_ARRAY] = "GL_COLOR_ARRAY";
		enumToStringMap[GL_COLOR_ARRAY_BUFFER_BINDING] = "GL_COLOR_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_COLOR_ARRAY_SIZE] = "GL_COLOR_ARRAY_SIZE";
		enumToStringMap[GL_COLOR_ARRAY_STRIDE] = "GL_COLOR_ARRAY_STRIDE";
		enumToStringMap[GL_COLOR_ARRAY_TYPE] = "GL_COLOR_ARRAY_TYPE";
		enumToStringMap[GL_COLOR_CLEAR_VALUE] = "GL_COLOR_CLEAR_VALUE";
		enumToStringMap[GL_COLOR_LOGIC_OP] = "GL_COLOR_LOGIC_OP";
		enumToStringMap[GL_COLOR_MATERIAL] = "GL_COLOR_MATERIAL";
		enumToStringMap[GL_COLOR_MATERIAL_FACE] = "GL_COLOR_MATERIAL_FACE";
		enumToStringMap[GL_COLOR_MATERIAL_PARAMETER] = "GL_COLOR_MATERIAL_PARAMETER";
		enumToStringMap[GL_COLOR_MATRIX] = "GL_COLOR_MATRIX";
		enumToStringMap[GL_COLOR_MATRIX_STACK_DEPTH] = "GL_COLOR_MATRIX_STACK_DEPTH";
		enumToStringMap[GL_COLOR_SUM] = "GL_COLOR_SUM";
		enumToStringMap[GL_COLOR_TABLE] = "GL_COLOR_TABLE";
		enumToStringMap[GL_COLOR_WRITEMASK] = "GL_COLOR_WRITEMASK";
		enumToStringMap[GL_COMPRESSED_TEXTURE_FORMATS] = "GL_COMPRESSED_TEXTURE_FORMATS";
		enumToStringMap[GL_CONVOLUTION_1D] = "GL_CONVOLUTION_1D";
		enumToStringMap[GL_CONVOLUTION_2D] = "GL_CONVOLUTION_2D";
		enumToStringMap[GL_CULL_FACE] = "GL_CULL_FACE";
		enumToStringMap[GL_CULL_FACE_MODE] = "GL_CULL_FACE_MODE";
		enumToStringMap[GL_CURRENT_COLOR] = "GL_CURRENT_COLOR";
		enumToStringMap[GL_CURRENT_FOG_COORD] = "GL_CURRENT_FOG_COORD";
		enumToStringMap[GL_CURRENT_INDEX] = "GL_CURRENT_INDEX";
		enumToStringMap[GL_CURRENT_NORMAL] = "GL_CURRENT_NORMAL";
		enumToStringMap[GL_CURRENT_PROGRAM] = "GL_CURRENT_PROGRAM";
		enumToStringMap[GL_CURRENT_RASTER_COLOR] = "GL_CURRENT_RASTER_COLOR";
		enumToStringMap[GL_CURRENT_RASTER_DISTANCE] = "GL_CURRENT_RASTER_DISTANCE";
		enumToStringMap[GL_CURRENT_RASTER_INDEX] = "GL_CURRENT_RASTER_INDEX";
		enumToStringMap[GL_CURRENT_RASTER_POSITION] = "GL_CURRENT_RASTER_POSITION";
		enumToStringMap[GL_CURRENT_RASTER_POSITION_VALID] = "GL_CURRENT_RASTER_POSITION_VALID";
		enumToStringMap[GL_CURRENT_RASTER_SECONDARY_COLOR] = "GL_CURRENT_RASTER_SECONDARY_COLOR";
		enumToStringMap[GL_CURRENT_RASTER_TEXTURE_COORDS] = "GL_CURRENT_RASTER_TEXTURE_COORDS";
		enumToStringMap[GL_CURRENT_SECONDARY_COLOR] = "GL_CURRENT_SECONDARY_COLOR";
		enumToStringMap[GL_CURRENT_TEXTURE_COORDS] = "GL_CURRENT_TEXTURE_COORDS";
		enumToStringMap[GL_DEPTH_BIAS] = "GL_DEPTH_BIAS";
		enumToStringMap[GL_DEPTH_BITS] = "GL_DEPTH_BITS";
		enumToStringMap[GL_DEPTH_CLEAR_VALUE] = "GL_DEPTH_CLEAR_VALUE";
		enumToStringMap[GL_DEPTH_FUNC] = "GL_DEPTH_FUNC";
		enumToStringMap[GL_DEPTH_RANGE] = "GL_DEPTH_RANGE";
		enumToStringMap[GL_DEPTH_SCALE] = "GL_DEPTH_SCALE";
		enumToStringMap[GL_DEPTH_TEST] = "GL_DEPTH_TEST";
		enumToStringMap[GL_DEPTH_WRITEMASK] = "GL_DEPTH_WRITEMASK";
		enumToStringMap[GL_DITHER] = "GL_DITHER";
		enumToStringMap[GL_DOUBLEBUFFER] = "GL_DOUBLEBUFFER";
		enumToStringMap[GL_DRAW_BUFFER] = "GL_DRAW_BUFFER";
		enumToStringMap[GL_DRAW_BUFFER0] = "GL_DRAW_BUFFER0";
		enumToStringMap[GL_DRAW_BUFFER1] = "GL_DRAW_BUFFER1";
		enumToStringMap[GL_DRAW_BUFFER2] = "GL_DRAW_BUFFER2";
		enumToStringMap[GL_DRAW_BUFFER3] = "GL_DRAW_BUFFER3";
		enumToStringMap[GL_DRAW_BUFFER4] = "GL_DRAW_BUFFER4";
		enumToStringMap[GL_DRAW_BUFFER5] = "GL_DRAW_BUFFER5";
		enumToStringMap[GL_DRAW_BUFFER6] = "GL_DRAW_BUFFER6";
		enumToStringMap[GL_DRAW_BUFFER7] = "GL_DRAW_BUFFER7";
		enumToStringMap[GL_DRAW_BUFFER8] = "GL_DRAW_BUFFER8";
		enumToStringMap[GL_DRAW_BUFFER9] = "GL_DRAW_BUFFER9";
		enumToStringMap[GL_DRAW_BUFFER10] = "GL_DRAW_BUFFER10";
		enumToStringMap[GL_DRAW_BUFFER11] = "GL_DRAW_BUFFER11";
		enumToStringMap[GL_DRAW_BUFFER12] = "GL_DRAW_BUFFER12";
		enumToStringMap[GL_DRAW_BUFFER13] = "GL_DRAW_BUFFER13";
		enumToStringMap[GL_DRAW_BUFFER14] = "GL_DRAW_BUFFER14";
		enumToStringMap[GL_DRAW_BUFFER15] = "GL_DRAW_BUFFER15";
		enumToStringMap[GL_EDGE_FLAG] = "GL_EDGE_FLAG";
		enumToStringMap[GL_EDGE_FLAG_ARRAY] = "GL_EDGE_FLAG_ARRAY";
		enumToStringMap[GL_EDGE_FLAG_ARRAY_BUFFER_BINDING] = "GL_EDGE_FLAG_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_EDGE_FLAG_ARRAY_STRIDE] = "GL_EDGE_FLAG_ARRAY_STRIDE";
		enumToStringMap[GL_ELEMENT_ARRAY_BUFFER_BINDING] = "GL_ELEMENT_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_FEEDBACK_BUFFER_SIZE] = "GL_FEEDBACK_BUFFER_SIZE";
		enumToStringMap[GL_FEEDBACK_BUFFER_TYPE] = "GL_FEEDBACK_BUFFER_TYPE";
		enumToStringMap[GL_FOG] = "GL_FOG";
		enumToStringMap[GL_FOG_COORD_ARRAY] = "GL_FOG_COORD_ARRAY";
		enumToStringMap[GL_FOG_COORD_ARRAY_BUFFER_BINDING] = "GL_FOG_COORD_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_FOG_COORD_ARRAY_STRIDE] = "GL_FOG_COORD_ARRAY_STRIDE";
		enumToStringMap[GL_FOG_COORD_ARRAY_TYPE] = "GL_FOG_COORD_ARRAY_TYPE";
		enumToStringMap[GL_FOG_COORD_SRC] = "GL_FOG_COORD_SRC";
		enumToStringMap[GL_FOG_COLOR] = "GL_FOG_COLOR";
		enumToStringMap[GL_FOG_DENSITY] = "GL_FOG_DENSITY";
		enumToStringMap[GL_FOG_END] = "GL_FOG_END";
		enumToStringMap[GL_FOG_HINT] = "GL_FOG_HINT";
		enumToStringMap[GL_FOG_INDEX] = "GL_FOG_INDEX";
		enumToStringMap[GL_FOG_MODE] = "GL_FOG_MODE";
		enumToStringMap[GL_FOG_START] = "GL_FOG_START";
		enumToStringMap[GL_FRAGMENT_SHADER_DERIVATIVE_HINT] = "GL_FRAGMENT_SHADER_DERIVATIVE_HINT";
		enumToStringMap[GL_FRONT_FACE] = "GL_FRONT_FACE";
		enumToStringMap[GL_GENERATE_MIPMAP_HINT] = "GL_GENERATE_MIPMAP_HINT";
		enumToStringMap[GL_GREEN_BIAS] = "GL_GREEN_BIAS";
		enumToStringMap[GL_GREEN_BITS] = "GL_GREEN_BITS";
		enumToStringMap[GL_GREEN_SCALE] = "GL_GREEN_SCALE";
		enumToStringMap[GL_HISTOGRAM] = "GL_HISTOGRAM";
		enumToStringMap[GL_INDEX_ARRAY] = "GL_INDEX_ARRAY";
		enumToStringMap[GL_INDEX_ARRAY_BUFFER_BINDING] = "GL_INDEX_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_INDEX_ARRAY_STRIDE] = "GL_INDEX_ARRAY_STRIDE";
		enumToStringMap[GL_INDEX_ARRAY_TYPE] = "GL_INDEX_ARRAY_TYPE";
		enumToStringMap[GL_INDEX_BITS] = "GL_INDEX_BITS";
		enumToStringMap[GL_INDEX_CLEAR_VALUE] = "GL_INDEX_CLEAR_VALUE";
		enumToStringMap[GL_INDEX_LOGIC_OP] = "GL_INDEX_LOGIC_OP";
		enumToStringMap[GL_INDEX_MODE] = "GL_INDEX_MODE";
		enumToStringMap[GL_INDEX_OFFSET] = "GL_INDEX_OFFSET";
		enumToStringMap[GL_INDEX_SHIFT] = "GL_INDEX_SHIFT";
		enumToStringMap[GL_INDEX_WRITEMASK] = "GL_INDEX_WRITEMASK";
		enumToStringMap[GL_LIGHTING] = "GL_LIGHTING";
		enumToStringMap[GL_LIGHT_MODEL_AMBIENT] = "GL_LIGHT_MODEL_AMBIENT";
		enumToStringMap[GL_LIGHT_MODEL_COLOR_CONTROL] = "GL_LIGHT_MODEL_COLOR_CONTROL";
		enumToStringMap[GL_LIGHT_MODEL_LOCAL_VIEWER] = "GL_LIGHT_MODEL_LOCAL_VIEWER";
		enumToStringMap[GL_LIGHT_MODEL_TWO_SIDE] = "GL_LIGHT_MODEL_TWO_SIDE";
		enumToStringMap[GL_LINE_SMOOTH] = "GL_LINE_SMOOTH";
		enumToStringMap[GL_LINE_SMOOTH_HINT] = "GL_LINE_SMOOTH_HINT";
		enumToStringMap[GL_LINE_STIPPLE] = "GL_LINE_STIPPLE";
		enumToStringMap[GL_LINE_STIPPLE_PATTERN] = "GL_LINE_STIPPLE_PATTERN";
		enumToStringMap[GL_LINE_STIPPLE_REPEAT] = "GL_LINE_STIPPLE_REPEAT";
		enumToStringMap[GL_LINE_WIDTH] = "GL_LINE_WIDTH";
		enumToStringMap[GL_LINE_WIDTH_GRANULARITY] = "GL_LINE_WIDTH_GRANULARITY";
		enumToStringMap[GL_LINE_WIDTH_RANGE] = "GL_LINE_WIDTH_RANGE";
		enumToStringMap[GL_LIST_BASE] = "GL_LIST_BASE";
		enumToStringMap[GL_LIST_INDEX] = "GL_LIST_INDEX";
		enumToStringMap[GL_LIST_MODE] = "GL_LIST_MODE";
		enumToStringMap[GL_LOGIC_OP_MODE] = "GL_LOGIC_OP_MODE";
		enumToStringMap[GL_MAP1_COLOR_4] = "GL_MAP1_COLOR_4";
		enumToStringMap[GL_MAP1_GRID_DOMAIN] = "GL_MAP1_GRID_DOMAIN";
		enumToStringMap[GL_MAP1_GRID_SEGMENTS] = "GL_MAP1_GRID_SEGMENTS";
		enumToStringMap[GL_MAP1_INDEX] = "GL_MAP1_INDEX";
		enumToStringMap[GL_MAP1_NORMAL] = "GL_MAP1_NORMAL";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_1] = "GL_MAP1_TEXTURE_COORD_1";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_2] = "GL_MAP1_TEXTURE_COORD_2";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_3] = "GL_MAP1_TEXTURE_COORD_3";
		enumToStringMap[GL_MAP1_TEXTURE_COORD_4] = "GL_MAP1_TEXTURE_COORD_4";
		enumToStringMap[GL_MAP1_VERTEX_3] = "GL_MAP1_VERTEX_3";
		enumToStringMap[GL_MAP1_VERTEX_4] = "GL_MAP1_VERTEX_4";
		enumToStringMap[GL_MAP2_COLOR_4] = "GL_MAP2_COLOR_4";
		enumToStringMap[GL_MAP2_GRID_DOMAIN] = "GL_MAP2_GRID_DOMAIN";
		enumToStringMap[GL_MAP2_GRID_SEGMENTS] = "GL_MAP2_GRID_SEGMENTS";
		enumToStringMap[GL_MAP2_INDEX] = "GL_MAP2_INDEX";
		enumToStringMap[GL_MAP2_NORMAL] = "GL_MAP2_NORMAL";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_1] = "GL_MAP2_TEXTURE_COORD_1";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_2] = "GL_MAP2_TEXTURE_COORD_2";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_3] = "GL_MAP2_TEXTURE_COORD_3";
		enumToStringMap[GL_MAP2_TEXTURE_COORD_4] = "GL_MAP2_TEXTURE_COORD_4";
		enumToStringMap[GL_MAP2_VERTEX_3] = "GL_MAP2_VERTEX_3";
		enumToStringMap[GL_MAP2_VERTEX_4] = "GL_MAP2_VERTEX_4";
		enumToStringMap[GL_MAP_COLOR] = "GL_MAP_COLOR";
		enumToStringMap[GL_MAP_STENCIL] = "GL_MAP_STENCIL";
		enumToStringMap[GL_MATRIX_MODE] = "GL_MATRIX_MODE";
		enumToStringMap[GL_MAX_3D_TEXTURE_SIZE] = "GL_MAX_3D_TEXTURE_SIZE";
		enumToStringMap[GL_MAX_CLIENT_ATTRIB_STACK_DEPTH] = "GL_MAX_CLIENT_ATTRIB_STACK_DEPTH";
		enumToStringMap[GL_MAX_ATTRIB_STACK_DEPTH] = "GL_MAX_ATTRIB_STACK_DEPTH";
		enumToStringMap[GL_MAX_CLIP_PLANES] = "GL_MAX_CLIP_PLANES";
		enumToStringMap[GL_MAX_COLOR_MATRIX_STACK_DEPTH] = "GL_MAX_COLOR_MATRIX_STACK_DEPTH";
		enumToStringMap[GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS] = "GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS";
		enumToStringMap[GL_MAX_CUBE_MAP_TEXTURE_SIZE] = "GL_MAX_CUBE_MAP_TEXTURE_SIZE";
		enumToStringMap[GL_MAX_DRAW_BUFFERS] = "GL_MAX_DRAW_BUFFERS";
		enumToStringMap[GL_MAX_ELEMENTS_INDICES] = "GL_MAX_ELEMENTS_INDICES";
		enumToStringMap[GL_MAX_ELEMENTS_VERTICES] = "GL_MAX_ELEMENTS_VERTICES";
		enumToStringMap[GL_MAX_EVAL_ORDER] = "GL_MAX_EVAL_ORDER";
		enumToStringMap[GL_MAX_FRAGMENT_UNIFORM_COMPONENTS] = "GL_MAX_FRAGMENT_UNIFORM_COMPONENTS";
		enumToStringMap[GL_MAX_LIGHTS] = "GL_MAX_LIGHTS";
		enumToStringMap[GL_MAX_LIST_NESTING] = "GL_MAX_LIST_NESTING";
		enumToStringMap[GL_MAX_MODELVIEW_STACK_DEPTH] = "GL_MAX_MODELVIEW_STACK_DEPTH";
		enumToStringMap[GL_MAX_NAME_STACK_DEPTH] = "GL_MAX_NAME_STACK_DEPTH";
		enumToStringMap[GL_MAX_PIXEL_MAP_TABLE] = "GL_MAX_PIXEL_MAP_TABLE";
		enumToStringMap[GL_MAX_PROJECTION_STACK_DEPTH] = "GL_MAX_PROJECTION_STACK_DEPTH";
		enumToStringMap[GL_MAX_TEXTURE_COORDS] = "GL_MAX_TEXTURE_COORDS";
		enumToStringMap[GL_MAX_TEXTURE_IMAGE_UNITS] = "GL_MAX_TEXTURE_IMAGE_UNITS";
		enumToStringMap[GL_MAX_TEXTURE_LOD_BIAS] = "GL_MAX_TEXTURE_LOD_BIAS";
		enumToStringMap[GL_MAX_TEXTURE_SIZE] = "GL_MAX_TEXTURE_SIZE";
		enumToStringMap[GL_MAX_TEXTURE_STACK_DEPTH] = "GL_MAX_TEXTURE_STACK_DEPTH";
		enumToStringMap[GL_MAX_TEXTURE_UNITS] = "GL_MAX_TEXTURE_UNITS";
		enumToStringMap[GL_MAX_VARYING_FLOATS] = "GL_MAX_VARYING_FLOATS";
		enumToStringMap[GL_MAX_VERTEX_ATTRIBS] = "GL_MAX_VERTEX_ATTRIBS";
		enumToStringMap[GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS] = "GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS";
		enumToStringMap[GL_MAX_VERTEX_UNIFORM_COMPONENTS] = "GL_MAX_VERTEX_UNIFORM_COMPONENTS";
		enumToStringMap[GL_MAX_VIEWPORT_DIMS] = "GL_MAX_VIEWPORT_DIMS";
		enumToStringMap[GL_MINMAX] = "GL_MINMAX";
		enumToStringMap[GL_MODELVIEW_MATRIX] = "GL_MODELVIEW_MATRIX";
		enumToStringMap[GL_MODELVIEW_STACK_DEPTH] = "GL_MODELVIEW_STACK_DEPTH";
		enumToStringMap[GL_NAME_STACK_DEPTH] = "GL_NAME_STACK_DEPTH";
		enumToStringMap[GL_NORMAL_ARRAY] = "GL_NORMAL_ARRAY";
		enumToStringMap[GL_NORMAL_ARRAY_BUFFER_BINDING] = "GL_NORMAL_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_NORMAL_ARRAY_STRIDE] = "GL_NORMAL_ARRAY_STRIDE";
		enumToStringMap[GL_NORMAL_ARRAY_TYPE] = "GL_NORMAL_ARRAY_TYPE";
		enumToStringMap[GL_NORMALIZE] = "GL_NORMALIZE";
		enumToStringMap[GL_NUM_COMPRESSED_TEXTURE_FORMATS] = "GL_NUM_COMPRESSED_TEXTURE_FORMATS";
		enumToStringMap[GL_PACK_ALIGNMENT] = "GL_PACK_ALIGNMENT";
		enumToStringMap[GL_PACK_IMAGE_HEIGHT] = "GL_PACK_IMAGE_HEIGHT";
		enumToStringMap[GL_PACK_LSB_FIRST] = "GL_PACK_LSB_FIRST";
		enumToStringMap[GL_PACK_ROW_LENGTH] = "GL_PACK_ROW_LENGTH";
		enumToStringMap[GL_PACK_SKIP_IMAGES] = "GL_PACK_SKIP_IMAGES";
		enumToStringMap[GL_PACK_SKIP_PIXELS] = "GL_PACK_SKIP_PIXELS";
		enumToStringMap[GL_PACK_SKIP_ROWS] = "GL_PACK_SKIP_ROWS";
		enumToStringMap[GL_PACK_SWAP_BYTES] = "GL_PACK_SWAP_BYTES";
		enumToStringMap[GL_PERSPECTIVE_CORRECTION_HINT] = "GL_PERSPECTIVE_CORRECTION_HINT";
		enumToStringMap[GL_PIXEL_MAP_A_TO_A_SIZE] = "GL_PIXEL_MAP_A_TO_A_SIZE";
		enumToStringMap[GL_PIXEL_MAP_B_TO_B_SIZE] = "GL_PIXEL_MAP_B_TO_B_SIZE";
		enumToStringMap[GL_PIXEL_MAP_G_TO_G_SIZE] = "GL_PIXEL_MAP_G_TO_G_SIZE";
		enumToStringMap[GL_PIXEL_MAP_I_TO_A_SIZE] = "GL_PIXEL_MAP_I_TO_A_SIZE";
		enumToStringMap[GL_PIXEL_MAP_I_TO_B_SIZE] = "GL_PIXEL_MAP_I_TO_B_SIZE";
		enumToStringMap[GL_PIXEL_MAP_I_TO_G_SIZE] = "GL_PIXEL_MAP_I_TO_G_SIZE";
		enumToStringMap[GL_PIXEL_MAP_I_TO_I_SIZE] = "GL_PIXEL_MAP_I_TO_I_SIZE";
		enumToStringMap[GL_PIXEL_MAP_I_TO_R_SIZE] = "GL_PIXEL_MAP_I_TO_R_SIZE";
		enumToStringMap[GL_PIXEL_MAP_R_TO_R_SIZE] = "GL_PIXEL_MAP_R_TO_R_SIZE";
		enumToStringMap[GL_PIXEL_MAP_S_TO_S_SIZE] = "GL_PIXEL_MAP_S_TO_S_SIZE";
		enumToStringMap[GL_PIXEL_PACK_BUFFER_BINDING] = "GL_PIXEL_PACK_BUFFER_BINDING";
		enumToStringMap[GL_PIXEL_UNPACK_BUFFER_BINDING] = "GL_PIXEL_UNPACK_BUFFER_BINDING";
		enumToStringMap[GL_POINT_DISTANCE_ATTENUATION] = "GL_POINT_DISTANCE_ATTENUATION";
		enumToStringMap[GL_POINT_FADE_THRESHOLD_SIZE] = "GL_POINT_FADE_THRESHOLD_SIZE";
		enumToStringMap[GL_POINT_SIZE] = "GL_POINT_SIZE";
		enumToStringMap[GL_POINT_SIZE_GRANULARITY] = "GL_POINT_SIZE_GRANULARITY";
		enumToStringMap[GL_POINT_SIZE_MAX] = "GL_POINT_SIZE_MAX";
		enumToStringMap[GL_POINT_SIZE_MIN] = "GL_POINT_SIZE_MIN";
		enumToStringMap[GL_POINT_SIZE_RANGE] = "GL_POINT_SIZE_RANGE";
		enumToStringMap[GL_POINT_SMOOTH] = "GL_POINT_SMOOTH";
		enumToStringMap[GL_POINT_SMOOTH_HINT] = "GL_POINT_SMOOTH_HINT";
		enumToStringMap[GL_POINT_SPRITE] = "GL_POINT_SPRITE";
		enumToStringMap[GL_POLYGON_MODE] = "GL_POLYGON_MODE";
		enumToStringMap[GL_POLYGON_OFFSET_FACTOR] = "GL_POLYGON_OFFSET_FACTOR";
		enumToStringMap[GL_POLYGON_OFFSET_UNITS] = "GL_POLYGON_OFFSET_UNITS";
		enumToStringMap[GL_POLYGON_OFFSET_FILL] = "GL_POLYGON_OFFSET_FILL";
		enumToStringMap[GL_POLYGON_OFFSET_LINE] = "GL_POLYGON_OFFSET_LINE";
		enumToStringMap[GL_POLYGON_OFFSET_POINT] = "GL_POLYGON_OFFSET_POINT";
		enumToStringMap[GL_POLYGON_SMOOTH] = "GL_POLYGON_SMOOTH";
		enumToStringMap[GL_POLYGON_SMOOTH_HINT] = "GL_POLYGON_SMOOTH_HINT";
		enumToStringMap[GL_POLYGON_STIPPLE] = "GL_POLYGON_STIPPLE";
		enumToStringMap[GL_POST_COLOR_MATRIX_COLOR_TABLE] = "GL_POST_COLOR_MATRIX_COLOR_TABLE";
		enumToStringMap[GL_POST_COLOR_MATRIX_RED_BIAS] = "GL_POST_COLOR_MATRIX_RED_BIAS";
		enumToStringMap[GL_POST_COLOR_MATRIX_GREEN_BIAS] = "GL_POST_COLOR_MATRIX_GREEN_BIAS";
		enumToStringMap[GL_POST_COLOR_MATRIX_BLUE_BIAS] = "GL_POST_COLOR_MATRIX_BLUE_BIAS";
		enumToStringMap[GL_POST_COLOR_MATRIX_ALPHA_BIAS] = "GL_POST_COLOR_MATRIX_ALPHA_BIAS";
		enumToStringMap[GL_POST_COLOR_MATRIX_RED_SCALE] = "GL_POST_COLOR_MATRIX_RED_SCALE";
		enumToStringMap[GL_POST_COLOR_MATRIX_GREEN_SCALE] = "GL_POST_COLOR_MATRIX_GREEN_SCALE";
		enumToStringMap[GL_POST_COLOR_MATRIX_BLUE_SCALE] = "GL_POST_COLOR_MATRIX_BLUE_SCALE";
		enumToStringMap[GL_POST_COLOR_MATRIX_ALPHA_SCALE] = "GL_POST_COLOR_MATRIX_ALPHA_SCALE";
		enumToStringMap[GL_POST_CONVOLUTION_COLOR_TABLE] = "GL_POST_CONVOLUTION_COLOR_TABLE";
		enumToStringMap[GL_POST_CONVOLUTION_RED_BIAS] = "GL_POST_CONVOLUTION_RED_BIAS";
		enumToStringMap[GL_POST_CONVOLUTION_GREEN_BIAS] = "GL_POST_CONVOLUTION_GREEN_BIAS";
		enumToStringMap[GL_POST_CONVOLUTION_BLUE_BIAS] = "GL_POST_CONVOLUTION_BLUE_BIAS";
		enumToStringMap[GL_POST_CONVOLUTION_ALPHA_BIAS] = "GL_POST_CONVOLUTION_ALPHA_BIAS";
		enumToStringMap[GL_POST_CONVOLUTION_RED_SCALE] = "GL_POST_CONVOLUTION_RED_SCALE";
		enumToStringMap[GL_POST_CONVOLUTION_GREEN_SCALE] = "GL_POST_CONVOLUTION_GREEN_SCALE";
		enumToStringMap[GL_POST_CONVOLUTION_BLUE_SCALE] = "GL_POST_CONVOLUTION_BLUE_SCALE";
		enumToStringMap[GL_POST_CONVOLUTION_ALPHA_SCALE] = "GL_POST_CONVOLUTION_ALPHA_SCALE";
		enumToStringMap[GL_PROJECTION_MATRIX] = "GL_PROJECTION_MATRIX";
		enumToStringMap[GL_PROJECTION_STACK_DEPTH] = "GL_PROJECTION_STACK_DEPTH";
		enumToStringMap[GL_READ_BUFFER] = "GL_READ_BUFFER";
		enumToStringMap[GL_RED_BIAS] = "GL_RED_BIAS";
		enumToStringMap[GL_RED_BITS] = "GL_RED_BITS";
		enumToStringMap[GL_RED_SCALE] = "GL_RED_SCALE";
		enumToStringMap[GL_RENDER_MODE] = "GL_RENDER_MODE";
		enumToStringMap[GL_RESCALE_NORMAL] = "GL_RESCALE_NORMAL";
		enumToStringMap[GL_RGBA_MODE] = "GL_RGBA_MODE";
		enumToStringMap[GL_SAMPLE_BUFFERS] = "GL_SAMPLE_BUFFERS";
		enumToStringMap[GL_SAMPLE_COVERAGE_VALUE] = "GL_SAMPLE_COVERAGE_VALUE";
		enumToStringMap[GL_SAMPLE_COVERAGE_INVERT] = "GL_SAMPLE_COVERAGE_INVERT";
		enumToStringMap[GL_SAMPLES] = "GL_SAMPLES";
		enumToStringMap[GL_SCISSOR_BOX] = "GL_SCISSOR_BOX";
		enumToStringMap[GL_SCISSOR_TEST] = "GL_SCISSOR_TEST";
		enumToStringMap[GL_SECONDARY_COLOR_ARRAY] = "GL_SECONDARY_COLOR_ARRAY";
		enumToStringMap[GL_SECONDARY_COLOR_ARRAY_BUFFER_BINDING] = "GL_SECONDARY_COLOR_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_SECONDARY_COLOR_ARRAY_SIZE] = "GL_SECONDARY_COLOR_ARRAY_SIZE";
		enumToStringMap[GL_SECONDARY_COLOR_ARRAY_STRIDE] = "GL_SECONDARY_COLOR_ARRAY_STRIDE";
		enumToStringMap[GL_SECONDARY_COLOR_ARRAY_TYPE] = "GL_SECONDARY_COLOR_ARRAY_TYPE";
		enumToStringMap[GL_SELECTION_BUFFER_SIZE] = "GL_SELECTION_BUFFER_SIZE";
		enumToStringMap[GL_SEPARABLE_2D] = "GL_SEPARABLE_2D";
		enumToStringMap[GL_SHADE_MODEL] = "GL_SHADE_MODEL";
		enumToStringMap[GL_SMOOTH_LINE_WIDTH_RANGE] = "GL_SMOOTH_LINE_WIDTH_RANGE";
		enumToStringMap[GL_SMOOTH_LINE_WIDTH_GRANULARITY] = "GL_SMOOTH_LINE_WIDTH_GRANULARITY";
		enumToStringMap[GL_SMOOTH_POINT_SIZE_RANGE] = "GL_SMOOTH_POINT_SIZE_RANGE";
		enumToStringMap[GL_SMOOTH_POINT_SIZE_GRANULARITY] = "GL_SMOOTH_POINT_SIZE_GRANULARITY";
		enumToStringMap[GL_STENCIL_BACK_FAIL] = "GL_STENCIL_BACK_FAIL";
		enumToStringMap[GL_STENCIL_BACK_FUNC] = "GL_STENCIL_BACK_FUNC";
		enumToStringMap[GL_STENCIL_BACK_PASS_DEPTH_FAIL] = "GL_STENCIL_BACK_PASS_DEPTH_FAIL";
		enumToStringMap[GL_STENCIL_BACK_PASS_DEPTH_PASS] = "GL_STENCIL_BACK_PASS_DEPTH_PASS";
		enumToStringMap[GL_STENCIL_BACK_REF] = "GL_STENCIL_BACK_REF";
		enumToStringMap[GL_STENCIL_BACK_VALUE_MASK] = "GL_STENCIL_BACK_VALUE_MASK";
		enumToStringMap[GL_STENCIL_BACK_WRITEMASK] = "GL_STENCIL_BACK_WRITEMASK";
		enumToStringMap[GL_STENCIL_BITS] = "GL_STENCIL_BITS";
		enumToStringMap[GL_STENCIL_CLEAR_VALUE] = "GL_STENCIL_CLEAR_VALUE";
		enumToStringMap[GL_STENCIL_FAIL] = "GL_STENCIL_FAIL";
		enumToStringMap[GL_STENCIL_FUNC] = "GL_STENCIL_FUNC";
		enumToStringMap[GL_STENCIL_PASS_DEPTH_FAIL] = "GL_STENCIL_PASS_DEPTH_FAIL";
		enumToStringMap[GL_STENCIL_PASS_DEPTH_PASS] = "GL_STENCIL_PASS_DEPTH_PASS";
		enumToStringMap[GL_STENCIL_REF] = "GL_STENCIL_REF";
		enumToStringMap[GL_STENCIL_TEST] = "GL_STENCIL_TEST";
		enumToStringMap[GL_STENCIL_VALUE_MASK] = "GL_STENCIL_VALUE_MASK";
		enumToStringMap[GL_STENCIL_WRITEMASK] = "GL_STENCIL_WRITEMASK";
		enumToStringMap[GL_STEREO] = "GL_STEREO";
		enumToStringMap[GL_SUBPIXEL_BITS] = "GL_SUBPIXEL_BITS";
		enumToStringMap[GL_TEXTURE_1D] = "GL_TEXTURE_1D";
		enumToStringMap[GL_TEXTURE_BINDING_1D] = "GL_TEXTURE_BINDING_1D";
		enumToStringMap[GL_TEXTURE_2D] = "GL_TEXTURE_2D";
		enumToStringMap[GL_TEXTURE_BINDING_2D] = "GL_TEXTURE_BINDING_2D";
		enumToStringMap[GL_TEXTURE_3D] = "GL_TEXTURE_3D";
		enumToStringMap[GL_TEXTURE_BINDING_3D] = "GL_TEXTURE_BINDING_3D";
		enumToStringMap[GL_TEXTURE_BINDING_CUBE_MAP] = "GL_TEXTURE_BINDING_CUBE_MAP";
		enumToStringMap[GL_TEXTURE_COMPRESSION_HINT] = "GL_TEXTURE_COMPRESSION_HINT";
		enumToStringMap[GL_TEXTURE_COORD_ARRAY] = "GL_TEXTURE_COORD_ARRAY";
		enumToStringMap[GL_TEXTURE_COORD_ARRAY_BUFFER_BINDING] = "GL_TEXTURE_COORD_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_TEXTURE_COORD_ARRAY_SIZE] = "GL_TEXTURE_COORD_ARRAY_SIZE";
		enumToStringMap[GL_TEXTURE_COORD_ARRAY_STRIDE] = "GL_TEXTURE_COORD_ARRAY_STRIDE";
		enumToStringMap[GL_TEXTURE_COORD_ARRAY_TYPE] = "GL_TEXTURE_COORD_ARRAY_TYPE";
		enumToStringMap[GL_TEXTURE_CUBE_MAP] = "GL_TEXTURE_CUBE_MAP";
		enumToStringMap[GL_TEXTURE_GEN_Q] = "GL_TEXTURE_GEN_Q";
		enumToStringMap[GL_TEXTURE_GEN_R] = "GL_TEXTURE_GEN_R";
		enumToStringMap[GL_TEXTURE_GEN_S] = "GL_TEXTURE_GEN_S";
		enumToStringMap[GL_TEXTURE_GEN_T] = "GL_TEXTURE_GEN_T";
		enumToStringMap[GL_TEXTURE_MATRIX] = "GL_TEXTURE_MATRIX";
		enumToStringMap[GL_TEXTURE_STACK_DEPTH] = "GL_TEXTURE_STACK_DEPTH";
		enumToStringMap[GL_TRANSPOSE_COLOR_MATRIX] = "GL_TRANSPOSE_COLOR_MATRIX";
		enumToStringMap[GL_TRANSPOSE_MODELVIEW_MATRIX] = "GL_TRANSPOSE_MODELVIEW_MATRIX";
		enumToStringMap[GL_TRANSPOSE_PROJECTION_MATRIX] = "GL_TRANSPOSE_PROJECTION_MATRIX";
		enumToStringMap[GL_TRANSPOSE_TEXTURE_MATRIX] = "GL_TRANSPOSE_TEXTURE_MATRIX";
		enumToStringMap[GL_UNPACK_ALIGNMENT] = "GL_UNPACK_ALIGNMENT";
		enumToStringMap[GL_UNPACK_IMAGE_HEIGHT] = "GL_UNPACK_IMAGE_HEIGHT";
		enumToStringMap[GL_UNPACK_LSB_FIRST] = "GL_UNPACK_LSB_FIRST";
		enumToStringMap[GL_UNPACK_SKIP_PIXELS] = "GL_UNPACK_SKIP_PIXELS";
		enumToStringMap[GL_UNPACK_SKIP_ROWS] = "GL_UNPACK_SKIP_ROWS";
		enumToStringMap[GL_UNPACK_SWAP_BYTES] = "GL_UNPACK_SWAP_BYTES";
		enumToStringMap[GL_VERTEX_ARRAY] = "GL_VERTEX_ARRAY";
		enumToStringMap[GL_VERTEX_ARRAY_BUFFER_BINDING] = "GL_VERTEX_ARRAY_BUFFER_BINDING";
		enumToStringMap[GL_VERTEX_ARRAY_SIZE] = "GL_VERTEX_ARRAY_SIZE";
		enumToStringMap[GL_VERTEX_ARRAY_STRIDE] = "GL_VERTEX_ARRAY_STRIDE";
		enumToStringMap[GL_VERTEX_ARRAY_TYPE] = "GL_VERTEX_ARRAY_TYPE";
		enumToStringMap[GL_VERTEX_PROGRAM_POINT_SIZE] = "GL_VERTEX_PROGRAM_POINT_SIZE";
		enumToStringMap[GL_VERTEX_PROGRAM_TWO_SIDE] = "GL_VERTEX_PROGRAM_TWO_SIDE";
		enumToStringMap[GL_VIEWPORT] = "GL_VIEWPORT";
		enumToStringMap[GL_ZOOM_X] = "GL_ZOOM_X";
		enumToStringMap[GL_ZOOM_Y] = "GL_ZOOM_Y";

		enumToStringMap[GL_NEAREST] = "GL_NEAREST";
		enumToStringMap[GL_LINEAR] = "GL_LINEAR";

		enumToStringMap[GL_ZERO] = "GL_ZERO";
		enumToStringMap[GL_ONE] = "GL_ONE";
		enumToStringMap[GL_SRC_COLOR] = "GL_SRC_COLOR";
		enumToStringMap[GL_ONE_MINUS_SRC_COLOR] = "GL_ONE_MINUS_SRC_COLOR";
		enumToStringMap[GL_DST_COLOR] = "GL_DST_COLOR";
		enumToStringMap[GL_ONE_MINUS_DST_COLOR] = "GL_ONE_MINUS_DST_COLOR";
		enumToStringMap[GL_SRC_ALPHA] = "GL_SRC_ALPHA";
		enumToStringMap[GL_ONE_MINUS_SRC_ALPHA] = "GL_ONE_MINUS_SRC_ALPHA";
		enumToStringMap[GL_DST_ALPHA] = "GL_DST_ALPHA";
		enumToStringMap[GL_ONE_MINUS_DST_ALPHA] = "GL_ONE_MINUS_DST_ALPHA";
		enumToStringMap[GL_CONSTANT_COLOR] = "GL_CONSTANT_COLOR";
		enumToStringMap[GL_ONE_MINUS_CONSTANT_COLOR] = "GL_ONE_MINUS_CONSTANT_COLOR";
		enumToStringMap[GL_CONSTANT_ALPHA] = "GL_CONSTANT_ALPHA";
		enumToStringMap[GL_ONE_MINUS_CONSTANT_ALPHA] = "GL_ONE_MINUS_CONSTANT_ALPHA";
		enumToStringMap[GL_SRC_ALPHA_SATURATE] = "GL_SRC_ALPHA_SATURATE";

		enumToStringMap[GL_POINTS] = "GL_POINTS";
		enumToStringMap[GL_LINES] = "GL_LINES";
		enumToStringMap[GL_LINE_STRIP] = "GL_LINE_STRIP";
		enumToStringMap[GL_LINE_LOOP] = "GL_LINE_LOOP";
		enumToStringMap[GL_TRIANGLES] = "GL_TRIANGLES";
		enumToStringMap[GL_TRIANGLE_STRIP] = "GL_TRIANGLE_STRIP";
		enumToStringMap[GL_TRIANGLE_FAN] = "GL_TRIANGLE_FAN";
		enumToStringMap[GL_QUADS] = "GL_QUADS";
		enumToStringMap[GL_QUAD_STRIP] = "GL_QUAD_STRIP";
		enumToStringMap[GL_POLYGON] = "GL_POLYGON";

	}

	map<int, string>::iterator it = enumToStringMap.find(inEnum);
	static char text[256] = { 0 };

	if( it != enumToStringMap.end() ) {
		strcpy(text, it->second.c_str());
	} else {
		sprintf(text, "%d", inEnum);
	}

	return text;
}

void _glPixelStorei(GLenum inEnum, int inValue) {
	glGetError();
	LogGLFunction2("glPixielStorei", GL_EnumToString(inEnum), inValue);
	
	::glPixelStorei(inEnum, inValue);

	GL_ReportError("glPixelStorei");
}

void _glDrawBuffer(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glPixielStorei", GL_EnumToString(inEnum));

	glDrawBuffer(inEnum);

	GL_ReportError("glDrawBuffer");
}

void _glClearColor(float inR, float inG, float inB, float inA) {
	glGetError();
	LogGLFunction4("glClearColor", inR, inG, inB, inA);

	glClearColor(inR, inG, inB, inA);

	GL_ReportError("glClearColor");
}

void _glClear(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glClear", GL_EnumToString(inEnum));

	glClear(inEnum);
	GL_ReportError("glClear");
}

void _glReadBuffer(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glReadBuffer", inEnum);

	glReadBuffer(inEnum);
	GL_ReportError("glReadBuffer");
}

void _glColorMask(GLboolean inR, GLboolean inG, GLboolean inB, GLboolean inA) {
	glGetError();
	LogGLFunction4("glColorMask", GL_BooleanToString(inR), GL_BooleanToString(inG), GL_BooleanToString(inB), GL_BooleanToString(inA));

	glColorMask(inR, inG, inB, inA);
	GL_ReportError("glColorMask");
}

void _glShadeModel(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glShadeModel", GL_EnumToString(inEnum));

	glShadeModel(inEnum);
	GL_ReportError("glShadeModel");
}

void _glPolygonOffset(float inFactor, float inBias) {
	glGetError();
	LogGLFunction2("glPolygonOffset", inFactor, inBias);

	glPolygonOffset(inFactor, inBias);
	GL_ReportError("glPolygonOffset");
}

void _glMatrixMode(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glMatrixMode", GL_EnumToString(inEnum));
	
	glMatrixMode(inEnum);
	GL_ReportError("glMatrixMode");
}

void _glLoadMatrixf(GLfloat * inMat) {
	glGetError();
	LogGLFunction1("glLoadMatrixf", inMat);
	glLoadMatrixf(inMat);
	GL_ReportError("glLoadMatrixf");
}

void _glEnable(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glEnable", GL_EnumToString(inEnum));

	glEnable(inEnum);
	GL_ReportError("glEnable");
}

void _glDisable(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glDisable", GL_EnumToString(inEnum));

	glDisable(inEnum);
	GL_ReportError("glDisable");
}

void _glBindBuffer(GLenum target, GLuint buffer) {
	glGetError();
	LogGLFunction2("glBindBuffer", GL_EnumToString(target), buffer);

	glBindBuffer(target, buffer);
	GL_ReportError("glBindBuffer");
}

void _glUnmapBuffer(GLenum target) {
	glGetError();
	LogGLFunction1("glUnmapBuffer", GL_EnumToString(target));

	glUnmapBuffer(target);
	GL_ReportError("glUnmapBuffer");
}

void _glBindFramebuffer(GLenum inEnum, GLuint inFBO) {
	glGetError();
	LogGLFunction2("glBindFramebuffer", GL_EnumToString(inEnum), inFBO);
	glBindFramebuffer(inEnum, inFBO);
	GL_ReportError("glBindFramebuffer");
}

void _glFramebufferRenderbuffer(GLenum inEnum1, GLenum inEnum2, GLenum inEnum3, GLuint inRB) {
	glGetError();
	LogGLFunction4("glFramebufferRenderbuffer", GL_EnumToString(inEnum1), GL_EnumToString(inEnum2), GL_EnumToString(inEnum3), inRB);

	glFramebufferRenderbuffer(inEnum1, inEnum2, inEnum3, inRB);
	GL_ReportError("glFramebufferRenderbuffer");
}

void _glBlitFramebuffer(GLint inX0, GLint inY0, GLint inW0, GLint inH0, GLint inX1, GLint inY1, GLint inW1, GLint inH1, GLbitfield inBitfield, GLenum inEnum) {
	glGetError();
	LogGLFunction10("glBlitFramebuffer", inX0, inY0, inW0, inH0, inX1, inY1, inW1, inH1, inBitfield, GL_EnumToString(inEnum));

	glBlitFramebuffer(inX0, inY0, inW0, inH0, inX1, inY1, inW1, inH1, inBitfield, inEnum);
	GL_ReportError("glBlitFramebuffer");
}
void _glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels) {
	glGetError();
	LogGLFunction7("glReadPixels", x, y, width, height, GL_EnumToString(format), GL_EnumToString(type), pixels);

	glReadPixels(x, y, width, height, format, type, pixels);
	GL_ReportError("glReadPixels");
}

void _glFlush() {
	glGetError();
	LogGLFunction0("glFlush");

	glFlush();
	GL_ReportError("glFlush");
}

void _glFinish() {
	glGetError();
	LogGLFunction0("glFinish");

	glFinish();
	GL_ReportError("glFinish");
}

void _glLoadIdentity() {
	glGetError();
	LogGLFunction0("glLoadIdentity");

	glLoadIdentity();
	GL_ReportError("glLoadIdentity");
}

void _glViewport(GLint x, GLint y, GLsizei width, GLsizei height) {
	glGetError();
	LogGLFunction4("glViewport", x, y, width, height);

	glViewport(x, y, width, height);
	GL_ReportError("glViewport");
}

void _glFramebufferTexture2D(GLenum inEnum1, GLenum inEnum2, GLenum inEnum3, GLuint inUint, GLint inInt) {
	glGetError();
	LogGLFunction5("glFramebufferTexture2D", GL_EnumToString(inEnum1), GL_EnumToString(inEnum2), GL_EnumToString(inEnum3), inUint, inInt);

	glFramebufferTexture2D(inEnum1, inEnum2, inEnum3, inUint, inInt);
	GL_ReportError("glFramebufferTexture2D");
}

void _glBindTexture(GLenum target, GLuint texture) {
	glGetError();
	LogGLFunction2("glBindTexture", GL_EnumToString(target), texture);

	glBindTexture(target, texture);
	GL_ReportError("glBindTexture");
}

void _glGetTexImage(GLenum target, GLint level, GLenum format, GLenum type, GLvoid *pixels) {
	glGetError();
	LogGLFunction5("glGetTexImage", GL_EnumToString(target), level, GL_EnumToString(format), GL_EnumToString(type), pixels);

	glGetTexImage(target, level, format, type, pixels);
	GL_ReportError("glGetTexImage");
}

void _glActiveTexture(GLenum texture) {
	glGetError();
	LogGLFunction1("glActiveTexture", GL_EnumToString(texture));

	glActiveTexture(texture);
	GL_ReportError("glActiveTexture");
}

void _glTexSubImage2D(GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum format, GLenum type, const GLvoid *pixels) {
	glGetError();
	LogGLFunction9("glTexSubImage2D", GL_EnumToString(target), level, xoffset, yoffset, width, height, GL_EnumToString(format), GL_EnumToString(type), pixels);

	glTexSubImage2D(target, level, xoffset, yoffset, width, height, format, type, pixels);
	GL_ReportError("glTexSubImage2D");
}

void _glDrawBuffers(GLsizei n, const GLenum *bufs) {
	glGetError();
	LogGLFunction2("glDrawBuffers", n, bufs);

	glDrawBuffers(n, bufs);
	GL_ReportError("glDrawBuffers");
}

void _glUseProgram(GLuint program) {
	glGetError();
	LogGLFunction1("glUseProgram", program);

	glUseProgram(program);
	GL_ReportError("glUseProgram");
}

void _glUniform3fv(GLint location, GLsizei count, const GLfloat *value) {
	glGetError();
	LogGLFunction3("glUniform3fv", location, count, value);

	glUniform3fv(location, count, value);
	GL_ReportError("glUniform3fv");
}

void _glUniform1i(GLint location, GLint v0) {
	glGetError();
	LogGLFunction2("glUniform1i", location, v0);

	glUniform1i(location, v0);
	GL_ReportError("glUniform1i");
}

void _glUniform1fv(GLint location, GLsizei count, const GLfloat *value) {
	glGetError();
	LogGLFunction3("glUniform1fv", location, count, value);

	glUniform1fv(location, count, value);
	GL_ReportError("glUniform1fv");
}

void _glUniform2fv(GLint location, GLsizei count, const GLfloat *value) {
	glGetError();
	LogGLFunction3("glUniform2fv", location, count, value);

	glUniform2fv(location, count, value);
	GL_ReportError("glUniform2fv");
}

void _glUniform1iv(GLint location, GLsizei count, const GLint *value) {
	glGetError();
	LogGLFunction3("glUniform1iv", location, count, value);

	glUniform1iv(location, count, value);
	GL_ReportError("glUniform1iv");
}

void _glBlendFunc(GLenum sfactor, GLenum dfactor) {
	glGetError();
	LogGLFunction2("glBlendFunc", GL_EnumToString(sfactor), GL_EnumToString(dfactor));

	glBlendFunc(sfactor, dfactor);
	GL_ReportError("glBlendFunc");
}

void _glDepthMask(GLboolean flag) {
	glGetError();
	LogGLFunction1("glDepthMask", GL_BooleanToString(flag));

	glDepthMask(flag);
	GL_ReportError("glDepthMask");
}

GLenum _glCheckFramebufferStatus(GLenum inEnum) {
	glGetError();
	LogGLFunction1("glCheckFramebufferStatus", GL_EnumToString(inEnum));

	GLenum status = glCheckFramebufferStatus(inEnum);
	GL_ReportError("glCheckFramebufferStatus");
	
	static map<GLenum, string> errorMap;
	if( errorMap.empty() ) {
		errorMap[GL_FRAMEBUFFER_UNDEFINED] = "GL_FRAMEBUFFER_UNDEFINED";
		errorMap[GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT] = "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT";
		errorMap[GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT] = "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT";
		errorMap[GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER] = "GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER";
		errorMap[GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER] = "GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER";
		errorMap[GL_FRAMEBUFFER_UNSUPPORTED] = "GL_FRAMEBUFFER_UNSUPPORTED";
		errorMap[GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE] = "GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE";
	}

	char text[512] = { 0 };
	map<GLenum, string>::iterator it = errorMap.find(status);
	if( it != errorMap.end() ) {
		sprintf(text, "OpenGL Error: glCheckFramebufferStatus returned %s\n", it->second.c_str());
		DMSG((kEveryone, text));
	} else if( status != GL_FRAMEBUFFER_COMPLETE ) {
		sprintf(text, "OpenGL Error: glCheckFramebufferStatus returned %d\n", status);
		DMSG((kEveryone, text));
	}
	
	if( text[0] != 0 && gGLDebugOut.is_open() ) {
		gGLDebugOut << text << endl;
		gGLDebugOut.flush();
	}

	return status;
}

void _glGetFramebufferAttachmentParameterivEXT(GLenum target, GLenum attachment, GLenum pname, GLint *params) {
	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT, params);
	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT, params);

	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT, params);
	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT, params);

	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT, params);
	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT, params);
	glGetError();
	
	glGetError();
	LogGLFunction4("glGetFramebufferAttachmentParameterivEXT", GL_EnumToString(target), GL_EnumToString(attachment), pname, params);
	
	glGetFramebufferAttachmentParameterivEXT(target, attachment, pname, params);
	GL_ReportError("glGetFramebufferAttachmentParameterivEXT");
}

void _glGetBooleanv(GLenum pname, GLboolean * params) {
	glGetError();
	LogGLFunction2("glGetBooleanv", GL_EnumToString(pname), params);

	glGetBooleanv(pname, params);
	GL_ReportError("glGetBooleanv");
}

void _glGetDoublev(GLenum pname, GLdouble * params) {
	glGetError();
	LogGLFunction2("glGetDoublev", GL_EnumToString(pname), params);

	glGetDoublev(pname, params);
	GL_ReportError("glGetDoublev");
}

void _glGetFloatv(GLenum pname, GLfloat * params) {
	glGetError();
	LogGLFunction2("glGetFloatv", GL_EnumToString(pname), params);

	glGetFloatv(pname, params);
	GL_ReportError("glGetFloatv");
}

void _glGetIntegerv(GLenum pname, GLint * params) {
	glGetError();
	LogGLFunction2("glGetIntegerv", GL_EnumToString(pname), params);

	glGetIntegerv(pname, params);
	GL_ReportError("glGetIntegerv");
}

void _glCallList(GLuint list) {
	glGetError();
	LogGLFunction1("glCallList", list);

	glCallList(list);
	GL_ReportError("glCallList");
}

void _glBegin(GLenum type) {
	glGetError();
	LogGLFunction1("glBegin", GL_EnumToString(type));

	glBegin(type);
	//GL_ReportError("glBegin"); // can't call glGetError between glBegin and glEnd

#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glBegin(type);
#endif
}

void _glEnd() {	
	// glGetError(); // can't call glGetError between glBegin and glEnd
	LogGLFunction0("glEnd");

	glEnd();
	GL_ReportError("glEnd");

#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glEnd();
#endif
}

void _glNormal3d(GLdouble x, GLdouble y, GLdouble z) {
	LogGLFunction3("glNormal3d", x, y, z);

	glNormal3d(x, y, z);
	
#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glNormal3d(x, y, z);
#endif
}

void _glNormal3dv(GLdouble * verts) {
	LogGLFunction1("glNormal3dv", verts);

	glNormal3dv(verts);
	
#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glNormal3dv(verts);
#endif
}

void _glVertex3d(GLdouble x, GLdouble y, GLdouble z) {
	LogGLFunction3("glVertex3d", x, y, z);

	glVertex3d(x,  y, z);
	
#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glVertex3d(x, y, z);
#endif
}

void _glVertex3dv(GLdouble * verts) {
	LogGLFunction1("glVertex3dv", verts);

	glVertex3dv(verts);
	
#if EXPORT_GEO_FOR_CULLING_DEMO
	CULL_glVertex3dv(verts);
#endif
}

void _glEdgeFlag(GLboolean flag) {
	LogGLFunction1("glEdgeFlag", GL_BooleanToString(flag));

	glEdgeFlag(flag);	
}

void _glEdgeFlagv(GLboolean * flag) {
	LogGLFunction1("glEdgeFlagv", flag);

	glEdgeFlagv(flag);
}


#undef LogGLFunction0
#undef LogGLFunction1
#undef LogGLFunction2
#undef LogGLFunction3
#undef LogGLFunction4
#undef LogGLFunction5
#undef LogGLFunction6
#undef LogGLFunction7
#undef LogGLFunction8
#undef LogGLFunction9
#undef LogGLFunction10

#if EXPORT_GEO_FOR_CULLING_DEMO

SymGeo * gSymGeo = 0;
map<BlankHandle, SymGeo *> gSymGeoMap;
bool gWritingMeshData;
bool gWritingVertexData;
ofstream geometryStream;
PipelineEnv * gEnv = 0;

double gNormal[3*3];
double gVertices[3*2];
int gFanState = 0;

void CULL_glBegin(GLenum type) {
	if( gWritingMeshData ) {
		gWritingVertexData = true;
		gFanState = (type==GL_QUADS || type==GL_POLYGON) ? 1 : 0;
	}
}

void CULL_glEnd() {
	if( gWritingVertexData ) {
		if( !gSymGeo ) {
			geometryStream << endl;
		}
		gWritingVertexData = false;
	}
}

void CULL_glNormal3d(GLdouble x, GLdouble y, GLdouble z) {
	if( gEnv ) {
		TransformMatrix gTransform = gEnv->GetTopMatrix();
		gTransform.P() /= 25.4;
		WorldPt3 p = VectorTransformN(WorldPt3(x, y, z), gTransform);
		x = p.x;
		y = p.y;
		z = p.z;
		
		if( gFanState == 0 || gFanState == 1 ) {		
			gNormal[0] = gNormal[3] = gNormal[6] = x;
			gNormal[1] = gNormal[4] = gNormal[7] = y;
			gNormal[2] = gNormal[5] = gNormal[8] = z;	
		} else if( gFanState == 2 ) {
			gNormal[6] = x;
			gNormal[7] = y;
			gNormal[8] = z;
		} else if( gFanState == 3 ) {
			gNormal[3] = gNormal[6];
			gNormal[4] = gNormal[7];
			gNormal[5] = gNormal[8];

			gNormal[6] = x;
			gNormal[7] = y;
			gNormal[8] = z;
		}
	}
}

void CULL_glNormal3dv(GLdouble * verts) {
	CULL_glNormal3d(verts[0], verts[1], verts[2]);
}

void CULL_glVertex3d(GLdouble x, GLdouble y, GLdouble z) {
	if( gEnv ) {
		TransformMatrix gTransform = gEnv->GetTopMatrix();
		gTransform.P() /= 25.4;
		WorldPt3 p = PointTransformN(WorldPt3(x, y, z), gTransform);
		x = p.x;
		y = p.y;
		z = p.z;
	
		if( gWritingVertexData ) {
			if( gFanState == 0 ) {
				if( gSymGeo ) {
					vector<double> & n = gSymGeo->normals;
					vector<double> & v = gSymGeo->vertices;
					n.push_back(gNormal[0]); n.push_back(gNormal[1]); n.push_back(gNormal[2]);
					v.push_back(x); v.push_back(y); v.push_back(z);
				} else {
					geometryStream << gNormal[0] << " " << gNormal[1] << " " << gNormal[2] << " ";
					geometryStream << x << " " << y << " " << z << " ";
				}
			} else if( gFanState == 1 ) {
				gVertices[0] = x;
				gVertices[1] = y;
				gVertices[2] = z;
				gFanState++;
			} else if( gFanState == 2 ) {
				gVertices[3] = x;
				gVertices[4] = y;
				gVertices[5] = z;
				gFanState++;
			} else if( gFanState == 3 ) {
				if( gSymGeo ) {
					vector<double> & n = gSymGeo->normals;
					vector<double> & v = gSymGeo->vertices;
				
					n.push_back(gNormal[0]); n.push_back(gNormal[1]); n.push_back(gNormal[2]);
					v.push_back(gVertices[0]); v.push_back(gVertices[1]); v.push_back(gVertices[2]);

					n.push_back(gNormal[3]); n.push_back(gNormal[4]); n.push_back(gNormal[5]);
					v.push_back(gVertices[3]); v.push_back(gVertices[4]); v.push_back(gVertices[5]);
			
					n.push_back(gNormal[6]); n.push_back(gNormal[7]); n.push_back(gNormal[8]);
					v.push_back(x); v.push_back(y); v.push_back(z);

				} else {
					geometryStream << gNormal[0] << " " << gNormal[1] << " " << gNormal[2] << " ";
					geometryStream << gVertices[0] << " " << gVertices[1] << " " << gVertices[2] << " ";
			
					geometryStream << gNormal[3] << " " << gNormal[4] << " " << gNormal[5] << " ";
					geometryStream << gVertices[3] << " " << gVertices[4] << " " << gVertices[5] << " ";
			
					geometryStream << gNormal[6] << " " << gNormal[7] << " " << gNormal[8] << " ";
					geometryStream << x << " " << y << " " << z << " ";
				}
			
				gVertices[3] = x;
				gVertices[4] = y;
				gVertices[5] = z;
			}
		}
	}
}

void CULL_glVertex3dv(GLdouble * verts) {
	CULL_glVertex3d(verts[0], verts[1], verts[2]);
}

#endif

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

#endif
