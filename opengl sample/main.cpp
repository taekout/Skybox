#define GS_WIN 1

#include "Resource.h"
#include <Windows.h>
#include "utilities.h"
#include <math.h>
#include <list>
#include <stack>
#include <iomanip>
#include <strstream>
#include <GdiPlus.h>
#include "sample.h"

using namespace Gdiplus;
using namespace std;

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;								// current instance
TCHAR szTitle[MAX_LOADSTRING];					// The title bar text
TCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name

// Forward declarations of functions included in this code module:
ATOM				MyRegisterClass(HINSTANCE hInstance);
BOOL				InitInstance(HINSTANCE, int);
LRESULT CALLBACK	WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	About(HWND, UINT, WPARAM, LPARAM);

HWND gMainWindow = 0; // main window

int gWindowWidth = 0;
int gWindowHeight = 0;

int WIN_GetWidth() {
	return gWindowWidth;
}

int WIN_GetHeight() {
	return gWindowHeight;
}

int gInitialTick = 0;
float TIME_GetTick() {
	return float(GetTickCount()-gInitialTick) / 1000.0f;
}

int gCurrentDrawTest = 0;
int APP_GetCurrentDrawTest() {
	return gCurrentDrawTest;
}

void CheckGLError() {
	GLenum err = glGetError();
	if( err != GL_NO_ERROR ) {
		//DebugBreak();
	}
}

/* Class for attaching an OpenGL context to an HWND (Windows) or LPane (Mac).
example:
void MyDrawingFunc(HWND myWindow) {

	LightweightOGL ogl;
	ogl.wnd = myWindow;
	if( ogl.BeginFrame() ) {
		glClear(GL_COLOR_BUFFER_BIT);
		ogl.SwapBuffers();
		ogl.EndFrame();
	}
	ogl.EndOGL();
}
*/
class LightweightOGL {
public:
	LightweightOGL();
	~LightweightOGL();
	
	void EndOGL();
	bool BeginFrame();
	void EndFrame();
	void SwapBuffers();
	
#if GS_WIN
public:
	HWND wnd;

private:
	HDC dc;
	HGLRC rc;
	
#else // mac
public:
	// set lPane xor wnd
	void * lPane;
	GrafPtr graf;

	GLint bufferRect[4];

private:
	AGLContext rc;
#endif
};





LightweightOGL::LightweightOGL() 
/* Initialization.  This class provides a simple interface for rendering OpenGL graphics
		on a Mac or Windows native window.
	return -
*/
{
#if GS_WIN
	wnd = NULL;
	dc = NULL;
	rc = NULL;
#else // mac
	lPane = NULL;
	graf = NULL;
	memset(bufferRect, 0, sizeof(GLint)*4);
	rc = NULL;
#endif
}

LightweightOGL::~LightweightOGL() 
/* Detaches from the window.
	return -
*/
{
	EndOGL();
}

void LightweightOGL::EndOGL() 
/* Detaches from the native window.  Should be called when OpenGL drawing is no longer needed by the app.
		The Destructor automatically calls this.  multiple calls to this will do nothing.
	return -
*/
{
#if GS_WIN
	if( rc != NULL ) {
		wglDeleteContext(rc);
		rc = NULL;
	}

	if( dc != NULL ) {
		// assert wnd != NULL
		ReleaseDC(wnd, dc);
		dc = NULL;
	}

	wnd = NULL;
	
#else // mac
	if ( rc != NULL ) {
		QDFlushPortBuffer(graf, 0);

		aglDestroyContext((AGLContext)rc);  
		rc = NULL;
	}
	
	lPane = NULL;
	graf = NULL;	
	memset(bufferRect, 0, sizeof(GLint)*4);

#endif

}

bool LightweightOGL::BeginFrame() 
/* Associates OpenGL with 'this->wnd'.  After a successful return, OpenGL calls will be directed
		into the window.
	return - true on success; otherwise false
*/
{
#if GS_WIN

	enum LocalError {
		InvalidHWND,
		DCRetrievalFailure,
		PixelFormatError,
		RCCreationFailure
	};

	try {
		if( wnd == NULL ) throw InvalidHWND;

		if( rc==NULL && dc==NULL ) {
			dc = GetDC(wnd);
			if( dc==NULL ) throw DCRetrievalFailure;
		
			PIXELFORMATDESCRIPTOR pfd = {
				sizeof(PIXELFORMATDESCRIPTOR),	// size of this pfd
				1,								// version number
				PFD_DRAW_TO_WINDOW |			// Draw to  window
				PFD_SUPPORT_OPENGL | 			// Support OpenGL calls in window
				PFD_DOUBLEBUFFER,				// Support double buffering 
				PFD_TYPE_RGBA,					// RGBA type
				32,								// color depth 
				0, 0, 0, 0, 0, 0,				// color bits ignored
				8,								// no alpha buffer
				0,								// shift bit ignored
				0,								// no accumulation buffer
				0, 0, 0, 0, 					// accum bits ignored
				0,								// 32-bit z-buffer	
				32,								// no stencil buffer
				0,								// no auxiliary buffer
				PFD_MAIN_PLANE,					// main layer
				0,								// reserved
				0, 0, 0							// layer masks ignored
			};

			// Choose a pixel format that best matches that described in pfd
			int pixelFormat = ChoosePixelFormat(dc, &pfd);

			DescribePixelFormat(dc, pixelFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

			// Set the pixel format for the device context
			BOOL success = SetPixelFormat(dc, pixelFormat, &pfd);
			if( !success ) throw PixelFormatError;

			rc = wglCreateContext(dc);
			if( rc == NULL ) throw RCCreationFailure;
		}

		wglMakeCurrent(dc, rc);

	} catch( LocalError /*localError*/ ) {
		EndOGL();
		return false;
	}

	return true;

#else // mac

	enum LocalError {
		InvalidWindowIdentifier,
		PixelFormatError,
		CreateContextFailed,
		InvalidGrafPtr
	};
	
	bool success = true;
	
	AGLPixelFormat pixelFormat = NULL;
	
	try {	
		
		if( lPane != NULL ) {
			graf = GrafPtrFromLPane(lPane);
			if( graf == NULL ) throw InvalidGrafPtr;

			// compute the bufferRect if it's not specified already
			if( bufferRect[2]==0 ) {
				Rect windowRect;
				GetPortBounds(graf, &windowRect);
			
				Rect rect = PaneRect(lPane);

				bufferRect[0] = rect.left;
				bufferRect[1] = windowRect.bottom-rect.bottom;
				bufferRect[2] = rect.right-rect.left;
				bufferRect[3] = rect.bottom-rect.top;
			}
		} else if( graf != NULL ) {

			// compute the bufferRect if it's not specified already
			if( bufferRect[2]==0 ) {
				Rect windowRect;
				GetPortBounds(graf, &windowRect);
			
				bufferRect[0] = 0;
				bufferRect[1] = 0;
				bufferRect[2] = windowRect.right-windowRect.left;
				bufferRect[3] = windowRect.bottom-windowRect.top;
			}
		} else {
			throw InvalidWindowIdentifier;
		}

		if( rc == NULL ) {
			// attirbutes
			unsigned long bitDepth = 32;
			GLint attrib[] = { AGL_RGBA, AGL_PIXEL_SIZE, bitDepth, AGL_RED_SIZE, 8, AGL_ALPHA_SIZE, 8, AGL_DEPTH_SIZE, 32, AGL_DOUBLEBUFFER,AGL_CLOSEST_POLICY, AGL_NO_RECOVERY, AGL_NONE };
			
			pixelFormat = aglChoosePixelFormat(NULL, 0, attrib); 
			if( pixelFormat == NULL ) throw PixelFormatError;
				
			rc = aglCreateContext(pixelFormat, NULL);	
			if( rc == NULL ) throw CreateContextFailed;
			
			if( aglEnable((AGLContext)rc, AGL_BUFFER_RECT) ) {
				aglSetInteger((AGLContext)rc, AGL_BUFFER_RECT, bufferRect);				
			}
			
			aglSetDrawable((AGLContext)rc, (AGLDrawable)graf);
					
		}
			
		// Set the current RC - the one that receives gl calls
		aglSetCurrentContext( (AGLContext)rc );
		
	} catch( LocalError /*localError*/ ) {
		success = false;
	}

	if( pixelFormat != NULL ) {
		aglDestroyPixelFormat(pixelFormat);
	}
		
	return success;	
#endif
}

void LightweightOGL::SwapBuffers() 
/* Swaps the back and front buffers so as to make the contents of the backbuffer visible.
		The OpenGL and window association should already be created - i.e. a successful call 
		to BeginFrame( ) should preceed this.
	return -
*/
{
#if GS_WIN
	::SwapBuffers(dc);
#else // mac
	aglSwapBuffers((AGLContext)rc);
#endif
}

void LightweightOGL::EndFrame() {
#if GS_WIN
	wglMakeCurrent(0, NULL);
#endif
}




const int WM_IDLE = WM_USER + 1;


int APIENTRY WinMain(HINSTANCE hInstance,
                     HINSTANCE /*hPrevInstance*/,
                     LPSTR    /*lpCmdLine*/,
                     int       nCmdShow)
{
	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	MSG msg;
	msg.wParam = 0;
	HACCEL hAccelTable;

	// Initialize global strings
	LoadString(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadString(hInstance, IDC_OPENGLSAMPLE, szWindowClass, MAX_LOADSTRING);
	MyRegisterClass(hInstance);

	// Perform application initialization:
	if (!InitInstance (hInstance, nCmdShow))
	{
		return FALSE;
	}

	hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_OPENGLSAMPLE));

	// Main message loop:
	bool quit = false;
	while( !quit ) {
		bool moreIdle = true;

		while( moreIdle ) {
			SendMessage(gMainWindow, WM_IDLE, (WPARAM)&moreIdle, 0);
			
			if( PeekMessage(&msg, 0, 0, 0, PM_NOREMOVE) ) {
				break;
			}
		};

		if( GetMessage(&msg, NULL, 0, 0) <= 0 ) {
			quit = true;
		} else {
			if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}

		int messageCount = 1;
		while( PeekMessage(&msg, 0, 0, 0, PM_NOREMOVE) ) {
			if( GetMessage(&msg, NULL, 0, 0) <= 0 ) {
				quit = true;
			} else {
				if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
				{
					TranslateMessage(&msg);
					DispatchMessage(&msg);
				}
			}
			if( messageCount-- == 0 ) { // if many messages have been processed; don't starve the idle loop
				break;
			}
		}

	}

	GdiplusShutdown(gdiplusToken);
   
	return (int) msg.wParam;
}



ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style			= CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS;
	wcex.lpfnWndProc	= WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= LoadIcon(hInstance, MAKEINTRESOURCE(IDI_OPENGLSAMPLE));
	wcex.hCursor		= LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground	= (HBRUSH)(COLOR_WINDOW+1);
	wcex.lpszMenuName	= 0;
	wcex.lpszClassName	= szWindowClass;
	wcex.hIconSm		= LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}

BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   HWND & hWnd = gMainWindow;

   hInst = hInstance; // Store instance handle in our global variable

   hWnd = CreateWindow(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW ,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, NULL, NULL, hInstance, NULL);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}


LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	static LightweightOGL ogl;
	static bool lButtonDown = false;
	static int gMouseX, gMouseY, gLastMouseX, gLastMouseY;
	static bool rButtonDown = false;
	static bool mButtonDown = false;
	
	RECT rect;
	GetClientRect(hWnd, &rect);

	PAINTSTRUCT ps;
	HDC hdc;

	static double theta = 0.7f;
	static double phi = 1.2;
			
	switch( message ) {
		case WM_KEYDOWN: {
			if( wParam == VK_UP ) {
				gCurrentDrawTest++;
			}
			else if( wParam == VK_DOWN ) {
				gCurrentDrawTest = max(0, gCurrentDrawTest-1);
			}
			else if ( wParam==VK_ESCAPE ) {
				PostQuitMessage(0);
			}
			break;
		}
		case WM_LBUTTONDOWN: {
			gCurrentDrawTest++;
			break;
		} 
		case WM_RBUTTONDOWN: {
			gCurrentDrawTest = max(0, gCurrentDrawTest-1);
			break;
		}
		case WM_CREATE: {
			gInitialTick = GetTickCount();
			break;
		}

		case WM_SIZE: {
			gWindowWidth = LOWORD(lParam);
			gWindowHeight = HIWORD(lParam);
			break;
		}

		case WM_PAINT: {
			hdc = BeginPaint(hWnd, &ps);

			if( ogl.wnd == NULL ) { // if uninitialized
				ogl.wnd = hWnd;
				ogl.BeginFrame();

				glewInit();
			}

			if( ogl.BeginFrame() ) { // if ready to draw in OpenGL

				GL_Draw();

				ogl.SwapBuffers();
				ogl.EndFrame();
			}

			EndPaint(hWnd, &ps);
			InvalidateRect(hWnd, 0, false);
			break;
		}
		case WM_DESTROY:
			if( ogl.BeginFrame() ) {
				ogl.EndFrame();
			}

			PostQuitMessage(0);
			break;
		case WM_ERASEBKGND:
			return 1;
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
	}

	return 0;
}
