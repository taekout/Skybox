// This file shows you how to do OpenGL 3.2 Core drawing.  
// Press the up or down arrow on the keyboard (or left click, right click) to cycle through the draw functions
// You can change the shader file (shader0.vs, shader1.fs, etc) while the program is running to see your changes.
#include "sample.h"
#include <math.h>
#include "utilities.h"
#include "VWMath.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// some drawing functions that do different things
void MyDraw0();
void MyDraw1();
void MyDraw2();


void GL_Draw() 
// press the up or down arrow on the keyboard to cycle through the draw functions
{
	typedef void (* DrawingFunction)();
	DrawingFunction drawingFunctions[] = { MyDraw2, MyDraw1, MyDraw0 };

	int numDrawingFunctions = sizeof(drawingFunctions) / sizeof(DrawingFunction);
	int drawingFunction = min(numDrawingFunctions-1, max(APP_GetCurrentDrawTest(), 0));  // don't go out of bounds

	drawingFunctions[drawingFunction]();  // execute a draw function
}


void MyDraw0() 
// This drawing function is very simple.  It draws a red rectangle.
{
	GLuint program = 0;
	GLuint inPosition = 0;
	GLuint kColor = 0;
	if( !GL_LoadGLSLProgramFromFile("shader0.vs", "shader0.fs", &program,
			eGLSLBindingAttribute, "inPosition", &inPosition,
			eGLSLBindingUniform, "kColor", &kColor,
			eGLSLBindingEnd) )
	// if the glsl shader successfully compiled and linked then we will use it
	{

		glViewport(0, 0, WIN_GetWidth(), WIN_GetHeight() );  // where are we drawing in the window (x, y, width, height) of pixel units.  (0, 0) is lower left corner.
		glClear(GL_COLOR_BUFFER_BIT);  // the color buffer will contain the output from the fragment shader

		// OpenGL's 'native' coordinate system, called normalized device coordinates, ranges from -1 on the left, bottom and near to +1 on the right, top and far.
		// Ant geometry outside of this range will not be visible.
		// Usually your vertices are in a local 3D space that you transform into global world space by multiplying your vertex by a local-to-world matrix, then into 
		// eye space (the eye of the viewer or camera) by multiplying by a world-to-eye matrix, next into clip coordinates by multiplying by a projection matrix (this
		// is the space the vertex shader output gl_Position should be in) and finally opengl will transform that coordinate into normalized device coordinates by
		// dividing by the coordinate's w.

		float vertices[] = {  // 4 vertices of a rectangle (z component is 0, w component is 1).
			-1, -1, 0, 1,
			0,  -1, 0, 1,
			-1,  0, 0, 1,
			0,   0, 0, 1,
		};

		// allocate a block of graphics memory (Vertex Buffer Object)
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		int numBytes = sizeof(vertices);
		glBufferData(GL_ARRAY_BUFFER, numBytes, vertices, GL_STATIC_DRAW);  // set the size of the memory block, also upload some data

		// make a structure (Vertex Array) that sets up the rendering state
		GLuint va = 0;
		glGenVertexArrays(1, &va);
		glBindVertexArray(va);
				
		glEnableVertexAttribArray(inPosition);  // this is an identifier to the vertex shader variable 'inPosition'
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(inPosition, 4, GL_FLOAT, GL_FALSE, 0, 0);  // specify that the data for 'inPosition' comes from offset 0 in 'vbo' and that it is 4 tightly packed floats
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);


		// Now we're ready to use the shader program
		glUseProgram(program);
		glBindVertexArray(va);
		glUniform4f(kColor, 1, 0, 0, 1);  // this sets the fragment shader variable 'kColor' to red (red=1, green=0, blue=0, alpha=1)
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);  // this draws a 'triangle strip' (search google images) that ends up drawing a red rectangle


		// Done.  Free resources.
		glUseProgram(0);  
		glDeleteProgram(program);
		glDeleteBuffers(1, &vbo);
		glDeleteVertexArrays(1, &va);
	}
}

void MyDraw1()
	// This function draws a textured image.
{

	GLuint program = 0;
	GLuint inPosition = 0, inUVs = 0;
	GLuint kImage = 0;
	if( !GL_LoadGLSLProgramFromFile("shader1.vs", "shader1.fs", &program,
		eGLSLBindingAttribute, "inPosition", &inPosition,
		eGLSLBindingAttribute, "inUVs", &inUVs,
		eGLSLBindingUniform, "kImage", &kImage,
		eGLSLBindingEnd) )
		// if the glsl shader successfully compiled and linked then we will use it
	{
		GL_LoadTextureImage("image1.jpg", kImage);

		glViewport(0, 0, WIN_GetWidth(), WIN_GetHeight() );  // where are we drawing in the window (x, y, width, height) of pixel units.  (0, 0) is lower left corner.
		glClear(GL_COLOR_BUFFER_BIT);  // the color buffer will contain the output from the fragment shader

		// OpenGL's 'native' coordinate system, called normalized device coordinates, ranges from -1 on the left, bottom and near to +1 on the right, top and far.
		// Ant geometry outside of this range will not be visible.
		// Usually your vertices are in a local 3D space that you transform into global world space by multiplying your vertex by a local-to-world matrix, then into 
		// eye space (the eye of the viewer or camera) by multiplying by a world-to-eye matrix, next into clip coordinates by multiplying by a projection matrix (this
		// is the space the vertex shader output gl_Position should be in) and finally opengl will transform that coordinate into normalized device coordinates by
		// dividing by the coordinate's w.

		float vertices[] = {  // 4 vertices of a rectangle (z component is 0, w component is 1).
			-1, -1, 0, 1,
			0,  -1, 0, 1,
			-1,  0, 0, 1,
			0,   0, 0, 1,
		};

		// allocate a block of graphics memory (Vertex Buffer Object)
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		int numBytes = sizeof(vertices);
		glBufferData(GL_ARRAY_BUFFER, numBytes, vertices, GL_STATIC_DRAW);  // set the size of the memory block, also upload some data
		glBindBuffer(GL_ARRAY_BUFFER, 0);


		// make a structure (Vertex Array) that sets up the rendering state
		GLuint va = 0;
		glGenVertexArrays(1, &va);
		glBindVertexArray(va);

		glEnableVertexAttribArray(inPosition);  // this is an identifier to the vertex shader variable 'inPosition'
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(inPosition, 4, GL_FLOAT, GL_FALSE, 0, 0);  // specify that the data for 'inPosition' comes from offset 0 in 'vbo' and that it is 4 tightly packed floats
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);


		// Now we're ready to use the shader program
		glUseProgram(program);
		glDisable(GL_CULL_FACE);
		glBindVertexArray(va);

		int textureUnit = 0;
		glActiveTexture(GL_TEXTURE0+textureUnit);  // this line and the one below associate kImage with GL_TEXTURE0
		glBindTexture(GL_TEXTURE_2D, kImage);
		glUniform1i(kImage, textureUnit);  // this sets the fragment shader variable 'kImage' to the value 0 - this means use the texture in GL_TEXTURE0

		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);  // this draws a 'triangle strip' (search google images) that ends up drawing a textured rectangle


		// Done.  Free resources.
		glUseProgram(0);  
		glDeleteProgram(program);
		glDeleteBuffers(1, &vbo);
		glDeleteVertexArrays(1, &va);
		glDeleteTextures(1, &kImage);
	}
}


void MyDraw2()
// This function draws a textured image.
{

	GLuint program = 0;
	GLuint inPosition = 0;
	GLuint inUV = 0;
	GLuint uProj = 0, uView = 0;
	GLuint kImage = 0;
	GLuint kTexID = 0;
	if( !GL_LoadGLSLProgramFromFile("Skybox.vert", "Skybox.frag", &program,
		eGLSLBindingAttribute, "inPosition", &inPosition,
		eGLSLBindingAttribute, "inUVs", &inUV,
		eGLSLBindingUniform, "kImage", &kImage,
		eGLSLBindingUniform, "Proj", &uProj,
		eGLSLBindingUniform, "View", &uView,
		eGLSLBindingEnd) )
		// if the glsl shader successfully compiled and linked then we will use it
	{
		GL_LoadTextureImage("image1.jpg", kTexID);

		glViewport(0, 0, WIN_GetWidth(), WIN_GetHeight() );  // where are we drawing in the window (x, y, width, height) of pixel units.  (0, 0) is lower left corner.
		glClear(GL_COLOR_BUFFER_BIT);  // the color buffer will contain the output from the fragment shader

		// OpenGL's 'native' coordinate system, called normalized device coordinates, ranges from -1 on the left, bottom and near to +1 on the right, top and far.
		// Ant geometry outside of this range will not be visible.
		// Usually your vertices are in a local 3D space that you transform into global world space by multiplying your vertex by a local-to-world matrix, then into 
		// eye space (the eye of the viewer or camera) by multiplying by a world-to-eye matrix, next into clip coordinates by multiplying by a projection matrix (this
		// is the space the vertex shader output gl_Position should be in) and finally opengl will transform that coordinate into normalized device coordinates by
		// dividing by the coordinate's w.

		float vertices[] = { //interleaved vertex(3) , uv(2)
			-1.0f,  1.0f, -1.0f, 0, 1,
			-1.0f, -1.0f, -1.0f, 0, 0,
			1.0f, -1.0f, -1.0f,  1, 0,
			1.0f, -1.0f, -1.0f,  1, 0,
			1.0f,  1.0f, -1.0f,  1, 1,
			-1.0f,  1.0f, -1.0f, 0, 1,

			-1.0f, -1.0f,  1.0f, 0, 1,
			-1.0f, -1.0f, -1.0f, 0, 0,
			-1.0f,  1.0f, -1.0f, 1, 0,
			-1.0f,  1.0f, -1.0f, 1, 0,
			-1.0f,  1.0f,  1.0f, 1, 1,
			-1.0f, -1.0f,  1.0f, 0, 1,

			1.0f, -1.0f, -1.0f, 0, 0,
			1.0f, -1.0f,  1.0f, 0, 1,
			1.0f,  1.0f,  1.0f, 1, 1,
			1.0f,  1.0f,  1.0f, 1, 1,
			1.0f,  1.0f, -1.0f, 1, 0,
			1.0f, -1.0f, -1.0f, 0, 0,

			-1.0f, -1.0f,  1.0f, 0, 0,
			-1.0f,  1.0f,  1.0f, 0, 1,
			1.0f,  1.0f,  1.0f,  1, 1,
			1.0f,  1.0f,  1.0f,  1, 1,
			1.0f, -1.0f,  1.0f,  1, 0,
			-1.0f, -1.0f,  1.0f, 0, 0,

			-1.0f,  1.0f, -1.0f, 0, 0,
			1.0f,  1.0f, -1.0f,  1, 0,
			1.0f,  1.0f,  1.0f,  1, 1,
			1.0f,  1.0f,  1.0f,  1, 1,
			-1.0f,  1.0f,  1.0f, 0, 1,
			-1.0f,  1.0f, -1.0f, 0, 0,

			-1.0f, -1.0f, -1.0f, 0, 0,
			-1.0f, -1.0f,  1.0f, 0, 1,
			1.0f, -1.0f, -1.0f,  1, 0,
			1.0f, -1.0f, -1.0f,  1, 0,
			-1.0f, -1.0f,  1.0f, 0, 1,
			1.0f, -1.0f,  1.0f,  1, 1
		};


		// allocate a block of graphics memory (Vertex Buffer Object)
		GLuint vertVBO;
		glGenBuffers(1, &vertVBO);
		glBindBuffer(GL_ARRAY_BUFFER, vertVBO);
		int numBytes = sizeof(vertices);
		glBufferData(GL_ARRAY_BUFFER, numBytes, vertices, GL_STATIC_DRAW);  // set the size of the memory block, also upload some data
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// make a structure (Vertex Array) that sets up the rendering state
		GLuint va = 0;
		glGenVertexArrays(1, &va);
		glBindVertexArray(va);

		glBindBuffer(GL_ARRAY_BUFFER, vertVBO);
		glEnableVertexAttribArray(inPosition);  // this is an identifier to the vertex shader variable 'inPosition'
		glVertexAttribPointer(inPosition, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 5, 0);  // specify that the data for 'inPosition' comes from offset 0 in 'vbo' and that it is 4 tightly packed floats

		glEnableVertexAttribArray(inUV);  // this is an identifier to the vertex shader variable 'inPosition'
		glVertexAttribPointer(inUV, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 5, (void *) (sizeof(float) * 3) );  // specify that the data for 'inPosition' comes from offset 0 in 'vbo' and that it is 4 tightly packed floats
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);


		// Now we're ready to use the shader program
		glUseProgram(program);
		glBindVertexArray(va);

		int textureUnit = 0;
		glActiveTexture(GL_TEXTURE0+textureUnit);  // this line and the one below associate kImage with GL_TEXTURE0
		glBindTexture(GL_TEXTURE_2D, kTexID);
		glUniform1i(kImage, textureUnit);  // this sets the fragment shader variable 'kImage' to the value 0 - this means use the texture in GL_TEXTURE0

		// set view / proj matrix.
		glm::mat4 projMat = glm::perspective(40.0f, 1.0f, 0.0f, 10.0f);
		glm::mat4 viewMat = glm::lookAt(glm::vec3(0,0,5), glm::vec3(0,0,0), glm::vec3(0,1,0));
		glUniformMatrix4fv(uProj, 1, false, &projMat[0][0]);
		glUniformMatrix4fv(uView, 1, false, &viewMat[0][0]);

		glDrawArrays(GL_TRIANGLES, 0, 36);  // this draws a 'triangle strip' (search google images) that ends up drawing a textured rectangle


		// Done.  Free resources.
		glUseProgram(0);  
		glDeleteProgram(program);
		glDeleteBuffers(1, &vertVBO);
		glDeleteVertexArrays(1, &va);
		glDeleteTextures(1, &kImage);
	}
}

