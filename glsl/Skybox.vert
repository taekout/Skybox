#version 330


layout(location=0) in vec4 inPosition;  // in for input
//attribute vec3 a_Color;
layout(location=1) in vec2 inUVs;
//attribute vec2 inUV;

uniform mat4 Proj;
uniform mat4 View;

out vec2 gUV;  // g for global.
//varying vec3 v_Color;

void main () {
	gUV = inUVs;
	gl_Position = Proj * View * inPosition;
}


/*
#version 120

attribute vec4 a_Vertex;
attribute vec4 a_Color;

uniform mat4 u_Proj;
uniform mat4 u_View;

varying vec4 v_Color;


void main()
{
    vec4 vert = u_View*a_Vertex;
    gl_Position = u_Proj*vert;
    v_Color = a_Color;
}

*/
