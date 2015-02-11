#version 330

layout(location=0) in vec4 inPosition;  // in for input
layout(location=1) in vec2 inUVs;  // in for input

out vec2 gUV;  // g for global.

void main() {
	gUV = inPosition.xy + vec2(1, 1);  // gUV's coordinates should range from 0 to 1
	gUV.y = 1 - gUV.y;  // invert the y component to flip the image vertically

	gl_Position = inPosition;
}