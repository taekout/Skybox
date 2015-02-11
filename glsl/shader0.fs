#version 330

uniform vec4 kColor;

layout(location=0) out vec4 outColor;

void main() {
	outColor.rgba = kColor;
    // outColor.g = 1;  // Run the demo and while its running uncomment this line to have the rectangle become yellow
}
