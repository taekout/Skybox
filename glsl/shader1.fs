#version 330

uniform sampler2D kImage;  // k for constant

in vec2 gUV;  // g for global.

layout(location=0) out vec4 outColor;  // out for output

void main() {
	outColor.rgba = texture(kImage, gUV);
}
