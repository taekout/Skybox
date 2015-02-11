#version 330

in vec2 gUV;
uniform sampler2D kImage;

layout(location=0) out vec4 outColor;  // out for output

void main () {
	//gl_FragColor = texture2D(cube_texture, v_UV);
	outColor.rgba = texture2D(kImage, gUV);
	//gl_FragColor.xyz = v_UV;
	//gl_FragColor.w = 1;
	//gl_FragColor.xyz = vec3(0,1,0);
}

/*

#version 120

varying vec4 v_Color;


void main()
{
	gl_FragColor = v_Color;
}


*/