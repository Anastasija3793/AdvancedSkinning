#version 330                                                                                                                                                      

//----------------------------------------------------------------------------------------------
// @file lbs_skinning_vs.glsl
// @brief lbs_skinning_vs is a vertex shader used for implementing Linear Blend Skinning (LBS)
// @author Anastasija Belaka
// @date 18/01/2019 Updated to NCCA Coding standard
//----------------------------------------------------------------------------------------------
layout (location = 0) in vec3 Position;                                             
layout (location = 1) in vec2 TexCoord;                                             
layout (location = 2) in vec3 Normal;                                               
layout (location = 3) in ivec4 BoneIDs;
layout (location = 4) in vec4 Weights;

out vec2 TexCoord0;
out vec3 Normal0;                                                                   
out vec3 WorldPos0;                                                                 

const int MAX_BONES = 100;

uniform mat4 gWVP;
uniform mat4 gWorld;
uniform mat4 gBones[MAX_BONES];

vec4 defV = vec4(0.0,0.0,0.0,0.0);
vec4 defN = vec4(0.0,0.0,0.0,0.0);

void main()
{
	for(int i = 0; i<3; i++)
	{
		mat4 bone = gBones[BoneIDs[i]];
		defV = defV + Weights[i] * (bone*vec4(Position,1.0));
		defN = defN + Weights[i] * (bone*vec4(Normal,1.0));
	}

	gl_Position  = gWVP * defV;
    TexCoord0    = TexCoord;
    vec4 NormalL = defN;
    Normal0      = (gWorld * NormalL).xyz;
    WorldPos0    = (gWorld * gl_Position).xyz;
}