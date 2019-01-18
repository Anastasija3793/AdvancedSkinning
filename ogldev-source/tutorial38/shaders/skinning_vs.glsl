#version 330                                                                        
                                                                                    
//----------------------------------------------------------------------------------------------
// @file skinning_vs.glsl
// @brief skinning_vs is a vertex shader used for implementing the proposed method
// @author Anastasija Belaka
// @date 18/01/2019 Updated to NCCA Coding standard
//----------------------------------------------------------------------------------------------
layout (location = 0) in vec3 Position;                                             
layout (location = 1) in vec2 TexCoord;                                             
layout (location = 2) in vec3 Normal;                                               
layout (location = 3) in ivec4 BoneIDs;
layout (location = 4) in vec4 Weights;
layout (location = 5) in vec3 Pi;

out vec2 TexCoord0;
out vec3 Normal0;                                                                   
out vec3 WorldPos0;                                                                 

const int MAX_BONES = 100;

uniform mat4 gWVP;
uniform mat4 gWorld;
uniform mat4 gBones[MAX_BONES];

mat3 b = mat3(gBones[BoneIDs[0]]);
vec4 pivotQ = vec4(b[0][0],b[0][1],b[0][2],b[1][0]);

vec3 defV = vec3(0.0,0.0,0.0);
vec3 defN = vec3(0.0,0.0,0.0);
vec4 Q = vec4(0.0,0.0,0.0,0.0);
vec3 transV = vec3(0.0,0.0,0.0);

void main()
{
	for(int i = 0; i<3; i++)
	{
		mat3 bone = mat3(gBones[BoneIDs[i]]);
		vec4 quat = vec4(bone[0][0],bone[0][1],bone[0][2],bone[1][0]);

		// normalize before conversion
		float n = sqrt(quat.x*quat.x + quat.y*quat.y + quat.z*quat.z + quat.w*quat.w);
		quat = quat/n;

		// quatToMat
		mat3 rot = mat3(
			1.0 - (2.0 * quat.y * quat.y) - (2.0 * quat.z * quat.z),
		  (2.0 * quat.x * quat.y) + (2.0 * quat.z * quat.w),
		  (2.0 * quat.x * quat.z) - (2.0 * quat.y * quat.w),

		  (2.0 * quat.x * quat.y) - (2.0 * quat.z * quat.w),
		  1.0 - (2.0 * quat.x * quat.x) - (2.0 * quat.z * quat.z),
		  (2.0 * quat.y * quat.z) + (2.0 * quat.x * quat.w),

		  (2.0 * quat.x * quat.z) + (2.0 * quat.y * quat.w),
		  (2.0 * quat.y * quat.z) - (2.0 * quat.x * quat.w),
		  1.0 - (2.0 * quat.x * quat.x) - (2.0 * quat.y * quat.y)
		  );
		  vec3 tra = vec3(bone[2][0],bone[2][1],bone[2][2]);

		  transV = transV + Weights[i] * (rot*Pi + tra);
		  defN = defN + Weights[i] * (rot*Normal);

		  if(i!=0)
		  {
			if(dot(pivotQ,quat) < 0)
			{
				quat=(-1)*quat;
			}
		  }
		  Q = Q + Weights[i] * quat;
	}

	// normalize Q
	float Qnorm = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z + Q.w*Q.w);
	Q = Q/Qnorm;

	mat3 QnMat = mat3(
			1.0 - (2.0 * Q.y * Q.y) - (2.0 * Q.z * Q.z),
		  (2.0 * Q.x * Q.y) + (2.0 * Q.z * Q.w),
		  (2.0 * Q.x * Q.z) - (2.0 * Q.y * Q.w),

		  (2.0 * Q.x * Q.y) - (2.0 * Q.z * Q.w),
		  1.0 - (2.0 * Q.x * Q.x) - (2.0 * Q.z * Q.z),
		  (2.0 * Q.y * Q.z) + (2.0 * Q.x * Q.w),

		  (2.0 * Q.x * Q.z) + (2.0 * Q.y * Q.w),
		  (2.0 * Q.y * Q.z) - (2.0 * Q.x * Q.w),
		  1.0 - (2.0 * Q.x * Q.x) - (2.0 * Q.y * Q.y)
		  );

	vec3 rotV = QnMat * (Position-Pi);
	defV = rotV + transV;
	

    gl_Position  = gWVP * vec4(defV,1.0);
    TexCoord0    = TexCoord;
    vec4 NormalL = mat4(QnMat) * vec4(Normal, 0.0);
    Normal0      = (gWorld * NormalL).xyz;
    WorldPos0    = (gWorld * gl_Position).xyz;

}
