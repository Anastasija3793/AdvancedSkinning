//----------------------------------------------------------------------------------------------
// @file ogldev_skinned_mesh.cpp
// @brief ogldev_skinned_mesh used for finding CoR, skinning, model animation
// @author of initial code - Etay Meiri
// @author Modified by Anastasija Belaka
// @date 18/01/2019 Updated to NCCA Coding standard
// note: the code belonging to Anastasija Belaka outlined/cited bellow in the code
//----------------------------------------------------------------------------------------------



#include "ogldev_skinned_mesh.h"
#include <iostream>
#include <cmath>

#include <limits.h>
#include <string>

#include "ogldev_util.h"

#include "technique.h"
#include "ogldev_math_3d.h"

#define POSITION_LOCATION    0
#define TEX_COORD_LOCATION   1
#define NORMAL_LOCATION      2
#define BONE_ID_LOCATION     3
#define BONE_WEIGHT_LOCATION 4
#define PI_LOCATION			 5

void SkinnedMesh::VertexBoneData::AddBoneData(uint BoneID, float Weight)
{
    for (uint i = 0 ; i < ARRAY_SIZE_IN_ELEMENTS(IDs) ; i++) {
        if (Weights[i] == 0.0) {
            IDs[i]     = BoneID;
            Weights[i] = Weight;
            return;
        }        
    }
    
    // should never get here - more bones than we have space for
    assert(0);
}

SkinnedMesh::SkinnedMesh()
{
    m_VAO = 0;
    ZERO_MEM(m_Buffers);
    m_NumBones = 0;
    m_pScene = NULL;
}


SkinnedMesh::~SkinnedMesh()
{
    Clear();
}


void SkinnedMesh::Clear()
{
    for (uint i = 0 ; i < m_Textures.size() ; i++) {
        SAFE_DELETE(m_Textures[i]);
    }

    if (m_Buffers[0] != 0) {
        glDeleteBuffers(ARRAY_SIZE_IN_ELEMENTS(m_Buffers), m_Buffers);
    }
       
    if (m_VAO != 0) {
        glDeleteVertexArrays(1, &m_VAO);
        m_VAO = 0;
    }
}


bool SkinnedMesh::LoadMesh(const string& Filename)
{
    // Release the previously loaded mesh (if it exists)
    Clear();
 
    // Create the VAO
    glGenVertexArrays(1, &m_VAO);   
    glBindVertexArray(m_VAO);
    
    // Create the buffers for the vertices attributes
    glGenBuffers(ARRAY_SIZE_IN_ELEMENTS(m_Buffers), m_Buffers);

    bool Ret = false;    
  
    m_pScene = m_Importer.ReadFile(Filename.c_str(), ASSIMP_LOAD_FLAGS);
    
    if (m_pScene) {  
        m_GlobalInverseTransform = m_pScene->mRootNode->mTransformation;
        m_GlobalInverseTransform.Inverse();
        Ret = InitFromScene(m_pScene, Filename);
    }
    else {
        printf("Error parsing '%s': '%s'\n", Filename.c_str(), m_Importer.GetErrorString());
    }

    // Make sure the VAO is not changed from the outside
    glBindVertexArray(0);	

    return Ret;
}


bool SkinnedMesh::InitFromScene(const aiScene* pScene, const string& Filename)
{  
    m_Entries.resize(pScene->mNumMeshes);
    m_Textures.resize(pScene->mNumMaterials);

    vector<Vector3f> Positions;
    vector<Vector3f> Normals;
    vector<Vector2f> TexCoords;
    vector<VertexBoneData> Bones;
    vector<uint> Indices;
	vector<Vector3f> Pi;
       
    uint NumVertices = 0;
    uint NumIndices = 0;
    
    // Count the number of vertices and indices
    for (uint i = 0 ; i < m_Entries.size() ; i++) {
        m_Entries[i].MaterialIndex = pScene->mMeshes[i]->mMaterialIndex;        
        m_Entries[i].NumIndices    = pScene->mMeshes[i]->mNumFaces * 3;
        m_Entries[i].BaseVertex    = NumVertices;
        m_Entries[i].BaseIndex     = NumIndices;
        
        NumVertices += pScene->mMeshes[i]->mNumVertices;
        NumIndices  += m_Entries[i].NumIndices;
    }
    
    // Reserve space in the vectors for the vertex attributes and indices
    Positions.reserve(NumVertices);
    Normals.reserve(NumVertices);
    TexCoords.reserve(NumVertices);
    Bones.resize(NumVertices);
    Indices.reserve(NumIndices);
	Pi.resize(NumVertices,Vector3f(0.0));
        
    // Initialize the meshes in the scene one by one
    for (uint i = 0 ; i < m_Entries.size() ; i++) {
        const aiMesh* paiMesh = pScene->mMeshes[i];
        InitMesh(i, paiMesh, Positions, Normals, TexCoords, Bones, Indices, Pi);
    }

    if (!InitMaterials(pScene, Filename)) {
        return false;
    }

    // Generate and populate the buffers with vertex attributes and the indices
  	glBindBuffer(GL_ARRAY_BUFFER, m_Buffers[POS_VB]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Positions[0]) * Positions.size(), &Positions[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(POSITION_LOCATION);
    glVertexAttribPointer(POSITION_LOCATION, 3, GL_FLOAT, GL_FALSE, 0, 0);    

    glBindBuffer(GL_ARRAY_BUFFER, m_Buffers[TEXCOORD_VB]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(TexCoords[0]) * TexCoords.size(), &TexCoords[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(TEX_COORD_LOCATION);
    glVertexAttribPointer(TEX_COORD_LOCATION, 2, GL_FLOAT, GL_FALSE, 0, 0);

   	glBindBuffer(GL_ARRAY_BUFFER, m_Buffers[NORMAL_VB]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Normals[0]) * Normals.size(), &Normals[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(NORMAL_LOCATION);
    glVertexAttribPointer(NORMAL_LOCATION, 3, GL_FLOAT, GL_FALSE, 0, 0);

   	glBindBuffer(GL_ARRAY_BUFFER, m_Buffers[BONE_VB]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Bones[0]) * Bones.size(), &Bones[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(BONE_ID_LOCATION);
    glVertexAttribIPointer(BONE_ID_LOCATION, 4, GL_INT, sizeof(VertexBoneData), (const GLvoid*)0);
    glEnableVertexAttribArray(BONE_WEIGHT_LOCATION);    
    glVertexAttribPointer(BONE_WEIGHT_LOCATION, 4, GL_FLOAT, GL_FALSE, sizeof(VertexBoneData), (const GLvoid*)16);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_Buffers[INDEX_BUFFER]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Indices[0]) * Indices.size(), &Indices[0], GL_STATIC_DRAW);

	//----------------------------------------------------------------------------------------------
    // @brief creating a buffer and giving a location for Pi (CoR)
    // @author Anastasija Belaka
    //---------------------------------------------------------------------------------------------- 
	glBindBuffer(GL_ARRAY_BUFFER, m_Buffers[PI_VB]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Pi[0]) * Pi.size(), &Pi[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(PI_LOCATION);
	glVertexAttribPointer(PI_LOCATION, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//----------------------------------------------------------------------------------------------

    return GLCheckError();
}


void SkinnedMesh::InitMesh(uint MeshIndex,
                    const aiMesh* paiMesh,
                    vector<Vector3f>& Positions,
                    vector<Vector3f>& Normals,
                    vector<Vector2f>& TexCoords,
                    vector<VertexBoneData>& Bones,
                    vector<uint>& Indices,
					vector<Vector3f>& Pi)
{    
    const aiVector3D Zero3D(0.0f, 0.0f, 0.0f);
    
    // Populate the vertex attribute vectors
    for (uint i = 0 ; i < paiMesh->mNumVertices ; i++) {
        const aiVector3D* pPos      = &(paiMesh->mVertices[i]);
        const aiVector3D* pNormal   = &(paiMesh->mNormals[i]);
        const aiVector3D* pTexCoord = paiMesh->HasTextureCoords(0) ? &(paiMesh->mTextureCoords[0][i]) : &Zero3D;

        Positions.push_back(Vector3f(pPos->x, pPos->y, pPos->z));
        Normals.push_back(Vector3f(pNormal->x, pNormal->y, pNormal->z));
        TexCoords.push_back(Vector2f(pTexCoord->x, pTexCoord->y));        
    }

    
    LoadBones(MeshIndex, paiMesh, Bones);
	//----------------------------------------------------------------------------------------------
	// @brief calling Function for finding Centre of Rotation (CoR)
	// @author Anastasija Belaka
	//----------------------------------------------------------------------------------------------
	FindPi(MeshIndex, paiMesh, Bones, Pi);
	//----------------------------------------------------------------------------------------------


    // Populate the index buffer
    for (uint i = 0 ; i < paiMesh->mNumFaces ; i++) {
        const aiFace& Face = paiMesh->mFaces[i];
        assert(Face.mNumIndices == 3);
        Indices.push_back(Face.mIndices[0]);
        Indices.push_back(Face.mIndices[1]);
        Indices.push_back(Face.mIndices[2]);
    }
}


void SkinnedMesh::LoadBones(uint MeshIndex, const aiMesh* pMesh, vector<VertexBoneData>& Bones)
{
	//auto res = 0.0;
    for (uint i = 0 ; i < pMesh->mNumBones ; i++) {                
        uint BoneIndex = 0;        
        string BoneName(pMesh->mBones[i]->mName.data);
        
        if (m_BoneMapping.find(BoneName) == m_BoneMapping.end()) {
            // Allocate an index for a new bone
            BoneIndex = m_NumBones;
            m_NumBones++;            
	        BoneInfo bi;			
			m_BoneInfo.push_back(bi);
            m_BoneInfo[BoneIndex].BoneOffset = pMesh->mBones[i]->mOffsetMatrix;            
            m_BoneMapping[BoneName] = BoneIndex;
        }
        else {
            BoneIndex = m_BoneMapping[BoneName];
        }                      
        
        for (uint j = 0 ; j < pMesh->mBones[i]->mNumWeights ; j++) {
            uint VertexID = m_Entries[MeshIndex].BaseVertex + pMesh->mBones[i]->mWeights[j].mVertexId;
            float Weight  = pMesh->mBones[i]->mWeights[j].mWeight;                   
            Bones[VertexID].AddBoneData(BoneIndex, Weight);
        }
    }    
}

//----------------------------------------------------------------------------------------------
// @brief Similarity Function
// @author Anastasija Belaka
//---------------------------------------------------------------------------------------------- 
typedef std::vector<float> weight_t;

float SimilarityFunctionW(weight_t& wp, weight_t& wv, float sigma)
{
	float ret = 0.0f;
	for (size_t j = 0; j < wp.size(); j++)
	{
		for (size_t k = 0; k < wp.size(); k++)
		{
			if (j == k) continue;
			float cr = wp[j] * wv[k] - wp[k] * wv[j];
			ret += wp[j] * wp[k] * wv[j] * wv[k] * exp(-cr * cr / (sigma*sigma));
		}
	}
	return ret;
}
//----------------------------------------------------------------------------------------------
// end of function
//----------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------
// @brief Function for Computation of CoR on Triangles Meshes
// @author Anastasija Belaka
//----------------------------------------------------------------------------------------------
void SkinnedMesh::FindPi(uint MeshIndex, const aiMesh* pMesh, vector<VertexBoneData>& Bones, vector<Vector3f>& Pi)
{
	float sum = 0;
	Vector3f avrPos;

	std::vector<weight_t> weights;
	uint max_bone = 0;
	for (uint i = 0; i < pMesh->mNumVertices; i++)
	{
		for (uint b = 0; b < NUM_BONES_PER_VEREX; b++)
		{
			uint pid = Bones[i].IDs[b];
			if (pid > max_bone)
				max_bone = pid;
		}
	}
	max_bone++;

	weights.resize(pMesh->mNumVertices);

	for (uint i = 0; i < pMesh->mNumVertices; i++)
	{
		weight_t w_i;
		w_i.resize(max_bone, 0.0);
		for (uint b = 0; b < NUM_BONES_PER_VEREX; b++)
		{
			uint pid = Bones[i].IDs[b];
			if (pid != 0 || Bones[i].Weights[b] != 0.0)
				w_i[pid] = Bones[i].Weights[b];
		}
		weights[i] = w_i;
	}

	for (uint i = 0; i < pMesh->mNumVertices; i++)
	{
		std::vector<float>& wI = weights[i];
		float pidet = 0.0;
		Vector3f pinum(0.0);

		for (uint n = 0; n < pMesh->mNumFaces; n++)
		{
			auto tA = pMesh->mFaces[n].mIndices[0];
			auto tB = pMesh->mFaces[n].mIndices[1];
			auto tC = pMesh->mFaces[n].mIndices[2];
			std::vector<float>& wA = weights[tA];
			std::vector<float>& wB = weights[tA];
			std::vector<float>& wC = weights[tA];
			auto vA = pMesh->mVertices[tA];
			auto vB = pMesh->mVertices[tB];
			auto vC = pMesh->mVertices[tC];
			auto crossPr = (vA.y*vB.z - vA.z*vB.y, vA.z*vB.x - vA.x*vB.z, vA.x*vB.y - vA.y*vB.x);
			auto area = fabs(crossPr) / 2.0;
			weight_t averageWeight;
			averageWeight.resize(max_bone);
			for (size_t b = 0; b < max_bone; b++)
			{
				averageWeight[b] = (wA[b] + wB[b] + wC[b]) / 3.0;
			}
			Vector3f averagePos;
			averagePos.x = (vA.x + vB.x + vC.x) / 3;
			averagePos.y = (vA.y + vB.y + vC.y) / 3;
			averagePos.z = (vA.z + vB.z + vC.z) / 3;

			// Similarity Function
			float s = SimilarityFunctionW(wI, averageWeight, 0.1);
			pinum += averagePos * s*area;
			pidet += s * area;
		}
		Pi[i].x = pinum.x / pidet;
		Pi[i].y = pinum.y / pidet;
		Pi[i].z = pinum.z / pidet;
	}
}
//----------------------------------------------------------------------------------------------
// end of function
//----------------------------------------------------------------------------------------------

bool SkinnedMesh::InitMaterials(const aiScene* pScene, const string& Filename)
{
    // Extract the directory part from the file name
    string::size_type SlashIndex = Filename.find_last_of("/");
    string Dir;

    if (SlashIndex == string::npos) {
        Dir = ".";
    }
    else if (SlashIndex == 0) {
        Dir = "/";
    }
    else {
        Dir = Filename.substr(0, SlashIndex);
    }

    bool Ret = true;

    // Initialize the materials
    for (uint i = 0 ; i < pScene->mNumMaterials ; i++) {
        const aiMaterial* pMaterial = pScene->mMaterials[i];

        m_Textures[i] = NULL;

        if (pMaterial->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
            aiString Path;

            if (pMaterial->GetTexture(aiTextureType_DIFFUSE, 0, &Path, NULL, NULL, NULL, NULL, NULL) == AI_SUCCESS) {
                string p(Path.data);
                
                if (p.substr(0, 2) == ".\\") {                    
                    p = p.substr(2, p.size() - 2);
                }
                               
                string FullPath = Dir + "/" + p;
                    
                m_Textures[i] = new Texture(GL_TEXTURE_2D, FullPath.c_str());

                if (!m_Textures[i]->Load()) {
                    printf("Error loading texture '%s'\n", FullPath.c_str());
                    delete m_Textures[i];
                    m_Textures[i] = NULL;
                    Ret = false;
                }
                else {
                    printf("%d - loaded texture '%s'\n", i, FullPath.c_str());
                }
            }
        }
    }

    return Ret;
}


void SkinnedMesh::Render()
{
    glBindVertexArray(m_VAO);
    
    for (uint i = 0 ; i < m_Entries.size() ; i++) {
        const uint MaterialIndex = m_Entries[i].MaterialIndex;

        assert(MaterialIndex < m_Textures.size());
        
        if (m_Textures[MaterialIndex]) {
            m_Textures[MaterialIndex]->Bind(GL_TEXTURE0);
        }

		glDrawElementsBaseVertex(GL_TRIANGLES, 
                                 m_Entries[i].NumIndices, 
                                 GL_UNSIGNED_INT, 
                                 (void*)(sizeof(uint) * m_Entries[i].BaseIndex), 
                                 m_Entries[i].BaseVertex);
    }

    // Make sure the VAO is not changed from the outside    
    glBindVertexArray(0);
}


uint SkinnedMesh::FindPosition(float AnimationTime, const aiNodeAnim* pNodeAnim)
{    
    for (uint i = 0 ; i < pNodeAnim->mNumPositionKeys - 1 ; i++) {
        if (AnimationTime < (float)pNodeAnim->mPositionKeys[i + 1].mTime) {
            return i;
        }
    }
    
    assert(0);

    return 0;
}


uint SkinnedMesh::FindRotation(float AnimationTime, const aiNodeAnim* pNodeAnim)
{
    assert(pNodeAnim->mNumRotationKeys > 0);

    for (uint i = 0 ; i < pNodeAnim->mNumRotationKeys - 1 ; i++) {
        if (AnimationTime < (float)pNodeAnim->mRotationKeys[i + 1].mTime) {
            return i;
        }
    }
    
    assert(0);

    return 0;
}


uint SkinnedMesh::FindScaling(float AnimationTime, const aiNodeAnim* pNodeAnim)
{
    assert(pNodeAnim->mNumScalingKeys > 0);
    
    for (uint i = 0 ; i < pNodeAnim->mNumScalingKeys - 1 ; i++) {
        if (AnimationTime < (float)pNodeAnim->mScalingKeys[i + 1].mTime) {
            return i;
        }
    }
    
    assert(0);

    return 0;
}


void SkinnedMesh::CalcInterpolatedPosition(aiVector3D& Out, float AnimationTime, const aiNodeAnim* pNodeAnim)
{
    if (pNodeAnim->mNumPositionKeys == 1) {
        Out = pNodeAnim->mPositionKeys[0].mValue;
        return;
    }
            
    uint PositionIndex = FindPosition(AnimationTime, pNodeAnim);
    uint NextPositionIndex = (PositionIndex + 1);
    assert(NextPositionIndex < pNodeAnim->mNumPositionKeys);
    float DeltaTime = (float)(pNodeAnim->mPositionKeys[NextPositionIndex].mTime - pNodeAnim->mPositionKeys[PositionIndex].mTime);
    float Factor = (AnimationTime - (float)pNodeAnim->mPositionKeys[PositionIndex].mTime) / DeltaTime;
    //assert(Factor >= 0.0f && Factor <= 1.0f);
    const aiVector3D& Start = pNodeAnim->mPositionKeys[PositionIndex].mValue;
    const aiVector3D& End = pNodeAnim->mPositionKeys[NextPositionIndex].mValue;
    aiVector3D Delta = End - Start;
    Out = Start + Factor * Delta;
}


void SkinnedMesh::CalcInterpolatedRotation(aiQuaternion& Out, float AnimationTime, const aiNodeAnim* pNodeAnim)
{
	// we need at least two values to interpolate...
    if (pNodeAnim->mNumRotationKeys == 1) {
        Out = pNodeAnim->mRotationKeys[0].mValue;
        return;
    }
    
    uint RotationIndex = FindRotation(AnimationTime, pNodeAnim);
    uint NextRotationIndex = (RotationIndex + 1);
    assert(NextRotationIndex < pNodeAnim->mNumRotationKeys);
    float DeltaTime = (float)(pNodeAnim->mRotationKeys[NextRotationIndex].mTime - pNodeAnim->mRotationKeys[RotationIndex].mTime);
    float Factor = (AnimationTime - (float)pNodeAnim->mRotationKeys[RotationIndex].mTime) / DeltaTime;
    //assert(Factor >= 0.0f && Factor <= 1.0f);
    const aiQuaternion& StartRotationQ = pNodeAnim->mRotationKeys[RotationIndex].mValue;
    const aiQuaternion& EndRotationQ   = pNodeAnim->mRotationKeys[NextRotationIndex].mValue;    
    aiQuaternion::Interpolate(Out, StartRotationQ, EndRotationQ, Factor);
    Out = Out.Normalize();
}


void SkinnedMesh::CalcInterpolatedScaling(aiVector3D& Out, float AnimationTime, const aiNodeAnim* pNodeAnim)
{
    if (pNodeAnim->mNumScalingKeys == 1) {
        Out = pNodeAnim->mScalingKeys[0].mValue;
        return;
    }

    uint ScalingIndex = FindScaling(AnimationTime, pNodeAnim);
    uint NextScalingIndex = (ScalingIndex + 1);
    assert(NextScalingIndex < pNodeAnim->mNumScalingKeys);
    float DeltaTime = (float)(pNodeAnim->mScalingKeys[NextScalingIndex].mTime - pNodeAnim->mScalingKeys[ScalingIndex].mTime);
    float Factor = (AnimationTime - (float)pNodeAnim->mScalingKeys[ScalingIndex].mTime) / DeltaTime;
    //assert(Factor >= 0.0f && Factor <= 1.0f);
    const aiVector3D& Start = pNodeAnim->mScalingKeys[ScalingIndex].mValue;
    const aiVector3D& End   = pNodeAnim->mScalingKeys[NextScalingIndex].mValue;
    aiVector3D Delta = End - Start;
    Out = Start + Factor * Delta;
}


void SkinnedMesh::ReadNodeHeirarchy(float AnimationTime, const aiNode* pNode, const Matrix4f& ParentTransform)
{
	string NodeName(pNode->mName.data);

	const aiAnimation* pAnimation = m_pScene->mAnimations[0];

	Matrix4f NodeTransformation(pNode->mTransformation);

	const aiNodeAnim* pNodeAnim = FindNodeAnim(pAnimation, NodeName);

	if (pNodeAnim) {
		// Interpolate scaling and generate scaling transformation matrix
		aiVector3D Scaling;
		CalcInterpolatedScaling(Scaling, AnimationTime, pNodeAnim);
		Matrix4f ScalingM;
		ScalingM.InitScaleTransform(Scaling.x, Scaling.y, Scaling.z);
		
			// Interpolate rotation and generate rotation transformation matrix
			aiQuaternion RotationQ;
			CalcInterpolatedRotation(RotationQ, AnimationTime, pNodeAnim);
			Matrix4f RotationM = Matrix4f(RotationQ.GetMatrix());

			// Interpolate translation and generate translation transformation matrix
			aiVector3D Translation;
			CalcInterpolatedPosition(Translation, AnimationTime, pNodeAnim);
			Matrix4f TranslationM;
			TranslationM.InitTranslationTransform(Translation.x, Translation.y, Translation.z);

			// Combine the above transformations
			NodeTransformation = TranslationM * RotationM * ScalingM;
	}
	Matrix4f GlobalTransformation = ParentTransform * NodeTransformation;

	if (m_BoneMapping.find(NodeName) != m_BoneMapping.end()) {
		uint BoneIndex = m_BoneMapping[NodeName];
		m_BoneInfo[BoneIndex].FinalTransformation = m_GlobalInverseTransform * GlobalTransformation * m_BoneInfo[BoneIndex].BoneOffset;
	}

	for (uint i = 0; i < pNode->mNumChildren; i++) {
		ReadNodeHeirarchy(AnimationTime, pNode->mChildren[i], GlobalTransformation);
	}
}


void SkinnedMesh::BoneTransform(float TimeInSeconds, vector<Matrix4f>& Transforms)
{
    Matrix4f Identity;
    Identity.InitIdentity();
    
    float TicksPerSecond = (float)(m_pScene->mAnimations[0]->mTicksPerSecond != 0 ? m_pScene->mAnimations[0]->mTicksPerSecond : 25.0f);
    float TimeInTicks = TimeInSeconds * TicksPerSecond;
    float AnimationTime = fmod(TimeInTicks, (float)m_pScene->mAnimations[0]->mDuration);

    ReadNodeHeirarchy(AnimationTime, m_pScene->mRootNode, Identity);

    Transforms.resize(m_NumBones);

    for (uint i = 0 ; i < m_NumBones ; i++) {
        Transforms[i] = m_BoneInfo[i].FinalTransformation;
    }
}


const aiNodeAnim* SkinnedMesh::FindNodeAnim(const aiAnimation* pAnimation, const string NodeName)
{
    for (uint i = 0 ; i < pAnimation->mNumChannels ; i++) {
        const aiNodeAnim* pNodeAnim = pAnimation->mChannels[i];
        
        if (string(pNodeAnim->mNodeName.data) == NodeName) {
            return pNodeAnim;
        }
    }
    
    return NULL;
}