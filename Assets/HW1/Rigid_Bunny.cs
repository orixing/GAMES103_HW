using System;
using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.4f;					// for collision
	float friction 	= 0.5f;	

	float g = 9.8f;

	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		int pointNum = 0;
		Vector3 avgPointPos = new Vector3();
		foreach (Vector3 point in vertices)
		{
			//算出点的世界坐标
			Vector3 pointWorldPos = transform.position + transform.rotation * point;
			//如果点和平面相交
			if (Vector3.Dot(P - pointWorldPos, N) < 0)
			{
				continue;
			}
			//算出点的总速度
			Vector3 pointTotalV = v + Vector3.Cross(w, transform.rotation * point);
			//如果是朝外运动，就不算了
			if (Vector3.Dot(pointTotalV, N) > 0)
			{
				continue;
			}
			//朝内运动，计数，记录位置
			pointNum++;
			avgPointPos += point;
		}
		//没有碰撞就返回
		if(pointNum == 0) return;
		//算出平均点
		avgPointPos = avgPointPos / (float)pointNum;
		Vector3 avgPointTotalV = v + Vector3.Cross(w, transform.rotation * avgPointPos);//总速度
		Vector3 avgPointTotalVN = Vector3.Dot(avgPointTotalV, N) * N;
		Vector3 avgPointTotalVT = avgPointTotalV - avgPointTotalVN;
		Vector3 new_VN = -restitution * avgPointTotalVN;
		Vector3 new_VT = avgPointTotalVT.magnitude<0.0001?Vector3.zero:avgPointTotalVT * Math.Max(0, 1 - friction * (1 + restitution) * avgPointTotalVN.magnitude / avgPointTotalVT.magnitude);
		Vector3 new_V = new_VN + new_VT;
		//根据新旧速度算冲量
		Matrix4x4 rotationMatrix = Matrix4x4.Rotate(transform.rotation);
		Matrix4x4 I_inverse = rotationMatrix * I_ref * Matrix4x4.Transpose(rotationMatrix);//转动惯量的逆
		Matrix4x4 Rr_Matrix = Vec2CrossMat(transform.rotation * avgPointPos);
		Matrix4x4 m1 = new Matrix4x4();
		m1[0, 0] = m1[1, 1] = m1[2, 2] = m1[3, 3] = 1.0f / mass;
		Matrix4x4 m2 = Rr_Matrix * I_inverse *Rr_Matrix;
		Matrix4x4 K = MatrixSub(m1, m2);
		Vector3 J = K.inverse.MultiplyVector(new_V-avgPointTotalV);//冲量
		//根据冲量更新线速度和角速度
		v += J / mass;
		w += I_inverse.MultiplyVector(Vector3.Cross(transform.rotation * avgPointPos, J));
	}

	private Matrix4x4 MatrixSub(Matrix4x4 m1, Matrix4x4 m2)
	{
		Matrix4x4 ret = new Matrix4x4();
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret[i, j] = m1[i, j] - m2[i, j];
			}
		}

		return ret;
	}

	private Matrix4x4 Vec2CrossMat(Vector3 v)
	{
		Matrix4x4 ret = Matrix4x4.zero;
		ret[0,1] = -v[2];
		ret[0,2] = v[1];
		ret[1,0] = v[2];
		ret[1,2] = -v[0];
		ret[2, 0] = -v[1];
		ret[2, 1] = v[0];
		ret[3, 3] = 1;
		return ret;
	}

	// Update is called once per frame
	void Update () 
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
			restitution = 0.5f;
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			launched=true;
		}

		if (!launched)
		{
			return;
		}

		// Part I: Update velocities

		v.y -= dt * g;
		v *= linear_decay;
		w *= angular_decay;

		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x = transform.position;
		x += v * dt;
		//Update angular status
		Quaternion q = transform.rotation;
		Vector3 dw = w * 0.5f *dt;
		q = QuaternionAdd(q, new Quaternion(dw.x, dw.y, dw.z, 0.0f) * q);
		
		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
		Debug.Log($"{v.x},{v.y},{v.z}");
	}

	private Quaternion QuaternionAdd(Quaternion q1, Quaternion q2)
	{
		return new Quaternion(q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w);
	}
}
