using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class IsotropicFVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;
    private float g = 9.8f;
    float restitution 	= 0.4f;			
    float friction 	= 0.5f;	

	int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	//For Laplacian smoothing.
	Vector3[]   V_sum;
	int[]		V_num;
	private float LaplacianBlend = 0.25f;

	SVD svd = new SVD();

    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/HW3/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);//第一项存储了四面体个数
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;//存储顶点索引
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/HW3/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;//存储顶点位置
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;//重新计算模型空间的顶点位置
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}
        // tet_number=1;
        // Tet = new int[tet_number*4];
        // Tet[0]=0;
        // Tet[1]=1;
        // Tet[2]=2;
        // Tet[3]=3;
        //
        // number=4;
        // X = new Vector3[number];
        // V = new Vector3[number];
        // Force = new Vector3[number];
        // X[0]= new Vector3(0, 0, 0);
        // X[1]= new Vector3(1, 0, 0);
        // X[2]= new Vector3(0, 1, 0);
        // X[3]= new Vector3(0, 0, 1);


        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];//所有三角面，1个四面体4个面，一共12个顶点
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

		//Need to allocate and assign inv_Dm
		inv_Dm = new Matrix4x4[tet_number];
		for (int i = 0; i < tet_number; i++)
		{
			inv_Dm[i] = Build_Edge_Matrix(i).inverse;
		}

    }

    //获取一个四棱锥 01 02 03三个边的矩阵
    Matrix4x4 Build_Edge_Matrix(int tet)
    {
    	Matrix4x4 ret=Matrix4x4.zero;

        ret[0, 0] = X[Tet[tet * 4 + 1]].x - X[Tet[tet * 4 + 0]].x;
        ret[1, 0] = X[Tet[tet * 4 + 1]].y - X[Tet[tet * 4 + 0]].y;
        ret[2, 0] = X[Tet[tet * 4 + 1]].z - X[Tet[tet * 4 + 0]].z;
        ret[3, 0] = 0;
        ret[0, 1] = X[Tet[tet * 4 + 2]].x - X[Tet[tet * 4 + 0]].x;
        ret[1, 1] = X[Tet[tet * 4 + 2]].y - X[Tet[tet * 4 + 0]].y;
        ret[2, 1] = X[Tet[tet * 4 + 2]].z - X[Tet[tet * 4 + 0]].z;
        ret[3, 1] = 0;
        ret[0, 2] = X[Tet[tet * 4 + 3]].x - X[Tet[tet * 4 + 0]].x;
        ret[1, 2] = X[Tet[tet * 4 + 3]].y - X[Tet[tet * 4 + 0]].y;
        ret[2, 2] = X[Tet[tet * 4 + 3]].z - X[Tet[tet * 4 + 0]].z;
        ret[3, 2] = 0;
        ret[0, 3] = 0;
        ret[1, 3] = 0;
        ret[2, 3] = 0;
        ret[3, 3] = 1;
		return ret;
    }
    private Matrix4x4 MatrixAdd(Matrix4x4 m1, Matrix4x4 m2)
    {
	    Matrix4x4 ret = new Matrix4x4();
	    for (int i = 0; i < 4; i++)
	    {
		    for (int j = 0; j < 4; j++)
		    {
			    ret[i, j] = m1[i, j] + m2[i, j];
		    }
	    }

	    return ret;
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
    
    private Matrix4x4 MatrixMultFloat(Matrix4x4 m, float f)
    {
	    Matrix4x4 ret = new Matrix4x4();
	    for (int i = 0; i < 4; i++)
	    {
		    for (int j = 0; j < 4; j++)
		    {
			    ret[i, j] = m[i, j] * f;
		    }
	    }

	    return ret;
    }
    
    private float Trace(Matrix4x4 m)
    {
	    return m[0,0]+m[1,1]+m[2,2];
    }


    void _Update()
    {
	    // Jump up.
	    if (Input.GetKeyDown(KeyCode.Space))
	    {
		    for (int i = 0; i < number; i++)
			    V[i].y += 0.2f;
	    }

	    for (int i = 0; i < number; i++)
	    {
		    //Add gravity to Force.
		    Force[i] = new Vector3(0, -g * mass, 0);
	    }

	    for (int tet = 0; tet < tet_number; tet++)
	    {
		    //Deformation Gradient
		    //F=当前边矩阵*参考边矩阵的逆
		    Matrix4x4 F = Build_Edge_Matrix(tet) * inv_Dm[tet];
		    
		    
		    //各向异性材料，可以用另一种方式求解
		    //不利用P=FS，P是F的函数P=f(F)，F分解为UDV后，P只跟D有关，P=Uf(D)V
		    //其中D是对角矩阵，f(D)实际上是W对d1 d2 d3的某种导数

		    Matrix4x4 U = new Matrix4x4(), D = new Matrix4x4(), V = new Matrix4x4();
		    svd.svd(F,ref U,ref D, ref V);
		    float I = Mathf.Pow(D[0, 0], 2) + Mathf.Pow(D[1, 1], 2) + Mathf.Pow(D[2, 2], 2);
		    float III = Mathf.Pow(D[0, 0], 2) * Mathf.Pow(D[1, 1], 2) +
		               Mathf.Pow(D[0, 0], 2) * Mathf.Pow(D[2, 2], 2) +
		               Mathf.Pow(D[2, 2], 2) * Mathf.Pow(D[1, 1], 2);
		    float II = Mathf.Pow(D[0, 0], 4) + Mathf.Pow(D[1, 1], 4) + Mathf.Pow(D[2, 2], 4);
		    //W对d1 d2 d3求导，先对I II III求导
		    float dWdI = 0.25f * stiffness_0 * (I - 3) - 0.5f * stiffness_1;
		    float dWdII = 0.25f * stiffness_1;
		    float dWdIII = 0;
		    //I II III再对d1 d2 d3求导
		    float dIdd1 = 2 * D[0, 0];
		    float dIdd2 = 2 * D[1, 1];
		    float dIdd3 = 2 * D[2, 2];
		    float dIIIdd1 = 2 * Mathf.Pow(D[1, 1], 2) * D[0, 0] + 2 * Mathf.Pow(D[2, 2], 2) * D[0, 0];
		    float dIIIdd2 = 2 * Mathf.Pow(D[0, 0], 2) * D[1, 1] + 2 * Mathf.Pow(D[2, 2], 2) * D[1, 1];
		    float dIIIdd3 = 2 * Mathf.Pow(D[1, 1], 2) * D[2, 2] + 2 * Mathf.Pow(D[0, 0], 2) * D[2, 2];
		    float dIIdd1 = 4 * Mathf.Pow(D[0, 0], 3);
		    float dIIdd2 = 4 * Mathf.Pow(D[1, 1], 3);
		    float dIIdd3 = 4 * Mathf.Pow(D[2, 2], 3);

		    Matrix4x4 P_ = Matrix4x4.identity;
		    P_[0, 0] = dWdI * dIdd1 + dWdII * dIIdd1 + dWdIII * dIIIdd1;
		    P_[1, 1] = dWdI * dIdd2 + dWdII * dIIdd2 + dWdIII * dIIIdd2;
		    P_[2, 2] = dWdI * dIdd3 + dWdII * dIIdd3 + dWdIII * dIIIdd3;
		    Matrix4x4 P = U * P_ * V.transpose;

		    Matrix4x4 elasticF = MatrixMultFloat(P * inv_Dm[tet].transpose, -1/(6 * inv_Dm[tet].determinant));
		    //Green Strain
		    //G = 1/2(Ft * F -I)
		    //Matrix4x4 G = MatrixMultFloat(MatrixSub(F.transpose * F, Matrix4x4.identity), 0.5f);
		    //Second PK Stress
		    //S是W对G求导，W跟模型有关，这里用一种简单模型
		    //S = 2*s1*G+s0*G的迹*I
		    //Matrix4x4 S = MatrixAdd(MatrixMultFloat(G, 2 * stiffness_1),
			    //MatrixMultFloat(Matrix4x4.identity, Trace(G) * stiffness_0));
		    //Elastic Force
		    //根据体积和FS求解
		    //Matrix4x4 elasticF = MatrixMultFloat(F * S * inv_Dm[tet].transpose, -1/(6 * inv_Dm[tet].determinant));
		    //F的三个列向量分别是点1 2 3受的力
		    Force[Tet[tet * 4 + 1]].x += elasticF[0, 0];
		    Force[Tet[tet * 4 + 1]].y += elasticF[1, 0];
		    Force[Tet[tet * 4 + 1]].z += elasticF[2, 0];
		    Force[Tet[tet * 4 + 2]].x += elasticF[0, 1];
		    Force[Tet[tet * 4 + 2]].y += elasticF[1, 1];
		    Force[Tet[tet * 4 + 2]].z += elasticF[2, 1];
		    Force[Tet[tet * 4 + 3]].x += elasticF[0, 2];
		    Force[Tet[tet * 4 + 3]].y += elasticF[1, 2];
		    Force[Tet[tet * 4 + 3]].z += elasticF[2, 2];
		    Force[Tet[tet * 4 + 0]].x += -elasticF[0, 0] - elasticF[0, 1] - elasticF[0, 2];
		    Force[Tet[tet * 4 + 0]].y += -elasticF[1, 0] - elasticF[1, 1] - elasticF[1, 2];
		    Force[Tet[tet * 4 + 0]].z += -elasticF[2, 0] - elasticF[2, 1] - elasticF[2, 2];
	    }

	    for (int i = 0; i < number; i++)
	    {
		    //Update X and V here.
		    V[i] = (V[i] + dt * Force[i] / mass) * damp;
		    V_sum[i] = new Vector3(0, 0, 0);
		    V_num[i] = 0;
	    }
	    //拉普拉斯平滑，每个点速度收到周围点速度的影响
	    for (int tet = 0; tet < tet_number; tet++)
	    {
		    //Update X and V here.
		    V_sum[Tet[tet * 4 + 0]] += V[Tet[tet * 4 + 1]] + V[Tet[tet * 4 + 2]] + V[Tet[tet * 4 + 3]];
		    V_sum[Tet[tet * 4 + 1]] += V[Tet[tet * 4 + 0]] + V[Tet[tet * 4 + 2]] + V[Tet[tet * 4 + 3]];
		    V_sum[Tet[tet * 4 + 2]] += V[Tet[tet * 4 + 1]] + V[Tet[tet * 4 + 0]] + V[Tet[tet * 4 + 3]];
		    V_sum[Tet[tet * 4 + 3]] += V[Tet[tet * 4 + 1]] + V[Tet[tet * 4 + 2]] + V[Tet[tet * 4 + 0]];
		    V_num[Tet[tet * 4 + 0]] += 3;
		    V_num[Tet[tet * 4 + 1]] += 3;
		    V_num[Tet[tet * 4 + 2]] += 3;
		    V_num[Tet[tet * 4 + 3]] += 3;
	    }

		for (int i = 0; i < number; i++) 
		{
			//平滑后的速度，更新位置
			V[i] = V[i] * (1 - LaplacianBlend) + LaplacianBlend * V_sum[i] / V_num[i];
			X[i] += V[i] * dt;
			//计算和地面的碰撞
		    Vector3 P = new Vector3(0, -3, 0);
		    Vector3 N = new Vector3(0, 1, 0);
		    if (Vector3.Dot(P - X[i], N) < 0)
		    {
			    continue;
		    }
		    
		    //如果是朝外运动，忽略
		    if (Vector3.Dot(V[i], N) > 0)
		    {
			    //TMP_X[i] += Vector3.Dot(V[i], N) * N;
			    continue;
		    }
		    
		    X[i] += Vector3.Dot(P - X[i], N) * N;
		    Vector3 VN = Vector3.Dot(V[i], N) * N;
		    Vector3 VT = V[i] - VN;
		    V[i] = (VT.magnitude<0.0001) ? -restitution *VN :
			    -restitution *VN + VT * Math.Max(0, 1 - friction * (1 + restitution) * VN.magnitude / VT.magnitude);
		    
	    }
    }

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}
