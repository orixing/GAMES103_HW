using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Rigid_Bunny_by_Shape_Matching : MonoBehaviour
{
	public bool launched = false;
	Vector3[] X;
	Vector3[] TMP_X;
	Vector3[] Q;
	Vector3[] V;
	Matrix4x4 QQt = Matrix4x4.zero;


    // Start is called before the first frame update
    void Start()
    {
    	Mesh mesh = GetComponent<MeshFilter>().mesh;
        V = new Vector3[mesh.vertices.Length];
        X = mesh.vertices;
        TMP_X = new Vector3[mesh.vertices.Length];
        Q = mesh.vertices;

        //Centerizing Q.
        Vector3 c=Vector3.zero;
        for(int i=0; i<Q.Length; i++)
        	c+=Q[i];
        c/=Q.Length;
        for(int i=0; i<Q.Length; i++)
        	Q[i]-=c;

        //Get QQ^t ready.
		for(int i=0; i<Q.Length; i++)
		{
			QQt[0, 0]+=Q[i][0]*Q[i][0];
			QQt[0, 1]+=Q[i][0]*Q[i][1];
			QQt[0, 2]+=Q[i][0]*Q[i][2];
			QQt[1, 0]+=Q[i][1]*Q[i][0];
			QQt[1, 1]+=Q[i][1]*Q[i][1];
			QQt[1, 2]+=Q[i][1]*Q[i][2];
			QQt[2, 0]+=Q[i][2]*Q[i][0];
			QQt[2, 1]+=Q[i][2]*Q[i][1];
			QQt[2, 2]+=Q[i][2]*Q[i][2];
		}
		QQt[3, 3]=1;

		for(int i=0; i<X.Length; i++)
			V[i][0]=4.0f;

		Update_Mesh(transform.position, Matrix4x4.Rotate(transform.rotation), 0);
		transform.position=Vector3.zero;
		transform.rotation=Quaternion.identity;
   	}

   	// Polar Decomposition that returns the rotation from F.
   	Matrix4x4 Get_Rotation(Matrix4x4 F)
	{
		Matrix4x4 C = Matrix4x4.zero;
	    for(int ii=0; ii<3; ii++)
	    for(int jj=0; jj<3; jj++)
	    for(int kk=0; kk<3; kk++)
	        C[ii,jj]+=F[kk,ii]*F[kk,jj];
	   
	   	Matrix4x4 C2 = Matrix4x4.zero;
		for(int ii=0; ii<3; ii++)
	    for(int jj=0; jj<3; jj++)
	    for(int kk=0; kk<3; kk++)
	        C2[ii,jj]+=C[ii,kk]*C[jj,kk];
	    
	    float det    =  F[0,0]*F[1,1]*F[2,2]+
	                    F[0,1]*F[1,2]*F[2,0]+
	                    F[1,0]*F[2,1]*F[0,2]-
	                    F[0,2]*F[1,1]*F[2,0]-
	                    F[0,1]*F[1,0]*F[2,2]-
	                    F[0,0]*F[1,2]*F[2,1];
	    
	    float I_c    =   C[0,0]+C[1,1]+C[2,2];
	    float I_c2   =   I_c*I_c;
	    float II_c   =   0.5f*(I_c2-C2[0,0]-C2[1,1]-C2[2,2]);
	    float III_c  =   det*det;
	    float k      =   I_c2-3*II_c;
	    
	    Matrix4x4 inv_U = Matrix4x4.zero;
	    if(k<1e-10f)
	    {
	        float inv_lambda=1/Mathf.Sqrt(I_c/3);
	        inv_U[0,0]=inv_lambda;
	        inv_U[1,1]=inv_lambda;
	        inv_U[2,2]=inv_lambda;
	    }
	    else
	    {
	        float l = I_c*(I_c*I_c-4.5f*II_c)+13.5f*III_c;
	        float k_root = Mathf.Sqrt(k);
	        float value=l/(k*k_root);
	        if(value<-1.0f) value=-1.0f;
	        if(value> 1.0f) value= 1.0f;
	        float phi = Mathf.Acos(value);
	        float lambda2=(I_c+2*k_root*Mathf.Cos(phi/3))/3.0f;
	        float lambda=Mathf.Sqrt(lambda2);
	        
	        float III_u = Mathf.Sqrt(III_c);
	        if(det<0)   III_u=-III_u;
	        float I_u = lambda + Mathf.Sqrt(-lambda2 + I_c + 2*III_u/lambda);
	        float II_u=(I_u*I_u-I_c)*0.5f;
	        
	        
	        float inv_rate, factor;
	        inv_rate=1/(I_u*II_u-III_u);
	        factor=I_u*III_u*inv_rate;
	        
	       	Matrix4x4 U = Matrix4x4.zero;
			U[0,0]=factor;
	        U[1,1]=factor;
	        U[2,2]=factor;
	        
	        factor=(I_u*I_u-II_u)*inv_rate;
	        for(int i=0; i<3; i++)
	        for(int j=0; j<3; j++)
	            U[i,j]+=factor*C[i,j]-inv_rate*C2[i,j];
	        
	        inv_rate=1/III_u;
	        factor=II_u*inv_rate;
	        inv_U[0,0]=factor;
	        inv_U[1,1]=factor;
	        inv_U[2,2]=factor;
	        
	        factor=-I_u*inv_rate;
	        for(int i=0; i<3; i++)
	        for(int j=0; j<3; j++)
	            inv_U[i,j]+=factor*U[i,j]+inv_rate*C[i,j];
	    }
	    
	    Matrix4x4 R=Matrix4x4.zero;
	    for(int ii=0; ii<3; ii++)
	    for(int jj=0; jj<3; jj++)
	    for(int kk=0; kk<3; kk++)
	        R[ii,jj]+=F[ii,kk]*inv_U[kk,jj];
	    R[3,3]=1;
	    return R;
	}

	// Update the mesh vertices according to translation c and rotation R.
	// It also updates the velocity.
	void Update_Mesh(Vector3 c, Matrix4x4 R, float inv_dt)
   	{
   		for(int i=0; i<Q.Length; i++)
		{
			Vector3 x=(Vector3)(R*Q[i])+c;
			
			// if(i==0 && float.IsNaN(x.x))
			// 	Debug.Log($"{x.x} {x.y} {x.z} {c} ");

			V[i]=(x-X[i])*inv_dt;
			X[i]=x;
		}
        Mesh mesh = GetComponent<MeshFilter>().mesh;
		mesh.vertices=X;
   	}

	//float restitution 	= 5f;		
	//float friction 	= 0.3f;	
	// void Collision(Vector3 P, Vector3 N)
	// {
	// 	//对每个点独立计算碰撞
	// 	for (int i = 0; i < Q.Length; i++)
	// 	{
	// 		//如果点和平面不相交
	// 		if (Vector3.Dot(P - X[i], N) < 0)
	// 		{
	// 			continue;
	// 		}
	// 		//如果是朝外运动，就不算了
	// 		if (Vector3.Dot(V[i], N) > 0)
	// 		{
	// 			continue;
	// 		}
	// 		//否则收到冲量，改变速度
	// 		Vector3 VN = Vector3.Dot(V[i], N)*N;
	// 		Vector3 VT = V[i] - VN;
	// 		Vector3 newVN = -restitution * VN;
	// 		Vector3 newVT = VT * Math.Max(0, 1 - friction * (1 + restitution) * VN.magnitude / VT.magnitude);
	// 		Vector3 new_V = newVN + newVT;
	// 		V[i] = new_V;
	// 	}
	// 	
	// }
	
	float restitution 	= 0.5f;		
	float friction 	= 0.5f;	
	
	void Collision(Vector3 P, Vector3 N)
	{
		//对每个点独立计算碰撞
		for (int i = 0; i < Q.Length; i++)
		{
			//如果点和平面不相交
			if (Vector3.Dot(P - TMP_X[i], N) < 0)
			{
				continue;
			}
			
			//如果是朝外运动，就移动到表面，速度如果小于某个阈值，也移动到表面
			//这样反弹会太强
			if (Vector3.Dot(V[i], N) > 0)
			{
				//TMP_X[i] += Vector3.Dot(V[i], N) * N;
				continue;
			}
			//速度如果小于某个阈值，也移动到表面
			if (V[i].magnitude < 0.001)
			{
				TMP_X[i] += Vector3.Dot(P - TMP_X[i], N) * N;
				if (float.IsNaN(TMP_X[i].x))
				{
	        
				}
				continue;
			}
			//应当模拟一个点在dt时间内撞墙后反弹的效果
			//实际进入物体的位移应当转化成在反弹方向上的位移
			//位移也分为N和T方向，两个方向位移转化比例等同于速度
			Vector3 XN_in = Vector3.Dot(TMP_X[i]-P, N) * N; //N方向上进入的位移
			Vector3 XN = Vector3.Dot(TMP_X[i]-X[i], N) * N;//dt时间内N方向上总位移
			Vector3 _X = TMP_X[i]- X[i];//dt时间内位移
			Vector3 XT = _X - XN;//dt时间内T方向上总位移
			//如果N方向位移过小，表示已经贴近边缘，也移动到表面
			if (XN.magnitude < 0.0001)
			{
				TMP_X[i] += Vector3.Dot(P - TMP_X[i], N) * N;
				continue;
			}
			Vector3 XT_in = XN_in.magnitude/XN.magnitude*XT; //T方向上进入的位移
			
			//弹射后N T两方向上的位移
			Vector3 XN_out = -restitution * XN_in;
			Vector3 XT_out = XT_in.magnitude < 0.0001 ? Vector3.zero : XT_in * Math.Max(0, 1 - friction * (1 + restitution) * XN_in.magnitude / XT_in.magnitude);
			
			//位移修正
			TMP_X[i] = TMP_X[i] - XN_in - XT_in + XN_out + XT_out;
			if (float.IsNaN(TMP_X[i].x))
			{
	        
			}
			// Vector3 fix = -XN_in - XT_in + XN_out + XT_out;
			// if(i==0)
			// 	Debug.Log($"{fix.x} {fix.y} {fix.z}");
		}
		
	}
	
	float linear_decay	= 0.99f;
	float g = 9.8f;
	float dt = 0.015f;
    // Update is called once per frame
    void Update()
    {
	    if(Input.GetKey("r"))
        {
	        for(int i=0; i<V.Length; i++)
	        {
		        V[i] = Vector3.zero;
	        }
	        Update_Mesh(new Vector3 (0, 0.6f, 0), Matrix4x4.Rotate(Quaternion.identity), 0);
	        launched=false;
        }
        if(Input.GetKey("l"))
        {
	        for(int i=0; i<V.Length; i++)
	        {
		        V[i] = new Vector3 (5, 2, 0);
	        }
	        launched=true;
        }

        if (!launched)
        {
	        return;
        }
  		//Step 1: run a simple particle system.
        for(int i=0; i<V.Length; i++)
        {
	        V[i].y -= dt * g;
	        V[i] *= linear_decay;
	        TMP_X[i] = X[i] + dt * V[i];
        }

        //Step 2: Perform simple particle collision.
		//碰撞改变速度，墙里的点获得向外的速度
		//Collision(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		// Collision(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));
		// //根据新速度计算位置
		// for(int i=0; i<V.Length; i++)
		// {
		// 	TMP_X[i] = X[i] + dt * V[i];
		// }
		//将点移动至平面外
		Collision(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));
		
		// Step 3: Use shape matching to get new translation c and 
		// new rotation R. Update the mesh by c and R.
        //Shape Matching (translation)
        Vector3 c = Vector3.zero;
        for (int i = 0; i < X.Length; i++)
        {
	        if (float.IsNaN(TMP_X[i].x))
	        {
	        
	        }
	        c += TMP_X[i];
        }
        c /= X.Length;
        //Shape Matching (rotation)
        Matrix4x4 A = Matrix4x4.zero;
        for (int i = 0; i < X.Length; i++)
        {
	        A = MatrixAdd(A,Vec3MulVec3(TMP_X[i]-c,Q[i]));
        }
        A[3, 3] = 1;
        A *= QQt.inverse;
        Matrix4x4 R = Get_Rotation(A);
        if (float.IsNaN(c.x))
        {
	        
        }
        Update_Mesh(c, R, 1/dt);
    }

    private Matrix4x4 Vec3MulVec3(Vector3 v1, Vector3 v2)
    {
	    Matrix4x4 ret = Matrix4x4.zero;
	    for (int i = 0; i < 3; i++)
	    {
		    for (int j = 0; j < 3; j++)
		    {
			    ret[i, j] = v1[i] * v2[j];
		    }
	    }

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
}
