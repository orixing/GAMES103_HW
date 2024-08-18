using System;
using UnityEngine;
using System.Collections;
using UnityEngine.Assertions;
using Random = UnityEngine.Random;

public class wave_motion : MonoBehaviour
{
    int size = 100;
    float rate = 0.005f;
    float gamma = 0.004f;
    float damping = 0.996f;
    float[,] old_h;
    float[,] low_h;
    float[,] vh;
    float[,] b;

    bool[,] cg_mask;
    float[,] cg_p;
    float[,] cg_r;
    float[,] cg_Ap;
    bool tag = true;

    Vector3 cube_v = Vector3.zero;
    Vector3 cube_w = Vector3.zero;

    public GameObject cube1;
    public GameObject cube2;


    // Use this for initialization
    void Start()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.Clear();

        Vector3[] X = new Vector3[size * size];

        for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            X[i * size + j].x = i * 0.1f - size * 0.05f;
            X[i * size + j].y = 0;
            X[i * size + j].z = j * 0.1f - size * 0.05f;
        }

        int[] T = new int[(size - 1) * (size - 1) * 6];
        int index = 0;
        for (int i = 0; i < size - 1; i++)
        for (int j = 0; j < size - 1; j++)
        {
            T[index * 6 + 0] = (i + 0) * size + (j + 0);
            T[index * 6 + 1] = (i + 0) * size + (j + 1);
            T[index * 6 + 2] = (i + 1) * size + (j + 1);
            T[index * 6 + 3] = (i + 0) * size + (j + 0);
            T[index * 6 + 4] = (i + 1) * size + (j + 1);
            T[index * 6 + 5] = (i + 1) * size + (j + 0);
            index++;
        }

        mesh.vertices = X;
        mesh.triangles = T;
        mesh.RecalculateNormals();

        low_h = new float[size, size];
        old_h = new float[size, size];
        vh = new float[size, size];
        b = new float[size, size];

        cg_mask = new bool [size, size];
        cg_p = new float[size, size];
        cg_r = new float[size, size];
        cg_Ap = new float[size, size];

        for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            low_h[i, j] = 99999;
            old_h[i, j] = 0;
            vh[i, j] = 0;
        }
    }

    void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
    {
        for (int i = li; i <= ui; i++)
        for (int j = lj; j <= uj; j++)
            if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j])
            {
                Ax[i, j] = 0;
                if (i != 0) Ax[i, j] -= x[i - 1, j] - x[i, j];
                if (i != size - 1) Ax[i, j] -= x[i + 1, j] - x[i, j];
                if (j != 0) Ax[i, j] -= x[i, j - 1] - x[i, j];
                if (j != size - 1) Ax[i, j] -= x[i, j + 1] - x[i, j];
            }
    }

    float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
    {
        float ret = 0;
        for (int i = li; i <= ui; i++)
        for (int j = lj; j <= uj; j++)
            if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j])
            {
                ret += x[i, j] * y[i, j];
            }

        return ret;
    }

    void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
    {
        //Solve the Laplacian problem by CG.
        A_Times(mask, x, cg_r, li, ui, lj, uj);

        for (int i = li; i <= ui; i++)
        for (int j = lj; j <= uj; j++)
            if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j])
            {
                cg_p[i, j] = cg_r[i, j] = b[i, j] - cg_r[i, j];
            }

        float rk_norm = Dot(mask, cg_r, cg_r, li, ui, lj, uj);

        for (int k = 0; k < 128; k++)
        {
            if (rk_norm < 1e-10f) break;
            A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
            float alpha = rk_norm / Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

            for (int i = li; i <= ui; i++)
            for (int j = lj; j <= uj; j++)
                if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j])
                {
                    x[i, j] += alpha * cg_p[i, j];
                    cg_r[i, j] -= alpha * cg_Ap[i, j];
                }

            float _rk_norm = Dot(mask, cg_r, cg_r, li, ui, lj, uj);
            float beta = _rk_norm / rk_norm;
            rk_norm = _rk_norm;

            for (int i = li; i <= ui; i++)
            for (int j = lj; j <= uj; j++)
                if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j])
                {
                    cg_p[i, j] = cg_r[i, j] + beta * cg_p[i, j];
                }
        }
    }

    void Shallow_Wave(float[,] old_h, float[,] h, float[,] new_h)
    {
        //Step 1:
        //Compute new_h based on the shallow wave model.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;
                if (i > 0)
                    new_h[i, j] += (h[i - 1, j] - h[i, j]) * rate;
                if (i < size - 1)
                    new_h[i, j] += (h[i + 1, j] - h[i, j]) * rate;
                if (j > 0)
                    new_h[i, j] += (h[i, j - 1] - h[i, j]) * rate;
                if (j < size - 1)
                    new_h[i, j] += (h[i, j + 1] - h[i, j]) * rate;
            }
        }

        //Step 2: Block->Water coupling
        //for block 1, calculate low_h.

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                vh[i, j] = 0;
                low_h[i, j] = 10;
            }
        }

        //vh是虚拟高度 b是线性方程组AX=B中的b，X是vh，A是拉普拉斯算子，用Conjugate_Gradient求解
        //这种代码在两个物体重叠的时候会有问题，但这里就不考虑物体重叠的情况了
        GetVirtualHeight(old_h, h, new_h, cube1);
        GetVirtualHeight(old_h, h, new_h, cube2);

        //Diminish vh. vh添加一个衰减，非物理，因为模拟本质上是显示积分
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (cg_mask[i, j])
                    vh[i, j] *= gamma;
            }
        }

        //Update new_h by vh.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (i != 0)
                    new_h[i, j] += (vh[i - 1, j] - vh[i, j]) * rate;
                if (i != size - 1)
                    new_h[i, j] += (vh[i + 1, j] - vh[i, j]) * rate;
                if (j != 0)
                    new_h[i, j] += (vh[i, j - 1] - vh[i, j]) * rate;
                if (j != size - 1)
                    new_h[i, j] += (vh[i, j + 1] - vh[i, j]) * rate;
            }
        }

        //Step 3
        //old_h <- h; h <- new_h;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                old_h[i, j] = h[i, j];
                h[i, j] = new_h[i, j];
            }
        }

        //Step 4: Water->Block coupling.
        Water2BolckCoupling(cube1);
        Water2BolckCoupling(cube2);
    }

    private void Water2BolckCoupling(GameObject cube)
    {
        float dt = 0.005f;
        float m = 200f;
        Vector3 force = new Vector3(0, -m * 9.8f, 0);
        Vector3 torque = new Vector3(0, 0, 0);

        // if (cube.TryGetComponent(out cube_motion cubeMotion) && cubeMotion.cube_move)
        // {
        //     cube.transform.rotation = Quaternion.identity;
        //     return;
        // }
        //
        // if (cube.TryGetComponent(out block_motion blockMotion) && blockMotion.block_move)
        // {
        //     cube.transform.rotation = Quaternion.identity;
        //     return;
        // }

        Vector3 position = cube.transform.position;
        GetBoundIndex(cube, out int left, out int right, out int top, out int bottom, out Vector3 worldMin, out Vector3 worldMax);
        Mesh cubeMesh = cube.GetComponent<MeshFilter>().mesh;
        Bounds bounds = cubeMesh.bounds;

        for (int i = left; i <= right; i++)
        {
            for (int j = top; j <= bottom; j++)
            {
                if (i < 0 || j < 0 || i >= size || j >= size) continue;
                if (vh[i, j] == 0) continue;

                Vector3 a1 = new Vector3(i * 0.1f - size * 0.05f, worldMin.y - 2, j * 0.1f - size * 0.05f);
                Vector3 a2 = new Vector3(i * 0.1f - size * 0.05f, worldMin.y - 1, j * 0.1f - size * 0.05f);
                a1 = cube.transform.InverseTransformPoint(a1);
                a2 = cube.transform.InverseTransformPoint(a2);

                float high = 10;
                Ray ray = new Ray(a1, a2 - a1);
                bounds.IntersectRay(ray, out high);

                float v = 0 - (worldMin.y - 2 + high); //排水量

                //根据排水量计算力
                Vector3 r = a1 + high * (a2 - a1);
                Vector3 f = new Vector3(0, v, 0) * 9.8f * 0.01f * 1000;
                force += f;
                torque += Vector3.Cross(r, f);
            }
        }

        cube_v *= 0.99f;
        cube_w *= 0.98f;
        cube_v += force * dt / m;
        position += cube_v * dt;
        cube_w += torque * dt / (m / 6);
        Quaternion rotation = cube.transform.rotation;
        Vector3 dw = cube_w * 0.5f * dt;
        rotation = QuaternionAdd(rotation, new Quaternion(dw.x, dw.y, dw.z, 0.0f) * rotation);
        cube.transform.position = position;
        cube.transform.rotation = rotation;
    }

    private void GetVirtualHeight(float[,] old_h, float[,] h, float[,] new_h, GameObject cube)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                b[i, j] = 0;
                cg_mask[i, j] = false;
            }
        }

        //提高，物体可能旋转，使用碰撞盒+射线来检测排水情况
        GetBoundIndex(cube, out int left, out int right, out int top, out int bottom, out Vector3 worldMin, out Vector3 worldMax);
        Mesh cubeMesh = cube.GetComponent<MeshFilter>().mesh;
        Bounds bounds = cubeMesh.bounds;

        for (int i = left; i <= right; i++)
        {
            for (int j = top; j <= bottom; j++)
            {
                if (i < 0 || j < 0 || i >= size || j >= size) continue;

                Vector3 a1 = new Vector3(i * 0.1f - size * 0.05f, worldMin.y - 2, j * 0.1f - size * 0.05f);
                Vector3 a2 = new Vector3(i * 0.1f - size * 0.05f, worldMin.y - 1, j * 0.1f - size * 0.05f);
                a1 = cube.transform.InverseTransformPoint(a1);
                a2 = cube.transform.InverseTransformPoint(a2);

                float high = 10;
                Ray ray = new Ray(a1, a2 - a1);
                bounds.IntersectRay(ray, out high);

                if (Math.Abs(high - (10)) > 0.001)
                {
                    low_h[i, j] = worldMin.y - 2 + high;
                }
            }
        }

        //Solve the Poisson equation to obtain vh (virtual height).
        //then set up b and cg_mask for conjugate gradient.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (low_h[i, j] < h[i, j])
                {
                    cg_mask[i, j] = true;
                    b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
                }
            }
        }

        //求解时扩大一圈，方便使用固定的拉普拉斯算子，视频中有提到
        Conjugate_Gradient(cg_mask, b, vh, left - 1, right + 1, top - 1, bottom + 1);
    }

    private void GetBoundIndex(GameObject cube, out int left, out int right, out int top, out int bottom, out Vector3 worldMin, out Vector3 worldMax)
    {
        Mesh cubeMesh = cube.GetComponent<MeshFilter>().mesh;
        Bounds bounds = cubeMesh.bounds;

        Vector3[] corners = new Vector3[8];
        corners[0] = bounds.min;
        corners[1] = new Vector3(bounds.min.x, bounds.min.y, bounds.max.z);
        corners[2] = new Vector3(bounds.min.x, bounds.max.y, bounds.min.z);
        corners[3] = new Vector3(bounds.min.x, bounds.max.y, bounds.max.z);
        corners[4] = new Vector3(bounds.max.x, bounds.min.y, bounds.min.z);
        corners[5] = new Vector3(bounds.max.x, bounds.min.y, bounds.max.z);
        corners[6] = new Vector3(bounds.max.x, bounds.max.y, bounds.min.z);
        corners[7] = bounds.max;

        worldMin = cube.transform.TransformPoint(corners[0]);
        worldMax = worldMin;

        for (int i = 1; i < corners.Length; i++)
        {
            Vector3 worldCorner = cube.transform.TransformPoint(corners[i]);
            worldMin = Vector3.Min(worldMin, worldCorner);
            worldMax = Vector3.Max(worldMax, worldCorner);
        }

        // Vector3 worldMin = cube1.transform.TransformPoint(bounds.min);
        // Vector3 worldMax = cube1.transform.TransformPoint(bounds.max);

        left = Mathf.RoundToInt((worldMin.x + 5.0f) * 10);
        right = Mathf.RoundToInt((worldMax.x + 5.0f) * 10);
        top = Mathf.RoundToInt((worldMin.z + 5.0f) * 10);
        bottom = Mathf.RoundToInt((worldMax.z + 5.0f) * 10);
    }

    private Quaternion QuaternionAdd(Quaternion q1, Quaternion q2)
    {
        return new Quaternion(q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w);
    }

    // Update is called once per frame
    void Update()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        float[,] new_h = new float[size, size];
        float[,] h = new float[size, size];

        //Load X.y into h.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                h[i, j] = X[i * size + j].y;
            }
        }

        if (Input.GetKeyDown("r"))
        {
            //Add random water.
            int i = Random.Range(0, size);
            int j = Random.Range(0, size);
            float addV = Random.value;

            int count = 4;
            if (i == 0) count--;
            if (i == size - 1) count--;
            if (j == 0) count--;
            if (j == size - 1) count--;

            h[i, j] += addV;

            if (i > 0) h[i - 1, j] -= addV / count;
            if (i < size - 1) h[i + 1, j] -= addV / count;
            if (j > 0) h[i, j - 1] -= addV / count;
            if (j < size - 1) h[i, j + 1] -= addV / count;
        }

        for (int l = 0; l < 8; l++)
        {
            Shallow_Wave(old_h, h, new_h);
        }

        //Store h back into X.y and recalculate normal.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                X[i * size + j].y = h[i, j];
            }
        }

        mesh.vertices = X;
        mesh.RecalculateNormals();
    }
}