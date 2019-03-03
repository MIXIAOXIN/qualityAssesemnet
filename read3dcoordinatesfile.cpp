// read3dcoordinatesfile.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"iostream"
#include"string.h"
#include"stdlib.h"
#include"stdio.h"
#include"math.h"

using namespace std;


//定义点的结构体
typedef struct POINT
{
	double x,y,z;//点的三维坐标
};
//定义边的结构体
typedef struct EDGE
{
	POINT p_start,p_end;//边的起点与终点
};
//定义面的结构
typedef struct FACE
{
	EDGE e[4];//假设构成模型的面均为四边形，且四边形的四条边按逆时针存贮，此时面的法向量垂直平面向外
	double f_normal[3];//平面法向量的三维坐标数组
	//normal_average(f_normal,e[4]);
};

//定义读取txt文件中的点云的三维坐标的函数，并按坐标次序存贮到特征点的结构体中
//void read_coordinates()


int _tmain(int argc, _TCHAR* argv[])
{
	FILE *fp;
	fp=fopen("points.txt","r");
	POINT p[48];
	if(fp==NULL)
	{
		printf("file open error.\n");
		exit(0);
	}
	for(int i=0;i<48;i++)
		//{fread(&p[i],sizeof(struct POINT),1,fp);
		{fscanf(fp,"%lf",&p[i].x);
		fscanf(fp,"%lf",&p[i].y);fscanf(fp,"%lf",&p[i].z);
		printf("%lf %lf %lf\n",p[i].x,p[i].y,p[i].z);}

	fclose(fp);

	//读入边索引文件
	FILE *fedge;
	fedge=fopen("edges_index.txt","r");
	EDGE ed[80];
	if(fedge==NULL)
	{
		printf("edges_index file open error.\n");
		exit(0);
	}
	int m[80],n[80];
	for(int j=0;j<80;j++)
	{
		fscanf(fedge,"%d",&m[j]);
		fscanf(fedge,"%d",&n[j]);
		ed[j].p_start=p[m[j]];ed[j].p_end=p[n[j]];
		printf("%lf %lf %lf %lf %lf %lf\n",ed[j].p_start.x,ed[j].p_start.y,ed[j].p_start.z,ed[j].p_end.x,ed[j].p_end.y,ed[j].p_end.z);
	}

	fclose(fedge);

	//读入面索引文件
	FILE *fface;
	fface=fopen("faces_index.txt","r");
	FACE f[31];
	if(fface==NULL)
	{
		printf("faces_index file open error.\n");
		exit(0);
	}
	int a,b,c,d;
	for(int k=0;k<31;k++)
	{
		fscanf(fface,"%d",&a);fscanf(fface,"%d",&b);fscanf(fface,"%d",&c);fscanf(fface,"%d",&d);
		f[k].e[0]=ed[a];f[k].e[1]=ed[b];f[k].e[2]=ed[c];f[k].e[3]=ed[d];
		normal_average(f[k].f_normal,f[k].e)

	}

	//waitKey(0);
	
	return 0;
}

