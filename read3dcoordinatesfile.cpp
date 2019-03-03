// read3dcoordinatesfile.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include"iostream"
#include"string.h"
#include"stdlib.h"
#include"stdio.h"
#include"math.h"

using namespace std;


//�����Ľṹ��
typedef struct POINT
{
	double x,y,z;//�����ά����
};
//����ߵĽṹ��
typedef struct EDGE
{
	POINT p_start,p_end;//�ߵ�������յ�
};
//������Ľṹ
typedef struct FACE
{
	EDGE e[4];//���蹹��ģ�͵����Ϊ�ı��Σ����ı��ε������߰���ʱ���������ʱ��ķ�������ֱƽ������
	double f_normal[3];//ƽ�淨��������ά��������
	//normal_average(f_normal,e[4]);
};

//�����ȡtxt�ļ��еĵ��Ƶ���ά����ĺ���������������������������Ľṹ����
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

	//����������ļ�
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

	//�����������ļ�
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

