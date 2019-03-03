// readobjfile.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"stdio.h"
#include"stdlib.h"
#include"iostream"
#include"string.h"

using namespace std;
#define NUMBER 131737
//定义点的结构体
typedef struct POINT
{
	double x,y,z;//点的三维坐标
};

int _tmain(int argc, _TCHAR* argv[])
{
	FILE *fobj;
	fobj=fopen("1516-1020-50mm - Cloud021.obj","r");
	if(fobj==NULL)
	{
		printf("the obj file open error.");
		exit(0);
	}
	POINT original_points[NUMBER];
	char v;
	int l=0;
	while(!feof(fobj))
	{
		v=fgetc(fobj);
		if(v=='v'&&l<NUMBER)			
		{
			fscanf(fobj,"%lf",&original_points[l].x);
			fscanf(fobj,"%lf",&original_points[l].y);
			fscanf(fobj,"%lf",&original_points[l].z);	
			printf("%3.7f %3.7f %3.7f\n",original_points[l].x,original_points[l].y,original_points[l].z);
			l++;
		}	
	}
	fclose(fobj);

	return 0;
}

