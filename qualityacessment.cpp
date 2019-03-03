// qualityacessment.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"iostream"
#include"string.h"
#include"stdlib.h"
#include"stdio.h"
#include"math.h"
//原始点云中点的数目
#define NUMBER 131737
#define thre_point 0.5  //原始点云中判断点是否在面附近的坐标阈值
#define thre_range 0.3  //原始点云中点是否是面上的距离阈值
using namespace std;


//void function_normal(double N[3],EDGE e1,EDGE e2);
//void normal_average(double normal_aver[3],EDGE e[4]);

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

//计算法向量,四边形单位法向量函数
void function_normal(double N[3],EDGE e1,EDGE e2)
{
	double x1,y1,z1,x2,y2,z2,x,y,z;
	x1=e1.p_end .x -e1.p_start.x;
	y1=e1.p_end .y -e1.p_start .y ;
	z1=e1.p_end .z -e1.p_start .z ;
	x2=e2.p_end .x -e2.p_start.x;
	y2=e2.p_end .y -e2.p_start .y ;
	z2=e2.p_end .z -e2.p_start .z ;
	x=y1*z2-z1*y2;
	y=z1*x2-x1*z2;
	z=x1*y2-x2*y1;
	double normal_length=sqrt(x*x+y*y+z*z);
	N[0]=x/normal_length;
	N[1]=y/normal_length;
	N[2]=z/normal_length;
}
//计算四边形平面的平均法向量
void normal_average(double normal_aver[3],EDGE e[4]) 
{
	double nor1[3],nor2[3];
	function_normal(nor1,e[0],e[1]);
	function_normal(nor2,e[2],e[3]);
	double xx,yy,zz,nor_l;
	xx=nor1[0]+nor2[0];yy=nor1[1]+nor2[1];zz=nor1[2]+nor2[2];
	nor_l=sqrt(xx*xx+yy*yy+zz*zz);
	normal_aver[0]=xx/nor_l;
	normal_aver[1]=yy/nor_l;
	normal_aver[2]=zz/nor_l;
}
//求8个数中的最大值
double function_max(double arr[8])
{
	double max;
	max=arr[0];
	for(int mx=1;mx<8;mx++)
	{	{if(max<arr[mx])
		max=arr[mx];}
	}
	return max;
}

//求8个数中的最小值
double function_min(double arr[8])
{
	double min;
	min=arr[0];
	for(int mn=1;mn<8;mn++)
	{
		if(min>arr[mn]) min=arr[mn];
	}
	return min;
}

//求数组平均值
//double function_aver(double arr[])

//计算点到面的距离
double function_range(POINT ptest,FACE ftest)
{
	double xv,yv,zv;//定义平面外一点与平面内一点构成的向量的三个分量
	xv=ptest.x-ftest.e[0].p_start.x;
	yv=ptest.y-ftest.e[0].p_start.y;
	zv=ptest.z-ftest.e[0].p_start.z;
	double ran,range;//ran为有符号的距离的过渡值，range为向量计算取绝对值的距离值
	ran=xv*ftest.f_normal[0]+yv*ftest.f_normal[1]+zv*ftest.f_normal[2];
	range=fabs(ran);
	return range;//返回计算出的点到面的距离值
}


int _tmain(int argc, _TCHAR* argv[])
{
	//初始化模型点
	FILE *fp;
	fp=fopen("points.txt","r");
	POINT p[48];
	if(fp==NULL){printf("file open error.\n");exit(0);}//如果点文件未被读入，则直接返回
	for(int i=0;i<48;i++)
	{
		fscanf(fp,"%lf",&p[i].x);fscanf(fp,"%lf",&p[i].y);fscanf(fp,"%lf",&p[i].z);//将文件中的坐标值依次存入点结构体中
		printf("%lf %lf %lf\n",p[i].x,p[i].y,p[i].z);//输出点的坐标，检验是否读入正确
	}
	fclose(fp);//关闭文件指针

	//初始化模型边线
	//读入边索引文件
	FILE *fedge;
	fedge=fopen("edges_index.txt","r");
	EDGE ed[80];
	if(fedge==NULL)
	{
		printf("edges_index file open error.\n");
		exit(0);
	}
	int m[80],n[80];//边的数目
	for(int j=0;j<80;j++)
	{
		fscanf(fedge,"%d",&m[j]);
		fscanf(fedge,"%d",&n[j]);
		ed[j].p_start=p[m[j]];ed[j].p_end=p[n[j]];
		printf("%lf %lf %lf %lf %lf %lf\n",ed[j].p_start.x,ed[j].p_start.y,ed[j].p_start.z,ed[j].p_end.x,ed[j].p_end.y,ed[j].p_end.z);
	}

	fclose(fedge);

	//初始化模型面
	FACE f[31];
	FILE *fface;
	fface=fopen("faces_index.txt","r");
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
		normal_average(f[k].f_normal,f[k].e);//计算每个面的平均单位法向量
		printf("此%d面的法向量为：%f %f %f\n",k,f[k].f_normal[0],f[k].f_normal[1],f[k].f_normal[2]);
	}	

	fclose(fface);

	//读入原始点云的坐标点文件
	FILE *fobj;
	fobj=fopen("1516-1020-50mm - Cloud0.21.obj","r");
	if(fobj==NULL)
	{
		printf("the obj file open error.");
		exit(0);
	}
	POINT original_points[NUMBER];//原始点云中点的数目
	char v;
	int l=0;
	while(!feof(fobj))
	{
		v=fgetc(fobj);
		if(v=='v'&&l<NUMBER)	//判断obj文件中是数字还是字母v，字母v后则为点的坐标		
		{
			fscanf(fobj,"%lf",&original_points[l].x);//将obj文件中的点云坐标依次存入程序中的点云结构体数组中
			fscanf(fobj,"%lf",&original_points[l].y);
			fscanf(fobj,"%lf",&original_points[l].z);	
			printf("%3.7f %3.7f %3.7f\n",original_points[l].x,original_points[l].y,original_points[l].z);
			l++;
		}	
	}
	fclose(fobj);

	//计算点到面的距离并统计一定间隔之间点的数量,并存到指定数组中
	FILE *stastics;
	stastics=fopen("stastic_output.txt","a");
	if(stastics==NULL){printf("cannot open stastic_output.txt.");exit(0);}
	double range_aver[31]={0};
	double sum[31]={0};
	int sum_num[31]={0};
	int statistic[31][30]={0};//给每个面片中点到面的距离间隔中点的数目统计结果分配内存
	for(int r=0;r<31;r++)
	{
		double x[8],y[8],z[8];//四边形四个顶点坐标数组，用于后续求对应坐标轴的最大、最小值
		for(int s=0;s<4;s++)
		{
			x[s]=f[r].e[s].p_start.x;
			y[s]=f[r].e[s].p_start.y;
			z[s]=f[r].e[s].p_start.z;
			x[s+4]=f[r].e[s].p_end.x;
			y[s+4]=f[r].e[s].p_end.y;
			z[s+4]=f[r].e[s].p_end.z;
		}
		double x_max=function_max(x);double x_min=function_min(x);
		double y_max=function_max(y);double y_min=function_min(y);
		double z_max=function_max(z);double z_min=function_min(z);

		for(int t=0;t<NUMBER;t++)//对原始点云中的每个点进行遍历
		{
			if(original_points[t].x>(x_min-thre_point)&&original_points[t].x<(x_max+thre_point)&&
				original_points[t].y>(y_min-thre_point)&&original_points[t].y<(y_max+thre_point)&&
				original_points[t].z>(z_min-thre_point)&&original_points[t].z<(z_max+thre_point))
			{
				double range=function_range(original_points[t],f[r]);//求出满足条件的点到当前面的距离
				if(range<thre_range)  //判断点到面的距离是否满足距离阈值条件
				{
					sum[r]=sum[r]+range;
					sum_num[r]++;
					for(int u=0;u<30;u++)
					{
						if(range<0.01*(u+1))
							{statistic[r][u]++;
							break;}
					}
				}

			}
		}
		printf("第%d个面的距离差分布：\n",r);
		for(int v=0;v<30;v++)
			{
				printf("%d ",statistic[r][v]);
				fprintf(stastics,"%d",statistic[r][v]);//将每个距离差内点的数量写到输出文件stastic_output.txt中
				fputc(' ',stastics);
			}
		fputc('\n',stastics);
		for(int u=0;u<30;u++)
			{
				printf("%lf ",((double)statistic[r][u])/sum_num[r]);
				fprintf(stastics,"%lf",((double)statistic[r][u])/sum_num[r]);//将每个距离差内点的数量写到输出文件stastic_output.txt中
				fputc(' ',stastics);
			}
		fputc('\n',stastics);
		fputs("range_aver:",stastics);
		range_aver[r]=sum[r]/sum_num[r];
		fprintf(stastics,"%lf",range_aver[r]);
		fputc('\n',stastics);
	}
	
	fclose(stastics);
		
	return 0;

}

