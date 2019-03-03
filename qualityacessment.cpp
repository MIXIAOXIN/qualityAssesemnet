// qualityacessment.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include"iostream"
#include"string.h"
#include"stdlib.h"
#include"stdio.h"
#include"math.h"
//ԭʼ�����е����Ŀ
#define NUMBER 131737
#define thre_point 0.5  //ԭʼ�������жϵ��Ƿ����渽����������ֵ
#define thre_range 0.3  //ԭʼ�����е��Ƿ������ϵľ�����ֵ
using namespace std;


//void function_normal(double N[3],EDGE e1,EDGE e2);
//void normal_average(double normal_aver[3],EDGE e[4]);

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

//���㷨����,�ı��ε�λ����������
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
//�����ı���ƽ���ƽ��������
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
//��8�����е����ֵ
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

//��8�����е���Сֵ
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

//������ƽ��ֵ
//double function_aver(double arr[])

//����㵽��ľ���
double function_range(POINT ptest,FACE ftest)
{
	double xv,yv,zv;//����ƽ����һ����ƽ����һ�㹹�ɵ���������������
	xv=ptest.x-ftest.e[0].p_start.x;
	yv=ptest.y-ftest.e[0].p_start.y;
	zv=ptest.z-ftest.e[0].p_start.z;
	double ran,range;//ranΪ�з��ŵľ���Ĺ���ֵ��rangeΪ��������ȡ����ֵ�ľ���ֵ
	ran=xv*ftest.f_normal[0]+yv*ftest.f_normal[1]+zv*ftest.f_normal[2];
	range=fabs(ran);
	return range;//���ؼ�����ĵ㵽��ľ���ֵ
}


int _tmain(int argc, _TCHAR* argv[])
{
	//��ʼ��ģ�͵�
	FILE *fp;
	fp=fopen("points.txt","r");
	POINT p[48];
	if(fp==NULL){printf("file open error.\n");exit(0);}//������ļ�δ�����룬��ֱ�ӷ���
	for(int i=0;i<48;i++)
	{
		fscanf(fp,"%lf",&p[i].x);fscanf(fp,"%lf",&p[i].y);fscanf(fp,"%lf",&p[i].z);//���ļ��е�����ֵ���δ����ṹ����
		printf("%lf %lf %lf\n",p[i].x,p[i].y,p[i].z);//���������꣬�����Ƿ������ȷ
	}
	fclose(fp);//�ر��ļ�ָ��

	//��ʼ��ģ�ͱ���
	//����������ļ�
	FILE *fedge;
	fedge=fopen("edges_index.txt","r");
	EDGE ed[80];
	if(fedge==NULL)
	{
		printf("edges_index file open error.\n");
		exit(0);
	}
	int m[80],n[80];//�ߵ���Ŀ
	for(int j=0;j<80;j++)
	{
		fscanf(fedge,"%d",&m[j]);
		fscanf(fedge,"%d",&n[j]);
		ed[j].p_start=p[m[j]];ed[j].p_end=p[n[j]];
		printf("%lf %lf %lf %lf %lf %lf\n",ed[j].p_start.x,ed[j].p_start.y,ed[j].p_start.z,ed[j].p_end.x,ed[j].p_end.y,ed[j].p_end.z);
	}

	fclose(fedge);

	//��ʼ��ģ����
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
		normal_average(f[k].f_normal,f[k].e);//����ÿ�����ƽ����λ������
		printf("��%d��ķ�����Ϊ��%f %f %f\n",k,f[k].f_normal[0],f[k].f_normal[1],f[k].f_normal[2]);
	}	

	fclose(fface);

	//����ԭʼ���Ƶ�������ļ�
	FILE *fobj;
	fobj=fopen("1516-1020-50mm - Cloud0.21.obj","r");
	if(fobj==NULL)
	{
		printf("the obj file open error.");
		exit(0);
	}
	POINT original_points[NUMBER];//ԭʼ�����е����Ŀ
	char v;
	int l=0;
	while(!feof(fobj))
	{
		v=fgetc(fobj);
		if(v=='v'&&l<NUMBER)	//�ж�obj�ļ��������ֻ�����ĸv����ĸv����Ϊ�������		
		{
			fscanf(fobj,"%lf",&original_points[l].x);//��obj�ļ��еĵ����������δ�������еĵ��ƽṹ��������
			fscanf(fobj,"%lf",&original_points[l].y);
			fscanf(fobj,"%lf",&original_points[l].z);	
			printf("%3.7f %3.7f %3.7f\n",original_points[l].x,original_points[l].y,original_points[l].z);
			l++;
		}	
	}
	fclose(fobj);

	//����㵽��ľ��벢ͳ��һ�����֮��������,���浽ָ��������
	FILE *stastics;
	stastics=fopen("stastic_output.txt","a");
	if(stastics==NULL){printf("cannot open stastic_output.txt.");exit(0);}
	double range_aver[31]={0};
	double sum[31]={0};
	int sum_num[31]={0};
	int statistic[31][30]={0};//��ÿ����Ƭ�е㵽��ľ������е����Ŀͳ�ƽ�������ڴ�
	for(int r=0;r<31;r++)
	{
		double x[8],y[8],z[8];//�ı����ĸ������������飬���ں������Ӧ������������Сֵ
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

		for(int t=0;t<NUMBER;t++)//��ԭʼ�����е�ÿ������б���
		{
			if(original_points[t].x>(x_min-thre_point)&&original_points[t].x<(x_max+thre_point)&&
				original_points[t].y>(y_min-thre_point)&&original_points[t].y<(y_max+thre_point)&&
				original_points[t].z>(z_min-thre_point)&&original_points[t].z<(z_max+thre_point))
			{
				double range=function_range(original_points[t],f[r]);//������������ĵ㵽��ǰ��ľ���
				if(range<thre_range)  //�жϵ㵽��ľ����Ƿ����������ֵ����
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
		printf("��%d����ľ����ֲ���\n",r);
		for(int v=0;v<30;v++)
			{
				printf("%d ",statistic[r][v]);
				fprintf(stastics,"%d",statistic[r][v]);//��ÿ��������ڵ������д������ļ�stastic_output.txt��
				fputc(' ',stastics);
			}
		fputc('\n',stastics);
		for(int u=0;u<30;u++)
			{
				printf("%lf ",((double)statistic[r][u])/sum_num[r]);
				fprintf(stastics,"%lf",((double)statistic[r][u])/sum_num[r]);//��ÿ��������ڵ������д������ļ�stastic_output.txt��
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

