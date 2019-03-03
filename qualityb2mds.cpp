// qualityb2mds.cpp : 定义控制台应用程序的入口点。
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

//空间两点之间的距离
double distance_points(POINT rp,POINT tp)
{
	double dist,d;
	d=(rp.x-tp.x)*(rp.x-tp.x)+(rp.y-tp.y)*(rp.y-tp.y)+(rp.z-tp.z)*(rp.z-tp.z);
	return dist=sqrt(d);
}

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

//计算两个平面的夹角
double function_includedangle(FACE f1,FACE f2)
{
	double cd,angle;
	cd=f1.f_normal[0]*f2.f_normal[0]+f1.f_normal[1]*f2.f_normal[1]+f1.f_normal[2]*f2.f_normal[2];
	angle=acos(fabs(cd));
	return angle;
}


int _tmain(int argc, _TCHAR* argv[])
{
	//初始化标准参考模型
	FILE *fp;
	fp=fopen("reference_points.txt","r");
	POINT rp[74];//参考模型中共有74个特征点
	if(fp==NULL){printf("file open error.\n");exit(0);}//如果点文件未被读入，则直接返回
	for(int i=0;i<74;i++)
	{
		fscanf(fp,"%lf",&rp[i].x);fscanf(fp,"%lf",&rp[i].y);fscanf(fp,"%lf",&rp[i].z);//将文件中的坐标值依次存入点结构体中
		printf("%lf %lf %lf\n",rp[i].x,rp[i].y,rp[i].z);//输出点的坐标，检验是否读入正确
	}
	fclose(fp);//关闭文件指针

	//初始化模型边线
	//读入边索引文件
	FILE *fedge;
	fedge=fopen("reference_edges.txt","r");
	EDGE ed[103];//参考模型中共有103条边
	if(fedge==NULL)
	{
		printf("edges_index file open error.\n");
		exit(0);
	}
	int m[103],n[103];//边的数目
	for(int j=0;j<103;j++)
	{
		fscanf(fedge,"%d",&m[j]);
		fscanf(fedge,"%d",&n[j]);
		ed[j].p_start=rp[m[j]];ed[j].p_end=rp[n[j]];
		printf("%lf %lf %lf %lf %lf %lf\n",ed[j].p_start.x,ed[j].p_start.y,ed[j].p_start.z,ed[j].p_end.x,ed[j].p_end.y,ed[j].p_end.z);
	}

	fclose(fedge);

	//初始化模型面
	FACE f[39];//参考模型中共有39个四边形平面
	FILE *fface;
	fface=fopen("reference_faces.txt","r");
	if(fface==NULL)
	{
		printf("faces_index file open error.\n");
		exit(0);
	}
	int a,b,c,d;
	for(int k=0;k<39;k++)
	{
		fscanf(fface,"%d",&a);fscanf(fface,"%d",&b);fscanf(fface,"%d",&c);fscanf(fface,"%d",&d);
		f[k].e[0]=ed[a];f[k].e[1]=ed[b];f[k].e[2]=ed[c];f[k].e[3]=ed[d];
		normal_average(f[k].f_normal,f[k].e);//计算每个面的平均单位法向量
		printf("此%d面的法向量为：%f %f %f\n",k,f[k].f_normal[0],f[k].f_normal[1],f[k].f_normal[2]);
	}	

	fclose(fface);

	//读入目标评价模型的点线面参数
	FILE *tpp;
	tpp=fopen("target_points.txt","r");
	POINT tp[14];//待评价目标模型的特征点数目为14
	if(tpp==NULL){printf("file open error.\n");exit(0);}//如果点文件未被读入，则直接返回
	for(int i=0;i<14;i++)
	{
		fscanf(tpp,"%lf",&tp[i].x);fscanf(tpp,"%lf",&tp[i].y);fscanf(tpp,"%lf",&tp[i].z);//将文件中的坐标值依次存入点结构体中
		printf("%lf %lf %lf\n",tp[i].x,tp[i].y,tp[i].z);//输出点的坐标，检验是否读入正确
	}
	fclose(tpp);//关闭文件指针

	//初始化模型边线
	//读入边索引文件
	FILE *ftedge;
	ftedge=fopen("target_edges.txt","r");
	EDGE ted[19];
	if(ftedge==NULL)
	{
		printf("edges_index file open error.\n");
		exit(0);
	}
	int mt[19],nt[19];//边的数目
	for(int j=0;j<19;j++)
	{
		fscanf(ftedge,"%d",&mt[j]);
		fscanf(ftedge,"%d",&nt[j]);
		ted[j].p_start=tp[mt[j]];ted[j].p_end=tp[nt[j]];
		printf("%lf %lf %lf %lf %lf %lf\n",ted[j].p_start.x,ted[j].p_start.y,ted[j].p_start.z,ted[j].p_end.x,ted[j].p_end.y,ted[j].p_end.z);
	}

	fclose(ftedge);

	//初始化模型面
	FACE tf[6];
	FILE *ftface;
	ftface=fopen("target_faces.txt","r");
	if(ftface==NULL)
	{
		printf("faces_index file open error.\n");
		exit(0);
	}
	int at,bt,ct,dt;
	for(int k=0;k<6;k++)
	{
		fscanf(ftface,"%d",&at);fscanf(ftface,"%d",&bt);fscanf(ftface,"%d",&ct);fscanf(ftface,"%d",&dt);
		tf[k].e[0]=ted[at];tf[k].e[1]=ted[bt];tf[k].e[2]=ted[ct];tf[k].e[3]=ted[dt];
		normal_average(tf[k].f_normal,tf[k].e);//计算每个面的平均单位法向量
		printf("此%d面的法向量为：%f %f %f\n",k,tf[k].f_normal[0],tf[k].f_normal[1],tf[k].f_normal[2]);
	}	

	fclose(ftface);

	//用读入的数据进行对应面的相关参数计算,当前只计算建筑物的大框架下的几何精度，忽略建筑物的细节，如：阳台、门窗和烟囱;并将计算结果输出到difference.txt文件中
	FILE *differ;
	differ=fopen("difference.txt","a");
	if(differ==NULL){printf("cannot open difference.txt file");exit(0);}
	//对应特征基元点的相关计算：角点距离偏移真值及平均值
	double p_offset[14],p_offset_aver,sum_dist=0;
	for(int ap=0;ap<14;ap++)
	{
		p_offset[ap]=distance_points(rp[ap],tp[ap]);
		sum_dist=sum_dist+p_offset[ap];
		fprintf(differ,"%lf",p_offset[ap]);
		fputc(' ',differ);
	}
	fputc('\n',differ);
	fputs("average distance",differ);
	p_offset_aver=sum_dist/14;
	fprintf(differ,"%lf",p_offset_aver);
	fputc('\n',differ);

	//计算对应基元线的差异：只考虑直线段
	//计算长度差值
	double dist_differ[19],dist_differ_aver,dist_differ_sum=0;
	for(int ad=0;ad<19;ad++)
	{
		dist_differ[ad]=fabs(distance_points(ed[ad].p_start,ed[ad].p_end)-distance_points(ted[ad].p_start,ted[ad].p_end));
		dist_differ_sum=dist_differ_sum+dist_differ[ad];
		fprintf(differ,"%lf",dist_differ[ad]);
		fputc(' ',differ);
	}
	fputc('\n',differ);
	fputs("average distance difference of edge",differ);
	dist_differ_aver=dist_differ_sum/19;
	fprintf(differ,"%lf",dist_differ_aver);
	fputc('\n',differ);
	//计算距离偏移
	POINT rmidpoint[19],tmidpoint[19];//分别为参考真值模型中边的中点和目标评价模型中边的中点
	double dist_offset[19],dist_offset_aver,dist_offset_sum=0;
	for(int am=0;am<19;am++)
	{
		rmidpoint[am].x=(ed[am].p_start.x+ed[am].p_end.x)/2;rmidpoint[am].y=(ed[am].p_start.y+ed[am].p_end.y)/2;rmidpoint[am].z=(ed[am].p_start.z+ed[am].p_end.z)/2;
		tmidpoint[am].x=(ted[am].p_start.x+ted[am].p_end.x)/2;tmidpoint[am].y=(ted[am].p_start.y+ted[am].p_end.y)/2;tmidpoint[am].z=(ted[am].p_start.z+ted[am].p_end.z)/2;
		dist_offset[am]=fabs((distance_points(ed[am].p_start,ted[am].p_start)+distance_points(ed[am].p_end,ted[am].p_end)+distance_points(rmidpoint[am],tmidpoint[am]))/3);
		dist_offset_sum=dist_offset_sum+dist_offset[am];
		fprintf(differ,"%lf",dist_offset[am]);
		fputc(' ',differ);
	}
	fputc('\n',differ);
	fputs("distance offset of edge",differ);
	dist_offset_aver=dist_offset_sum/19;
	fprintf(differ,"%lf",dist_offset_aver);
	fputc('\n',differ);
	
	//计算对应边的角度偏移
	double fenzi=0.,fenmu=0.,degree,sum=0,ct_aver;
	for(int ac=0;ac<19;ac++)
	{
		fenzi=(ed[ac].p_start.x-ed[ac].p_end.x)*(ted[ac].p_start.x-ted[ac].p_end.x)+(ed[ac].p_start.y-ed[ac].p_end.y)*(ted[ac].p_start.y-ted[ac].p_end.y)+(ed[ac].p_start.z-ed[ac].p_end.z)*(ted[ac].p_start.z-ted[ac].p_end.z);
		fenmu=distance_points(ed[ac].p_start,ed[ac].p_end)*distance_points(ted[ac].p_start,ted[ac].p_end);
		degree=acos(fabs(fenzi/fenmu));
		sum=sum+degree;
		fprintf(differ,"%lf",degree);
		fputc(' ',differ);
	}
	fputc('\n',differ);
	fputs("angle offset of edge",differ);
	ct_aver=sum/19;
	fprintf(differ,"%lf",ct_aver);
	fputc('\n',differ);

	//基元面评价
	//对应面的面积差值
	double r_area[6],t_area[6],differ_area,aver_darea,sum_darea=0;
	double adist_face[6],dist_pf1,dist_pf2,sum1=0,sum2=0,aver_all;//计算两面距离的相应参数
	double angle[6],sum_a=0,aver_a;//计算两平面的夹角的相关参数
	for(int aa=0;aa<6;aa++)
	{
		double x[4],y[4],z[4],xt[4],yt[4],zt[4];//每个四边形面的边的向量，为计算面积时使用
		for(int bb=0;bb<4;bb++)
		{
			x[bb]=f[aa].e[bb].p_end.x-f[aa].e[bb].p_start.x;
			y[bb]=f[aa].e[bb].p_end.y-f[aa].e[bb].p_start.y;
			z[bb]=f[aa].e[bb].p_end.z-f[aa].e[bb].p_start.z;
			xt[bb]=tf[aa].e[bb].p_end.x-tf[aa].e[bb].p_start.x;
			yt[bb]=tf[aa].e[bb].p_end.y-tf[aa].e[bb].p_start.y;
			zt[bb]=tf[aa].e[bb].p_end.z-tf[aa].e[bb].p_start.z;
		}
		r_area[aa]=(std::sqrt((y[0]*z[1]-z[0]*y[1])*(y[0]*z[1]-z[0]*y[1])+(z[0]*x[1]-x[0]*z[1])*(z[0]*x[1]-x[0]*z[1])+(x[0]*y[1]-y[0]*x[1])*(x[0]*y[1]-y[0]*x[1]))+
					std::sqrt((y[2]*z[3]-z[2]*y[3])*(y[2]*z[3]-z[2]*y[3])+(z[2]*x[3]-x[2]*z[3])*(z[2]*x[3]-x[2]*z[3])+(x[2]*y[3]-y[2]*x[3])*(x[2]*y[3]-y[2]*x[3])))/2;
		t_area[aa]=(std::sqrt((yt[0]*zt[1]-zt[0]*yt[1])*(yt[0]*zt[1]-zt[0]*yt[1])+(zt[0]*xt[1]-xt[0]*zt[1])*(zt[0]*xt[1]-xt[0]*zt[1])+(xt[0]*yt[1]-yt[0]*xt[1])*(xt[0]*yt[1]-yt[0]*xt[1]))+
					std::sqrt((yt[2]*zt[3]-zt[2]*yt[3])*(yt[2]*zt[3]-zt[2]*yt[3])+(zt[2]*xt[3]-xt[2]*zt[3])*(zt[2]*xt[3]-xt[2]*zt[3])+(xt[2]*yt[3]-yt[2]*xt[3])*(xt[2]*yt[3]-yt[2]*xt[3])))/2;
		differ_area=fabs(r_area[aa]-t_area[aa]);
		sum_darea=sum_darea+differ_area;
		fprintf(differ,"%lf",differ_area);
		fputc(' ',differ);

		//两面距离
		for(int cc=0;cc<4;cc++)
		{
			dist_pf1=function_range(tf[aa].e[cc].p_start,f[aa]);
			dist_pf2=function_range(tf[aa].e[cc].p_end,f[aa]);
			sum1=sum1+dist_pf1+dist_pf2;
		}
		adist_face[aa]=sum1/8;
		sum1=0;
		sum2=sum2+adist_face[aa];

		//计算两面夹角
		angle[aa]=function_includedangle(f[aa],tf[aa]);
		sum_a=sum_a+angle[aa];
	}
	aver_darea=sum_darea/6;
	fputs("average difference of face area",differ);
	fprintf(differ,"%lf",aver_darea);
	fputc('\n',differ);//以上为计算输出对应面的面积的差值

	//以下为计算输出两面距离的数据
	aver_all=sum2/6;
	for(int dd=0;dd<6;dd++)
	{fprintf(differ,"%lf",adist_face[dd]);fputc(' ',differ);}
	fputc('\n',differ);
	fputs("average distance among all faces",differ);
	fprintf(differ,"%lf",aver_all);
	fputc('\n',differ);
		
	//以下为输出两面夹角的相关数据
	aver_a=sum_a/6;
	for(int ee=0;ee<6;ee++)
	{fprintf(differ,"%lf",angle[ee]);fputc(' ',differ);}
	fputc('\n',differ);
	fputs("average angle between allfaces",differ);
	fprintf(differ,"%lf",aver_a);
	fputc('\n',differ);

	//计算房屋规则细节的体积或面积，根据数据的一致性判断重建的细节质量
	//对于阳台，有体的特征，则计算体积
	double length,width,depth,volume[3],av_volume_balcony,variance_balcony;
	length=(distance_points(rp[38],rp[44])+distance_points(rp[40],rp[42]))/2;
	width=(distance_points(rp[38],rp[40])+distance_points(rp[44],rp[42]))/2;
	depth=(distance_points(rp[38],rp[39])+distance_points(rp[40],rp[41])+distance_points(rp[42],rp[43])+distance_points(rp[44],rp[45]))/4;
	volume[0]=length*width*depth;
	length=(distance_points(rp[46],rp[52])+distance_points(rp[48],rp[50]))/2;
	width=(distance_points(rp[46],rp[48])+distance_points(rp[50],rp[52]))/2;
	depth=(distance_points(rp[46],rp[47])+distance_points(rp[48],rp[49])+distance_points(rp[50],rp[51])+distance_points(rp[52],rp[53]))/4;
	volume[1]=length*width*depth;
	length=(distance_points(rp[54],rp[60])+distance_points(rp[56],rp[58]))/2;
	width=(distance_points(rp[54],rp[56])+distance_points(rp[58],rp[60]))/2;
	depth=(distance_points(rp[54],rp[55])+distance_points(rp[56],rp[57])+distance_points(rp[58],rp[59])+distance_points(rp[60],rp[61]))/4;
	volume[2]=length*width*depth;
	av_volume_balcony=(volume[0]+volume[1]+volume[2])/3;//体积的均值
	//体积的方差
	variance_balcony=sqrt((volume[0]-av_volume_balcony)*(volume[0]-av_volume_balcony)+(volume[1]-av_volume_balcony)*(volume[1]-av_volume_balcony)+(volume[2]-av_volume_balcony)*(volume[2]-av_volume_balcony));
	//输出三个阳台的各自的体积、均值、方差
	fputs("the volume of balcony",differ);
	for(int ii=0;ii<3;ii++)
	{fprintf(differ,"%lf",volume[ii]);fputc(' ',differ);}
	fputc('\n',differ);fputs("the average volume of balcony",differ);fprintf(differ,"%lf",av_volume_balcony);
	fputc('\n',differ);fputs("the variance volume of balcony",differ);fprintf(differ,"%lf",variance_balcony);
	fputc('\n',differ);

	//对于烟囱，有体的特征，计算体积
	double len_ch,wid_ch,dep_ch1,dep_ch2,volume_ch[3],aver_volume_ch,variance_ch;
	len_ch=(distance_points(rp[14],rp[20])+distance_points(rp[16],rp[18]))/2;
	wid_ch=(distance_points(rp[14],rp[16])+distance_points(rp[20],rp[18]))/2;
	dep_ch1=(distance_points(rp[16],rp[17])+distance_points(rp[18],rp[19]))/2;
	dep_ch2=(distance_points(rp[14],rp[15])+distance_points(rp[20],rp[21]))/2;
	volume_ch[0]=len_ch*wid_ch*(dep_ch1+(dep_ch2-dep_ch1)/2);//第一个烟囱的体积计算完毕
	len_ch=(distance_points(rp[22],rp[28])+distance_points(rp[24],rp[26]))/2;
	wid_ch=(distance_points(rp[22],rp[24])+distance_points(rp[26],rp[28]))/2;
	dep_ch1=(distance_points(rp[24],rp[25])+distance_points(rp[26],rp[27]))/2;
	dep_ch2=(distance_points(rp[22],rp[23])+distance_points(rp[28],rp[29]))/2;
	volume_ch[1]=len_ch*wid_ch*(dep_ch1+(dep_ch2-dep_ch1)/2);//第二个烟囱的体积计算完毕
	len_ch=(distance_points(rp[30],rp[36])+distance_points(rp[32],rp[34]))/2;
	wid_ch=(distance_points(rp[30],rp[32])+distance_points(rp[34],rp[36]))/2;
	dep_ch1=(distance_points(rp[32],rp[33])+distance_points(rp[34],rp[35]))/2;
	dep_ch2=(distance_points(rp[30],rp[31])+distance_points(rp[36],rp[37]))/2;
	volume_ch[2]=len_ch*wid_ch*(dep_ch1+(dep_ch2-dep_ch1)/2);//第三个烟囱的体积计算完毕
	aver_volume_ch=(volume_ch[0]+volume_ch[1]+volume_ch[2])/3;//三个烟囱体积的均值
	//三个烟囱体积的均方差根
	variance_ch=sqrt((volume_ch[0]-aver_volume_ch)*(volume_ch[0]-aver_volume_ch)+(volume_ch[1]-aver_volume_ch)*(volume_ch[1]-aver_volume_ch)+(volume_ch[2]-aver_volume_ch)*(volume_ch[2]-aver_volume_ch));
	//输出三个烟囱的各自的体积、均值、方差
	fputs("the volume of chimny",differ);
	for(int ij=0;ij<3;ij++)
	{fprintf(differ,"%lf",volume_ch[ij]);fputc(' ',differ);}
	fputc('\n',differ);fputs("the average volume of chimny",differ);fprintf(differ,"%lf",aver_volume_ch);
	fputc('\n',differ);fputs("the variance volume of chimny",differ);fprintf(differ,"%lf",variance_ch);
	fputc('\n',differ);

	//计算三个门的表面积，没有体的特征
	double wid_door,height_door,area_door[3],av_area_door,variance_door;
	wid_door=(distance_points(rp[62],rp[65])+distance_points(rp[63],rp[64]))/2;
	height_door=(distance_points(rp[62],rp[63])+distance_points(rp[64],rp[65]))/2;
	area_door[0]=wid_door*height_door;
	wid_door=(distance_points(rp[66],rp[69])+distance_points(rp[67],rp[68]))/2;
	height_door=(distance_points(rp[66],rp[67])+distance_points(rp[68],rp[69]))/2;
	area_door[1]=wid_door*height_door;
	wid_door=(distance_points(rp[70],rp[73])+distance_points(rp[71],rp[72]))/2;
	height_door=(distance_points(rp[70],rp[71])+distance_points(rp[72],rp[73]))/2;
	area_door[2]=wid_door*height_door;
	av_area_door=(area_door[0]+area_door[1]+area_door[2])/3;
	variance_door=sqrt((area_door[0]-av_area_door)*(area_door[0]-av_area_door)+(area_door[1]-av_area_door)*(area_door[1]-av_area_door)+(area_door[2]-av_area_door)*(area_door[2]-av_area_door));
	//输出三个门的各自的面积、均值、方差
	fputs("the area of door",differ);
	for(int ik=0;ik<3;ik++)
	{fprintf(differ,"%lf",area_door[ik]);fputc(' ',differ);}
	fputc('\n',differ);fputs("the average area of door",differ);fprintf(differ,"%lf",av_area_door);
	fputc('\n',differ);fputs("the variance area of door",differ);fprintf(differ,"%lf",variance_door);
	fputc('\n',differ);


	return 0;
}

