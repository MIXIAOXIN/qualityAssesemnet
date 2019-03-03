// reconstruction.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"string.h"
#include"iostream"
#include"stdlib.h"
#include"stdio.h"

typedef int ElemType;
using namespace std;
#define List_Init_Size 10//线性表存储空间的初始分配量
#define List_Increment 5//线性表存储空间增量

struct POINT
{
	double x,y,z;//点的三维坐标
};

struct EDGE
{
	struct POINT p_start,p_end;
};


//构成面的边的顺序链表
struct FACE
{
	int *elem;//空间存储基址
	int length;//当前存贮长度即当前面的边数
	int listsize;//当前分配的存贮容量即最大多边形边数
};
void InitList(FACE &L)//构造一个空的顺序链表
{
	L.elem =(ElemType*)malloc(List_Init_Size*sizeof(ElemType));
	if(!L.elem )
		exit(0);//内存分配失败
	else 
		L.length =0;//空表长为0
		L.listsize =List_Init_Size;//初始存贮容量

}
int ListInsert(FACE &L,ElemType e)//在顺序链表
{
	int *p;
	p=L.elem ;
	*p=e;
	L.length ++;
	L.elem ++;
	return 1;
}

int _tmain(int argc, _TCHAR* argv[])
{
	struct POINT p0={0,0,0},p1={1,0,0},p2={1,1,0},p3={0,1,0},p4={0,0,1},p5={1,0,1},p6={1,1,1},p7={0,1,1};
	struct EDGE l0={p0,p1},l1={p1,p2},l2={p2,p3},l3={p3,p0},l4={p0,p4},l5={p1,p5},l6={p2,p6},l7={p3,p7},l8={p4,p5},l9={p5,p6},l10={p6,p7},l11={p7,p4};
	struct FACE f1;
	InitList(f1);
	printf("初始化后：f1.length=%d f1.listsize=%d",f1.length ,f1.listsize );
	printf("请输入当前面的边数：");
	int num_poly,n;
	//scanf("%d",&num_poly);
	cin>>num_poly;
	for(int i=0;i<num_poly;i++)
	{
		printf("请输入当前多边形边号：");
		//scanf("%d",&n);
		cin>>n;
		ListInsert(f1,n);
	}

	for(int j=0;j<num_poly;j++)
	{
		printf("第%d条边为：%d",j,*(f1.elem -f1.length +j));
	}
	


	return 0;
}

