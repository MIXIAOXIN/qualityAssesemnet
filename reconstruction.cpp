// reconstruction.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include"string.h"
#include"iostream"
#include"stdlib.h"
#include"stdio.h"

typedef int ElemType;
using namespace std;
#define List_Init_Size 10//���Ա�洢�ռ�ĳ�ʼ������
#define List_Increment 5//���Ա�洢�ռ�����

struct POINT
{
	double x,y,z;//�����ά����
};

struct EDGE
{
	struct POINT p_start,p_end;
};


//������ıߵ�˳������
struct FACE
{
	int *elem;//�ռ�洢��ַ
	int length;//��ǰ�������ȼ���ǰ��ı���
	int listsize;//��ǰ����Ĵ���������������α���
};
void InitList(FACE &L)//����һ���յ�˳������
{
	L.elem =(ElemType*)malloc(List_Init_Size*sizeof(ElemType));
	if(!L.elem )
		exit(0);//�ڴ����ʧ��
	else 
		L.length =0;//�ձ�Ϊ0
		L.listsize =List_Init_Size;//��ʼ��������

}
int ListInsert(FACE &L,ElemType e)//��˳������
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
	printf("��ʼ����f1.length=%d f1.listsize=%d",f1.length ,f1.listsize );
	printf("�����뵱ǰ��ı�����");
	int num_poly,n;
	//scanf("%d",&num_poly);
	cin>>num_poly;
	for(int i=0;i<num_poly;i++)
	{
		printf("�����뵱ǰ����αߺţ�");
		//scanf("%d",&n);
		cin>>n;
		ListInsert(f1,n);
	}

	for(int j=0;j<num_poly;j++)
	{
		printf("��%d����Ϊ��%d",j,*(f1.elem -f1.length +j));
	}
	


	return 0;
}

