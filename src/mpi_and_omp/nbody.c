#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"mpi.h"
#include<omp.h>

typedef struct
{
	double x;
	double y;
}Point;
typedef struct
{
	double m;
	double fx;
	double fy;
	double px;
	double py;
	double vx;
	double vy;
}Bodies;
double dt;
int n;
double c;
MPI_Datatype MPI_BODY;
MPI_Status status;
int steps=1000;
//int steps;
int *start;
int *end;
int numprocs;
int myid;
double G;
Point size;
Point delv;
Point delp;
Bodies *body;
Point *temPoint;

void init()
{
	int k=0;
	int j=0;
	FILE *file;
	/*���ļ�*/
	file=fopen("//home//200820110945//mpiandopenmp//sample_input.in","r");
	if(!file)
	{
		printf("cannot open file\n");
	}
	else
	{
		/*���ļ��е�����n,G,size.x,size.y,dt,c*/
		fscanf(file,"%d%lf%lf%lf%lf%lf",&n,&G,&size.x,&size.y,&dt,&c);
		body=(Bodies*)malloc(sizeof(Bodies)*n);
		for(j=0;j<n;j++)
		{
			/*���ļ�������m,px,py,vx,vy*/
			fscanf(file,"%lf%lf%lf%lf%lf",&body[j].m,&body[j].px,&body[j].py,&body[j].vx,&body[j].vy);
			body[j].fx=body[j].fy=0;
		}
	}
	fclose(file);
}
void calForces()
{
	double distance,magnitude;
	double temp;
	int i,j;
	Point direction;	
	for(i=start[myid];i<end[myid];i++){
		for(j=0;i<n;i++)
		{
			if(i==j) continue;
			/*�����������*/
			distance=sqrt((body[i].px-body[j].px)*(body[i].px-body[j].px)+(body[i].py-body[j].py)*(body[i].py-body[j].py));
			temp=(distance*distance)+(c*c);
			/*����������������*/
			magnitude=(G*body[i].m*body[j].m)/temp;
			/*x��y�����ϵľ���*/
			direction.x=body[j].px-body[i].px;
			direction.y=body[j].py-body[i].py;
			/*x��y�����ϵ���������*/
			body[i].fx=body[i].fx+magnitude*direction.x/distance;
			body[j].fx=body[j].fx+magnitude*direction.x/distance;
			body[i].fy=body[i].fy+magnitude*direction.y/distance;
			body[j].fy=body[j].fy+magnitude*direction.y/distance;
		}	
	}
}

void calnbody()
{
	int i;
	temPoint=(Point*)malloc(sizeof(Point)*n);
	for(i=start[myid];i<end[myid];i++)
	{
		/*x��y�����ϵ��ٶ�*/
		delv.x=(body[i].fx/body[i].m)*dt;
		delv.y=(body[i].fy/body[i].m)*dt;
		delp.x=(body[i].vx+delv.x/2)*dt;
		delp.y=(body[i].vy+delv.y/2)*dt;
		/*��ʱ����*/
		temPoint[i].x=body[i].px+delp.x;
		temPoint[i].y=body[i].py+delp.y;
		/*Խ����*/
		if(temPoint[i].x<0||temPoint[i].x>size.x||temPoint[i].y<0||temPoint[i].y>size.y)
		{
			body[i].px=body[i].px;
			body[i].py=body[i].py;
			body[i].vx=-body[i].vx;
			body[i].vy=-body[i].vy;
		}
		else
		{
			body[i].px=temPoint[i].x;
			body[i].py=temPoint[i].y;
			body[i].vx=body[i].vx+delv.x;
			body[i].vy=body[i].vy+delv.y;
		}
		body[i].fx=body[i].fy=0.0;
	}
	
}


int main(int argc,char *argv[])
{
	int i,q,j,step,k,r,p;
	double stime,etime;
	
	MPI_Init(&argc,&argv);//��ʼ��
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//�ҵ��Լ���id
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//�ҵ�������
	steps=atoi(argv[1]);	
	stime=MPI_Wtime();//��¼��ʼʱ��
	init();
	omp_set_num_threads(numprocs);
	#pragma omp parrallel
	start=(int *)malloc(numprocs*sizeof(int));
	end=(int*)malloc(numprocs*sizeof(int));
	start[0]=0;
	end[0]=n/numprocs;
	/*��1000������ƽ�������numprocs������*/
	for(i=1;i<numprocs;i++){
		start[i]=end[i-1];
		end[i]=start[i]+n/numprocs;
	}
	end[numprocs-1]=n;

	#pragma omp for
	for (step = 0; step < steps; step++) {	
		/*���м���*/
		calForces();
		calnbody();
		/*һ��ͬ����*/
	//	MPI_Barrier(MPI_COMM_WORLD);
	}		
	etime=MPI_Wtime();//��¼����ʱ��
	/*������ݵ��ı��ļ���*/
	if (myid == 0) {
		FILE *fp;
		fp=fopen("//home//200820110945//mpiandopenmp//result1000.data","w");
		if(!fp)
		{
			printf("open file error");
		}else{			
				fprintf(fp,"%d\r\n%lf\r\n%lf %lf\r\n%lf\r\n%lf\r\n",n,G,size.x,size.y,dt,c);			
				for(r=0;r<n;r++){	
					fprintf(fp,"%f %f %f %f %f\r\n",body[r].m,body[r].px,body[r].py,body[r].vx,body[r].vy);
				}		
			}
		printf("\nthis process cost %f s\n", etime - stime);
		fclose(fp);
		}
	MPI_Finalize();
	free(start);
	free(end);
	free(body);
	free(temPoint);
	return 0;
}
