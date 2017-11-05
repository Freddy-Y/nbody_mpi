#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"mpi.h"

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
typedef struct
{
	double G;
	double sizex;
	double sizey;
	double dt;
	double c;
}Para;
Para para;

int n;

MPI_Datatype MPI_PARA;
MPI_Datatype MPI_BODY;
MPI_Status status;
int steps=1000;
int *start;
int *end;
int numprocs=100000;
int myid;
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
	/*打开文件*/
	file=fopen("//home//200820110945//MPI/sample_input.in","r");
	if(!file)
	{
		printf("cannot open file\n");
	}
	else
	{
		/*读文件中的数据n,G,size.x,size.y,dt,c*/
		fscanf(file,"%d%lf%lf%lf%lf%lf",&n,&para.G,&para.sizex,&para.sizey,&para.dt,&para.c);
		body=(Bodies*)malloc(sizeof(Bodies)*n);
		for(j=0;j<n;j++)
		{
			/*读文件的数据m,px,py,vx,vy*/
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
			/*计算星体距离*/
			distance=sqrt((body[i].px-body[j].px)*(body[i].px-body[j].px)+(body[i].py-body[j].py)*(body[i].py-body[j].py));
			temp=(distance*distance)+(para.c*para.c);
			/*计算星体万有引力*/
			magnitude=(para.G*body[i].m*body[j].m)/temp;
			/*x和y分量上的距离*/
			direction.x=body[j].px-body[i].px;
			direction.y=body[j].py-body[i].py;
			/*x和y分量上的万有引力*/
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
		/*x和y分量上的速度*/
		delv.x=(body[i].fx/body[i].m)*para.dt;
		delv.y=(body[i].fy/body[i].m)*para.dt;
		delp.x=(body[i].vx+delv.x/2)*para.dt;
		delp.y=(body[i].vy+delv.y/2)*para.dt;
		/*临时变量*/
		temPoint[i].x=body[i].px+delp.x;
		temPoint[i].y=body[i].py+delp.y;
		/*越界检查*/
		if(temPoint[i].x<0||temPoint[i].x>para.sizex||temPoint[i].y<0||temPoint[i].y>para.sizey)
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
	

	MPI_Init(&argc,&argv);//初始化
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//找到自己的id
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//找到进程数
	MPI_Type_contiguous(7, MPI_DOUBLE, &MPI_BODY);//把body里的7个数构造成一个数据结构
	MPI_Type_commit(&MPI_BODY);//提交MPI_BODY这个数据类型
	MPI_Type_contiguous(5,MPI_DOUBLE,&MPI_PARA);
	MPI_Type_commit(&MPI_PARA);
	 steps=atoi(argv[1]);
	if(myid==0)
	{
		stime=MPI_Wtime();//记录起始时间
		init();
		/*广播数据*/
		MPI_Bcast(&n,1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&para,1, MPI_PARA, 0, MPI_COMM_WORLD);		
	}else{
		body=(Bodies*)malloc(sizeof(Bodies)*(n/numprocs));
		memset(body,0,sizeof(Bodies)*(n/numprocs));
	}
	start=(int *)malloc(numprocs*sizeof(int));
	end=(int*)malloc(numprocs*sizeof(int));
	start[0]=0;
	end[0]=n/numprocs;
	/*把1000个星体平均分配给numprocs个进程*/
	for(i=1;i<numprocs;i++){
		start[i]=end[i-1];
		end[i]=start[i]+n/numprocs;
	}
	end[numprocs-1]=n;
	for (step = 0; step < steps; step++) {
		/*广播body里面的数据*/
			MPI_Scatter(body,end[0],MPI_BODY,body,end[0],MPI_BODY,0,MPI_COMM_WORLD);
	
		/*进行计算*/
		calForces();
		calnbody();
       
		if(myid==0){
			for(k=1;k<numprocs;k++){
		
				/*接收计算结果*/
				MPI_Recv(&body[start[k]],end[k]-start[k],MPI_BODY,k,12,MPI_COMM_WORLD,&status);
			}			
		}
		else{
			/*发送计算结果*/
			MPI_Send(&body[start[myid]],end[myid]-start[myid],MPI_BODY,0,12,MPI_COMM_WORLD);	
		}
	}		
	etime=MPI_Wtime();//记录结束时间
	/*输出数据到文本文件里*/
	if (myid == 0) {
		FILE *fp;
		fp=fopen("//home//200820110945//MPI//result1000.data","w");
		if(!fp)
		{
			printf("open file error");
		}else{			
				fprintf(fp,"%d\r\n%lf\r\n%lf %lf\r\n%lf\r\n%lf\r\n",n,para.G,para.sizex,para.sizey,para.dt,para.c);			
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
