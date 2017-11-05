#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
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
int steps=1000;
int *start;
int *end;
int numprocs=1;
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
	/*打开文件*/
	file=fopen("//home//200820110945//OpenMP/sample_input.in","r");
	if(!file)
	{
		printf("cannot open file\n");
	}
	else
	{
		/*读文件中的数据n,G,size.x,size.y,dt,c*/
		fscanf(file,"%d%lf%lf%lf%lf%lf",&n,&G,&size.x,&size.y,&dt,&c);
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
	#pragma omp parallel	
	for(i=start[myid];i<end[myid];i++){
		#pragma omp parallel
		for(j=0;i<n;i++)
		{
			if(i==j) continue;
			/*计算星体距离*/
			distance=sqrt((body[i].px-body[j].px)*(body[i].px-body[j].px)+(body[i].py-body[j].py)*(body[i].py-body[j].py));
			temp=(distance*distance)+(c*c);
			/*计算星体万有引力*/
			magnitude=(G*body[i].m*body[j].m)/temp;
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
	#pragma omp parallel
	for(i=start[myid];i<end[myid];i++)
	{
		/*x和y分量上的速度*/
		delv.x=(body[i].fx/body[i].m)*dt;
		delv.y=(body[i].fy/body[i].m)*dt;
		delp.x=(body[i].vx+delv.x/2)*dt;
		delp.y=(body[i].vy+delv.y/2)*dt;
		/*临时变量*/
		temPoint[i].x=body[i].px+delp.x;
		temPoint[i].y=body[i].py+delp.y;
		/*越界检查*/
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
	if(argc==1)
	{
		printf("please input steps and np\n");
	}	
	else if(argc==3){
		steps=atoi(argv[1]);
		numprocs=atoi(argv[2]);
	}
	else if(argc==2){
		numprocs=atoi(argv[2]);	
	}
	stime=omp_get_wtime();//记录起始时间
	init();
	omp_set_num_threads(numprocs);
 
	 #pragma omp parallel shared(numprocs) private(myid)
	{
	myid=omp_get_thread_num();
	
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
	
	#pragma omp for
	for (step = 0; step < steps; step++) {	
		/*进行计算*/
		calForces();
		calnbody(); 
	}	
	#pragma omp barrier
	
	etime=omp_get_wtime();//记录结束时间
	/*输出数据到文本文件里*/
	if (myid == 0) {
		FILE *fp;
		fp=fopen("//home//200820110945//OpenMP//result2.data","w");
		if(!fp)
		{
			printf("open file error");
		}else{			
				fprintf(fp,"%d\r\n%lf\r\n%lf %lf\r\n%lf\r\n%lf\r\n",n,G,size.x,size.y,dt,c);			
				for(r=0;r<n;r++){	
					fprintf(fp,"%f %f %f %f %f\r\n",body[r].m,body[r].px,body[r].py,body[r].vx,body[r].vy);
				}		
			}
		printf("\nProcessing took %f seconds\n", etime - stime);
		fclose(fp);
		printf("\nsucess\n");
		}
	}
	free(start);
	free(end);
	free(body);
	free(temPoint);
	return 0;
}
