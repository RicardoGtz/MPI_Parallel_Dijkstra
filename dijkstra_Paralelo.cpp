#include <iostream>
#include<stdio.h>
#include <cfloat>
#include<stdlib.h>
#include <time.h>
#include "mpi.h"
using namespace std;

#define SEED 5
#define MASTER 0

//Creamos un puntero de tipo FILE
FILE *fp;
double **graph=NULL;
double *minDist=NULL;
int *parents=NULL;
int *visited=NULL;
int num=0;
int	taskid,    // ID del procesador */
	numtasks,    // Numero de procesadores a usar
  start,       // Index de incio de su parte del arreglo
  size,        // Tamaño de elementos del procesador
	initNode=0,  //Nodo incial
	actNode,
	totalVisited=0;

void loadFile(char c[]);
void printGraph();
void splitWork();
void updateMinDist();
void collectFromWorkers();
void reportToMaster();

int main(int argc, char *argv[]) {
	/* Inicializa MPI */
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	char c[]="graphs/tamaulipas.txt";
	loadFile(c);
	MPI_Barrier(MPI_COMM_WORLD);
	splitWork();

	if(taskid==MASTER){
 		//printGraph();
		actNode=initNode;
	}
	//Broadcast actual Node
	MPI_Bcast(&actNode,				//Inforacion a difundir
						1,							//Numero de elementos
						MPI_INT,				//Tipo de dato
						0,							//Origen
						MPI_COMM_WORLD);//Comunicador
	//Parametros iniciales
	parents[actNode]=-1;
	minDist[actNode]=0;

	while(totalVisited<num){
		updateMinDist();

		/*if(taskid==2){
			for(int i=start;i<start+size;i++)
				printf("[%d]%2.2f ",taskid,minDist[i]==DBL_MAX ? -1:minDist[i]);
			printf("\n");

			printf("Mindist de ");
			for(int i=0;i<num;i++)
				printf("[%d]%2.2f ",taskid,minDist[i]==DBL_MAX ? -1:minDist[i]);
			printf("\n");
		}*/
		if(taskid==MASTER)
			collectFromWorkers();
		else
			reportToMaster();

		if(taskid==MASTER){
			visited[actNode]=1;
			totalVisited++;
	    double min=DBL_MAX;
	    int index=0;
	    for(int i=0;i<num;i++){
	      if(visited[i]!=1 && minDist[i]<min){
	        min=minDist[i];
	        index=i;
	      }
	    }
			actNode=index;
		}

		//Broadcast actual Node
		MPI_Bcast(&actNode,1,MPI_INT,0,MPI_COMM_WORLD);
		//Broadcast arreglo de distancias minimas
		MPI_Bcast(minDist,num,MPI_DOUBLE,0,MPI_COMM_WORLD);
		//Broadcast arreglo de distancias minimas
		MPI_Bcast(parents,num,MPI_INT,0,MPI_COMM_WORLD);
		//Broadcast totalVisited
		MPI_Bcast(&totalVisited,1,MPI_INT,0,MPI_COMM_WORLD);

	}
	if(taskid==MASTER){
	//Añade los resultados del nodo maestro al resultados
	for(int i=0;i<num;i++)
		printf("[%d]%2.2f ",taskid,minDist[i]==DBL_MAX ? -1:minDist[i]);

	printf("\n");
	}

	//printf ("tarea=%2d  Nodo actual %d minDist %f parents %d \n", taskid,actNode,minDist[actNode], parents[actNode]);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	//free(graph);
	//free(minDist);
	//free(parents);
	//free(visited);
  return 0;
}
void loadFile(char c[]){
	//Abrimos el archivo a leer
 	if((fp = fopen (c, "r" ))==NULL){
 		printf("No se pudo leer el archivo\n");
 	}
	//leemos los datos
 	fscanf(fp, "%d" ,&num);
    //Incializa los arreglos
 	graph = (double**)malloc(num*sizeof(double*));
  minDist=(double*)malloc(num*sizeof(double));
  parents=(int *)malloc(num*sizeof(int));
	if(taskid==MASTER)
  	visited=(int *)malloc(num*sizeof(int));
	// Crea las columnas
	for(int i=0;i<num;i++){
		graph[i] = (double*)malloc(num*sizeof(double));
		for(int j=0;j<num;j++){
      // Inicializa la matriz en infinito (Max double)
			graph[i][j]=DBL_MAX;
		}
    // Inicializa el vector en infinito (Max double)
    minDist[i]=DBL_MAX;
		parents[i]=i;
		if(taskid==MASTER)
    	visited[i]=0;
	}
	int a,b;
  double d;
 	//Leemos las aristas de cada vertice
 	while(feof(fp)==0){
 		fscanf(fp,"%d\t%d\t%lf",&a,&b,&d);
 		graph[a-1][b-1]=d;
    graph[b-1][a-1]=d;
 	}
}
void printGraph(){
	//Imprimimos el grafo
 	for(int i = 0; i < num; i++){
 		for (int j = 0; j < num; ++j){
      if(graph[i][j]!=DBL_MAX)
        printf("%.2f\t",graph[i][j]);
      else
        printf("INF\t");
    }
    printf("\n");
 	}
}
void splitWork(){
	int nmin, nleft, nnum;
  //Define el tamaño de particion
  nmin=num/numtasks;
  nleft=num%numtasks;
  int k=0;
  for (int i = 0; i < numtasks; i++) {
    nnum = (i < nleft) ? nmin + 1 : nmin;
    if(i==taskid){
       start=k;
       size=nnum;
       printf ("tarea=%2d  Inicio=%2d  tamaño=%2d \n", taskid,start, size);

    	}
  	k+=nnum;
  }
}
void updateMinDist(){
	for(int i=start;i<start+size;i++){
		if(graph[actNode][i]<DBL_MAX){
			if((graph[actNode][i]+minDist[actNode])<minDist[i]){
				minDist[i]=graph[actNode][i]+minDist[actNode];
				parents[i]=actNode;
			}
		}
	}
}
void collectFromWorkers(){
	MPI_Status status;

  int buffer[2];
  for(int i=1;i<numtasks;i++){
    //Recibe parametros de los trabajadores
    MPI_Recv(buffer,2,MPI_INT,i,0,MPI_COMM_WORLD,&status);
    int iStart=buffer[0];
    int iSize=buffer[1];
    //Recibe la parte del arreglo de distancias minimas de cada trabajador
    MPI_Recv(&minDist[iStart],iSize,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
		//Recibe la parte del arreglo de padres de cada trabajador
    MPI_Recv(&parents[iStart],iSize,MPI_INT,i,2,MPI_COMM_WORLD,&status);
  }
}
void reportToMaster(){
  int buffer[2];
  buffer[0]=start;
  buffer[1]=size;
  //Envia los parametros al proceso Maestro
  MPI_Send(buffer,2,MPI_INT,0,0,MPI_COMM_WORLD);
  //Envia su parte del arrgelo de distancias minimas al proceso Maestro
  MPI_Send(&minDist[start],size,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
	//Envia su parte del arrgelo de padres al proceso Maestro
  MPI_Send(&parents[start],size,MPI_INT,0,2,MPI_COMM_WORLD);
}
