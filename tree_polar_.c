// This code is to calculate the 2-dimensional 2-point correlation function with
// spherical bin. 
// It only calculates DR and DD pairs.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

#define leafsize 1000

typedef struct {
  int numdata;
  int leafpoint;
  int leftchild;
  int rightchild;
  int depth;
} tree_st;

typedef struct {
  float hrect[2][3];
} hrect_st;

typedef struct {
  float leafpoints[leafsize][4];
} leaf_st;

typedef struct {
  float **data;
  hrect_st hrect;
  int depth;
  int parent;
  int left;
  int ndata;
} stack_st;

double mysecond();

long **dualTreeCountAuto(float BIN_MIN, float MIN_MAX, float BIN_WIDTH, tree_st *tree1, int len_tree1, hrect_st *hrectree1, leaf_st *leaf_points1, int node_idx1, int int_idx2);
long **dualTreeCountCross(float BIN_MIN, float MIN_MAX, float BIN_WIDTH, tree_st *tree1, tree_st *tree2,int len_tree1, int len_tree2, hrect_st *hrectree1, hrect_st *hrectree2, leaf_st *leaf_points1, leaf_st *leaf_points2, int node_idx1, int node_idx2);
int kdtree(float **data, int ndata, tree_st **tree, hrect_st **hrect, leaf_st **leaf_points);
int ReadData(char *filename, float ***data);
int ReadbinData(char *filename, float ***data);
int *makeNodeList(tree_st *tree, int N, int len_tree1, int nproc, int *nodesize);
float min(const float *arr, size_t length);
float max(const float *arr, size_t length);
int cmp_x(const void *pa, const void *pb);
int cmp_y(const void *pa, const void *pb);
int cmp_z(const void *pa, const void *pb);
long **dualTreeParallel(tree_st *tree1, tree_st *tree2, int len_tree1, int len_tree2, hrect_st *hrectree1, hrect_st *hrectree2, leaf_st *leaf_points1, leaf_st *leaf_points2, int *nodelist1, int nodesize1, int *nodelist2, int nodesize2, int cross, float BIN_MIN, float BIN_MAX, float BIN_WIDTH, int nproc);

double mysecond()
{
   struct timeval tp;
   struct timezone tzp;
   int i;
   i = gettimeofday(&tp,&tzp);
   return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

typedef char* CString; /* CString is type 'pointer to char' ... */
	 
	/* returns NULL for EOF ... or a pointer to a NEW CString ... */

float min(const float *arr, size_t length) {
    // returns the minimum value of array
    size_t i;
    float minimum = arr[0];
    for (i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
        }
    }
    return minimum;
}

float max(const float *arr, size_t length) {
    // returns the minimum value of array
    size_t i;
    float maximum = arr[0];
    for (i = 1; i < length; ++i) {
        if (maximum < arr[i]) {
            maximum = arr[i];
        }
    }
    return maximum;
}

inline float float_max(float a, float b){
  return a > b ? a : b;
}

int cmp_z ( const void *pa, const void *pb ) {
  float* a = *(float**) pa;
  float* b = *(float**) pb;
  if ( a[2] < b[2] ) return -1;
  else if ( a[2] > b[2] ) return +1;
  else return 0;
}
int cmp_x ( const void *pa, const void *pb ) {
  float* a = *(float**) pa;
  float* b = *(float**) pb;
  if ( a[0] < b[0] ) return -1;
  else if ( a[0] > b[0] ) return +1;
  else return 0;
}
int cmp_y ( const void *pa, const void *pb ) {
  float* a = *(float**)pa;
  float* b = *(float**)pb;
  if ( a[1] < b[1] ) return -1;
  else if ( a[1] > b[1] ) return +1;
  else return 0;
}

CString readLine(FILE* f)
{
    static int c = 0; /* static ... to remember if a previous call set EOF */
    static int lineCount = 0; /* to handle the empty file ... */
    int bufSize = 10000, i = 0; /* adjust 255 to whatever is 'right' for your data */
	    if(c == EOF) {c=0; lineCount=0; return NULL;} /* reset so can rewind and read again */
    CString line = (CString) calloc(bufSize, sizeof(char));
    while ((c = fgetc(f)) != EOF && c != '\n')
    {
        if(i >= bufSize)
        {
            bufSize += 10001; /* adjust 256 to whatever is 'right' for your data */
            line = (CString) realloc(line, bufSize*sizeof(char));
        }
        line[i++] = c;
    }
    /* handle special case of empty file ...*/
    if(lineCount++ == 0 && c == EOF) {free(line); return NULL;}
 
    line[i] = '\0'; /* confirm terminal 0 */
    return realloc(line, i+1); /* total len =  last index i ... + 1 more */
}

int ReadData(char *filename, float ***da)
{

  FILE *fp;
  CString cs;
  int i,j;
  fp = fopen(filename,"r");
  float *dat, **data, num;
  int ndata;
  printf("start to read %s\n",filename);
  i = 0;
  while((cs = readLine(fp)))
    {
      i += 1;
      free(cs);
    }
  ndata = i-1;
  fclose(fp);

  fp = fopen(filename,"r");
  dat = (float*)malloc(ndata * 4 * sizeof(float));
  data = (float**)malloc(ndata * sizeof(float*) );
  for(j = 0;j<=ndata-1;j++){
    data[j] = &dat[j*4];
  }

  for(i=0; i<ndata     ;i++){
    for(j = 0;j < 3; j++){
      fscanf(fp,"%f", &num);
      data[i][j] = num;
    }
  }

  for(i=0; i<ndata     ;i++){
     data[i][3] = sqrt(data[i][0]*data[i][0]+data[i][1]*data[i][1]+data[i][2]*data[i][2]);
  }
  *da = *&data;
  fclose(fp);
  return ndata;
}
int ReadbinData(char *filename, float ***da)
{
struct rec
{
        float x,y,z;
};
  FILE *fp;
  CString cs;
  int i,j;
  fp = fopen(filename,"rb");
  struct rec my_record;
  float *dat, **data, num;
  int ndata;
  printf("start to read %s\n",filename);

  fp = fopen(filename,"rb");
  fread(&ndata,sizeof(int),1,fp);
  printf("%d\n",ndata);
  dat = (float*)malloc(ndata * 4 * sizeof(float));
  data = (float**)malloc(ndata * sizeof(float*) );
  for(i=0;i<=ndata-1;i++)
  {
    data[i] = &dat[4*i];
  }
  for( i=2;i<=ndata+1;i++)
  {
        fread(&my_record,sizeof(struct rec),1,fp);
        data[i-2][0] = my_record.x;
        data[i-2][1] = my_record.y;
        data[i-2][2] = my_record.z;
	data[i-2][3] = sqrt(my_record.x*my_record.x+my_record.y*my_record.y+my_record.z*my_record.z);
  }
  *da = *&data;
  printf("%f %f %f %f\n",data[10][0],data[10][1],data[10][2],data[10][3]);
  fclose(fp);
  return ndata;
}


int kdtree(float **data, int ndata, tree_st **tr, hrect_st **hr, leaf_st **lp)
{
  int right, left, leftbranch;
  hrect_st left_hrect, right_hrect, hrect, *hrectree;
  tree_st leaf, *tree;
  float *dat_tmp, **data_tmp, *left_dat, *right_dat, **left_data, **right_data, xvec[ndata],yvec[ndata],zvec[ndata],r[ndata];
  int ndata_l, ndata_r;
  leaf_st leaf_data, *leaf_points;
  stack_st *stack;
  int node_ptr;
  int splitdim ;
  int numpoints;
  int stacksize, depth, parent, num_leaf = 0;
  int i, j, l;
  FILE *hre;
  char *hrectfile;
  int num = 0;
  int n_left_tot, n_right_tot, mx_dpt;  
  for(i=0;i<=ndata-1;i++){
    xvec[i] = data[i][0];
    yvec[i] = data[i][1];
    zvec[i] = data[i][2];
    r[i] = data[i][3];
  }
  //split data into position vectors
  hrect.hrect[0][0] = min(xvec,ndata);
  hrect.hrect[0][1] = min(yvec,ndata);
  hrect.hrect[0][2] = min(zvec,ndata);
  hrect.hrect[1][0] = max(xvec,ndata);
  hrect.hrect[1][1] = max(yvec,ndata);
  hrect.hrect[1][2] = max(zvec,ndata);
  //hrect constructed
printf("Construction of tree started\n");
  //sorting...
  qsort(data, ndata, sizeof(float*), cmp_x);

  n_left_tot = 0;
  n_right_tot = 0;
  mx_dpt = (int)(log(2*ndata/leafsize)/log(2.));
  printf("maximum depth is %d\n",mx_dpt);
  //sorting done
  ndata_l = ndata/2;
  ndata_r = ndata - ndata_l;

  n_left_tot += ndata_l;
  n_right_tot += ndata_r;

  left_data = (float**) malloc(5*mx_dpt* ndata_l *    sizeof(float*));
  left_dat = (float*) malloc(5*mx_dpt*ndata_l * 4* sizeof(float));
  for(j = 0;j<=5*mx_dpt*ndata_l-1;j++) left_data[j] = &left_dat[j*4];

  right_data = (float**) malloc(5*mx_dpt*ndata_r *    sizeof(float*));
  right_dat = (float*) malloc(5*mx_dpt*ndata_r * 4* sizeof(float));
  for(j = 0;j<=5*mx_dpt*ndata_r-1;j++) right_data[j] = &right_dat[j*4];
  for(j = 0; j <= ndata_l-1; j++){
    memmove(left_data[j],data[j],4*sizeof(float*));
  }
   
  for(j = 0; j <= ndata_r-1; j++){   
    memmove(right_data[j],data[j+ndata_l],4*sizeof(float*));
  }
  

  for(i = 0; i <= 1; i++){
    for(j = 0; j <=2; j++){
      *&left_hrect.hrect[i][j] = *&hrect.hrect[i][j];
      *&right_hrect.hrect[i][j] = *&hrect.hrect[i][j];
    }
  }
/*
	left_hrect = hrect;
	right_hrect = hrect;
*/
  *&left_hrect.hrect[1][0] = *&left_data[ndata_l-1][0];
  *&right_hrect.hrect[0][0] = *&right_data[0][0];

  node_ptr = 1;

  tree = calloc(node_ptr , sizeof(tree_st));
  tree[0].numdata = ndata;
  tree[0].leafpoint = -1;
  tree[0].leftchild = -1;
  tree[0].rightchild = -1;
  tree[0].depth = 0;

  stacksize = 2;
  stack = calloc(stacksize,sizeof(stack_st));
  stack[1].data= left_data;
  stack[1].hrect = left_hrect;
  stack[1].depth = 1;
  stack[1].parent = 0;
  stack[1].left = 1;
  stack[1].ndata = ndata_l;
  stack[0].data= right_data;
  stack[0].hrect = right_hrect;
  stack[0].depth = 1;
  stack[0].parent = 0;
  stack[0].left = 0;
  stack[0].ndata = ndata_r;

  hrectree = calloc(node_ptr,sizeof(hrect_st));
  hrectree[node_ptr-1] = hrect;

  leaf_points = (leaf_st*)malloc(1 * sizeof(leaf_st));

   data_tmp = calloc(ndata,sizeof(float*));
   dat_tmp = calloc(4*ndata,sizeof(float));
  printf("%p %p %p %p %p %p\n",data_tmp,dat_tmp,left_dat,left_data,right_dat,right_data);
    for(i=0;i<ndata;i++) data_tmp[i] = &dat_tmp[4*i];

  while (stacksize != 0){
    l = stacksize - 1;
    ndata = stack[l].ndata;
    for(i=0;i<ndata;i++) for(j=0;j<4;j++){
         data_tmp[i][j] = stack[l].data[i][j];
     }
    hrect = stack[l].hrect;
    depth = stack[l].depth;
    parent = stack[l].parent;
    leftbranch = stack[l].left;
    stacksize -= 1;

    numpoints = tree[parent].numdata;
    left = tree[parent].leftchild;
    right = tree[parent].rightchild;
    if(leftbranch == 1){
      tree[parent].leftchild = node_ptr;
    }
    else{
      tree[parent].rightchild = node_ptr;
    }

    node_ptr += 1;
    hrectree = realloc(hrectree, node_ptr * sizeof(hrect_st));

    hrectree[node_ptr-1] = hrect; 
    if(ndata <= leafsize){
      num_leaf += 1;
      for (i = 0; i < 4; i++){
	for (j = 0; j < ndata; j++){
	  memmove(&leaf_data.leafpoints[j][i],&data_tmp[j][i],sizeof(float));
	}
	if(ndata < leafsize){
	  for(j = ndata; j < leafsize; j++){
	    leaf_data.leafpoints[j][i] = -1.;
	  }
	}
      }
      
      leaf_points = realloc(leaf_points, num_leaf * sizeof(leaf_st));
      leaf_points[num_leaf - 1] = leaf_data;
      leaf.leafpoint = num_leaf-1; 
      leaf.numdata = ndata;
      leaf.leftchild = 0;
      leaf.rightchild = 0;
      leaf.depth = depth;
      tree = realloc(tree, node_ptr * sizeof(tree_st));
      tree[node_ptr-1] = leaf; 
    }
    else{
      
      splitdim = depth % 3;

      if(splitdim == 0){
	qsort(data_tmp,ndata,sizeof(float*),cmp_x);
      }
     
      else if(splitdim == 1){
	qsort(data_tmp,ndata,sizeof(float*),cmp_y);
      }

      else{
	qsort(data_tmp,ndata,sizeof(float*),cmp_z);
      }
      ndata_l = ndata/2;
      ndata_r = ndata - ndata_l;

      n_left_tot += ndata_l;
      n_right_tot += ndata_r;
for(i=0;i<4;i++){
      for(j = 0; j <= ndata_l-1; j++){
	memmove(&left_data[n_left_tot-ndata_l+j][i],&data_tmp[j][i],sizeof(float));
      }

      for(j = 0; j <= ndata_r-1; j++){
	memmove(&right_data[n_right_tot-ndata_r+j][i],&data_tmp[j+ndata_l][i],sizeof(float));
      }
}
   /*   for(j = 0; j <= ndata_l-1; j++){
	memmove(left_data[n_left_tot-ndata_l+j],data_tmp[j],3*sizeof(float*));
      }

      for(j = 0; j <= ndata_r-1; j++){
	memmove(right_data[n_right_tot-ndata_r+j],data_tmp[j+ndata_l],3*sizeof(float*));
      }
     */ for(i = 0; i <= 1; i++){
	for(j = 0; j<=2; j++){
	  *&left_hrect.hrect[i][j] = *&hrect.hrect[i][j];
	  *&right_hrect.hrect[i][j] = *&hrect.hrect[i][j];
	}
      }

      *&left_hrect.hrect[1][splitdim] = *&left_data[n_left_tot-1][splitdim];
      *&right_hrect.hrect[0][splitdim] = *&right_data[n_right_tot-ndata_r][splitdim];
      tree = realloc(tree, node_ptr*sizeof(tree_st));

      tree[node_ptr-1].numdata = ndata;
      tree[node_ptr-1].leafpoint = -1;
      tree[node_ptr-1].leftchild = -1;
      tree[node_ptr-1].rightchild = -1;
      tree[node_ptr-1].depth = depth;
      stacksize += 2;

      stack = realloc(stack, stacksize*sizeof(stack_st));
      
      stack[stacksize-1].data = left_data+(n_left_tot-ndata_l);
      stack[stacksize-1].hrect = left_hrect;
      stack[stacksize-1].depth = depth+1;
      stack[stacksize-1].parent = node_ptr-1;
      stack[stacksize-1].left = 1;
      stack[stacksize-1].ndata = ndata_l;

      stack[stacksize-2].data = right_data+(n_right_tot-ndata_r);
      stack[stacksize-2].hrect = right_hrect;
      stack[stacksize-2].depth = depth+1;
      stack[stacksize-2].parent = node_ptr-1;
      stack[stacksize-2].left = 0;
      stack[stacksize-2].ndata = ndata_r;
    }
  }
  *tr = *&tree;
  *hr = *&hrectree;  
  *lp = *&leaf_points;

  printf("size of the tree is %d \n",node_ptr);
  printf("number of leaves is %d \n",num_leaf);
  printf("%p %p %p %p %p %p\n",data_tmp,dat_tmp,left_dat,left_data,right_dat,right_data);
/*  for(i=0;i<node_ptr;i++){
     if(tree[i].depth == 12){
	     num+=tree[i].numdata;
     }
  }
  printf("ndata for depth 12 is : %d\n",num);
*/
//  hrectfile = "hrectree.txt";
//  hre = fopen(hrectfile,"w");
//  for(i=0;i<node_ptr;i++){
//  fprintf(hre,"%f %f %f %f %f %f\n",hrectree[i].hrect[0][0],hrectree[i].hrect[1][0],hrectree[i].hrect[0][1],hrectree[i].hrect[1][1],hrectree[i].hrect[0][2],hrectree[i].hrect[1][2]);
//  }
//  fclose(hre);
  free(stack);
  free(left_data);
  free(right_data);
  free(right_dat);
  free(left_dat);
  free(data_tmp);
  free(dat_tmp);
/*
*/
  return node_ptr;

  printf("construction of tree done\n");
}

long **dualTreeParallel(tree_st *tree1, tree_st *tree2,  int len_tree1, int len_tree2, hrect_st *hrectree1, hrect_st *hrectree2, leaf_st *leaf_points1, leaf_st *leaf_points2, int *nodelist1, int nodesize1, int *nodelist2, int nodesize2, int cross, float BIN_MIN, float BIN_MAX, float BIN_WIDTH, int nproc)
{
  int numbins = (int)((BIN_MAX-BIN_MIN)/(BIN_WIDTH));
  int node1, node2 , i, j, k, l,n;
  long **bins, **counts,**temps;
  long *bin,*count,*temp;
  bins = calloc(20,sizeof(long*));
  bin = calloc(numbins*20,sizeof(long));
  for(i=0;i<20;i++) bins[i] = &bin[i*numbins];
  if(cross == 1){
  printf("parallel region started, Cross\n"); 
  #pragma omp parallel private (n,temp,temps,counts,count, node1, node2, i,j,k,l) shared(bins,tree1,tree2,hrectree1,hrectree2,leaf_points1,leaf_points2,len_tree1,len_tree2)
  {
    n = omp_get_thread_num();
    temps = calloc(20,sizeof(long*));
    temp = calloc(numbins*20,sizeof(long));
    count = calloc(numbins*20,sizeof(long));
    for(k=0;k<20;k++){
      temps[k] = &temp[k*numbins];
    }
//    printf("temporary array allocation done\n");
    #pragma omp for schedule(dynamic,30)
    for(i=0;i<nodesize1;i++){
      for(j = 0; j<nodesize2; j++){
         counts = calloc(20,sizeof(long*));
//    printf("I am getting into parallel for region %d\n",n);
         for(k=0;k<20;k++){
           counts[k] = &count[k*numbins];
         }
   	 node1 = nodelist1[i];
	 node2 = nodelist2[j];
	 counts =  dualTreeCountCross(BIN_MIN, BIN_MAX, BIN_WIDTH, tree1, tree2, len_tree1, len_tree2, hrectree1, hrectree2, leaf_points1, leaf_points2, node1, node2);
	 for(l=0;l<20;l++){
	   for(k = 0; k<numbins; k++){
	     temps[l][k] += counts[l][k];
 	   }
	 }
         free(counts);
       }
    }
    free(count);
    #pragma omp critical
    {
    printf("I am in the critical region %d\n",n);
      for(l=0;l<20;l++){
        for(k = 0; k<numbins; k++){
          bins[l][k] += temps[l][k];
        }
      }
    printf("I am out of critical region %d\n",n);
    }
    free(temps);
    free(temp);
//    printf("freeing temporary array allocation done\n");
    }
    printf("pair count complete\n");
  }
  else{
    printf("parallel region started, Auto\n"); 
    #pragma omp parallel private (n,temp,temps,counts,count, node1, node2, i,j,k,l) shared(bins,tree1,hrectree1,leaf_points1,len_tree1)
    {
      temps = calloc(20,sizeof(long*));
      temp = calloc(numbins*20,sizeof(long));
      count = calloc(numbins*20,sizeof(long));
      for(k=0;k<20;k++){
        temps[k] = &temp[k*numbins];
      }
//      printf("temporary array allocation done\n");
      n = omp_get_thread_num();
      #pragma omp for schedule(dynamic,30) 
      for(i=0;i<nodesize1;i++){
        for(j = i; j<nodesize1; j++){
          counts = calloc(20,sizeof(long*));
          for(k=0;k<20;k++){
            counts[k] = &count[k*numbins];
          }
	node1 = nodelist1[i];
	node2 = nodelist1[j];
	counts =  dualTreeCountAuto(BIN_MIN, BIN_MAX, BIN_WIDTH, tree1, len_tree1, hrectree1, leaf_points1, node2, node1);
        for(l=0;l<20;l++){
	  for(k = 0; k<numbins; k++){
	    temps[l][k] += counts[l][k];
 	  }
	}
      free(counts);
      }
    }
    free(count);
    #pragma omp critical
    {
    printf("I am in the critical region %d\n",n);
      for(l=0;l<20;l++){
        for(k = 0; k<numbins; k++){
          bins[l][k] += temps[l][k];
        }
      }
    printf("I am out of critical region %d\n",n);
    }
    free(temps);
    free(temp);
//    printf("freeing temporary array allocation done\n");
    }
  printf("pair count complete\n");
  }
  printf("returning bins\n");
  return bins;
  printf("bins returned\n");
}
  
long **dualTreeCountAuto(float BIN_MIN, float BIN_MAX, float BIN_WIDTH, tree_st *tree1,  int len_tree1,  hrect_st *hrectree1, \
		       leaf_st *leaf_points1, int node_idx1, int node_idx2){
  int numbins = (int)((BIN_MAX-BIN_MIN)/BIN_WIDTH);
  int stacksize = 1;
  int idx1, idx2, numpoints1, numpoints2, left1, right1, left2, right2, i, j, r_index;
  int leaf_ptr1, leaf_ptr2,mu_idx;
  float dmin, dmax, dist,z,mu;
  hrect_st hrect1, hrect2;
  long **bins, *bin;
  int stack[32]; 
  //printf("length of the tree %d \n",len_tree1);
  bins = calloc(20, sizeof(long*));
  bin = calloc(numbins * 20, sizeof(long));
  for(i = 0;i<20;i++){
	bins[i] = &bin[i*numbins];
  }
  stack[0] = node_idx2;
  stack[1] = node_idx1;

  while(stacksize >= 0){
    idx1 = stack[stacksize];
    stacksize -= 1;
    idx2 = stack[stacksize];
    stacksize -= 1;
    if (idx2 <= idx1){
      numpoints1 = tree1[idx1].numdata;
      leaf_ptr1 = tree1[idx1].leafpoint;
      left1 = tree1[idx1].leftchild;
      right1 = tree1[idx1].rightchild;
      numpoints2 = tree1[idx2].numdata;
      leaf_ptr2 = tree1[idx2].leafpoint;
      left2 = tree1[idx2].leftchild;
      right2 = tree1[idx2].rightchild;
      hrect1 = hrectree1[idx1];
      hrect2 = hrectree1[idx2];
      dmin = powf((powf(float_max(0.,float_max(hrect1.hrect[0][0]-hrect2.hrect[1][0],hrect2.hrect[0][0]-hrect1.hrect[1][0])),2.) + \
		  powf(float_max(0.,float_max(hrect1.hrect[0][1]-hrect2.hrect[1][1],hrect2.hrect[0][1]-hrect1.hrect[1][1])),2.) + \
		   powf(float_max(0.,float_max(hrect1.hrect[0][2]-hrect2.hrect[1][2],hrect2.hrect[0][2]-hrect1.hrect[1][2])),2.)),0.5);
      if (dmin <= BIN_MAX){
	dmax = powf((powf(float_max(hrect1.hrect[1][0]-hrect2.hrect[0][0],hrect2.hrect[1][0]-hrect1.hrect[0][0]),2.) + \
		     powf(float_max(hrect1.hrect[1][1]-hrect2.hrect[0][1],hrect2.hrect[1][1]-hrect1.hrect[0][1]),2.) + \
		     powf(float_max(hrect1.hrect[1][2]-hrect2.hrect[0][2],hrect2.hrect[1][2]-hrect1.hrect[0][2]),2.)),0.5);
	if(dmax >= BIN_MIN){
	  if(idx1 == idx2){
//	    if((int)((dmin-BIN_MIN)/BIN_WIDTH) == (int)((dmax-BIN_MIN)/BIN_WIDTH)){
//	      bins[(int)((dmin-BIN_MIN)/BIN_WIDTH)] += (numpoints1 * (numpoints1-1))/2;
//              printf("multiplied in a leaf\n");
//	    }
//	    else if(leaf_ptr1 != -1){
	    if(leaf_ptr1 != -1){
	      //printf("leaf!\n");
	      for(i = 0; i <= numpoints1-2; i++){
		for(j=i+1; j <= numpoints1-1; j++){
		  dist = powf((powf((leaf_points1[leaf_ptr1].leafpoints[i][0]-leaf_points1[leaf_ptr1].leafpoints[j][0]),2.) + \
			       powf((leaf_points1[leaf_ptr1].leafpoints[i][1]-leaf_points1[leaf_ptr1].leafpoints[j][1]),2.) + \
			       powf((leaf_points1[leaf_ptr1].leafpoints[i][2]-leaf_points1[leaf_ptr1].leafpoints[j][2]),2.)),0.5);
		  if(dist > BIN_MIN){
		    r_index = (int)((dist-BIN_MIN)/BIN_WIDTH);
		    if(r_index < numbins){
		      mu = fabsf(leaf_points1[leaf_ptr1].leafpoints[i][3]-leaf_points1[leaf_ptr1].leafpoints[j][3])*sqrt(0.5+0.5*(leaf_points1[leaf_ptr1].leafpoints[i][0]*leaf_points1[leaf_ptr1].leafpoints[j][0]+leaf_points1[leaf_ptr1].leafpoints[i][1]*leaf_points1[leaf_ptr1].leafpoints[j][1]+leaf_points1[leaf_ptr1].leafpoints[i][2]*leaf_points1[leaf_ptr1].leafpoints[j][2])/(leaf_points1[leaf_ptr1].leafpoints[i][3]*leaf_points1[leaf_ptr1].leafpoints[j][3]))/dist;
		      mu_idx = (int)(20.*mu);
//			mu_idx+=10;
		      if(mu_idx == 20) mu_idx -= 1;
		      //printf("mu_idx : %d\n",mu_idx);
		      bins[mu_idx][r_index] += 1;
		    }
		  }
		}
	      }
	    }
	    else{
		stacksize += 1;
		stack[stacksize] = idx2;
		stacksize += 1;
		stack[stacksize] = right1;
		stacksize += 1;
		stack[stacksize] = idx2;
		stacksize += 1;
		stack[stacksize] = left1;
	    }
	  }
	  else{
//idx2 < idx1
	    //if((int)((dmin-BIN_MIN)/BIN_WIDTH) == (int)((dmax-BIN_MIN)/BIN_WIDTH)){
	      //bins[(int)((dmin-BIN_MIN)/BIN_WIDTH)] += numpoints1 * numpoints2;
               //printf("multiplied\n");
	    //}
	    //else if((leaf_ptr1 != -1) &&( leaf_ptr2 != -1)){
	    if((leaf_ptr1 != -1) &&( leaf_ptr2 != -1)){
	      //printf("leaf!\n");
	      for(i = 0; i <= numpoints1-1; i++){
		for(j=0; j <= numpoints2-1; j++){
		  dist = powf((powf((leaf_points1[leaf_ptr1].leafpoints[i][0]-leaf_points1[leaf_ptr2].leafpoints[j][0]),2.) + \
			       powf((leaf_points1[leaf_ptr1].leafpoints[i][1]-leaf_points1[leaf_ptr2].leafpoints[j][1]),2.) + \
			       powf((leaf_points1[leaf_ptr1].leafpoints[i][2]-leaf_points1[leaf_ptr2].leafpoints[j][2]),2.)),0.5);
		 
		  if(dist > BIN_MIN){
		    r_index = (int)((dist-BIN_MIN)/BIN_WIDTH);
		    if(r_index < numbins){
		      mu = fabsf(leaf_points1[leaf_ptr1].leafpoints[i][3]-leaf_points1[leaf_ptr2].leafpoints[j][3])*sqrt(0.5+0.5*(leaf_points1[leaf_ptr1].leafpoints[i][0]*leaf_points1[leaf_ptr2].leafpoints[j][0]+leaf_points1[leaf_ptr1].leafpoints[i][1]*leaf_points1[leaf_ptr2].leafpoints[j][1]+leaf_points1[leaf_ptr1].leafpoints[i][2]*leaf_points1[leaf_ptr2].leafpoints[j][2])/(leaf_points1[leaf_ptr1].leafpoints[i][3]*leaf_points1[leaf_ptr2].leafpoints[j][3]))/dist;
		      mu_idx = (int) (20.*mu);
//			mu_idx+=10;
		      if(mu_idx == 20) mu_idx -= 1;
		      //printf("mu_idx : %d\n",mu_idx);
		      bins[mu_idx][r_index] += 1;
		    }
		  }
		}
	      }
	    }
	    else{
	      if (numpoints2 > numpoints1){
		stacksize += 1;
		stack[stacksize] = right2;
		stacksize += 1;
		stack[stacksize] = idx1;
		stacksize += 1;
		stack[stacksize] = left2;
		stacksize += 1;
		stack[stacksize] = idx1;
	      }
	      else{	      
		stacksize += 1;
		stack[stacksize] = idx2;
		stacksize += 1;
		stack[stacksize] = right1;
		stacksize += 1;
		stack[stacksize] = idx2;
		stacksize += 1;
		stack[stacksize] = left1;
	      } 
	    }
	  }
	}     
      }
    }
  }
  //free(stack);
  return bins;
  //free(bins);
}
long **dualTreeCountCross(float BIN_MIN, float BIN_MAX, float BIN_WIDTH, tree_st *tree1, tree_st *tree2, int len_tree1, int len_tree2, hrect_st *hrectree1, \
			hrect_st *hrectree2, leaf_st *leaf_points1, leaf_st *leaf_points2, int node_idx1, int node_idx2){
  int numbins = (int)((BIN_MAX-BIN_MIN)/BIN_WIDTH);
  int stacksize = 1;
  int idx1, idx2, numpoints1, numpoints2, left1, right1, left2, right2, i, j, r_index;
  int leaf_ptr1, leaf_ptr2,mu_idx;
  float dmin, dmax, dist, z, mu;
  hrect_st hrect1, hrect2;
  long **bins,*bin;
  //int stack[(len_tree1+len_tree2)];
  int stack[32];
  bins = calloc(20, sizeof(long*));
  bin = calloc(numbins * 20, sizeof(long));
  for(i = 0;i<20;i++){
	bins[i] = &bin[i*numbins];
  }
  
  stack[0] = node_idx2;
  stack[1] = node_idx1;
 
  while(stacksize >= 0){
    //printf("stacksize : %d\n",stacksize);
    idx1 = stack[stacksize];
    stacksize -= 1;
    idx2 = stack[stacksize];
    stacksize -= 1;
    numpoints1 = tree1[idx1].numdata;
    leaf_ptr1 = tree1[idx1].leafpoint;
    left1 = tree1[idx1].leftchild;
    right1 = tree1[idx1].rightchild;
    numpoints2 = tree2[idx2].numdata;
    leaf_ptr2 = tree2[idx2].leafpoint;
    left2 = tree2[idx2].leftchild;
    right2 = tree2[idx2].rightchild;
    hrect1 = hrectree1[idx1];
    hrect2 = hrectree2[idx2];

    dmin = powf((powf(float_max(0.,float_max(hrect1.hrect[0][0]-hrect2.hrect[1][0],hrect2.hrect[0][0]-hrect1.hrect[1][0])),2.) + \
		 powf(float_max(0.,float_max(hrect1.hrect[0][1]-hrect2.hrect[1][1],hrect2.hrect[0][1]-hrect1.hrect[1][1])),2.) + \
		 powf(float_max(0.,float_max(hrect1.hrect[0][2]-hrect2.hrect[1][2],hrect2.hrect[0][2]-hrect1.hrect[1][2])),2.)),0.5);
    if (dmin <= BIN_MAX){
      dmax = powf((powf(float_max(hrect1.hrect[1][0]-hrect2.hrect[0][0],hrect2.hrect[1][0]-hrect1.hrect[0][0]),2.) + \
		   powf(float_max(hrect1.hrect[1][1]-hrect2.hrect[0][1],hrect2.hrect[1][1]-hrect1.hrect[0][1]),2.) + \
		   powf(float_max(hrect1.hrect[1][2]-hrect2.hrect[0][2],hrect2.hrect[1][2]-hrect1.hrect[0][2]),2.)),0.5);
        if(dmax >= BIN_MIN){
//	if((int)((dmin-BIN_MIN)/BIN_WIDTH) == (int)((dmax-BIN_MIN)/BIN_WIDTH)){
//	  bins[(int)((dmin-BIN_MIN)/BIN_WIDTH)] += numpoints1 * numpoints2;
//               printf("multiplied\n");
//	}
	if((leaf_ptr1 != -1) && (leaf_ptr2 != -1)){
	      //printf("leaf!\n");
	  for(i = 0; i < numpoints1; i++){
	    for(j=0; j < numpoints2; j++){
//      printf("test %d %f %f %f %f\n",i,leaf_points1[leaf_ptr1].leafpoints[i][0],leaf_points1[leaf_ptr1].leafpoints[i][1],leaf_points1[leaf_ptr1].leafpoints[i][2],leaf_points1[leaf_ptr1].leafpoints[i][3]);
	      dist = powf((powf((leaf_points1[leaf_ptr1].leafpoints[i][0]-leaf_points2[leaf_ptr2].leafpoints[j][0]),2.) + \
			   powf((leaf_points1[leaf_ptr1].leafpoints[i][1]-leaf_points2[leaf_ptr2].leafpoints[j][1]),2.) + \
			   powf((leaf_points1[leaf_ptr1].leafpoints[i][2]-leaf_points2[leaf_ptr2].leafpoints[j][2]),2.)),0.5);
	      if(dist >= BIN_MIN){
		r_index = (int)((dist-BIN_MIN)/BIN_WIDTH);
		
		if(r_index < numbins){
		      mu = fabsf(leaf_points1[leaf_ptr1].leafpoints[i][3]-leaf_points2[leaf_ptr2].leafpoints[j][3])*sqrt(0.5+0.5*(leaf_points1[leaf_ptr1].leafpoints[i][0]*leaf_points2[leaf_ptr2].leafpoints[j][0]+leaf_points1[leaf_ptr1].leafpoints[i][1]*leaf_points2[leaf_ptr2].leafpoints[j][1]+leaf_points1[leaf_ptr1].leafpoints[i][2]*leaf_points2[leaf_ptr2].leafpoints[j][2])/(leaf_points1[leaf_ptr1].leafpoints[i][3]*leaf_points2[leaf_ptr2].leafpoints[j][3]))/dist;
			mu_idx = (int) (20. * mu);
//			mu_idx+=10;
		        if(mu_idx == 20) mu_idx -= 1;
//			if(mu_idx==0)printf("%d test %d %f %f %f %f\n",r_index,i,leaf_points1[leaf_ptr1].leafpoints[i][3],leaf_points2[leaf_ptr2].leafpoints[j][3],dist,leaf_points1[leaf_ptr1].leafpoints[i][3]-leaf_points2[leaf_ptr2].leafpoints[j][3]);
		  bins[mu_idx][r_index] += 1;
		}
	      }
	    }
	  }
	}
	else{
	  if (numpoints2 > numpoints1){
	    stacksize += 1;
	    stack[stacksize] = right2;
	    stacksize += 1;
	    stack[stacksize] = idx1;
	    stacksize += 1;
	    stack[stacksize] = left2;
	    stacksize += 1;
	    stack[stacksize] = idx1;
	  }
	  else{	      
	    stacksize += 1;
	    stack[stacksize] = idx2;
	    stacksize += 1;
	    stack[stacksize] = right1;
	    stacksize += 1;
	    stack[stacksize] = idx2;
	    stacksize += 1;
	    stack[stacksize] = left1;
	  } 
	}
      }
    }
  }
  //free(stack);
  return bins;
  //free(bins);
}     

int *makeNodeList(tree_st *tree, int N, int len_tree, int nprocs, int *nodesize){
  int *nodelist;
  double min_level, max_level;
  int level;
  int i, k, size ;
  printf("creating nodelist\n");
  min_level = log(nprocs)/log(2);
  max_level = log(2*N/leafsize)/log(2);
  if (min_level < max_level){
    level = (int)(0.25*min_level + 0.75*max_level);
  }
  else{
    level = max_level; 
  }
  size = pow(2,level);
  printf("nodelevel : %d \n",level);
  printf("nodesize : %d \n",size);
  nodelist = (int *)malloc(size*sizeof(int));
  k = 0;
  for(i =0; i < len_tree; i++){
    if(tree[i].depth == level-1){
      nodelist[k] = tree[i].leftchild;
      nodelist[k+1] = tree[i].rightchild;
      k += 2;
    }
  }
 // printf("Crosschecking the size of node %d \n",k);
  *nodesize = *&size;
  return nodelist;
  //free(nodelist);
}

 
int main(){
  int i,j;
  float BIN_MIN = 0.;
  float BIN_MAX = 150.;
  
  float BIN_WIDTH = 2.;
  int ndata1, ndata2;
  char filename1[400];
  //char *filename1;
  FILE *out1;
  printf("data filename : \n");
  scanf("%s",filename1);
  char tmp_filename[500],filename2[400], select[10], num_R[10];
  printf("selection : \n");
  scanf("%s",select);
  printf("Random data level : \n");
  scanf("%s",num_R);
  sprintf(filename2,"/ptmp/mpa/ctchiang/HETDEX_pkfinal/galdata/1000_realizations/74ifus_RAsep=0.2435/sel_%s/gpicc_my_tiling_74ifus_RAsep=0.2435_randoms_sel_%s_Ngal=%s.bin",select,select,num_R);
  float **data1, **data2;
  int nproc; 
  char *num = getenv("OMP_NUM_THREADS");
  nproc = atoi(num);
  int cross = 1;
  int bin_ = 0;
  long **bins1,*bin;
  int numbins = (int)((BIN_MAX-BIN_MIN)/BIN_WIDTH);
  bins1 = calloc(20,sizeof(long*));
  bin = calloc(20*numbins,sizeof(long));
  for(i=0;i<20;i++) bins1[i] = &bin[i*numbins];
  tree_st *tree1, *tree2;
  hrect_st *hrectree1, *hrectree2;
  leaf_st *leaf_points1, *leaf_points2;
  int  *nodelist1, *nodelist2, nodesize1, nodesize2;
  int len_tree1, len_tree2;
  printf("program started\n");
  printf("type of file : binary = 1, no =0\n");
  scanf("%d",&bin_);
  char *fileout,DD[500],DR[500],RR[500];
  if(bin_==0){
  ndata1 = ReadData(filename1, &data1);
  }
  else{
  ndata1 = ReadbinData(filename1, &data1);
  }
  len_tree1 = kdtree(data1, ndata1, &tree1, &hrectree1, &leaf_points1);
  free(data1); 
  if(bin_==0){
  ndata2 = ReadData(filename2, &data2);
  }
  else{
  ndata2 = ReadbinData(filename2, &data2);
  }
  len_tree2 = kdtree(data2, ndata2, &tree2, &hrectree2, &leaf_points2);
  free(data2);
  cross = 1;
  nodelist1 = makeNodeList(tree1, ndata1, len_tree1, nproc, &nodesize1);
  nodelist2 = makeNodeList(tree2, ndata2, len_tree2, nproc, &nodesize2);
  bins1 = dualTreeParallel(tree1, tree2, len_tree1, len_tree2, hrectree1, hrectree2, leaf_points1, leaf_points2, nodelist1, nodesize1, nodelist2, nodesize2, cross, BIN_MIN, BIN_MAX, BIN_WIDTH, nproc);
  printf("dualTreeParallel done\n");
  fileout = "DRpol";
  unsigned int length=strlen(filename1);
  strcpy(tmp_filename,filename1+length-8);
  sprintf(DR,"%s%s_%s_%s.txt",fileout,tmp_filename,select,num_R);
  out1 = fopen(DR,"w");
  fprintf(out1,"%f\n",BIN_MIN);
  fprintf(out1,"%f\n",BIN_MAX);
  fprintf(out1,"%f\n",BIN_WIDTH);
  fprintf(out1,"%d\n",-99);
  for(j=0;j<20;j++){
  for(i=0;i<numbins;i++){
    fprintf(out1,"\t%ld",bins1[j][i]);
  }
    fprintf(out1,"\n");
  }
  fprintf(out1,"\n");
  fclose(out1);
  
  free(tree2);
  free(hrectree2);
  free(leaf_points2);

  cross = 0;
  bins1 = dualTreeParallel(tree1, tree1, len_tree1, len_tree1, hrectree1, hrectree1, leaf_points1, leaf_points1, nodelist1, nodesize1, nodelist1, nodesize1, cross, BIN_MIN, BIN_MAX, BIN_WIDTH, nproc);
  printf("dualTreeParallel done");
  fileout = "DDpol";
  sprintf(DD,"%s%s_%s.txt",fileout,tmp_filename,select);
  out1 = fopen(DD,"w");
  fprintf(out1,"%f\n",BIN_MIN);
  fprintf(out1,"%f\n",BIN_MAX);
  fprintf(out1,"%f\n",BIN_WIDTH);
  fprintf(out1,"%d\n",ndata1);
  for(j=0;j<20;j++){
  for(i=0;i<numbins;i++){
    fprintf(out1,"\t%ld",bins1[j][i]);
  }
    fprintf(out1,"\n");
  }
  fprintf(out1,"\n");
  fclose(out1);

  free(tree1);
  free(hrectree1);
  free(leaf_points1);
/*
  cross = 0;
  bins1 = dualTreeParallel(tree2, tree2, len_tree2, len_tree2, hrectree2, hrectree2, leaf_points2, leaf_points2, nodelist2, nodesize2, nodelist2, nodesize2, cross, BIN_MIN, BIN_MAX, BIN_WIDTH, nproc);
  fileout = "RR";
  sprintf(RR,"polar%s%d_%s_%s.txt",fileout,nproc,select,num_R);
  out1 = fopen(RR,"w");

  fprintf(out1,"%f\n",BIN_MIN);
  fprintf(out1,"%f\n",BIN_MAX);
  fprintf(out1,"%f\n",BIN_WIDTH);
  fprintf(out1,"%d\n",ndata2);

  for(j=0;j<20;j++){
  for(i=0;i<numbins;i++){
    fprintf(out1,"\t%ld",bins1[j][i]);
  }
    fprintf(out1,"\n");
  }
  fprintf(out1,"\n");
  fclose(out1);
  
  free(tree2);
  free(hrectree2);
  free(leaf_points2);
*/  return 0;
}
