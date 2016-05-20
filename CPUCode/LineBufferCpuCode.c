/**
 * Document: MaxCompiler Tutorial (maxcompiler-tutorial.pdf)
 * Chapter: 8      Example: 3      Name: Two-dimensional average variable
 * MaxFile name: LineBuffer
 * Summary:
 *    Averages within an 8-point window.
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "Top.h"
#include <MaxSLiCInterface.h>
//#include "types.h"

void generateInputData(float *dataIn, int H, int W)
{
  for (int i = 0; i < H*W; i++) {
    dataIn[i] = 1;
  }
  //dataIn[nx * (ny / 2) + nx / 2] = 9;
}
inline bool neq(float a, float b) {
  if(a>b) 
    return (a-b) > 1e-5;
  return (b-a) > 1e-5;
}

int check(int size, float *dataOut, float *expectedOut)
{
  int status = 0;
  for (int i = 0; i < size ; i++) {
    if (neq(dataOut[i],expectedOut[i])) {
      //fprintf(stderr, "Error in output data at location %d, expected %g, received %g\n",i, expectedOut[i], dataOut[i]);
      status = 1;
    }
  }
  return status;
}

void print_2d_data_int(int *data, char* name, int H, int W)
{
  printf("\n%s\n", name);
  for (int y = 0; y < H; y++) {
    for (int x = 0; x < W; x++) {
      printf("%d ", data[y * W + x]);
    }
    printf("\n");
  }
}
void print_2d_data(float *data, char* name, int H, int W)
{
  printf("\n%s\n", name);
  for (int y = 0; y < H; y++) {
    for (int x = 0; x < W; x++) {
      printf("%5.2g", data[y * W + x]);
      //printf("%d ", (int)data[y * W + x]);
    }
    printf("\n");
  }
}

void pad(int inH, int inW, int padH, int padW, float *dataIn, float *dataOut){
  for (int h = 0; h < inH+2*padH; h++) {
    for (int w = 0; w < inW+2*padW; w++) {
      bool read = (w >= padW && w<inW+padW) & (h >= padH && h<inH+padH);
      float data =  read ? dataIn[(h-padH)*inW+w-padW] : 0.0;
      //printf("(%d,%d)=%f,",h,w,data);
      dataOut[h * (inW+2*padW) + w] = data;
    }
  }
}



void lineBuffer(int inH, int inW, int kH, int kW, int stride,float *dataIn, float **dataOut){
  int outH = (inH-kH)/stride + 1;
  int outW = (inW-kW)/stride + 1;
  for (int h=0; h<outH; h++) {
    for (int w=0; w<outW; w++) {
      for (int khh=0; khh<kH; khh++) {
        for (int kww=0; kww<kW; kww++) {
          int idxH = (h*stride+khh);
          int idxW = (w*stride+kww);
          dataOut[h*outW+w][khh*kW+kww] = dataIn[idxH*inW+idxW];
        }
      }
    }
  }

}

void LineBufferCPU(int inH, int inW, int padH, int padW, int kH, int kW, int stride, float *kernel, float *dataIn, float *dataOut){
  int paddedH = inH+2*padH;
  int paddedW = inW+2*padW;
  float *padded = malloc(paddedH*paddedW*sizeof(float));
  pad(inH,inW,padH,padW,dataIn,padded); 
  int outH = (paddedH-kH) + 1;
  int outW = (paddedW-kW) + 1;
  float **lb = malloc(outH*outW*sizeof(float*));
  for (int i=0; i<outH*outW; i++) lb[i] = malloc(kH*kW*sizeof(float));
  lineBuffer(paddedH,paddedW,kH,kW,stride,padded,lb);
  print_2d_data(padded, "PADDED",paddedH,paddedW);
  free(padded);
  for (int i=0; i<outH*outW; i++) {
    dataOut[i] = 0;
    for(int j=0; j<kH*kW;j++) {
      dataOut[i] += kernel[j] *lb[i][j];
    }
  }
  for (int i=0; i<outH*outW; i++) free(lb[i]);
  free(lb);
  
}

/*
void LineBufferCPU(int size, int nx, float *dataIn, float *dataOut)
{
  const int ny = size / nx;
  for (int i = 1; i < nx - 1; i++) {
    for (int j = 1; j < ny - 1; j++) {
      dataOut[j * nx + i] =
        (dataIn[(j-1)*nx + i - 1] +
        dataIn[(j-1)*nx + i ] +
        dataIn[(j-1)*nx + i + 1] +
        dataIn[j*nx + i - 1] +
        dataIn[j*nx + i] +
        dataIn[j*nx + i + 1]  +
        dataIn[(j+1)*nx + i - 1] +
        dataIn[(j+1)*nx + i ] +
        dataIn[(j+1)*nx + i + 1]) / 9.0;
    }
  }
}
*/

inline int get_min_in(int p,int k,int s) {
  float n_float = ((2*p-k)*1.0/s+1.0)/4.0;
  int n = 1;
  if (n_float > 1) {
    n = (int)(n_float+0.5);
  }
  return (4*n-1)*s+k-2*p;
}

//adds to x no make a multiple of n
inline int adjust(int x,int n){
  return x -1 + n - ((x-1)%n);
}

int main()
{
  const int padH = Top_pad;
  const int padW = Top_pad;
  const int kH = Top_k;
  const int kW = Top_k;
  const int stride = Top_stride;


  int inH = 7;//get_min_in(padH,kH,stride);
  int inW = 7;//get_min_in(padW,kW,stride);
  printf("inH,inW=%d,%d\n",inH,inW);
  const int NY = inH+2*padH;
  const int NX = inW+2*padW;


  const float kernel[kH*kW] = {0,1,0,1,2,1,0,1,0};
      
  int paddedH = inH+2*padH;
  int paddedW = inW+2*padW;
  int outH = 1+(paddedH-kH)/stride;
  int outW = 1+(paddedW-kW)/stride;
  uint32_t inSizeBytes = adjust(inH * inW*sizeof(float),384);
  uint32_t outSizeBytes = adjust(outH * outW*sizeof(float),384);
  uint32_t inAddr = 0;
  uint32_t outAddr = inAddr + inSizeBytes;

  if ((paddedH-kH)%stride || (paddedW-kW)%stride) {
    printf("BAD!\n");
    exit(1);
  }
  printf("inSize,outSize=%d,%d\n",inSizeBytes,outSizeBytes);
  float *dataIn = malloc(inSizeBytes);
  float *dataOut = malloc(outSizeBytes);
  float *expectedOut = malloc(outSizeBytes);
  uint32_t *debug = malloc(outSizeBytes);
  
  if ( (paddedH-kH)%stride !=0 ||
       (paddedW-kW)%stride !=0 ) {
    printf("Constants are not compatible!!\n");
    exit(1);
  }

  int nxMax = Top_nxMax;
  if (NX > nxMax) {
    printf("2D filter with maximum size nxMax=%d can not process data with nx=%d\n",nxMax,NX);
    exit(1);
  }
  
  generateInputData(dataIn, inH, inW);

  LineBufferCPU(inH,inW,padH,padW,kH,kW,stride,kernel,dataIn, expectedOut);
  print_2d_data(dataIn, "INPUT DATA", inH, inW);
  print_2d_data(expectedOut, "BEFORE: EXPECTED DATA", outH, outW);
  
  max_file_t *maxfile = Top_init();
  max_engine_t *engine = max_load(maxfile, "local:*");
  
  Top_writeLMem_actions_t write;
  write.param_size = inSizeBytes;
  write.param_start = inAddr;
  write.instream_fromcpu = dataIn;
  
  Top_writeLMem_run(engine,&write);

  
  Top_actions_t lb;
  lb.param_NX = NX;
  lb.param_inH = inH;
  lb.param_inW = inW;
  lb.param_inBytes = inSizeBytes;
  lb.param_inAddr = inAddr;
  lb.param_outBytes = outSizeBytes;
  lb.param_outAddr = outAddr;

  printf("BEFORE run\n"); 
  Top_run(engine,&lb);
  printf("AFTER run\n"); 
  Top_readLMem_actions_t read;
  read.param_size = outSizeBytes;
  read.param_start = outAddr;
  read.outstream_tocpu = dataOut;

  Top_readLMem_run(engine,&read);
  printf("AFTER read\n"); 

  print_2d_data(dataOut, "OUTPUT DATA", outH, outW);
  print_2d_data(expectedOut, "EXPECTED DATA", outH, outW);

  int status = check(outW*outH, dataOut, expectedOut);

  if (status)
    printf("Test failed.\n");
  else
    printf("Test passed OK!\n");



  //cleanup
  free(dataIn);
  free(dataOut);
  max_unload(engine);



  return 0;
}
