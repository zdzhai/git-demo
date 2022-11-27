

#ifndef __MYHEADFILE_H__
#define __MYHEADFILE_H__

#include "xdrfile.h"

#define EPSILON0 138.935485
#define     e    1.0
rvec *x;

typedef struct gro_info
{
    char name[10];
    char atom[10];
    int num;
    float x;
    float y;
    float z;
} stru_Gro;

stru_Gro * readGro(char* filename,int TOTAL,stru_Gro data[],
                        float *px,float *py,float *pz);
float pbcr(float a1,float b1,float c1,float a2,float b2,float c2);
float dist2(int o1,int o2,float xbox,float ybox,float zbox,rvec *x);
float calcVdw(float sigma,float epsilon,float r);
float calcCoul(float charge1,float charge2,float r);
float calcEnergy(int k,int i,int j,rvec *x,float xbox,float ybox,float zbox);

#endif
