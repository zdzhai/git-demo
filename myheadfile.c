
#include "myheadfile.h"

/**
 * @brief To read gro file which is modified
 * 
 * @param filename      The .gro filename 
 * @param TOTAL         The total rows of .gro file
 * @param data          The structure array which is used to accommodate .gro file
 * @param px            The pionter of the variable xbox
 * @param py            The pionter of the variable ybox 
 * @param pz            The pionter of the variable zbox 
 * @return stru_Gro*    To return the pointer of the structure array
 */
stru_Gro * readGro(char* filename,int TOTAL,stru_Gro data[],
                        float *px,float *py,float *pz){
    FILE *fp;
    fp = fopen(filename, "r");
    char read[10];
    int i = 0;
    for (i = 0; i <= 1; i++)
    {
        fscanf(fp, "%s\n", read);
        printf("%s\n", read);
    }
    for (i = 0; i < TOTAL; i++)
    {
        fscanf(fp, "%s%s%d%f%f%f", &data[i].name, &data[i].atom, &data[i].num,
               &data[i].x, &data[i].y, &data[i].z);
    }
    fscanf(fp, "%f%f%f\n", px, py, pz);
    // printf("%f %f %f\n", *px, *py, *pz);
    fclose(fp);
    
    return  data;
}

float pbcr(float a1,float b1,float c1,float a2,float b2,float c2)
{
	float X,Y,Z,R;
	X=a1-a2;
	Y=b1-b2;
	Z=c1-c2;
	R=sqrt(X*X+Y*Y+Z*Z);
	return R; 
}
/**
 * @brief calculate the distance between two atom
 * 
 * @param o1        index of first atom
 * @param o2        index of second atom
 * @param x          all cordination
 * @param xbox      box of x
 * @param ybox      box of y
 * @param zbox      box of z
 * @return float 
 */
float dist2(int o1,int o2,float xbox,float ybox,float zbox,rvec *x)
{
	float o1x,o1y,o1z,o2x,o2y,o2z;
	o1x=x[o1][0],o1y=x[o1][1],o1z=x[o1][2];
	o2x=x[o2][0],o2y=x[o2][1],o2z=x[o2][2];
	float X,Y,Z,r; 
	X=o2x-o1x,Y=o2y-o1y,Z=o2z-o1z;
	if(X>xbox/2){X=X-xbox;}	if(X<-xbox/2){X=X+xbox;}
	if(Y>ybox/2){Y=Y-ybox;}	if(Y<-ybox/2){Y=Y+ybox;}
	if(Z>zbox/2){Z=Z-zbox;}	if(Z<-zbox/2){Z=Z+zbox;}
	r=sqrt((X*X)+(Y*Y)+(Z*Z));
	return r;
}
float calcVdw(float sigma,float epsilon,float r){
    return 4.0*epsilon*((pow((sigma/r),12))-(pow((sigma/r),6)));
}
float calcCoul(float charge1,float charge2,float r){
    return charge1*charge2*e*e*EPSILON0/r;
}
/**
 * @brief calc HBonds energy
 * 
 * @param k idnex of SGP
 * @param i index of C=O(O)
 * @param j index of H2O(O)
 * @param r1 O...O distance
 * @param x all cordination
 * @param xbox box of x
 * @param ybox box of y
 * @param zbox box of z
 * @return float 
 */
float calcEnergy(int k,int i,int j,rvec *x,float xbox,float ybox,float zbox){
    //定义SGP CH2的sigma epsilon
    float sigC =3.75 ,sigO_3 =2.96 , sigOH = 3.0, sigOW =3.166,sigS = 3.56359,sigCT = 3.5,sigHC = 2.5;   //埃
    float epsC =0.43932 ,epsO_3= 0.87864, epsOH =0.71128,epsOW =0.6507,epsS = 1.046,epsCT = 0.276144,epsHC = 0.125520; // KJ/mol
    float chaC=0.52,chaO_3 = -0.44,chaOH = -0.53,chaHO = 0.45,chaOW=-0.8476,chaHW=0.4238;//e
    float chaS=0.0,chaCT = -0.12,chaHC = 0.06;
    float vdw_ener = 0.0,coul_ener = 0.0,all_ener = 0.0;

    float sigmaCO = sqrt(sigOW * sigC), epsilonCO = sqrt(epsOW * epsC);
    float sigmaO_3O = sqrt(sigOW * sigO_3), epsilonO_3O = sqrt(epsOW * epsO_3);
    float sigmaOHO = sqrt(sigOW * sigOH), epsilonOHO = sqrt(epsOW * epsOH);

    float sigmaSO = sqrt(sigOW * sigS), epsilonSO = sqrt(epsOW * epsS);
    float sigmaCTO = sqrt(sigOW * sigCT), epsilonCTO = sqrt(epsOW * epsCT);
    float sigmaHCO = sqrt(sigOW * sigHC), epsilonHCO = sqrt(epsOW * epsHC);

    float r1=0.0,r2=0.0,r3=0.0,r4=0.0,r5=0.0,r6=0.0;
    float r7=0.0,r8=0.0,r9=0.0,r10=0.0,r11=0.0,r12=0.0;
    r1 = dist2(i-1,j,xbox,ybox,zbox,x); //C OW
    r2 = dist2(i-1,j+1,xbox,ybox,zbox,x); //C HW1
    r3 = dist2(i-1,j+2,xbox,ybox,zbox,x); //C HW2

    r4 = dist2(i,j,xbox,ybox,zbox,x); //O_3 HW1
    r5 = dist2(i,j+1,xbox,ybox,zbox,x); //O_3 HW1
    r6 = dist2(i,j+2,xbox,ybox,zbox,x); //O_3 HW2
   
    r7 = dist2(i+1,j,xbox,ybox,zbox,x); //OH OW
    r8 = dist2(i+1,j+1,xbox,ybox,zbox,x); //OH HW1
    r9 = dist2(i+1,j+2,xbox,ybox,zbox,x); //OH HW2

    r10 = dist2(i+2,j,xbox,ybox,zbox,x); //HO OW
    r11 = dist2(i+2,j+1,xbox,ybox,zbox,x); //HO HW1
    r12 = dist2(i+2,j+2,xbox,ybox,zbox,x); //HO HW2

    vdw_ener=calcVdw(sigmaCO,epsilonCO,r1*10.0)+
            calcVdw(sigmaO_3O,epsilonO_3O,r4*10.0)+
            calcVdw(sigmaOHO,epsilonOHO,r7*10.0);
    
//	printf("vdw-%.1f %d %d %.8lf %.8lf\n",time,i,j,v,vdw_ener);
    coul_ener =  calcCoul(chaC,chaOW,r1)+calcCoul(chaC,chaHW,r2)+calcCoul(chaC,chaHW,r3)+
    calcCoul(chaO_3,chaOW,r4)+calcCoul(chaO_3,chaHW,r5)+calcCoul(chaO_3,chaHW,r6)+
    calcCoul(chaOH,chaOW,r7)+calcCoul(chaOH,chaHW,r8)+calcCoul(chaOH,chaHW,r9)+
    calcCoul(chaHO,chaOW,r10)+calcCoul(chaHO,chaHW,r11)+calcCoul(chaHO,chaHW,r12);

    //printf("coul-%.1f %d %d %.5f %.5f\n",time,i,j,c,coul_ener);
    r1 = dist2(k,j,xbox,ybox,zbox,x);//S OW
    r2 = dist2(k,j+1,xbox,ybox,zbox,x);//S HW1
    r3 = dist2(k,j+2,xbox,ybox,zbox,x);//S HW2

    vdw_ener += calcVdw(sigmaSO,epsilonSO,r1*10.0);
    int a = 0;
    for (a = k+1; a < i-1; a+=3)
    {
        r4 = dist2(a,j,xbox,ybox,zbox,x);//CT OW
        r5 = dist2(a,j+1,xbox,ybox,zbox,x);//CT HW1
        r6 = dist2(a,j+2,xbox,ybox,zbox,x);//CT HW2
        r7 = dist2(a+1,j,xbox,ybox,zbox,x);//HC OW
        r8 = dist2(a+1,j+1,xbox,ybox,zbox,x);//HC HW1
        r9 = dist2(a+1,j+2,xbox,ybox,zbox,x);//HC HW2
        r10 = dist2(a+2,j,xbox,ybox,zbox,x);//HC OW
        r11 = dist2(a+2,j+1,xbox,ybox,zbox,x);//HC HW1
        r12 = dist2(a+2,j+2,xbox,ybox,zbox,x);//HC HW2

        vdw_ener += (calcVdw(sigmaCTO,epsilonCTO,r4*10.0)+
        calcVdw(sigmaHCO,epsilonHCO,r7*10.0)+
        calcVdw(sigmaHCO,epsilonHCO,r10*10.0));

        coul_ener += (calcCoul(chaCT,chaOW,r4)+calcCoul(chaCT,chaHW,r5)+calcCoul(chaCT,chaHW,r6)+
                calcCoul(chaHC,chaOW,r7)+calcCoul(chaHC,chaHW,r8)+calcCoul(chaHC,chaHW,r9)+
                calcCoul(chaHC,chaOW,r10)+calcCoul(chaHC,chaHW,r11)+calcCoul(chaHC,chaHW,r12));
    }

    all_ener=vdw_ener+coul_ener;
    return   all_ener;
}
/**
 * @brief calculate the energy between two water molecules
 * 
 * @param i         index of first water(O)
 * @param j         index of second water(O)
 * @param x         all cordination
 * @param xbox      box of x
 * @param ybox      box of y
 * @param zbox      box of z
 * @return float    return the energy
 */
float calcOOEnergy(int i,int j,rvec *x,float xbox,float ybox,float zbox){
    float  epsOW=0.6507;//KJ/mol 
    float  sigOW=3.166;//
    float  chaOW=-0.8476,chaHW=0.4238;
    float r1=0.0,r2=0.0,r3=0.0,r4=0.0,r5=0.0,r6=0.0,r7=0.0,r8=0.0,r9=0.0;
    float vdw_ener = 0.0,coul_ener = 0.0,all_ener = 0.0;

    r1 = dist2(i,j,xbox,ybox,zbox,x);       //OW OW
    r2 = dist2(i,j+1,xbox,ybox,zbox,x);     //OW HW
    r3 = dist2(i,j+2,xbox,ybox,zbox,x);     //OW HW

    r4 = dist2(i+1,j,xbox,ybox,zbox,x);     //HW OW
    r5 = dist2(i+1,j+1,xbox,ybox,zbox,x);   //HW HW
    r6 = dist2(i+1,j+2,xbox,ybox,zbox,x);   //HW HW

    r7 = dist2(i+2,j,xbox,ybox,zbox,x);     //HW OW
    r8 = dist2(i+2,j+1,xbox,ybox,zbox,x);   //HW HW
    r9 = dist2(i+2,j+2,xbox,ybox,zbox,x);   //HW HW

    vdw_ener = calcVdw(sigOW,epsOW,r1*10);
        
    coul_ener = calcCoul(chaOW,chaOW,r1)+calcCoul(chaOW,chaHW,r2)+calcCoul(chaOW,chaHW,r3)+
                calcCoul(chaHW,chaOW,r4)+calcCoul(chaHW,chaHW,r5)+calcCoul(chaHW,chaHW,r6)+
                calcCoul(chaHW,chaOW,r7)+calcCoul(chaHW,chaHW,r8)+calcCoul(chaHW,chaHW,r9);

    all_ener=vdw_ener+coul_ener;

    return all_ener;
}