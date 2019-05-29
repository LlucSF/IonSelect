/***********************************************
 *   Copyright (C) 2016 by Esteban del Castillo   
 *   esteban.delcastillo@urv.cat   
 *
 *   Proyecto 
 *  	diciembre 2016
************************************************/
#include "ionsTreat.h"

IonsTreat::IonsTreat(int matrixRows, int matrixCols, float**matrixIn)
{
  m_matrixRows=matrixRows;
  m_matrixCols=matrixCols;
  m_matrixIn=matrixIn;
  m_pMatrix=NULL;
  m_fcMatrix=NULL;
  m_threshold=EPSILON_D;
}

IonsTreat::IonsTreat(double threshold)
{
  m_matrixRows=0;
  m_matrixCols=0;
  m_matrixIn=NULL;
  m_pMatrix=NULL;
  m_fcMatrix=NULL;
  /*
  if(threshold<=0.0)   
    m_threshold=EPSILON_D;
  else
  */
    //permite que no actue la gestión de zeros si threshold es negativo
    m_threshold=threshold;
}



IonsTreat::~IonsTreat()
{
//  printf("IonTreat destructor input\n");
  if(m_pMatrix)
    {
    for(int row=0; row<m_matrixRows; row++)
    if(m_pMatrix[row]) 
      delete []m_pMatrix[row];
    delete []m_pMatrix; 
    }
  
  if(m_fcMatrix)
    {
    for(int row=0; row<m_matrixRows; row++)
    if(m_fcMatrix[row]) 
      delete []m_fcMatrix[row];
    delete []m_fcMatrix;
    }  
//  printf("IonTreat destructor output\n");
}

float** IonsTreat::getPMatrix()
{
  float p;
  if((m_pMatrix=new float *[m_matrixRows])==NULL)
    {printf("ERROR new() in getPMatrix()\n"); return NULL;}
  for(int row=0; row< m_matrixRows; row++)
    if((m_pMatrix[row]=new float [m_matrixCols])==NULL)
      {printf("ERROR new() in getPMatrix()\n"); return NULL;}

  for(int rowA=0; rowA< m_matrixRows; rowA++)  
    for(int rowB=rowA; rowB< m_matrixRows; rowB++)
    {
     if(rowA==rowB) 
       p=0.5;
     else
      p=getStudentP(m_matrixIn[rowA], m_matrixIn[rowB], m_matrixCols, m_matrixCols);
     if(p!=-1)
      {
       m_pMatrix[rowA][rowB]=p;
       m_pMatrix[rowB][rowA]=p;       
      }
     else
      {
       m_pMatrix[rowA][rowB]=0;
       m_pMatrix[rowB][rowA]=0;       
      }
    }
    
    
//   FILE* fp;
//   fp=fopen("matrixP.txt", "w");
//   for(int rowA=0; rowA< m_matrixRows; rowA++)
//     {
//     fprintf(fp,"[%3d]", rowA);
//     for(int rowB=0; rowB< m_matrixRows; rowB++)
// //      fprintf(fp, "%5.3f ", m_pMatrix[rowA][rowB]);
//       fprintf(fp, "%+.3e ", m_pMatrix[rowA][rowB]);
//     fprintf(fp,"\n");
//     }
//   fclose(fp);
  
//   fp=fopen("matrixP_gr2.txt", "w");
//   for(int rowA=0; rowA< m_matrixRows; rowA++)
//     {
//     fprintf(fp,"[%3d]", rowA);
//     for(int rowB=0; rowB< m_matrixRows; rowB++)
//       if(m_pMatrix[rowA][rowB]>=0.25)
// 	fprintf(fp, "%d ",rowB);
//     fprintf(fp,"\n");
//     }
//   fclose(fp);
    
 return m_pMatrix;
}

double IonsTreat::getNormalP(double *dataA, double* dataB, int sizeA, int sizeB)
{
 double meanA, meanB, varA, varB, A, B, P;
  if(sizeA<=1 || sizeB<=1)
    return -1;
 
  meanA = gsl_stats_mean(dataA, 1, sizeA);
  meanB = gsl_stats_mean(dataB, 1, sizeB);
  varA = gsl_stats_variance(dataA, 1, sizeA);    
  varB = gsl_stats_variance(dataB, 1, sizeB);    
 
  A=varA/sizeA;
  B=varB/sizeB;
  
  double z;
  if(A+B<EPSILON_LD) A=EPSILON_LD;
  z=(meanA-meanB)/sqrt(A+B);  
  
  if(z>0)
    P=2*gsl_cdf_ugaussian_Q (z);	//p
  else
    P=2*gsl_cdf_ugaussian_P (z);	//p
  return (float)P;
}


double IonsTreat::getStudentP(float *dataA_f, float* dataB_f, int sizeA, int sizeB)
{
  double meanA, meanB, varA, varB, stdDevA, stdB, P, A, B, C, D, t, dof;
  if(sizeA<=1 || sizeB<=1)
    return -1;
  double dataA[sizeA];
  double dataB[sizeA];
  
  for(int i=0; i<sizeA; i++)
    dataA[i]=(double)dataA_f[i];
  for(int i=0; i<sizeB; i++)
    dataB[i]=(double)dataB_f[i];
  
  meanA = gsl_stats_mean(dataA, 1, sizeA);
  meanB = gsl_stats_mean(dataB, 1, sizeB);
  varA = gsl_stats_variance(dataA, 1, sizeA);    
  varB = gsl_stats_variance(dataB, 1, sizeB);    

 
  A=varA/sizeA;
  B=varB/sizeB;
  C=(A*A)/(sizeA-1);
  D=(B*B)/(sizeB-1);
  
   //Welch test
  if(A+B<EPSILON_D) 
    t=(meanA-meanB)/sqrt(EPSILON_LD); 	//t-test
  else
    t=(meanA-meanB)/sqrt(A+B); 	//t-test
    
  if(C+D<EPSILON_D)
    dof=((A+B)*(A+B))/EPSILON_LD;	//grados de libertad
  else
    dof=((A+B)*(A+B))/(C+D);	//grados de libertad
    
  if(t>0)
    P=2.0*gsl_cdf_tdist_Q(t, dof);	//p
  else
    P=2.0*gsl_cdf_tdist_P(t, dof);	//p
  return (float)P;
  
}

double IonsTreat::getStudentP(double *dataA, double* dataB, int sizeA, int sizeB)
{
  double meanA, meanB, varA, varB, stdDevA, stdB, P, A, B, C, D, t, dof;
  if(sizeA<=1 || sizeB<=1)
    return -1;
 
  meanA = gsl_stats_mean(dataA, 1, sizeA);
  meanB = gsl_stats_mean(dataB, 1, sizeB);
  varA = gsl_stats_variance(dataA, 1, sizeA);    
  varB = gsl_stats_variance(dataB, 1, sizeB);    

 
  A=varA/sizeA;
  B=varB/sizeB;
  C=(A*A)/(sizeA-1);
  D=(B*B)/(sizeB-1);
  
  //Welch test
  if(A+B<EPSILON_D) 
    t=(meanA-meanB)/sqrt(EPSILON_LD); 	//t-test
  else
    t=(meanA-meanB)/sqrt(A+B); 	//t-test
    
  if(C+D<EPSILON_D)
    dof=((A+B)*(A+B))/EPSILON_LD;	//grados de libertad
  else
    dof=((A+B)*(A+B))/(C+D);	//grados de libertad
  
  if(t>0)
    P=2.0*gsl_cdf_tdist_Q(t, dof);	//p
  else
    P=2.0*gsl_cdf_tdist_P(t, dof);	//p
  return (float)P;
  
// return (float)t; 
}

float ** IonsTreat::getFoldChangeMatrix()
{
  float fc;
  if((m_fcMatrix=new float *[m_matrixRows])==NULL)
    {printf("ERROR new() in getPMatrix()\n"); return NULL;}
  for(int row=0; row< m_matrixRows; row++)
    if((m_fcMatrix[row]=new float [m_matrixCols])==NULL)
      {printf("ERROR new() in getPMatrix()\n"); return NULL;}

  for(int rowA=0; rowA< m_matrixRows; rowA++)  
    for(int rowB=rowA; rowB< m_matrixRows; rowB++)
    {
     if(rowA==rowB) 
       fc=1;
     else
      fc=getFoldChange(m_matrixIn[rowA], m_matrixIn[rowB], m_matrixCols, m_matrixCols);
     if(fc!=-1)
      {
       m_fcMatrix[rowA][rowB]=fc;
       m_fcMatrix[rowB][rowA]=fc;       
      }
     else
      {
       m_fcMatrix[rowA][rowB]=0;
       m_fcMatrix[rowB][rowA]=0;       
      }
    }
    
  // FILE* fp;
  // fp=fopen("matrixFC.txt", "w");
  // for(int rowA=0; rowA< m_matrixRows; rowA++)
  //   {
  //   fprintf(fp,"[%3d]", rowA);
  //   for(int rowB=0; rowB< m_matrixRows; rowB++)
  //     fprintf(fp, "%5.1f ", m_fcMatrix[rowA][rowB]);
  //   fprintf(fp,"\n");
  //   }
 return m_fcMatrix;

}


//retorna el fold change (cociente de medianas) de dos distribuciones con elementos del tipo float 
//Los valores nulos no contribuyen
double IonsTreat::getFoldChange(float *dataA_f, float* dataB_f, int sizeA, int sizeB)
{
  double medianA, medianB, fc;
  int iA=0, iB=0, sizeA2=0, sizeB2=0;
//  double dataA[]={0, 0, 0, 0, 0, 0, 0, 15, 20, 5, 7, 22, 3, 22};
//  double dataB[]={0, 0, 0, 9, 10, 7, 1, 4, 22};
//  int sizeA=14, sizeB=9;
  double dataA2[sizeA];
  double dataB2[sizeB];
  double signo=1.0;
  
  if(sizeA<=1 || sizeB<=1)
    return -1;

  for(iA=0; iA<sizeA; iA++)
    {
    if(fabs(dataA_f[iA])<=m_threshold) continue;//se eliminan ceros
    dataA2[sizeA2++]=(double)dataA_f[iA];
    }
    
  for(iB=0; iB<sizeB; iB++)
    {
    if(fabs(dataB_f[iB])<=m_threshold) continue;//se eliminan ceros
    dataB2[sizeB2++]=(double)dataB_f[iB];
    }
    
  if(sizeA2/sizeA < 0.1 || sizeB2/sizeB < 0.1) //si el tamaño de datos es pequeño frente a los zeros
      signo= -1.0;
  
  gsl_sort (dataA2, 1, sizeA);  
  gsl_sort (dataB2, 1, sizeB);  
  
  medianA=gsl_stats_median_from_sorted_data (dataA2, 1, sizeA2); //se eliminan ceros
  medianB=gsl_stats_median_from_sorted_data (dataB2, 1, sizeB2);
  
  if (medianB==0) medianB=EPSILON_D;
  else if(fabs(medianB)<EPSILON_D && medianB<0) medianB=-EPSILON_D;
  else if(fabs(medianB)<EPSILON_D && medianB>0) medianB=EPSILON_D;

  fc=medianA/medianB;
  if(isnan(fc)) 
    //printf("%e %e %d %d\n", medianA, medianB, sizeA2, sizeB2);
  if(fc>=9.9999 && signo>0) 
    fc=9.9999;
  else if(fc>=9.9999 && signo<0) 
    fc=9.9990;
  else if(fc<0.1) fc=0.1;
  
  return (float)(fc*signo);
  }

//retorna el fold change (cociente de medianas) de dos distribuciones con elementos del tipo double
//Los valores nulos no contribuyen
double IonsTreat::getFoldChange(double *dataA, double* dataB, int sizeA, int sizeB)
{
  double medianA, medianB, fc;
  int iA=0, iB=0, sizeA2=0, sizeB2=0;
//  double dataA[]={0, 0, 0, 0, 0, 0, 0, 15, 20, 5, 7, 22, 3, 22};
//  double dataB[]={0, 0, 0, 9, 10, 7, 1, 4, 22};
//  int sizeA=14, sizeB=9;
  double dataA2[sizeA];
  double dataB2[sizeB];
  double signo=1.0;
  
  if(sizeA<=1 || sizeB<=1)
    return -1;

  gsl_sort (dataA, 1, sizeA);   //ordenación de los datos
  gsl_sort (dataB, 1, sizeB);  
  
  for(iA=0; iA<sizeA; iA++)
    {
    if(fabs(dataA[iA])<=m_threshold) continue;//se eliminan ceros
    dataA2[sizeA2++]=dataA[iA];
    }
    
  for(iB=0; iB<sizeB; iB++)
    {
    if(fabs(dataB[iB])<=m_threshold) continue;//se eliminan ceros
    dataB2[sizeB2++]=dataB[iB];
    }
    
  if((float)sizeA2/sizeA < 0.1 || (float)sizeB2/sizeB < 0.1) //si el tamaño de datos es pequeño frente a los zeros
      signo= -1.0;
  
  medianA=gsl_stats_median_from_sorted_data (dataA2, 1, sizeA2); //se eliminan ceros
  medianB=gsl_stats_median_from_sorted_data (dataB2, 1, sizeB2);
  
  if (medianB==0) medianB=EPSILON_LD;
  else if(fabs(medianB)<EPSILON_LD && medianB<0) medianB=-EPSILON_LD;
  else if(fabs(medianB)<EPSILON_LD && medianB>0) medianB=EPSILON_LD;

  fc=medianA/medianB;
  if(isnan(fc)) 
    //printf("%e %e %d %d\n", medianA, medianB, sizeA2, sizeB2);
  if(fc>=9.9999 && signo>0) 
    fc=9.9999;
  else if(fc>=9.9999 && signo<0) 
    fc=9.9990;
  else if(fc<0.1) fc=0.1;
  
  return (float)(fc*signo);
}

double IonsTreat::getHistoP(double *ionsA_p, double *ionsB_p, int pSize, int sSize, double ionMin, double ionMax, int nBins)
{
  double binValueA, binValueB, cumula=0, cumulaA=0, cumulaB=0, maxDen, P;
  
  if(pSize<=1 || sSize<=1)
    return 0;
  
  gsl_histogram * hA = gsl_histogram_alloc (nBins);
  gsl_histogram * hB = gsl_histogram_alloc (nBins);
  if(ionMin>=ionMax) 
    {
    printf("ERROR: xMin (%.4f) must be less than xMax(%.4f) in gsl_histogram_set_ranges_uniform()\n", ionMin, ionMax);
    gsl_histogram_free (hA);
    gsl_histogram_free (hB);
    return 00;
    }
  
  gsl_histogram_set_ranges_uniform (hA, ionMin, ionMax);
  gsl_histogram_set_ranges_uniform (hB, ionMin, ionMax);
  
  for(int pix=0; pix<pSize; pix++) //info a histograma
    if(gsl_histogram_increment(hA, ionsA_p[pix])!=0)
      printf("\n");
  
  for(int pix=0; pix<sSize; pix++) //info a histograma
    if(gsl_histogram_increment(hB, ionsB_p[pix])!=0)
      printf("\n");
double C=0;  
  for(int i=0; i<nBins; i++)
    {
    binValueA=gsl_histogram_get (hA, i)/pSize;
    binValueB=gsl_histogram_get (hB, i)/sSize;
    cumula+=binValueA*binValueB;
    cumulaA+=binValueA*binValueA;
    cumulaB+=binValueB*binValueB;
  C+=binValueA;
    }
  //Normalización
  maxDen=(cumulaA > cumulaB) ? cumulaA: cumulaB;
  if(maxDen<EPSILON_D) 
    P=0;
  else
    P=cumula/maxDen;

    
//  P=cumula;
//printf("cumula:%.4f cumulaA:%.4f cumulaB:%.4f C::%.4f\n", cumula, cumulaA, cumulaB, C);    
//printf("U:%+10.3f mean:%+10.3f  sigma:%10.3f z:%+10.4f p:%10.4f\n",cumula,  mean, sigma, z, P);    
  gsl_histogram_free (hA);
  gsl_histogram_free (hB);
  return P;
}

//Test U de Mann-Whitney
//Recibe dos array con muestras y determina si esas muestras provienen de la misma distribución 
//Retorna valores entre 0 y 1: 1 si la distribución original es la misma para ambos arrays muestras
//Retorna valor negativo si el resultado no se deriva de datos consistentes
//'sigma' determina la desviación estándar a utilizar
//en 'zeroRate' se retorna el cociente entre los zeros de ambos arrays. Limitado entre 0.1 y 9.999.
//'zeroRate' vale negativo si el cociente no se deriva de datos consistentes
double IonsTreat::getMannWhitneyUTest(double *dataA, double* dataB, int sizeA, int sizeB, double *zeroRate, double sigma)
{  
  bool indexZeroA[sizeA], indexZeroB[sizeB];
  double p1;

  int countA_0=0, countB_0=0;
  int countA_1, countB_1;
  double zeroFactor;
  
  countA_0=getZeros(dataA, sizeA, NULL);
  countB_0=getZeros(dataB, sizeB, NULL);

  countA_1=sizeA-countA_0;
  countB_1=sizeB-countB_0;
  double dataA_1[countA_1]; //aloja la parte de no ceros de A
  double dataB_1[countB_1];

  //Se separan los datos de forma que dataA_1 y dataB_1 mantengan las listas de nozeros 
  getNotZeros(dataA, sizeA, dataA_1);  
  getNotZeros(dataB, sizeB, dataB_1);  
  
  //tratamiento de los valores nulos
  double z1, rateA, rateB, signo=1.0;
  if(sizeA==0 && sizeB==0) *zeroRate=-1; //no se pueden tomar decisiones
  else if(countA_0==0 && countB_0==0) *zeroRate=-1; //no se pueden tomar decisiones
  else
    {
    if(sizeA==0) sizeA=EPSILON_LD;
    if(sizeB==0) sizeB=EPSILON_LD;
//    rateA=(double)countA_0/sizeA; //se consideran valores nulos
//    rateB=(double)countB_0/sizeB;
    rateA=1.0-(double)countA_1/sizeA; //se consideran valores no nulos
    rateB=1.0-(double)countB_1/sizeB;
    if(rateA<0.1 && rateB<0.1)
      signo=-1.0;
    if(rateB<EPSILON_LD)
      rateB=EPSILON_LD;
    *zeroRate=rateA/rateB;
    if(*zeroRate>9.999)  *zeroRate=9.999;
    else if(*zeroRate<0.1) *zeroRate=0.1;
    *zeroRate*=signo;
    }
    
  //Se marcan como negativos los resultados con pocos datos relativos > 0
  signo=1.0;
  if((double)countA_1/sizeA<0.1 || (double)countB_1/sizeB<0.1)
    signo=-1.0;
  z1=getMannWhitneyZ(dataA_1, dataB_1, countA_1, countB_1); //test 'Z' de Mann-Whitney
    
  //Se determina la probabilidad de que se dé ese valor de Z según la sigma pasada
    if(z1<=0) 
      p1=2.0*gsl_cdf_gaussian_P(z1, sigma); 
    else 
      p1=2.0*gsl_cdf_gaussian_Q(z1, sigma); 
    
  //Para determinar los niveles de corte, se excluyen los valores nulos (valor=0.0)
  //Posteriormente, esos valores sí se concideran en la selección
  //se observa que los valores de corte están más equilibrados (Z, FC y P) si se dejan como cero las p que sean cero
//    if(p1==0) p1=1e-300; 
    
  return (p1*signo);
} 
  
//Test Z de Mann-Whitney
//Recibe dos array con muestras y determina si esas muestras provienen de la misma distribución 
//Retorna un valor nulo si si la distribución original es la misma para ambos arrays muestras
//La varioable aleatoria retornada está normalizada según gausiana de media cero y sigma 1
double IonsTreat::getMannWhitneyZ(double *dataA, double* dataB, int sizeA, int sizeB)
{
  int iA=0, iB=0, iC=0;
  bool input;
  double Ua, Ub, P, Pa, Pb, acu;
//  double dataA[]={0, 0, 0, 0, 15, 20, 5, 7, 22, 3, 22};
//  double dataB[]={0, 0, 0, 9, 10, 7, 1, 4, 22};
//  int sizeA=11, sizeB=9;
  if(sizeA==sizeB && sizeA==0) return 0;
  if(sizeA==0 || sizeB==0) return 1000; //nada parecidos
  double array[sizeA+sizeB];
  bool   group[sizeA+sizeB];
  
  //ordenación de menor a mayor
  gsl_sort (dataA, 1, sizeA);  
  gsl_sort (dataB, 1, sizeB); 

  //intercadalo de elementos
  for(int i=0; i<sizeA+sizeB; i++)
    {
    if(dataA[iA]<=dataB[iB]) 
      {
      array[iC]=dataA[iA];
      group[iC]=false; //el elemento pertenece al grupo 0
      iA++;
      }
    else 
      {
      array[iC]=dataB[iB]; 
      group[iC]=true; //el elemento pertenece al grupo 1
      iB++;
      }
    iC++;
    if(iA>=sizeA)
      {
      while(iB<sizeB) {array[iC]=dataB[iB++]; group[iC++]=1;}
      break;
      }
    if(iB>=sizeB)
      {
      while(iA<sizeA) {array[iC]=dataA[iA++]; group[iC++]=0;}
      break;
      }
    }

  //trato de elementos redundantes
  iA=0; iB=0; 
  input=false;
  double nLinks, acuLinks=0;;
  for(int i=0; i<sizeA+sizeB; i++)
    {
      if(i<sizeA+sizeB-1 && array[i]==array[i+1])
        {
        if(input==true) iB++; //si dentro de ligaduras (elementos redundantes)
        else {iA=i; iB=i+1; input=true; acu=0;} //inicio de ligadura
        }
      else if(input==true) //fin de ligadura: se establecen sus rangos (valor promediado)
        {
        input=false;
        nLinks=iB-iA+1;
        acuLinks+=nLinks*(nLinks*nLinks-1); //para determinar la varianza
        
        for(int j=iA; j<=iB; j++)
        acu+=j+1; //se suman los rangos (el primero es '1' de ahí '+1')
        for(int j=iA; j<=iB; j++)
        array[j]=acu/(1.0+iB-iA); //se establecen los nuevos rangos promedio
        }
      else
        array[i]=i+1; //rango corregido (el primero es '1' de ahí '+1')
	  
    }
     
  double acuRankA=0, acuRankB=0;
  for(int i=0; i<sizeA+sizeB; i++)
    {
    if(group[i]==0) acuRankA+=array[i];//el elemento pertenece al grupo cero
    else acuRankB+=array[i];//el elemento pertenece al grupo uno
    }
  Ua=acuRankA-(sizeA*(sizeA+1)/2.0);
  Ub=acuRankB-(sizeB*(sizeB+1)/2.0);
    double A=sizeA+sizeB;
    double stdDev=sqrt(((sizeA*sizeB)/12.0)*(A+1-acuLinks/(A*(A-1))));
    if(stdDev<EPSILON_D) stdDev=EPSILON_D;
    
    //Normalización del estadístico
    //Se relativiza para que z=0 si Ua==media y z=1 si Ua se desvia 1 stdDev
    //ello permite aplicar la gausiana de media cero y stdDev=1
    double z=(Ua-(sizeA*sizeB)/2.0)/(stdDev);
//    double z1=(Ub-(sizeA*sizeB)/2.0)/(stdDev);
//    printf("sizeA:%3d sizeB:%3d Ua:%7.3f Ub:%7.3f stdDev:%7.3f, acuLinks:%7.3f z%7.3f, z1:%7.3f\n", sizeA, sizeB, Ua, Ub, stdDev, acuLinks, z, z1);
  return z;
}

//Retorna la cantidad de elementos con valor nulo en data
//Si zerosIndex!=NULL, establece los índices de los elementos nulos en data
int IonsTreat::getZeros(double *data, int size, int *zerosIndex)
{
int k=0;
  if(zerosIndex)
    {
    for(int i=0; i<size; i++)
      {
//printf("%e\n", fabs(data[i]));  
	
      if(fabs(data[i])<=m_threshold)
	zerosIndex[k++]=i;
      }
    }
  else
    {
    for(int i=0; i<size; i++)
      {
      if(fabs(data[i])<=m_threshold)
	k++;
//printf("%e\n", fabs(data[i]));  
      }
    }
  return k;
}

//Retorna la cantidad de elementos con valor no nulo en data
//Si notZeros!=NULL, copia los elementos no nulos de data
int IonsTreat::getNotZeros(double *data, int size, double *notZeros)
{
int k=0;
  if(notZeros)
    {
    for(int i=0; i<size; i++)
      if(fabs(data[i])>m_threshold)
	notZeros[k++]=data[i];
    }
  else
    {
    for(int i=0; i<size; i++)
      if(fabs(data[i])>m_threshold)
	k++;
    }
  return k;
}


