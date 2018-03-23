/***********************************************
 *   Copyright (C) 2017 by Esteban del Castillo
 *   esteban.delcastillo@urv.cat
 *
 *   Proyecto Ions_Select
 *  	julio 2017
************************************************/
#include "ionsSelect.h"

//constructor
//recibe información sobre las muestras a considerar (lsp_p & nSamples), sobre los grupos a analizar (pTestGroups & nPTestGroups)
//sobre los grupos de la etapa de segmentación (allGroups_p) y sobre las matriz de magnitudes pixels-ions (allSpectrums_p)
IonsSelect::IonsSelect(LOAD_SP *lsp_p, int nSamples, int *pTestGroups, int nPTestGroups, GROUP* allGroups_p, float** allSpectrums_p)
{
  m_error=false;
  m_nPixels=0;
  for(int i=0; i<nSamples; i++)
    m_nPixels+=lsp_p[i].nPixels;
  m_pTestRows=lsp_p[nSamples-1].mzWidth;
  m_mzAxis_p=lsp_p[nSamples-1].mzAxis_p;
  m_pTestCols=nPTestGroups*nPTestGroups;

  for(int i=0; i<nPTestGroups; i++) //copia de grupos totales
    m_pTestGroups[i]=pTestGroups[i];
  m_nPTestGroups=nPTestGroups;

  m_allGroups_p=allGroups_p;
  m_allSpectrums_p=allSpectrums_p;
  m_ionsData_p=NULL;
  m_ionsClusterV_p=NULL;
  m_ionsClusterZ_p=NULL;
  m_probability=-1;

    if((m_ionsClusterV_p=new int*[m_pTestRows])==NULL)
      {printf("ERROR new() in p-test\n"); m_error=true;}
    for(int i=0; i<m_pTestRows; i++)
      m_ionsClusterV_p[i]=NULL;
    for(int i=0; i<m_pTestRows; i++)
      if((m_ionsClusterV_p[i]=new int[nPTestGroups])==NULL)
	{printf("ERROR new() in p-test\n"); m_error=true;}


    if((m_ionsClusterZ_p=new int*[m_pTestRows])==NULL)
      {printf("ERROR new() in p-test\n"); m_error=true;}
    for(int i=0; i<m_pTestRows; i++)
      m_ionsClusterZ_p[i]=NULL;
    for(int i=0; i<m_pTestRows; i++)
      if((m_ionsClusterZ_p[i]=new int[nPTestGroups])==NULL)
	{printf("ERROR new() in p-test\n"); m_error=true;}

}

//Destruye memoria reservada
IonsSelect::~IonsSelect()
{
//printf("IonsSelect destructor start\n");
  //memoria para matriz IonsData
  if(m_ionsData_p)
    {
    for(int i=0; i<m_pTestRows; i++)
      if(m_ionsData_p[i]) delete m_ionsData_p[i];
    delete [] m_ionsData_p;
    m_ionsData_p=NULL;
    }

  if(m_ionsClusterV_p)
    {
    for(int i=0; i<m_pTestRows; i++)
      if(m_ionsClusterV_p[i]) delete []m_ionsClusterV_p[i];
    if(m_ionsClusterV_p) delete []m_ionsClusterV_p;
    }

  if(m_ionsClusterZ_p)
    {
    for(int i=0; i<m_pTestRows; i++)
      if(m_ionsClusterZ_p[i]) delete []m_ionsClusterZ_p[i];
    if(m_ionsClusterZ_p) delete []m_ionsClusterZ_p;
    }

//printf("IonsSelect destructor end\n");
}

//Extrae las medidas de z, P y FC
//Genera la matriz m_ionsData_p a partir de los datos en las matrices m_allSpectrums_p y de grupos m_allGroups_p
//m_ionsData_p: filas=iones; columnas=variables Z, P, FC para cada combinación de dos grupos de m_allGroups_p (0/0, 0/1, 0/2..0/n; 1/0, 1/1..)
//-1 en m_ionsData_p indica que el dato no es válido
int IonsSelect::getMeasures()
  {
  int maxGroupSize=0;
  double *ionsA_p=NULL;
  double *ionsB_p=NULL;
  int pGr, sGr, pSize, sSize;
  int pGrSize, sGrSize, pixel;
  double pixelValue;
  IonsTreat ionsTreat;

  //se detecta el tamaño del grupo mayor
  for(int i=0; i<m_nPTestGroups; i++)
    if(m_allGroups_p[m_pTestGroups[i]].size>maxGroupSize)
      maxGroupSize=m_allGroups_p[m_pTestGroups[i]].size;

  //memoria para array de iones a comparar
  if((ionsA_p=new double[maxGroupSize])==NULL)
    {printf("ERROR new() in getMeasures()"); return -1;}
  if((ionsB_p=new double[maxGroupSize])==NULL)
    {printf("ERROR new() in getMeasures()"); return -1;}

  //Memoria para la matriz ionsData
  if((m_ionsData_p=new double*[m_pTestRows])==NULL) //iones
    {printf("ERROR new() in getMeasures()"); return -1;}
  for(int i=0; i<m_pTestRows; i++)
    if((m_ionsData_p[i]=new double[m_pTestCols*3])==NULL)//combinaciones de todos los grupos
      {printf("ERROR new() in getMeasures()"); return -1;}

double p0, p1, zeroRate;
int col=0;
  //Se da valor a una matriz de filas=nº iones y columnas=nGrops*nGroups (todas las combinaciones posibles!!!!)
  for(int pGri=0; pGri<m_nPTestGroups; pGri++) //para todos los elementos de los grupos considerados
     {
    pGr=m_pTestGroups[pGri]; //grupo a tratar (primario)
    pGrSize=m_allGroups_p[pGr].size;

    for(int sGri=0; sGri<m_nPTestGroups; sGri++) //para todos los elementos de los grupos considerados
      {
      sGr=m_pTestGroups[sGri]; //grupo a tratar (secundario)
      sGrSize=m_allGroups_p[sGr].size;
      for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
	{
	pSize=0;
	for(int pix=0; pix<pGrSize; pix++) //para todos los pixeles del grupo primario
	  {
	   pixel=m_allGroups_p[pGr].set[pix]; //un pixel del grupo primario
	   pixelValue=m_allSpectrums_p[pixel][ion]; //magnitud de un ion de ese pixel
	   ionsA_p[pSize++]=pixelValue;
	  }
	sSize=0;
	for(int pix=0; pix<sGrSize; pix++) //para todos los pixeles del grupo secundario
	  {
	   pixel=m_allGroups_p[sGr].set[pix]; //un pixel del grupo secundario
	   pixelValue=m_allSpectrums_p[pixel][ion]; //magnitud de un ion de ese pixel
	   ionsB_p[sSize++]=pixelValue;
	  }

	m_ionsData_p[ion][col+1]=ionsTreat.getMannWhitneyUTest(ionsA_p, ionsB_p, pSize, sSize, &zeroRate, 1);
	m_ionsData_p[ion][col+0]=zeroRate;
	m_ionsData_p[ion][col+2]=ionsTreat.getFoldChange(ionsA_p, ionsB_p, pSize, sSize);
	}
      col+=3;
      }
    }
    if(ionsA_p) delete []ionsA_p;
    if(ionsB_p) delete []ionsB_p;
  }


  //Se determinan los limites en Z, P y FC para la selección de iones
  //Recibe la probabilidad (%) a la que debe limitarse los valores de corte
  //Se integra el histograma desde ambos extremos hasta que se alcanza el valor de probabilidad dado
  //Luego se determina el valor en X (magnitud) correspondiente. Sirve de valor de corte
  //La idea es hacer que cada variable considere un rango de valores cuya probabilidad sea menor o igual a una dada.
  //Como FC y P se consideran simultáneamente, la probabilida conjunta es el producto de ambas y se ajusta
  //Se hace uso de histogramas ante la falta de normalidad en las distribuciones
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int IonsSelect::setBoundaries(double probability)
{
  int pCol, pGr, sGr;
  double min[3], max[3], tmp[3];
//  int eCount[3];
  int nBins[3]={1000, 1000, 1000};
  gsl_histogram * h[3];
  int hSize[3]={0, 0, 0};
  double prob[3];
  double acu[6]={0,0,0, 0,0,0};
  int bin[6];
  double xValue[6];

  m_probability=probability;

  prob[0]=m_probability/100.0;
  prob[1]=sqrt(prob[0]); //se ajusta la prob de FC y P para que su producto coincida con el dado (que se asigna a Z)
  prob[2]=prob[1];
  min[0]=1e300; min[1]=1e300; min[2]=1e300;
  max[0]=-1e300; max[1]=-1e300; max[2]=-1e300;

 //Valores extremos
 for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
  {
  for(int pGri=0; pGri<m_nPTestGroups; pGri++) //para todos los grupos
    {
    pGr=pGri;
    pCol=m_nPTestGroups*3*pGr;
    for(int sGri=0; sGri<m_nPTestGroups; sGri++)//para todos los grupos
      {
      sGr=sGri;
      if(sGr==pGr) continue;
      tmp[0]=m_ionsData_p[ion][pCol+sGri*3+0];
      tmp[1]=m_ionsData_p[ion][pCol+sGri*3+1];
      tmp[2]=m_ionsData_p[ion][pCol+sGri*3+2];

      if(tmp[0]>0)
	{
	if(tmp[0]<min[0]) min[0]=tmp[0];
	if(tmp[0]>max[0]) max[0]=tmp[0];
	hSize[0]++;
	}
      if(tmp[1]>0)
	{
	  double a=log10(tmp[1]);
	if(a<min[1]) min[1]=a;
	if(a>max[1]) max[1]=a;
	hSize[1]++;
	}
      if(tmp[2]>0)
	{
	if(tmp[2]<min[2]) min[2]=tmp[2];
	if(tmp[2]>max[2]) max[2]=tmp[2];
	hSize[2]++;
	}
      }
    }
  }
  //Se evalua la posibilidad de que ningún valor sea válido
  for(int i=0; i<3; i++)
    {
    if(min[i]==1e300) min[i]=-1;
    if(max[i]==-1e300) max[i]=-0.9;
    if(min[i]>=max[i])
      {min[i]=-1; max[i]=-0.9;}//evita posibles errores
    }

  //Histogramas: creación
  h[0] = gsl_histogram_alloc (nBins[0]);
  h[1] = gsl_histogram_alloc (nBins[1]);
  h[2] = gsl_histogram_alloc (nBins[2]);

  gsl_histogram_set_ranges_uniform (h[0], min[0], max[0]);
  gsl_histogram_set_ranges_uniform (h[1], min[1], max[1]);
  gsl_histogram_set_ranges_uniform (h[2], min[2], max[2]);

 //Histogramas: inicialización
 for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
  {
  for(int pGri=0; pGri<m_nPTestGroups; pGri++) //para todos los grupos
    {
    pGr=pGri;
    pCol=m_nPTestGroups*3*pGr;
    for(int sGri=0; sGri<m_nPTestGroups; sGri++)//para todos los grupos
      {
      sGr=sGri;
      if(sGr==pGr) continue;
      tmp[0]=m_ionsData_p[ion][pCol+sGri*3+0]; //Z
      tmp[1]=m_ionsData_p[ion][pCol+sGri*3+1]; //P
      tmp[2]=m_ionsData_p[ion][pCol+sGri*3+2]; //FC
      if(tmp[0]>0)
	gsl_histogram_increment(h[0], tmp[0]);
      if(tmp[1]>0)
	gsl_histogram_increment(h[1], log10(tmp[1])); //mejora el resultado con log
      if(tmp[2]>0)
	gsl_histogram_increment(h[2], tmp[2]);
      }
    }
  }

  //Se acumula la probabilidad desde ambos extremos de la distribución
  //eliminando los bins extremos para Z y FC que están saturados
  //límites para Z
  for(int i=1; i<nBins[0]-1; i++) //se eliminan los extremos saturados
    if(acu[0]<=prob[0])
      {
      acu[0]+=gsl_histogram_get (h[0], i)/hSize[0];
      bin[0]=i; //índice a la posición en X que acumula la probabilidad dada por la izquierda
      }
    else break;
  for(int i=nBins[0]-2; i>0; i--) //se eliminan los extremos saturados
    if(acu[3]<=prob[0])
      {
      acu[3]+=gsl_histogram_get (h[0], i)/hSize[0];
      bin[3]=i; //índice a la posición en X que acumula la probabilidad dada por la derecha
      }
    else break;

  //límites para p
  for(int i=1; i<nBins[1]-1; i++) //se eliminan los extremos saturados
    if(acu[1]<=prob[1])
      {
      acu[1]+=gsl_histogram_get (h[1], i)/hSize[1];
      bin[1]=i;
      }
    else break;
  for(int i=nBins[1]-2; i>0; i--) //se eliminan los extremos saturados
    if(acu[4]<=prob[1])
      {
      acu[4]+=gsl_histogram_get (h[1], i)/hSize[1];
      bin[4]=i;
      }
    else break;
  //límites para FC
  for(int i=1; i<nBins[2]-1; i++) //se eliminan los extremos saturados
    if(acu[2]<=prob[2])//se eliminan los extremos saturados
      {
      acu[2]+=gsl_histogram_get (h[2], i)/hSize[2];
      bin[2]=i;
      }
    else break;
  for(int i=nBins[2]-2; i>0; i--) //se eliminan los extremos saturados
    if(acu[5]<=prob[2])
      {
      acu[5]+=gsl_histogram_get (h[2], i)/hSize[2];
      bin[5]=i;
      }
    else break;

  //valores de las variables Z, P y FC correspondientes a la probabilidad acumulada deseada
  double delta_Z=(max[0]-min[0])/(double)nBins[0];
  double delta_Vx=(max[2]-min[2])/(double)nBins[2];
  double delta_Vy=(max[1]-min[1])/(double)nBins[1];
  double b;

  xValue[0]=min[0]+((double)bin[0]+0.5)*delta_Z; //Z low

   b=min[1]+fabs(((double)bin[1]+0.5)*delta_Vy);
  xValue[1]=pow(10.0, b); //Vy low (estaba en escala logaritmica)

  xValue[2]=min[2]+((double)bin[2]+0.5)*delta_Vx; //Vx low

  xValue[3]=min[0]+((double)bin[3]+0.5)*delta_Z; //Z high

   b=min[1]+fabs(((double)bin[4]+0.5)*delta_Vy);
  xValue[4]=pow(10.0, b); //Vy high (estaba en escala logaritmica)

  xValue[5]=min[2]+((double)bin[5]+0.5)*delta_Vx; //Vx high

  m_xLowLimit_Z=xValue[0];
  m_xHighLimit_Z=xValue[3];
  m_yLowLimit_V=xValue[1];
  m_yHighLimit_V=xValue[4];
  m_xLowLimit_V=xValue[2];
  m_xHighLimit_V=xValue[5];

  if(m_xHighLimit_Z>10) m_xHighLimit_Z=10;
  if(m_xLowLimit_Z<0.1) m_xLowLimit_Z=0.1;
  if(m_xHighLimit_V>10) m_xHighLimit_V=10;
  if(m_xLowLimit_V<0.1) m_xLowLimit_V=0.1;
  if(m_yHighLimit_V>1) m_yHighLimit_V=1;
  if(m_yLowLimit_V<1e-320) m_yLowLimit_V=1e-320;

  printf("p-test constrains:\nZ  range:[0,1...%.3f] [10,0...%.3f] \nFC range:[0,1...%.3f] [10,0...%.3f]\nP  range:[0,0...%.2e] \n",
       m_xLowLimit_Z, m_xHighLimit_Z, m_xLowLimit_V, m_xHighLimit_V, m_yLowLimit_V);

  for(int i=0; i<3; i++)
    gsl_histogram_free(h[i]);

  return 0;
}



//A partir de los límites ya determinados y de la matriz m_ionsData_p se extraen los iones significativos
//las matrices m_ionsClusterZ_p y m_ionsClusterV_p mantienen los resultados (filas=iones; col=grupos focales)
//int IonsSelect::selection(int *grIndex, int nGr)
int IonsSelect::selection(int *gr, int nGr)
{
  int pCol, pGr, sGr;
  int grIndex[nGr];

  //determinamos los índices asociados a cada grupo pasado para acceso a la matriz
  for(int i=0; i<nGr; i++)
    for(int j=0; j<m_nPTestGroups; j++) //copia de grupos totales
      if(gr[i]==m_pTestGroups[j])
	grIndex[i]=j;

  for(int i=0; i<nGr; i++)  //Copia de los grupos a analizar
    m_analyzeGr_p[i]=grIndex[i];
  m_nAnalyzeGr=nGr;

  //Se extraen los iones significativos para cada grupo a partir del volcano test
  //Son significativos aquellos iones que están situados dentro de los límites antes determinados
  bool cumple=false;
  for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
    {
    cumple=true; //evita que un ion cumpla las condiciones límite en todos los grupos
    for(int pGri=0; pGri<nGr; pGri++) //para todos los grupos
      {
      pGr=grIndex[pGri];
      pCol=m_nPTestGroups*3*pGr; //Se mueve por la matrix principal donde se analizan m_nPTestGroups grupos
      m_ionsClusterV_p[ion][pGri]=2; //code 2
      for(int sGri=0; sGri<m_nAnalyzeGr; sGri++) //para todos los grupos
	{
	sGr=grIndex[sGri];
	if(sGr==pGr) continue;
	if(m_ionsData_p[ion][pCol+sGr*3+1]<0 || m_ionsData_p[ion][pCol+sGr*3+2]<0 ||
	  fabs(m_ionsData_p[ion][pCol+sGr*3+1])>=m_yLowLimit_V ||
	  fabs(m_ionsData_p[ion][pCol+sGr*3+2])>=m_xLowLimit_V)
	    {m_ionsClusterV_p[ion][pGri]=0; break;} //no cumple
	}

      if(m_ionsClusterV_p[ion][pGri]==0)
	{
	m_ionsClusterV_p[ion][pGri]=4; //code 4
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=grIndex[sGri];
	  if(sGr==pGr) continue;
	  if(m_ionsData_p[ion][pCol+sGr*3+1]<0 || m_ionsData_p[ion][pCol+sGr*3+2]<0 ||
	    fabs(m_ionsData_p[ion][pCol+sGr*3+1])>=m_yLowLimit_V ||
	    fabs(m_ionsData_p[ion][pCol+sGr*3+2])<=m_xHighLimit_V)
	      {m_ionsClusterV_p[ion][pGri]=0; break;} //no cumple
	  }
	}
      if(m_ionsClusterV_p[ion][pGri]==0) //falló en este grupo
	  cumple=false;
      }

      if(cumple==true && m_nAnalyzeGr>2) //cumple en todos los grupos primarios => no puede ser
	for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	  m_ionsClusterV_p[ion][pGri]=0;
    }

  //Se extraen los iones significativos para cada grupo a partir del Zero test
  //Son significativos aquellos iones que están situados dentro de los límites antes determinados
  for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
    {
    cumple=true; //evita que un ion cumpla las condiciones límite en todos los grupos
    for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
      {
      pGr=grIndex[pGri];
      pCol=m_nPTestGroups*3*pGr;
      m_ionsClusterZ_p[ion][pGri]=1; //code 1

      for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	{
	sGr=grIndex[sGri];
	if(sGr==pGr) continue;

	if(m_ionsData_p[ion][pCol+sGr*3+0]<0 ||  //porcentajes de ceros <10% en ambos grupos
	  fabs(m_ionsData_p[ion][pCol+sGr*3+0])>=m_xLowLimit_Z) //valor superior
	    {m_ionsClusterZ_p[ion][pGri]=0; break;} //no cumple
	}

      if(m_ionsClusterZ_p[ion][pGri]==0) //no cumple con code 1 y se prueba con code 3
	{
	m_ionsClusterZ_p[ion][pGri]=3; //code 3
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=grIndex[sGri];
	  if(sGr==pGr) continue;
	  if(m_ionsData_p[ion][pCol+sGr*3+0]<0 ||
	    fabs(m_ionsData_p[ion][pCol+sGr*3+0])<=m_xHighLimit_Z)
	      {m_ionsClusterZ_p[ion][pGri]=0; break;} //no cumple
	  }
	}
      if(m_ionsClusterZ_p[ion][pGri]==0) //falló en este grupo
	  cumple=false;
      }
    if(cumple==true && m_nAnalyzeGr>2) //cumple en todos los grupos primarios => no puede ser
      for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	m_ionsClusterZ_p[ion][pGri]=0;
    }
}

  //Ordena los iones pasados atendiendo a su posición. Los más alejados primero
  //Se consideran las variables Z, P y FC y se extrae el vector. Mayor módulo primero
  //Se considera el origen de los iones Z, V, ZV para conformar el vector en una, dos o tres dimensiones
  //En bestIonsIndex_p se retornan los índices ordenados segun 'direct' -> 'true'=de mayor a menor significación
  //En weight_p se retornan los pesos asociados (0..1) ('direct' = 'true => 1=más significativo)
int IonsSelect::getBestIonDistances(bool direct, int group, int *ions_p, int ionsCount, char *ionsCode_p, int *bestIonsIndex_p, float *weight_p)
  {
  int ion, *vLengthIndex_p;
  double Z, P, FC, minSVLength, tmp;
  double *vLength_p=NULL;
  int pGr, sGr;

  if(ionsCount==0)
    return 0;

  if((vLength_p=new double[ionsCount])==NULL)
    {printf("ERROR new() in getBestIon2()\n"); return -1;}

  for(int i=0; i<ionsCount; i++)
    {
    ion=ions_p[i];
    pGr=m_analyzeGr_p[group];
    int pCol=m_nPTestGroups*3*pGr;
    minSVLength=1e300;
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      {
      sGr=m_analyzeGr_p[sGri];
      if(pGr==sGr) continue;
      switch(ionsCode_p[i])
	{
	case 1: case 3: //source Z
	  Z= m_ionsData_p[ion][pCol+sGr*3+0];
	  P=0;
	  FC=0;
	  break;
	case 2: case 4: //source V
	  Z=0;
	  P= m_ionsData_p[ion][pCol+sGr*3+1];
	  FC=m_ionsData_p[ion][pCol+sGr*3+2];
	  break;
	case 14: case 23: //source Z+V
	  Z= m_ionsData_p[ion][pCol+sGr*3+0];
	  P= m_ionsData_p[ion][pCol+sGr*3+1];
	  FC=m_ionsData_p[ion][pCol+sGr*3+2];
	  break;
	default: //case 0
	  Z=0;
	  P=1.0;
	  FC=0;
	}
      if(Z>0)
	Z=log10(Z);
      else
	Z=0;
      P=1*(1-P);
      if(FC>0)
	FC=log10(FC);
      else
	FC=0;
      tmp=Z*Z+P*P+FC*FC; //longitud del vector al cuadrado
      if(tmp<minSVLength) minSVLength=tmp; //se considera el menor módulo (más desfavorable)
      }
    vLength_p[i]=minSVLength;
    }
    sortDown(vLength_p, bestIonsIndex_p, ionsCount, true); //ordenación de índices (mayor módulo primero)

  if(weight_p && vLength_p[bestIonsIndex_p[0]]>EPSILON_LD)
  for(int i=0; i<ionsCount; i++)
     weight_p[i]=vLength_p[i]/vLength_p[bestIonsIndex_p[0]]; //Normalización 0-1

  if(vLength_p) delete []vLength_p ;

  if(ionsCount && vLength_p[bestIonsIndex_p[0]]<EPSILON_LD)
    return -1;
  return 0;
  }

  //Se analizan los valores medios de cada ion con el propósito de ordenarlos
  //por significación. Es más significativo cuando el valor medio del grupo
  //pasado respecto al valor medio de los grupos considerados es elevado
  //En bestIonsIndex_p se retornan los índices ordenados segun 'direct' -> 'true'=de mayor a menor significación
  //En weight_p se retornan los pesos asociados (0..1) ('direct' = 'true => 1=más significativo)
int IonsSelect::getBestIonIntensity(bool direct, int group, int *grTest_p, int nGrTest, int *ions_p, int ionsCount, int *bestIonsIndex_p, float *weight_p)
  {
  int *ionsIndex_p=NULL, grSize;
  int gr, nPixels, px, pxCount;
  double *ionsMean_p=NULL,*ionsValue_p=NULL, *ionsValueGr_p=NULL;
  double ionMeanGr, ionMean;
//  int nPixels=m_nPixels;

      grSize=m_allGroups_p[group].size;
      if((ionsValue_p=new double[m_nPixels])==NULL)
	{printf("ERROR new() in makeFileFromFile()\n"); return -1;}
      if((ionsValueGr_p=new double[grSize])==NULL)
	{printf("ERROR new() in makeFileFromFile()\n"); return -1;}
      if((ionsMean_p=new double[ionsCount])==NULL)
	{printf("ERROR new() in makeFileFromFile()\n"); return -1;}

      for(int ion=0; ion<ionsCount; ion++)
	{
	for(int i=0; i<grSize; i++)
	  {
	  px=m_allGroups_p[group].set[i];
	  ionsValueGr_p[i]=m_allSpectrums_p[px][ions_p[ion]];
	  }
//	ionVarGr=gsl_stats_variance(ionsValue_p, 1, grSize);
	ionMeanGr=gsl_stats_mean(ionsValueGr_p, 1, grSize); //media del ion en el grupo focal

	//Se obtiene la media de todos los grupos considerados (que pueden no ser todos los de la muestra)
	pxCount=0;
	for(int gri=0; gri<nGrTest; gri++)
	  {
	  gr=grTest_p[gri]; //un grupos del conjunto de grupos considerados
	  nPixels=m_allGroups_p[gr].size; //tamaño del grupo
	  for(int i=0; i<nPixels; i++) //para todos los píxeles del grupo
	    {
	    px=m_allGroups_p[gr].set[i]; //un pixel del grupo
	    ionsValue_p[pxCount]=m_allSpectrums_p[px][ions_p[ion]]; //magnitud del ion para ese pixel
	    pxCount++; //contador
	    }
	  }
	//Se obtiene la media de todos los grupos de la muestra
	ionMean=gsl_stats_mean(ionsValue_p, 1, pxCount); //media del ion en todos los grupos considerados
	if(ionMean<EPSILON_LD && ionMeanGr<EPSILON_LD)
	  ionsMean_p[ion]=0;
	else
	  {
	  if(fabs(ionMean)<EPSILON_LD) ionMean=EPSILON_LD;
	  if(fabs(ionMeanGr)<EPSILON_LD) ionMeanGr=EPSILON_LD;
	  if(direct)
	    ionsMean_p[ion]=ionMeanGr/ionMean; //directo->valores elevados
	  else
	    ionsMean_p[ion]=ionMean/ionMeanGr; //inverso->valores pequeños
	  //Saturación: evita ordenaciones extrañas
	  if(ionsMean_p[ion]>100)ionsMean_p[ion]=100.0;
	  if(ionsMean_p[ion]<0.01)ionsMean_p[ion]=0.01;
	  }
	}
      sortDown(ionsMean_p, bestIonsIndex_p, ionsCount, true); //ordenación up->down

      if(weight_p)
	for(int ion=0; ion<ionsCount; ion++)
	  {
	  if(fabs(ionsMean_p[bestIonsIndex_p[0]])<EPSILON_LD)
	   {weight_p[ion]=0; continue;}
	  weight_p[ion]=ionsMean_p[ion]/ionsMean_p[bestIonsIndex_p[0]];//normalización
	  if(isnan(weight_p[ion])) weight_p[ion]=0;
	  }

      if(ionsValue_p) 	delete []ionsValue_p;
      if(ionsValueGr_p) delete []ionsValueGr_p;
      if(ionsMean_p) 	delete []ionsMean_p;
      return 0;
}

//Genera el fichero fileName con info extraida del test de selección de iones
//presenta los iones seleccionados segregados por grupos indicando si es Up o down regulated
//y con info particular de cada variable asociada (Z, P, FC)
//también se dan las restricciones del test
//Genera el fichero filename_b con otro formato
int IonsSelect::ionsSelectToFile(const char* fileName, int *grTest_p, int nGrTest, int ordination)
{
  FILE *fp, *fp2;
  int *ions_14=NULL, *ions_23=NULL;
  char *ions_14Code=NULL, *ions_23Code=NULL;
  int count_14, count_23;
  float *weight14_p=NULL;
  float *weight23_p=NULL;
  char source_0, source_1;
  int pGr, sGr, pCol;
  int **ionsGroups_p=NULL;

  if((ions_14=new int [m_pTestRows])==NULL) //mantiene índices a iones con códigos 1 y 4
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}
  if((ions_14Code=new char [m_pTestRows])==NULL) //códigos 1, 4 ó 14
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}
  if((ions_23=new int [m_pTestRows])==NULL) //mantiene índices a iones con códigos 2 y 3
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}
  if((ions_23Code=new char [m_pTestRows])==NULL) //códigos 2, 3 ó 23
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}

  //Mantiene info de ions-clusters sobre/sub expresados
  if((ionsGroups_p=new int *[m_pTestRows])==NULL)
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}
  for(int i=0; i<m_pTestRows; i++)
    if((ionsGroups_p[i]=new int [m_nAnalyzeGr])==NULL)
      {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}

  for(int i=0; i<m_pTestRows; i++)
    for(int j=0; j<m_nAnalyzeGr; j++)
      ionsGroups_p[i][j]=0;

  if((fp=fopen(fileName, "w"))==NULL)
    {perror("ERROR:"); return -1;}

  if((weight14_p=new float[m_pTestRows])==NULL) //memoria para pesos upregulated
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}
  if((weight23_p=new float[m_pTestRows])==NULL) //memoria para pesos downregulated
    {printf("ERROR new() in ionsSelectToFile()\n"); return -1;}

  //cabecera del fichero con las restricciones
  fprintf(fp, "Results for ions selection test\r\n");
  fprintf(fp, "  Constrains:\r\n");
  fprintf(fp, "\tprobability:%.2f %%\r\n", m_probability);
  fprintf(fp, "\t         Zero range left:(0,1...%.3f)   right:(10,0...%.3f)\r\n",
	 m_xLowLimit_Z, m_xHighLimit_Z);
  fprintf(fp, "\t  Fold change range left:(0,1...%.3f)   right:(10,0...%.3f)\r\n\t            p range high:0.0 low:%.2e \r\n",
	 m_xLowLimit_V, m_xHighLimit_V, m_yLowLimit_V);

  if(ordination==0)
    fprintf(fp, "  Ordination: distances\r\n");
  else
    fprintf(fp, "  Ordination: magnitudes\r\n");


  //para todos los grupos analizados
  for(int focalGroup=0; focalGroup<m_nAnalyzeGr; focalGroup++)
    {
    count_14=0; count_23=0;
    //Se extraen los iones upregulated (ions_14) y downregulated (ions_23)
    for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
      {
//if(ion==224)
//  printf("\n");
      if( m_ionsClusterZ_p[ion][focalGroup]==1 || m_ionsClusterV_p[ion][focalGroup]==4)
	{
	ions_14[count_14]=ion;
	if( m_ionsClusterZ_p[ion][focalGroup]==1 && m_ionsClusterV_p[ion][focalGroup]==4)
	  ions_14Code[count_14]=14;
	else if( m_ionsClusterZ_p[ion][focalGroup]==1)
	  ions_14Code[count_14]=1;
	else
	  ions_14Code[count_14]=4;
	ionsGroups_p[ion][focalGroup]=ions_14Code[count_14];
	count_14++;
	}
      else if( m_ionsClusterZ_p[ion][focalGroup]==3 || m_ionsClusterV_p[ion][focalGroup]==2)
	{
	ions_23[count_23]=ion;
	if( m_ionsClusterZ_p[ion][focalGroup]==3 && m_ionsClusterV_p[ion][focalGroup]==2)
	  ions_23Code[count_23]=23;
	else if( m_ionsClusterZ_p[ion][focalGroup]==3)
	  ions_23Code[count_23]=3;
	else
	  ions_23Code[count_23]=2;
	ionsGroups_p[ion][focalGroup]=ions_23Code[count_23];
	count_23++;
	}
      }

    //Se ordenan según los valores medios de sus magnitudes en el grupo frente al total de grupos
    int ionsIndex_14[count_14]; //códigos 1 y 4 -> sobre-expresados
    int ionsIndex_23[count_23]; //códigos 2 y 3 ->   sub-expresados
    bool direct=true;

    switch(ordination) //tipo de ordenación de las listas de iones
      {
      case 0: //por distancias
	direct=true;
	getBestIonDistances(direct, focalGroup, ions_14, count_14, ions_14Code, ionsIndex_14, weight14_p);
	direct=false;
	getBestIonDistances(direct, focalGroup, ions_23, count_23, ions_23Code, ionsIndex_23, weight23_p);
	break;
      case 1: //por intensidades
	direct=true;
	getBestIonIntensity(direct, focalGroup, grTest_p, nGrTest, ions_14, count_14, ionsIndex_14, weight14_p);
	direct=false;
	getBestIonIntensity(direct, focalGroup, grTest_p, nGrTest, ions_23, count_23, ionsIndex_23, weight23_p);
	break;
      }

  //info a fichero
    char cluster[50], tmpStr[10];
    cluster[0]=0;
    for(int i=0; i<m_nAnalyzeGr; i++)
      {
      sprintf(tmpStr, "%d,",m_pTestGroups[m_analyzeGr_p[i]]); //prepara la lista de grupos
      strcat(cluster, tmpStr);
      }
      cluster[strlen(cluster)-1]=0; //elimina la última coma
//    fprintf(fp, "\r\np-test results for cluster %d over clusters %d-%d\r\n", m_pTestGroups[focalGroup], m_pTestGroups[0], m_pTestGroups[m_nPTestGroups-1] );
    fprintf(fp, "\r\np-test results for cluster %d over clusters %s\r\n", m_pTestGroups[m_analyzeGr_p[focalGroup]], cluster);
    fprintf(fp, "\tUpregulated ions=%d; downregulated ions=%d\r\n", count_14, count_23);
    //caso de que existan iones upregulated en este grupo
    if(count_14>0)
      {
      fprintf(fp, "\r\n---Upregulated ions list--- \r\n");
      fprintf(fp,  "------------groups----------> ");

      for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	{
	pGr=m_analyzeGr_p[pGri];
	pCol=m_nPTestGroups*3*pGr;
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=m_analyzeGr_p[sGri];
	  if(pGr>=sGr) continue; //Para limitar la presentación se excluyen los complementarios aunque difieren en Z y FC
	  fprintf(fp, "-----------%d/%d----------- ", m_pTestGroups[pGr], m_pTestGroups[sGr]);
	  }
	}

      fprintf(fp, "\r\n------ions----- weight source ");
      for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	{
	pGr=m_analyzeGr_p[pGri];
	pCol=m_nPTestGroups*3*pGr;
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=m_analyzeGr_p[sGri];
	  if(pGr>=sGr) continue;//Para limitar la presentación se excluyen los complementarios aunque difieren en Z y FC
	  fprintf(fp, "---Z--- ----P---- ---FC-- ");
	  }
	}
      fprintf(fp, "\r\n");
      //iones (índice y masa) y origen (Zero/Volcano)
      for(int i=0; i<count_14; i++)
	{
	int ion=ions_14[ionsIndex_14[i]];
	fprintf(fp, "[%3d](%8.3f) ", ion, m_mzAxis_p[ion]);
	if(ions_14Code[ionsIndex_14[i]]==14)
	  {source_0='Z'; source_1='V';}
	else if(ions_14Code[ionsIndex_14[i]]==1)
	  {source_0=' '; source_1='Z';}
	else
	  {source_0=' '; source_1='V';}
	fprintf(fp, "%.3f   %c%c    ", weight14_p[ionsIndex_14[i]], source_0, source_1);

	//Info de variables Z, P y FC para este ion sobre todos los grupos (no se repite)
	double Z, P, FC;
	for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	  {
//if(ion==366 && pGri==3)
//  printf("\n");
	  pGr=m_analyzeGr_p[pGri];
	  pCol=m_nPTestGroups*3*pGr;
	  for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	    {
	    sGr=m_analyzeGr_p[sGri];
	    if(pGr>=sGr) continue;//Para limitar la presentación se excluyen los complementarios aunque difieren en Z y FC
	    Z= m_ionsData_p[ion][pCol+sGr*3+0];
	    P= m_ionsData_p[ion][pCol+sGr*3+1];
	    FC=m_ionsData_p[ion][pCol+sGr*3+2];
//	    fprintf(fp, "%+.4f %+.1e %+.4f ", Z, P, FC);
	    if(P>0 && P<1.0e-99)
	      fprintf(fp, "%+.4f %+.1e %+.4f ", Z, P, FC);
	    else
	      fprintf(fp, "%+.4f %+.1e  %+.4f ", Z, P, FC);
	    }
	  }
	fprintf(fp, "\r\n");
	}
      }

    //caso de que existan iones downregulated en este grupo
    if(count_23>0)
      {
      fprintf(fp, "\r\n");

      fprintf(fp, "\r\n---Downregulated ions list --- \r\n");
      fprintf(fp,  "------------groups----------> ");

      for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	{
	pGr=m_analyzeGr_p[pGri];
	pCol=m_nPTestGroups*3*pGr;
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=m_analyzeGr_p[sGri];
	  if(pGr>=sGr) continue;
	  fprintf(fp, "-----------%d/%d----------- ", m_pTestGroups[pGr], m_pTestGroups[sGr]);
	  }
	}

      fprintf(fp, "\r\n------ions----- weight source ");
      for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	{
	pGr=m_analyzeGr_p[pGri];
	pCol=m_nPTestGroups*3*pGr;
	for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	  {
	  sGr=m_analyzeGr_p[sGri];
	  if(pGr>=sGr) continue;//Para limitar la presentación se excluyen los complementarios aunque difieren en Z y FC
	  fprintf(fp, "---Z--- ----P---- ---FC-- ");
	  }
	}
      fprintf(fp, "\r\n");

      for(int i=0; i<count_23; i++)
	{
	int ion=ions_23[ionsIndex_23[i]];
	fprintf(fp, "[%3d](%8.3f) ", ion, m_mzAxis_p[ion]);
	if(ions_23Code[ionsIndex_23[i]]==23)
	  {source_0='Z'; source_1='V';}
	else if(ions_23Code[ionsIndex_23[i]]==3)
	  {source_0=' '; source_1='Z';}
	else
	  {source_0=' '; source_1='V';}
	fprintf(fp, "%.3f   %c%c    ", weight23_p[ionsIndex_23[i]], source_0, source_1);

	double Z, P, FC;
	for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
	  {
//if(ion==366 && pGri==3)
//  printf("\n");
	  pGr=m_analyzeGr_p[pGri];
	  pCol=m_nPTestGroups*3*pGr;
	  for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
	    {
	    sGr=m_analyzeGr_p[sGri];
	    if(pGr>=sGr) continue;
	    Z= m_ionsData_p[ion][pCol+sGr*3+0];
	    P= m_ionsData_p[ion][pCol+sGr*3+1];
	    FC=m_ionsData_p[ion][pCol+sGr*3+2];
//	    fprintf(fp, "%+.4f %+.1e %+.4f ", Z, P, FC);
	    if(P>0 && P<1.0e-99)
	      fprintf(fp, "%+.4f %+.1e %+.4f ", Z, P, FC);
	    else
	      fprintf(fp, "%+.4f %+.1e  %+.4f ", Z, P, FC);
	    }
	  }
	fprintf(fp, "\r\n");
	}
      fprintf(fp, "\r\n");
    }
  }

  fclose(fp);

  //lo que sigue genera el fichero filename_b (misma info en distinto formato)
  char str[300];
  int index;
  bool hit;
  for(index=strlen(fileName); index>0; index--)
    {
    if(fileName[index]=='.') { break;}
    }
  for(int i=0; i<index; i++)
      str[i]=fileName[i];
   str[index]=0;
  strcat(str, "_b.txt");
  if((fp2=fopen(str, "w"))==NULL)
    {perror("ERROR:"); return -1;}

  FILE *fp3;
   str[index]=0;
  strcat(str, "_c.txt");
  if((fp3=fopen(str, "w"))==NULL)
    {perror("ERROR:"); return -1;}

  //cabecera del fichero con las restricciones
  fprintf(fp2, "Results for ions selection test\r\n");
  fprintf(fp2, "  Constrains:\r\n");
  fprintf(fp2, "\tprobability:%.2f %%\r\n", m_probability);
  fprintf(fp2, "\t         Zero range left:(0,1...%.3f)   right:(10,0...%.3f)\r\n",
	 m_xLowLimit_Z, m_xHighLimit_Z);
  fprintf(fp2, "\t  Fold change range left:(0,1...%.3f)   right:(10,0...%.3f)\r\n\t            p range high:0.0 low:%.2e \r\n",
	 m_xLowLimit_V, m_xHighLimit_V, m_yLowLimit_V);

  fprintf(fp2, "\r\n  Codes list:\r\n\t- => not regulated\r\n\t1 =>   upregulated (Z)\r\n\t4 =>   upregulated (p & fc)\
  \r\n\t3 => downregulated (Z)\r\n\t2 => downregulated (p & fc)\r\n");
  fprintf(fp2, "\r\n------ions-----  ");
  for(int j=0; j<m_nAnalyzeGr; j++)
    fprintf(fp2, "-%d- ", m_pTestGroups[m_analyzeGr_p[j]]);
  fprintf(fp2, "\r\n");

  for(int i=0; i<m_pTestRows; i++)
    {
    hit=false;
    for(int j=0; j<m_nAnalyzeGr; j++)
      if(ionsGroups_p[i][j]!=0) //si existe el ion en la lista de seleccionados
	{hit=true; break;}
    if(!hit) continue;
    if(hit)
      {
      fprintf(fp2, "[%3d](%8.3f) ", i, m_mzAxis_p[i]);
      fprintf(fp3, "%9.4f\r\n", m_mzAxis_p[i]);
      for(int j=0; j<m_nAnalyzeGr; j++)
	{
	if(ionsGroups_p[i][j]!=0)
	  fprintf(fp2, " %2d ", ionsGroups_p[i][j]);
	else
	  fprintf(fp2," -- ");
	}
      }
    fprintf(fp2, "\r\n");
    }
  fclose(fp2);
  fclose(fp3);

  if(ions_14) delete [] ions_14;
  if(ions_23) delete [] ions_23;
  if(ions_14Code) delete [] ions_14Code;
  if(ions_23Code) delete [] ions_23Code;

  if(weight14_p) delete [] weight14_p;
  if(weight23_p) delete [] weight23_p;

  for(int i=0; i<m_pTestRows; i++)
    if(ionsGroups_p[i]) delete []ionsGroups_p[i];
  if(ionsGroups_p) delete []ionsGroups_p;
}

//Retorna la cantidad de iones sobre-expresados en el grupo dado respecto al resto de grupos
int IonsSelect::getUpRegulatedNumber(int focalGroup)
  {
  int count_14=0;
  for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
     if( m_ionsClusterZ_p[ion][focalGroup]==1 || m_ionsClusterV_p[ion][focalGroup]==4)
	count_14++;
  return count_14;
  }

//Retorna la lista de iones sobre-expresados ordenados según el criterio dado
//retorna la cantidad de iones en la lista (debe ser <= nIons)
//los punteros ionsCode y weight pueden ser nulos si no interesa su contenido
int IonsSelect::getUpRegulated(int focalGroup, int *grTest_p, int nGrTest, char ordination, int *ionList, char *ionsCode, float *weight, int nIons)
  {
  int count_14=0;
  int ions_14[nIons]; //códigos 1 y 4 -> sobre-expresados
  char ions_14Code[nIons];
  float weight14_p[nIons];

    for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
      {
      if(count_14>nIons) return -1;
      if( m_ionsClusterZ_p[ion][focalGroup]==1 || m_ionsClusterV_p[ion][focalGroup]==4)
	{
	ions_14[count_14]=ion;
	if( m_ionsClusterZ_p[ion][focalGroup]==1 && m_ionsClusterV_p[ion][focalGroup]==4)
	  ions_14Code[count_14]=14;
	else if( m_ionsClusterZ_p[ion][focalGroup]==1)
	  ions_14Code[count_14]=1;
	else
	  ions_14Code[count_14]=4;
	count_14++;
	}
      }

    //Se ordenan según los valores medios de sus magnitudes en el grupo frente al total de grupos
    int ionsIndex_14[nIons]; //códigos 1 y 4 -> sobre-expresados
    bool direct=true;

    switch(ordination) //tipo de ordenación de las listas de iones
      {
      case 0: //por distancias
	getBestIonDistances(direct, focalGroup, ions_14, count_14, ions_14Code, ionsIndex_14, weight14_p);
	break;
      case 1: //por intensidades
	getBestIonIntensity(direct, focalGroup, grTest_p, nGrTest, ions_14, count_14, ionsIndex_14, weight14_p);
	break;
      }
  for(int i=0; i<count_14; i++)
    ionList[i]=ions_14[ionsIndex_14[i]];
  if(ionsCode)
    for(int i=0; i<count_14; i++)
      ionsCode[i]=ions_14Code[ionsIndex_14[i]];
  if(weight)
    for(int i=0; i<count_14; i++)
      weight[i]=weight14_p[ionsIndex_14[i]];

  return count_14;
  }

//Retorna la cantidad de iones sub-expresados en el grupo dado respecto al resto de grupos
int IonsSelect::getDownRegulatedNumber(int focalGroup)
  {
  int count_23=0;
  for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
     if( m_ionsClusterZ_p[ion][focalGroup]==3 || m_ionsClusterV_p[ion][focalGroup]==2)
	count_23++;
  return count_23;
  }

//Retorna la lista de iones sub-expresados ordenados según el criterio dado
//retorna la cantidad de iones en la lista (debe ser <= nIons)
//los punteros ionsCode y weight pueden ser nulos si no interesa su contenido
int IonsSelect::getDownRegulated(int focalGroup, int *grTest_p, int nGrTest, char ordination, int *ionList, char *ionsCode, float *weight, int nIons)
  {
  int count_23=0;
  int ions_23[nIons]; //códigos 1 y 4 -> sobre-expresados
  char ions_23Code[nIons];
  float weight23_p[nIons];

    for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
      {
      if(count_23>nIons) return -1;
      if( m_ionsClusterZ_p[ion][focalGroup]==3 || m_ionsClusterV_p[ion][focalGroup]==2)
	{
	ions_23[count_23]=ion;
	if( m_ionsClusterZ_p[ion][focalGroup]==3 && m_ionsClusterV_p[ion][focalGroup]==2)
	  ions_23Code[count_23]=23;
	else if( m_ionsClusterZ_p[ion][focalGroup]==3)
	  ions_23Code[count_23]=3;
	else
	  ions_23Code[count_23]=2;
	count_23++;
	}
      }

    //Se ordenan según los valores medios de sus magnitudes en el grupo frente al total de grupos
    int ionsIndex_23[count_23]; //códigos 2 y 3 ->   sub-expresados
    bool direct=false;

    switch(ordination) //tipo de ordenación de las listas de iones
      {
      case 0: //por distancias
	getBestIonDistances(direct, focalGroup, ions_23, count_23, ions_23Code, ionsIndex_23, weight23_p);
	break;
      case 1: //por intensidades
	getBestIonIntensity(direct, focalGroup, grTest_p, nGrTest, ions_23, count_23, ionsIndex_23, weight23_p);
	break;
      }
  for(int i=0; i<count_23; i++)
    ionList[i]=ions_23[ionsIndex_23[i]];
  if(ionsCode)
    for(int i=0; i<count_23; i++)
      ionsCode[i]=ions_23Code[ionsIndex_23[i]];
  if(weight)
    for(int i=0; i<count_23; i++)
      weight[i]=weight23_p[ionsIndex_23[i]];
  return count_23;
  }

//Retorna los índices a un array de doubles de forma que quede ordenado de manera descendente
//'absolutas' es 'true' si no se debe considerar el signo en los valores del array
int IonsSelect::sortDown(double* bufferIn, int *sort, int size, bool absolutas)
{
int k=0, j;// end=m_m_distributionSize;
double  mayor;
double  list[size];
#ifndef DBL_MAX
#define DBL_MAX 1e-300
#endif

if(size==0) return 0;
for(j=0; j<size; j++)
{
  if(absolutas)
    list[j]=fabs(bufferIn[j]); //copia ya que se producen cambios en el contenido
  else
    list[j]=bufferIn[j]; //copia ya que se producen cambios en el contenido
}

if(size==1) sort[0]=0;//Si sólo hay uno ya está ordenado. Si no hay ninguno -> no hacer nada
do
	{
	mayor=-DBL_MAX; //minor<=minValue.
	for(j=0; j<size; j++)
		{
		if((list[j])>mayor)
			{
			mayor=list[j];
			sort[k]=j;
			}
		}

	list[sort[k]]=-DBL_MAX;
	k++;
	}
while(k<size);
return 0;
}

  //Se presentan tres ficheros con la info de m_ionsData_p
  ////////////////////////////////////////////////////////
  //ZMeasToFile->variable Z
int IonsSelect::ZMeasToFile(const char *fileName)
{
  FILE* fp;
  int pGr, sGr;

  fp=fopen(fileName, "w");
  fprintf(fp, "      ");
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      fprintf(fp, "  %d/%d  ", m_pTestGroups[m_analyzeGr_p[pGri]], m_pTestGroups[m_analyzeGr_p[sGri]]);
  fprintf(fp, "\r\n");
 for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
  {
  fprintf(fp, "[%3d] ",ion);
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
    {
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      {
      pGr=m_analyzeGr_p[pGri];
      sGr=m_analyzeGr_p[sGri];
      if(m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+0]<0)
	fprintf(fp, "%.3f ",  m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+0]);
      else
	fprintf(fp, "%6.3f ", m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+0]);
      }
    }
  fprintf(fp, "\r\n");
  }
  fclose(fp);
}

  //PMeasToFile->variable P
int IonsSelect::PMeasToFile(const char *fileName)
{
  FILE* fp;
  int pGr, sGr;

  fp=fopen(fileName, "w");
  fprintf(fp, "      ");
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      fprintf(fp, "   %d/%d    ", m_pTestGroups[m_analyzeGr_p[pGri]], m_pTestGroups[m_analyzeGr_p[sGri]]);
  fprintf(fp, "\r\n");
 for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
  {
  fprintf(fp, "[%3d] ",ion);
 double tmpData;
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
    {
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      {
      pGr=m_analyzeGr_p[pGri];
      sGr=m_analyzeGr_p[sGri];
      tmpData=m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+1];
      if(tmpData>0 && tmpData<1.0e-99)
	fprintf(fp, "%+.1e ", tmpData);
      else
	fprintf(fp, "%+.1e  ", tmpData);
      }
    }
  fprintf(fp, "\r\n");
  }
  fclose(fp);
}

  //FCMeasToFile->variable FC
int IonsSelect::FCMeasToFile(const char *fileName)
{
  FILE* fp;
  int pGr, sGr;

  fp=fopen(fileName, "w");
  fprintf(fp, "      ");
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
//    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
//      fprintf(fp, "   %d/%d    ", m_pTestGroups[m_analyzeGr_p[pGri]], m_pTestGroups[m_analyzeGr_p[sGri]]);
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      fprintf(fp, "  %d/%d  ", m_pTestGroups[m_analyzeGr_p[pGri]], m_pTestGroups[m_analyzeGr_p[sGri]]);
  fprintf(fp, "\r\n");
 for(int ion=0; ion<m_pTestRows; ion++) //para todos los iones
  {
  fprintf(fp, "[%3d] ",ion);
  for(int pGri=0; pGri<m_nAnalyzeGr; pGri++) //para todos los grupos
    {
    for(int sGri=0; sGri<m_nAnalyzeGr; sGri++)//para todos los grupos
      {
      pGr=m_analyzeGr_p[pGri];
      sGr=m_analyzeGr_p[sGri];
      if(m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+2]<0)
	fprintf(fp, "%+.3f ", m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+2]);
      else
	fprintf(fp, "%.4f ",  m_ionsData_p[ion][m_nPTestGroups*pGr*3+sGr*3+2]);
    }
    }
  fprintf(fp, "\r\n");
  }
  fclose(fp);
}

//Retorna 'z-value' desde la matriz m_ionsData_p
double IonsSelect::getZvalue(int ion, int pGr, int sGr)
{
  int pGri=-1, sGri=-1;

  if(ion<0 || ion>=m_pTestRows) //control de límites para ion
    return -20;

  for(int i=0; i<m_nPTestGroups; i++) //captura de índices de grupos
    {
    if(pGr==m_pTestGroups[i])
      pGri=i;
    if(sGr==m_pTestGroups[i])
      sGri=i;
    }
  if(pGri==-1 || sGri==-1) //control de límites para grupos
    return -20;

  return m_ionsData_p[ion][3*(m_nPTestGroups*pGri+sGri)+0];
}

//Retorna 'p-value' desde la matriz m_ionsData_p
double IonsSelect::getPvalue(int ion, int pGr, int sGr)
{
  int pGri=-1, sGri=-1;

  if(ion<0 || ion>=m_pTestRows) //control de límites para ion
    return -20;

  for(int i=0; i<m_nPTestGroups; i++) //captura de índices de grupos
    {
    if(pGr==m_pTestGroups[i])
      pGri=i;
    if(sGr==m_pTestGroups[i])
      sGri=i;
    }
  if(pGri==-1 || sGri==-1) //control de límites para grupos
    return -20;

  return m_ionsData_p[ion][3*(m_nPTestGroups*pGri+sGri)+1];
}

//Retorna 'fold_change-value' desde la matriz m_ionsData_p
double IonsSelect::getFCvalue(int ion, int pGr, int sGr)
{
  int pGri=-1, sGri=-1;

  if(ion<0 || ion>=m_pTestRows) //control de límites para ion
    return -20;

  for(int i=0; i<m_nPTestGroups; i++) //captura de índices de grupos
    {
    if(pGr==m_pTestGroups[i])
      pGri=i;
    if(sGr==m_pTestGroups[i])
      sGri=i;
    }
  if(pGri==-1 || sGri==-1) //control de límites para grupos
    return -20;

  return m_ionsData_p[ion][3*(m_nPTestGroups*pGri+sGri)+2];
}

//Retorna la lista con los grupos que conforman la combinación dada y su tamaño
//Retorna -1 si 'combination' está fuera de límites
//El tamaño de geList debe ser mayor o igual que el número de grupos máximo a evaluar (pow(2, m_nPTestGroups))
  //IDEA:
  //	 Los '1' en 'combination' indican los grupos de ' m_pTestGroups' que intervienen en cada combinación
  //por ejemplo, si m_pTestGroups={1, 3, 5}, 2^3=8 combinaciones. A partir de todas las combi binarias:
  //0:000->none
  //1:001->1
  //2:010->3
  //3:011->1, 3
  //4:100->5
  //5:101->1, 5
  //6:110->3, 5
  //7:111->1, 3, 5
int IonsSelect::getGroupsFromCombination(int combination, int *grList)
  {
  if(combination<0 || combination>=pow(2, m_nPTestGroups))
    return -1;

  int mask=1;
  int nGrTest=0;
  for(int i=0; i<m_nPTestGroups; i++)
    {
    if(mask & combination) //unos
      grList[nGrTest++]=m_pTestGroups[i];
    mask<<=1;
    }
  return nGrTest;
  }
