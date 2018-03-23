
#ifndef CUSTERING_IONTREAT_H
#define CUSTERING_IONTREAT_H

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include "types.h"

//Tratamiento estadístico de iones
class IonsTreat
{
public:
  //contructor A
  IonsTreat();

  //contructor B
  IonsTreat(int matrixRows, int matrixCols, float**matrixIn);

  //Destructor
  ~IonsTreat();

  double getStudentP(float *dataA, float* dataB, int sizeA, int sizeB);
  double getStudentP(double *dataA, double* dataB, int sizeA, int sizeB);
  double getHistoP(double *ionsA_p, double *ionsB_p, int pSize, int sSize, double ionMin, double ionMax, int nBins);
  float** getPMatrix();
  float** getFoldChangeMatrix();

  //retorna el fold change (cociente de medianas) de dos distribuciones con elementos del tipo float
  //Los valores nulos no contribuyen
  double getFoldChange(float *dataA_f, float* dataB_f, int sizeA, int sizeB);

  //retorna el fold change (cociente de medianas) de dos distribuciones con elementos del tipo double
  //Los valores nulos no contribuyen
  double getFoldChange(double *dataA, double* dataB, int sizeA, int sizeB);

  //Test U de Mann-Whitney
  //Recibe dos array con muestras y determina si esas muestras provienen de la misma distribución
  //Retorna valores entre 0 y 1: 1 si la distribución original es la misma para ambos arrays muestras
  //Retorna valor negatrivo si el resultado no se deriva de datos consistentes
  //'sigma' determina la desviación estándar a utilizar
  //en 'zeroRate' se retorna el cociente entre los zeros de ambos arrays. Limitado entre 0.1 y 9.999.
  //'zeroRate' vale negatrivo si el cociente no se deriva de datos consistentes
  double getMannWhitneyUTest(double *dataA, double* dataB, int sizeA, int sizeB, double *zeroRate, double sigma);

private:
  //Test Z de Mann-Whitney
  //Recibe dos array con muestras y determina si esas muestras provienen de la misma distribución
  //Retorna un valor nulo si si la distribución original es la misma para ambos arrays muestras
  //La varioable aleatoria retornada está normalizada según gausiana de media cero y sigma 1
  double getMannWhitneyZ(double *dataA, double* dataB, int sizeA, int sizeB);

  //Retorna la cantiad de elementos con valor nulo en data
  //Si zerosIndex!=NULL, establece los índices de los elementos nulos en data
  int getZeros(double *data, int size, int *zerosIndex);

  //Retorna la cantiad de elementos con valor no nulo en data
  //Si notZeros!=NULL, copia los elementos no nulos de data
  int getNotZeros(double *data, int size, double *notZeros);



private:
  int 	m_matrixRows,
	m_matrixCols;
  float	**m_matrixIn,
	**m_pMatrix,
	**m_fcMatrix;
};

#endif
