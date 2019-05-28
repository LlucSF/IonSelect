#ifndef IONS_SELECT_H
#define IONS_SELECT_H

#include "types.h"
#include "ionsTreat.h"
#include <string.h>
#include <RcppGSL.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

class IonsSelect
{
public:
  //constructor
  //recibe información sobre las muestras a considerar (lsp_p & nSamples), sobre los grupos a analizar (pTestGroups & nPTestGroups)
  //sobre los grupos de la etapa de segmentación (allGroups_p) y sobre las matriz de magnitudes pixels-ions (allSpectrums_p)
  IonsSelect(LOAD_SP *lsp, int nSamples, int *pTestGroups, int nPTestsGroups, GROUP* allGroups_p, float** allSpectrums_p);

  //Destruye memoria reservada
  ~IonsSelect();
  

  void setThreshold(double threshold);
  
  //Extrae las medidas de z, P y FC
  //Genera la matriz m_ionsData_p a partir de los datos en las matrices m_allSpectrums_p y de grupos m_allGroups_p
  //m_ionsData_p: filas=iones; columnas=variables Z, P, FC para cada combinación de dos grupos de m_allGroups_p (0/0, 0/1, 0/2..0/n; 1/0, 1/1..)
  //-1 en m_ionsData_p indica que el dato no es válido  
  int getMeasures();
  
  //A partir de los límites ya determinados y de la matriz m_ionsData_p se extraen los iones significativos
  //las matrices m_ionsClusterZ_p y m_ionsClusterV_p mantienen los resultados (filas=iones; col=grupos focales)
  //int IonsSelect::selection(int *grIndex, int nGr)  
  int selection(int *grIndex, int nGrIndex);
  
  //Genera el fichero fileName con info extraida del test de selección de iones
  //presenta los iones seleccionados segregados por grupos indicando si es Up o down regulated
  //y con info particular de cada variable asociada (Z, P, FC)
  //también se dan las restricciones del test
  //Genera el fichero filename_b con otro formato
  int ionsSelectToFile(const char* fileName, int *grTest, int nGrTest, int ordination);
  
  
  //Genera un fichero en el propio directorio del proyecto con los valores de 'zeros'
  int ZMeasToFile(const char *fileName);
  
  //Genera un fichero en el propio directorio del proyecto con los valores de 'p'
  int PMeasToFile(const char *fileName);
  
  //Genera un fichero en el propio directorio del proyecto con los valores de 'fold change'
  int FCMeasToFile(const char *fileName);
  
  //Retorna la cantidad de iones upRegulated en todos los grupos
  int getAllUpRegulatedNumber();
  
  //Retorna la cantidad de iones downRegulated en todos los grupos
  int getAllDownRegulatedNumber();
  
  //Retorna la cantidad de iones sobre-expresados en el grupo dado respecto al resto de grupos
  int getUpRegulatedNumber(int focalGroup);
  
  //Retorna la lista de iones sobre-expresados ordenados según el criterio dado
  //retorna la cantidad de iones en la lista (debe ser <= nIons)
  //los punteros ionsCode y weight pueden ser nulos si no interesa su contenido
  int getUpRegulated(int focalGroup, int *grTest, int nGrTest, char ordination, int *ionList, char *ionsCode, float *weight, int nIons);
  
  //Retorna la cantidad de iones sub-expresados en el grupo dado respecto al resto de grupos
  int getDownRegulatedNumber(int focalGroup);
  
  //Retorna la lista de iones sub-expresados ordenados según el criterio dado
  //retorna la cantidad de iones en la lista (debe ser <= nIons)
  //los punteros ionsCode y weight pueden ser nulos si no interesa su contenido
  int getDownRegulated(int focalGroup, int *grTest, int nGrTest, char ordination, int *ionList, char *ionsCode, float *weight, int nIons);
  
  //Se establecen los percentiles de los que se derivarán los limites en Z, P y FC para la selección de iones
  //Recibe un percentil único
  //Como FC y P se consideran simultáneamente, la probabilida conjunta es el producto de ambas y se ajusta
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int setBoundaries(double percentil);
  
  //Se establecen los percentiles de las que se derivarán los limites en Z, P y FC para la selección de iones
  //Recibe los percentiles asociada a cada límite
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int setBoundaries(double zPercentil, double pPercentil, double fcPercentil);
  
//Retorna 'z-value' correspondiente a dos cluster
  double getZvalue(int ion, int pGr, int sGr);
  
//Retorna 'P-value' correspondiente a dos cluster
  double getPvalue(int ion, int pGr, int sGr);
  
//Retorna 'FC-value' correspondiente a dos cluster
  double getFCvalue(int ion, int pGr, int sGr);
  
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
  int getGroupsFromCombination(int combination, int *grList);

  //Se analizan los valores medios de cada ion con el propósito de ordenarlos
  //por significación. Es más significativo cuando el valor medio del grupo 
  //pasado respecto al valor medio de los grupos considerados es elevado
  //En bestIonsIndex_p se retornan los índices ordenados segun 'direct' -> 'true'=de mayor a menor significación
  //En weight_p se retornan los pesos asociados (0..1) ('direct' = 'true => 1=más significativo)
  int getBestIonIntensity(bool direct, int group, int *grTest, int nGrTest,
                          int *ions_p, int ionsCount, int *bestIonsIndex_p,
                          float *weight_p);
  
private:
  //Se determinan los limites en Z, P y FC para la selección de iones
  //Recibe los percentiles para delimitar los valores de corte
  //Se generan tres histogramas con los datos de todas las combinaciones evaluadas.
  //Se integra cada histograma desde ambos extremos hasta que se alcanza el valor de probabilidad dado
  //Luego se determina el valor en X (magnitud) correspondiente. Sirve de valor de corte
  //La idea es hacer que cada variable considere un rango de valores cuya probabilidad sea menor o igual a una dada.
  //Se hace uso de histogramas ante la falta de normalidad en las distribuciones
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int setZPFcBoundaries();
  
  //Ordena los iones pasados atendiendo a su posición. Los más alejados primero
  //Se consideran las variables Z, P y FC y se extrae el vector. Mayor módulo primero
  //Se considera el origen de los iones Z, V, ZV para conformar el vector en una, dos o tres dimensiones
  //En bestIonsIndex_p se retornan los índices ordenados segun 'direct' -> 'true'=de mayor a menor significación
  //En weight_p se retornan los pesos asociados (0..1) ('direct' = 'true => 1=más significativo)
  int getBestIonDistances(bool direct, int group, int *ions_p, int ionsCount, char *ionsCode_p, int *bestIonsIndex_p, float *weight_p);
  

  //Retorna los índices a un array de doubles de forma que quede ordenado de manera descendente
  //'absolutas' es 'true' si no se debe considerar el signo en los valores del array
  int sortDown(double* bufferIn, int *sort, int size, bool absolutas);
  
  int m_pTestRows, //filas de la matriz de test (iones)
    m_pTestCols, //columnas de la matriz de test (combinaciones de grupos no repetitivos)
		m_nPixels,	//píxeles totales en la matriz de carga de espectros
		m_pTestGroups[50], //lista de grupos (clusters) a evaluar
		m_nPTestGroups; //número de grupos
  GROUP		*m_allGroups_p; //info de clustering
  float		**m_allSpectrums_p; //matriz principal (o de carga de espectros)

  double 	m_xLowLimit_Z,  //Límites asociados a cada variable para selección de iones
		m_xHighLimit_Z,
 		m_xLowLimit_V, 
		m_xHighLimit_V, 
		m_yLowLimit_V, 
		m_yHighLimit_V;
  float		*m_mzAxis_p; //puntero a array de masas
  

  double    m_zProb,    //probablidad acumulada (percentil/100) para los umbrales de las medidas Z, P y FC
            m_pProb,
            m_fcProb;
public:
  int 		m_analyzeGr_p[50], //mantiene los índices a los grupos a analizar(ej. si grupos=3,0,4->índices=1,0,2)
                     m_nAnalyzeGr; //lamaño de la lista
  double 	**m_ionsData_p; //matriz que alberga las variables Z, P, FC sobre todas las combinaciones
  int		**m_ionsClusterZ_p, //matriz con iones seleccionados, variable Z(filas=iones; col=grupos focales)
  **m_ionsClusterV_p; //matriz con iones seleccionados, variables P/FC(filas=iones; col=grupos focales)
  bool		m_error; //'true' si error en el constructor
  double 	m_threshold; //nivel máximo que delimita la zona de 'zeros'  
};

#endif
