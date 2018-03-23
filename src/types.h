#ifndef MALDI_TYPES_H
#define MALDI_TYPES_H

#include <sys/types.h>

#define EPSILON_F 	1e-6 
#define EPSILON_D 	1e-15 
#define EPSILON_LD 	1e-30
#define DOUBLE_MIN 		+1.0e-300
#define DOUBLE_MAX 		+1.0e+300
#define MAX_IMAGES 	50
#define PCA_FILE_COLUMS	20   
#define NEPER	 2.7182818 
#define MZ_AXIS_FROM_MATRIX	1
#define MZ_AXIS_FROM_FILE	2
#define MZ_AXIS_FROM_PARAMS	3
#define PIXELS_IONS		0
#define IONS_PIXELS		1
#define DIRECT			1
#define INVERSE			0
#define MAX_THREADS 		20
#define MATRIX_FILE_FORMAT	0
#define PERE_FILE_FORMAT	1
#define IMZML_FILE_FORMAT	2
#define GAUSS_PDF		0
#define MAXWELL_PDF 		1
#define MEAN_MAGNITUDE		0
#define MP_MAGNITUDE 		1
   

typedef struct
  {
  double
    probability,
    xLowLimit_Z,  
    xHighLimit_Z,  
    xLowLimit_V,  
    xHighLimit_V,  
    yLowLimit_V,  
    yHighLimit_V;  
  }CONSTRAINTS;

typedef struct {
  int	itemType, samePosSize, index;
	
  double x, y, width, height;
}INFO_MOUSE;

typedef struct {
double x, y;
}POINT_2D;

typedef struct {
double ionMass, x, y, xb, yb;
int itemType;
}ION_PLOT;

typedef struct {
  float r, g, b; 
}RGB_COLOR; 

typedef struct {
  float h, s, v;
}HSV_COLOR;

typedef struct {
  HSV_COLOR **image;
  int width,
      high;
}HSV_IMAGE;


typedef struct {
  float mzDelta,	//segmento del eje mz representado por cada entrada en el array mz_p
	*mz_p;		//array de masas válidas (con iones)
  int 	mzSize;		//tamaño del array
}ION_MZ;

typedef struct 
{
  float mean,
	sigma,
	magnitude;
}GAUSS;

typedef struct xMAXWELL_ION{
  float 	mzIniIndex,
		mzMax,
		maxY,
		xMaxRel,
		offset, //offset respecto a indexIni de rawIon
		stdDev,
		factor;
  struct xMAXWELL_ION	*nextIon_p; 
//			*previousIon_p; 
}MAXWELL_ION;

typedef struct{
  int 		indexIni,
		indexEnd,
		indexMax;
  float		maxY;
  float 	*magnitude_p,
		*mz_p;
  int 		magnitudeSize;
}RAW_ION;

typedef struct{
  int 		minMassPoints;
  double	
		minIniValue,
		minPeakValue;  
}ION_PARAM;

typedef struct xGROUP{
    int 	size;
    int 	*set;
    struct xGROUP 	*group;
  }GROUP;
   
typedef struct{ 
double maxValue,
      TIC, 
      Zr, Lr, Cr, Rr, Wr, Dr, Ed,  
      iLa, iLb, 
      iCa, iCb, 
      iRa, iRb;
bool  valid;
float	factor;
}DATA_VIEW;
  
typedef struct
  {
    double 	ML1, ML2, ML3; 	//Coeficientes de masas ajustados según el experimento 
    double 	delay, 		//Retardo temporal antes de iniciar las medidas
		dw;		//Tamaño de la ventana temporal de análisis (resolución para TOF)
    int 	td;		//Número de puntos de medida.
  }PARAM;

typedef struct
  {
    bool	
	      nil,
	      cte,
	      gauss,
	      rayleigh,
	      lineUp,
	      lineDown,
	      noise;
    float 
	      nilPosition,
	      ctePosition,
	      gaussPosition,
	      rayleighPosition,
	      lineUpPosition,
	      lineDownPosition,
	      noisePosition,
	      nilRange,
	      cteRange,
	      gaussRange,
	      rayleighRange,
	      lineUpRange,
	      lineDownRange,
	      noiseRange;
    float
	      cteValue_a,
	      gaussValue_a,
	      gaussValue_b,
	      rayleighValue_a,
	      rayleighValue_b,
	      lineUpValue_a,
	      lineDownValue_a,
	      noiseValue_a,
	      noiseValue_b;    
  }
  SP_PARAMS;
  
  typedef struct{
   int 	*set, 
	size;
  float *iLevel;
  }PIXEL_LEVEL;
  
  typedef struct
  {
    int x, y;
//    u_int32_t *spectrum;
    float *spectrum;
  }STR_SPECTRUM;
  
  typedef struct
  {
   int 	x, 
	y,
	nIons;
  float	area;		//area acumulada de la cadena de iones
   MAXWELL_ION  *ionsChain_p;
  }IONS_CHAIN;

  typedef struct{ //Masas representativas
    float	maxValue,
		var;
   int 		maxIndex; 
  }MASS_MAX;
  
  typedef struct{
    float x, y;
  }XY_f;
  
  typedef struct{
    int x, y;
  }PIXEL_XY;
    
  typedef struct{
    int mzLow, mzHigh;
  }MZ_VALID;
  
  typedef struct{
    int x1, y1,
	x2, y2,
	xRange, yRange;
  }REC;
  
typedef struct
{  
  int		
	xMinFisic,
	xMaxFisic,
	yMinFisic,
	yMaxFisic,
	xMax, 
	yMax,
	xRange,
	yRange,
	nPixels;
}SPECIMEN_COOR;  

typedef struct{
  char 	*dirName;		//nombre del directorio que contiene los espectros
  REC	rec;			//región de interés (ROI) en píxeles
  float partialMzLow, 		//masa del primer elemento del eje mz cargado
	partialMzHigh, 		//masa del último elemento del eje mz cargado
	totalMzLow, 		//masa del primer elemento del eje mz completo
	totalMzHigh, 		//masa del último elemento del eje mz completo
	cutLevel, 		//nivel de corte de ruido (valores inferiores se anulan)
	minTIC;			//pixeles con valores inferiores se eliminan
  int 	ionMinMassPoints,	//mínima cantidad de puntos másicos para que se considere como ion gapSize,
	nPixels,		//numero de pixeles del paquete
	firstPixelIndex, 	//índice al primer pixel de carga del paquete
	mzWidth; 		//numero de columnas en la matriz principal
  float* mzAxis_p;		//puntero al eje de masas
  int	partialMzAxisSize,	//tamaño de la parte de mz a cargar
	totalMzAxisSize;	//tamaño total del eje mz original
  int 	partialMzLowIndex,	//índice al primer elemento cargado de mz 
	partialMzHighIndex;	//índice al último elemento cargado de mz 
  int	matrixRows,
	matrixCols;
  SPECIMEN_COOR *spc_p;		//puntero a coordenadas del espécimen	
}LOAD_SP; 
   
typedef struct
{
  bool		newProject;
  bool		valid;		//válido cuando se ha acabado de crear un nuevo proyecto o se abre
  char 		*dir,  		//directorio donde se aloja el proyecto
		*name; 		//nombre del proyecto (nombre.prj
  int 		nSamples; 	//número de muestras
  char 		*sampleName[50];//pointer to absolute path of de samples
  REC  		sampleROI[50];	//ROI de cada muestra
  PIXEL_XY	
		sampleMaxXY[50];//dimensión de la muestra en píxeles
  REC		sampleMaxFisicROI[50]; //ROI fïsico máximo (puede que el pixel menor no sea cero)
  double	mzLow,		//valor mínimo del eje de masas cargado
		mzHigh;		//valor máximo del eje de masas cargado
  int		normalizationType, //tipo de normalización
		normalizationValue;//posible argumento para la normalización
  double	cutLevel,	//nivel de corte sobre la magnitud de cada masa para eliminación de ruido (ruido <= cutLevel)
		minTIC;		//Valor mínimo para que se conserve el pixel
  bool		matrixType,	//(filas=pixels; columnas=mz) o traspuesta
		ionsSelect,	//se usan iones particulares
		excludeAu;	//se excluyen los iones de oro
  char 		*ionsSelectFile;//nombre del fichero usado en la selección de iones
  double	deltaColor,	//colores de la imágenes
		H_range,
		H_low,
		H_factor,
		V_factor,
		S_level;
		
 //valores para CLEVEL
 double		MCSlevel;	//Clevel -> nivel de similaridad 
 int		minPxGroup,	//Clevel -> mínimo tamaño de un grupo
		zoom,		//zoom para imágenes
		nThreads;	//Clevel -> número de threads
  char 		groupsToSeeStr[30]; //grupos a visualizar
 //valores para K-MEANS
  int		kmeansGrNumber,	//kMeans -> número de grupos deseados
		kmeansMaxIterations;//kMeans -> iteraciones máximas
  int		fileFormatType;
		
		
}PROJECT;

typedef struct
  {
  float 	axisLineSize,
		itemSize,
		itemLineSize,
		XYtextSize,
		titleSize,
		itemColor_h,
		itemColor_s,
		itemColor_v;
  int		itemType;
  }PLOT_GRAPHICS_PARAM;
#endif
    