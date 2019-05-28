/*************************************************************************
 *
 *     Copyright (C) 2018 - Lluc Sementé Fernàndez
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/
#include <RcppGSL.h>
#include <Rcpp.h>
#include "types.h"
#include "ionsSelect.h"
using namespace Rcpp;

// [[Rcpp::export]]

List IonSelectC(double m_focalProb, int numPixels, NumericVector SP_Pixels, int numCols,
                NumericVector massAxis, int numSamples, int nPTestGroups,
                NumericVector R_pTestGroups, NumericVector ClustersSize,
                List ClustersPixels, NumericMatrix data, double zeroThreshold)
{
  //// Reserving memory for the data structures & passing variables ////

  // LOAD_SP //
  LOAD_SP* myLOAD_SP = new LOAD_SP[numSamples];
  float* p_massAxis = new float[numCols];
  for (int i = 0; i<numCols; i++)
  {
    p_massAxis[i] = (float)massAxis[i];
  }
  for (int i = 0; i<numSamples; i++)
  {
    myLOAD_SP[i].nPixels = SP_Pixels[i]; //total number of pixels
    myLOAD_SP[i].mzWidth = numCols;   //numb of columns of the principal matrix
    myLOAD_SP[i].mzAxis_p = p_massAxis; //pointer to the axis of masses
  }


  // GROUP //
  GROUP* myGROUP = new GROUP[ClustersSize.length()];
  int** ClusterId = new int*[ClustersSize.length()];
  for (int i = 0; i<(ClustersSize.length()); i++)
  {
    myGROUP[i].size = ClustersSize[i];  //number of pixels of each cluster
    ClusterId[i] = new int[myGROUP[i].size];
    SEXP ll = ClustersPixels[i];
    NumericVector y(ll);
    for (int j=0; j<myGROUP[i].size; j++)
    {
      ClusterId[i][j]= y[j];
    }
    myGROUP[i].set = ClusterId[i];  //pixels index of each cluster
  }


  // pTestGroups //
  int *pTestGroups = new int[R_pTestGroups.length()];
  for (int i=0; i<R_pTestGroups.length(); i++)
  {
    pTestGroups[i] = R_pTestGroups[i];  //Clusters to analyze
  }


  // allSpectrums_p //
  float** allSpectrums_p;
  allSpectrums_p = new float*[numPixels];
  for (int i = 0; i < numPixels; i++)
  {
    allSpectrums_p[i] = new float[numCols];
    for (int j=0; j<numCols; j++)
    {
      allSpectrums_p[i][j] = data(i,j);
    }
  }

  //// IonsSelect algorithm ////
  /**************************************************************************************************/
  
  IonsSelect myIonsSelect(myLOAD_SP, numSamples, pTestGroups, nPTestGroups, myGROUP, allSpectrums_p);
  myIonsSelect.setThreshold(zeroThreshold);
  myIonsSelect.getMeasures();
  myIonsSelect.setBoundaries(m_focalProb);
  
  int totalCombi = std::pow(2,nPTestGroups); //total number of combinations
  
  List Results_v(totalCombi);
  List Results_z(totalCombi);
  NumericMatrix IonsData(numCols, nPTestGroups*nPTestGroups*3); //totalCombi*3
  NumericMatrix tmp(numCols, nPTestGroups);
  // New variables starts here
  int totalUpRegulatedIons=0, upRegulatedIons;
  int totalDownRegulatedIons=0, downRegulatedIons, regulatedCount=0;
  bool direct;
  int *ions_14=NULL, *ions_23=NULL;
  float *weight14_p=NULL;
  float *weight23_p=NULL;
  ions_14 = new int[numCols];
  ions_23 = new int[numCols];
  weight14_p = new float[numCols];
  weight23_p = new float[numCols];
  // Variables end
  
  for(int gr=0; gr<totalCombi; gr++) //TODO meter la funcion del contraste
  {
    nPTestGroups = myIonsSelect.getGroupsFromCombination(gr, pTestGroups);
    if(nPTestGroups<2)
      continue;
    myIonsSelect.selection(pTestGroups, nPTestGroups);
    
    // New code starts here //////////
    // upRegulatedIons = myIonsSelect.getAllUpRegulatedNumber();
    // downRegulatedIons= myIonsSelect.getAllDownRegulatedNumber();
    // totalUpRegulatedIons += upRegulatedIons;
    // totalDownRegulatedIons+=downRegulatedIons;
    // 

    for(int focalGroup=0; focalGroup<myIonsSelect.m_nAnalyzeGr; focalGroup++)
    {
      Rcpp::Rcout << "focalGroup " << focalGroup << std::endl;
      int count_14=0, count_23=0;   
      
      //Se extraen los iones upregulated (ions_14) y downregulated (ions_23)
      for(int ion=0; ion<numCols; ion++) //para todos los iones
      {
        if( myIonsSelect.m_ionsClusterZ_p[ion][focalGroup]==1 || myIonsSelect.m_ionsClusterV_p[ion][focalGroup]==4)
        {
          ions_14[count_14]=ion;
          count_14++;
        }
        else if( myIonsSelect.m_ionsClusterZ_p[ion][focalGroup]==3 || myIonsSelect.m_ionsClusterV_p[ion][focalGroup]==2)
        {
          ions_23[count_23]=ion;
          count_23++;
        }
      }
      
      Rcpp::Rcout << "before int ionsindex_14" << std::endl;
      //Se ordenan según los valores medios de sus magnitudes en el grupo frente al total de grupos
      int ionsIndex_14[count_14]; //códigos 1 y 4 -> sobre-expresados
      int ionsIndex_23[count_23]; //códigos 2 y 3 ->   sub-expresados
      direct=true;
      Rcpp::Rcout << "getBestIonIntensity call" << std::endl;
      myIonsSelect.getBestIonIntensity(direct, focalGroup, pTestGroups, nPTestGroups, ions_14, count_14, ionsIndex_14, weight14_p);
      direct=false;
      myIonsSelect.getBestIonIntensity(direct, focalGroup, pTestGroups, nPTestGroups, ions_23, count_23, ionsIndex_23, weight23_p);
    }
    // New code ends here //////////////////
    
    for (int i = 0; i < numCols; i++)
    {
      for (int j = 0; j < nPTestGroups; j++)
      {
        tmp(i,j) = myIonsSelect.m_ionsClusterV_p[i][j];
      }
    }
    Results_v[gr] = tmp(Range(0,numCols-1),Range(0,nPTestGroups-1));

    for (int i = 0; i < numCols; i++)
    {
      for (int j = 0; j < nPTestGroups; j++)
      {
        tmp(i,j) = myIonsSelect.m_ionsClusterZ_p[i][j];
      }
    }
    Results_z[gr] = tmp(Range(0,numCols-1),Range(0,nPTestGroups-1));
  }

  for (int i = 0; i < numCols; i++)
  {
    for (int j = 0; j < (nPTestGroups*nPTestGroups*3); j++) //totalCombi
    {
      IonsData(i,j) = myIonsSelect.m_ionsData_p[i][j];
    }
  }

  List Results; //Output structure
  Results["ionsFromVolcano"] = Results_v;
  Results["ionsFromZeros"] = Results_z;
  Results["ionsData"] = IonsData;

  //// Free memory ////

  delete[] myLOAD_SP;
  delete[] p_massAxis;

  for (int i = 0; i<(ClustersSize.length()); i++)
  {
    delete[] ClusterId[i];
  }
  delete[] ClusterId;
  delete[] myGROUP;

  delete[] pTestGroups;

  for (int i = 0; i < numPixels; i++)
  {
    delete[] allSpectrums_p[i];
  }
  delete[] allSpectrums_p;

  return Results;
}
