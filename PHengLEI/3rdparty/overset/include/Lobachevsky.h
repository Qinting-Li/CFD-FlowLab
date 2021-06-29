#pragma once
#include "mpi.h"
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

void OGSM_SplitOversetGrid();

void OGSM_RegisterCommonWorld( MPI_Comm comm );

void OGSM_CreateOversetWorld( int * processorList, int * groupList, int numberOfZones );

void OGSM_ImportGridTopology( double * x, double * y, double * z, int ni, int nj, int nk, 
						      int * range, int * physicalType, int * normalSign, int * normalSeries, int * target, int numberOfBoundaryRegions, 
						      int * neighborZoneIndexList, int * neighborFaceIndexList, int * localFaceIndexMapping, int numberOfInterfaces, 
						      int iZone );
void OGSM_SewOversetGrid();
void OGSM_ExportOversetInfo( int * iblank, int numberOfCells, int iZone );

void OGSM_ImportPhysicalField( double * q, int dimension, int numberOfNodes );
void OGSM_TradePhysicalField( int dimension );
void OGSM_ExportPhysicalField( double * q, int dimension, int numberOfOversetCells );

void OGSM_DeleteOversetWorld();

#ifdef __cplusplus
}
#endif