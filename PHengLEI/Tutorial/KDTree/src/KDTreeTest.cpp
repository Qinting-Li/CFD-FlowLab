//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      ToolKitTest.cpp
//! @brief     It tests all public interfaces of ToolKit.
//! @author    Meng Liyuan.
#pragma once
#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "DataStruct_KDTree.h"
#include "Precision.h"
#include "TK_Exit.h"

#define DEF_NUM_PTS 200000
#define DEF_INNER_PTS 1000000
#define DEF_DIMS 3

using namespace PHSPACE;

//! Declare and define the functions of several methods that terminate the program .
RDouble DistSQ(RDouble *a1, RDouble *a2, int dim);
RDouble Random(void);
void GenerateDims(int num, RDouble *nums);

int main(int argc, char **argv) 
{
    int i, numberOfPoints = DEF_NUM_PTS;
    KDTree *ptree = new KDTree();
    char *data, *pch;
    KDRes *pResults;
    KDRes *minResults;
    RDouble pos[DEF_DIMS], dist;
    RDouble pt[DEF_DIMS];
    RDouble nums[DEF_DIMS];
    RDouble radius = 1.0;

    if ( argc > 1 && isdigit(argv[1][0]) )
    {
        numberOfPoints = atoi(argv[1]);
    }

    if ( !(data = new char[numberOfPoints]) )
    {
        cout << "New Failed" << endl;
    }

    srand(time(0));
    
    ptree = CreatKDTree(DEF_DIMS);
    cout << "Seeding " << numberOfPoints << " entries of " << DEF_DIMS << " dimensions." << endl;
    for (int i=0; i<numberOfPoints; ++i)
    {
        data[i] = 0;
        GenerateDims(DEF_DIMS, nums);
        KDInsert(ptree, nums, &data[i]);
    }

    dist = LARGE;
    for (int n=0; n<DEF_INNER_PTS; ++n)
    {
        for (i=0; i<DEF_DIMS; ++i)
        {
            pt[i] = Random();
        }

        //pResults = NearestRange(ptree, pt, radius);
        //
        //cout << "Found " << pResults->Size() << " results!" << endl;

        //while (!pResults->End())
        //{
        //    pch = (char*)pResults->Item(pos);
        //    RDouble dist1 = sqrt(DistSQ(pt, pos, DEF_DIMS));
        //    //cout << "Node at (" << pos[0] << " "<< pos[1] << " "<< pos[2] << ") to point (" 
        //    //    << pt[0] << " "<< pt[1] << " "<< pt[2] << ") is " << dist1 << " away and has data " << *pch << endl;
        //    pResults->Next();
        //}

        minResults = KDNearest(ptree, pt);

        pch = (char*)minResults->Item(pos);
        RDouble dist0 = sqrt(DistSQ(pt, pos, DEF_DIMS));
        RDouble dist1 = sqrt(minResults->rlist->distSQ);
        dist = min(dist, dist0);
    }

    cout << endl;
    cout << "Node at (" << pos[0] << " "<< pos[1] << " "<< pos[2] << ") to point (" 
    << pt[0] << " "<< pt[1] << " "<< pt[2] << ") is min " << dist << " away and has data " << *pch << endl;

    return 0;
}

RDouble DistSQ(RDouble *a1, RDouble *a2, int dim)
{
    RDouble distSQ = 0.0;
    while (--dim >=0)
    {
        RDouble diff = a1[dim] - a2[dim];
        distSQ += diff * diff;
    }
    return distSQ;
}

RDouble Random(void)
{
    return (RDouble)rand() / RAND_MAX * 20.0 - 10.0;
}

void GenerateDims(int num, RDouble *nums)
{
    for (int i=0; i<num; ++i)
    {
        nums[i] = Random();
    }
}
