#include "MBModualCenter.h"
#include "GlobalDataBase.h"
#include "MBGridFactory.h"
#include "MBGridGroup.h"

using namespace std;

namespace PHSPACE
{
MBModualCenter::MBModualCenter()
{
    this->factory      = NULL;
    this->groupCluster = NULL;

    this->codeOfDynamicMultiBodies = GetCodeOfDynamicMultiBodies();
    this->numberOfGroups = GlobalDataBase::GetIntParaFromDB("numberOfHelicopterParts");
    vector< int > span(numberOfGroups);
    int * procSpanForEachParts = reinterpret_cast< int * >(GlobalDataBase::GetDataPtr("procSpanForEachParts"));
    
    int counter = 0;
    for (int i = 0; i < numberOfGroups; ++ i)
    {
        span[i] = procSpanForEachParts[i];
        counter += span[i];
    }
    this->numberOfZones = counter;

    DetermineGroupRange(span);
}

MBModualCenter::~MBModualCenter()
{
    delete factory; factory = NULL;

    for (int iGroup = 0; iGroup < numberOfGroups; ++ iGroup)
    {
        delete (* groupCluster)[iGroup]; (* groupCluster)[iGroup] = NULL;
    }
    delete groupCluster; groupCluster = NULL;
}

void MBModualCenter::DetermineGroupRange(vector< int > & span)
{
    this->zone_st.resize(numberOfGroups, 0);
    this->zone_ed.resize(numberOfGroups, numberOfZones - 1);

    int counter = 0;
    for (int iGroup = 0; iGroup < numberOfGroups - 1; ++ iGroup)
    {
        counter += span[iGroup];

        zone_ed[iGroup    ] = counter - 1;
        zone_st[iGroup + 1] = counter;
    }
    return;
}

bool MBModualCenter::GetCodeOfDynamicMultiBodies()
{
    return GlobalDataBase::GetIntParaFromDB("codeOfDynamicMultiBodies") != 0;
}

}