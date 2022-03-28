#pragma once
#include <vector>
using namespace std;

namespace PHSPACE
{
class MBGridFactory;
class MBGridGroup;

class MBModualCenter
{
private:
    bool codeOfDynamicMultiBodies;
    int numberOfZones, numberOfGroups;
    vector< int > zone_st, zone_ed;
    MBGridFactory * factory;
    vector< MBGridGroup * > * groupCluster;
public:
    MBModualCenter();
    ~MBModualCenter();
public:
    void SetNumberOfGroups(int ngs) { this->numberOfGroups = ngs;  }
    int  GetNumberOfGroups()          { return this->numberOfGroups; }
private:
    void DetermineGroupRange(vector< int > & span);
    bool GetCodeOfDynamicMultiBodies();
};

}