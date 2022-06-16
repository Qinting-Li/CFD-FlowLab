#pragma once
#include <vector>
using namespace std;

namespace PHSPACE
{
class UnstructGrid;

class RotaryCoordinateSystem
{
private:
    vector< double > omiga;
    UnstructGrid * grid;
    int zone_st, zone_ed;
public:
    RotaryCoordinateSystem();
    ~RotaryCoordinateSystem();
public:
    void SetAngularVelocity(vector< double > & omiga) { this->omiga = omiga; }
    void SetZoneIndexRange(int zone_st, int zone_ed)  { this->zone_st = zone_st; this->zone_ed = zone_ed; }
private:
    void RecognizeGrid(UnstructGrid * grid);
};

}