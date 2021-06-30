#include "Pre_FluentToPHengLEI_Unstruct.h"
#include "OverLappingGrid.h"

using namespace std;

namespace PHSPACE
{

void Fluent2PHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");

    int codeOfOversetGrid = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    Pre_FluentToPHengLEI_Unstruct fluentGridConversion(fromGridFile);

    fluentGridConversion.Run();
    int numberOfBlocks = fluentGridConversion.GetNumberofBlocks();
    Grid **grids = fluentGridConversion.GetGrids();

    ConstructOverlappingInfo(grids, numberOfBlocks, codeOfOversetGrid);
    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);
    fluentGridConversion.WriteAdditionalInformation(outGridFile);
    DumpCharacteristicBoundary(outGridFile, grids, numberOfBlocks);

    fluentGridConversion.Clear();
}

}
