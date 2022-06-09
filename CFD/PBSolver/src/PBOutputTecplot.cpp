#include <iostream>
#include <fstream>

#include "PBOutputTecplot.h"
#include "Geo_UnstructGrid.h"
#include "Zone.h"
#include "Geometry.h"
#include "HyList.h"
#include "Post_WriteTecplot.h"

using namespace std;

namespace PHSPACE
{
    PBOutputTecplot::PBOutputTecplot(const string &name):PBOutput(name)
    {

    }

    PBOutputTecplot::~PBOutputTecplot()
    {

    }

    void PBOutputTecplot::Run(Zone *zone, vector<std::string> &varNameList)
    {
        UnstructGrid *grid = UnstructGridCast(zone->GetGeometry()->GetGrid());

        if (grid->GetDim() == 2)
        {
            std::string fileName;

            if (ObjGetIntParaFromDB("iunsteady", 0)) //! Transient.
            {
                int intervalStepPlot = ObjGetIntParaFromDB("intervalStepPlot", 1);
                int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

                if (outIterStep % intervalStepPlot != 0)
                {
                    return;
                }

                RDouble currentTime = 0.0;
                GlobalDataBase::GetData("currentTime", &currentTime, PHDOUBLE, 1);
                fileName = Name() + "_time" + toString(currentTime) + "_zone" + toString(zone->GetIndex()) + ".dat";
            }
            else
            {
                int innerIter = 0;
                GlobalDataBase::GetData("innerIter", &innerIter, PHINT, 1);
                fileName = Name() + "_iter" + toString(innerIter) + "_zone" + toString(zone->GetIndex()) + ".dat";
            }

            std::transform(fileName.begin(), fileName.end(), fileName.begin(), ::tolower);

            std::string outputFileName = "results/" + fileName;
            std::ofstream file(outputFileName.c_str());

            std::vector<RDouble*> phi;

            for (vector<string>::iterator iter = varNameList.begin(); iter != varNameList.end(); ++ iter)
            {
                phi.push_back((RDouble*)grid->GetDataPtr(*iter));
            }

            if (file.is_open())
            {
                int nTotalNode = grid->GetNTotalNode();
                int nTotalCell = grid->GetNTotalCell();
                int nBoundFace = grid->GetNBoundFace();
                int nTotalFace = grid->GetNTotalFace();

                int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
                int *face2node = grid->GetFace2Node();
                int *left_cell_of_face = grid->GetLeftCellOfFace();
                int *right_cell_of_face = grid->GetRightCellOfFace();

                //! TODO: Replace hylist and save it in advance.
                //! Find out cellToFace.
                HyList < int > *cell2face = new HyList < int >[nTotalCell];
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    cell2face[iCell].SetAverNnode(4);
                }
                for (int iFace = 0; iFace < nBoundFace; ++ iFace)
                {
                    int le = left_cell_of_face[iFace];
                    cell2face[le].insert(iFace);
                }
                for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
                {
                    int le = left_cell_of_face[iFace];
                    int re = right_cell_of_face[iFace];
                    cell2face[le].insert(iFace);
                    cell2face[re].insert(iFace);
                }

                //! Find out cellToNode.
                HyList < int > *cell2node = new HyList < int >[nTotalCell];
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    cell2node[iCell].SetAverNnode(4);
                }

                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    int face0 = cell2face[iCell][0];
                    int ip1 = face2node[2 * face0];
                    int ip2 = face2node[2 * face0 + 1];
                    cell2node[iCell].insert(ip1);
                    cell2node[iCell].insert(ip2);

                    int targetNode = ip2;
                    int checkNode = ip1;

                    int nFaceOfCell = cell2face[iCell].size();
                    int isClockWise = false;
                    while (!isClockWise)
                    {
                        for (int iface = 1; iface < nFaceOfCell; ++ iface)
                        {
                            int face = cell2face[iCell][iface];
                            int node1 = face2node[2 * face];
                            int node2 = face2node[2 * face + 1];

                            if (node1 == targetNode)
                            {
                                if (node2 != checkNode)
                                {
                                    cell2node[iCell].insert(node2);
                                    targetNode = node2;
                                }
                                else
                                {
                                    isClockWise = true;
                                    break;
                                }
                            }
                            else if (node2 == targetNode)
                            {
                                if (node1 != checkNode)
                                {
                                    cell2node[iCell].insert(node1);
                                    targetNode = node1;
                                }
                                else
                                {
                                    isClockWise = true;
                                    break;
                                }

                            }

                        }   //! Iterate over all faces of cell.

                    }       //! Judge whether the nodes of cell have searched closed.

                }           //! iCell.

                file << " VARIABLES = \"X\",\"Y\"";
                for (int i = 0; i < varNameList.size(); ++ i)
                {
                    file << ",\"" << varNameList[i] << "\"";
                }
                file << " \n";

                if (ObjGetIntParaFromDB("isunsteady", 0))
                {
                    file << "ZONE T=\"" << toString(ObjGetDoubleParaFromDB("currentTime", 0.0)) 
                        << "s\",SOLUTIONTIME=" << toString(ObjGetDoubleParaFromDB("currentTime", 0.0))
                        << ",N = " << grid->GetNTotalNode()
                        << ",E=" << grid->GetNTotalCell()
                        << ",DATAPACKING=BLOCK,ZONETYPE=FEQUADRILATERAL,VARLOCATION=([";
                }
                else
                {
                    file << " ZONE T=\"Rampant\",N=" << grid->GetNTotalNode()
                        << ",E=" << grid->GetNTotalCell()
                        << ",DATAPACKING=BLOCK,ZONETYPE=FEQUADRILATERAL,VARLOCATION=([";
                }

                for (int i = 0; i < varNameList.size(); ++ i)
                {
                    file << (i + 3);

                    if (i != (varNameList.size() - 1))
                    {
                        file << ",";
                    }
                }
                file << "]=CELLCENTERED)\n";

                for (int iNode = 0; iNode < grid->GetNTotalNode(); ++ iNode)
                {
                    file << grid->GetX()[iNode] << "\n";
                }

                for (int iNode = 0; iNode < grid->GetNTotalNode(); ++ iNode)
                {
                    file << grid->GetY()[iNode] << "\n";
                }

                for (int i = 0; i < varNameList.size(); ++i)
                {
                    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
                    {
                        file << phi[i][iCell] << "\n";
                    }
                }

                for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
                {
                    if (cell2node[iCell].size() == 3)
                    {
                        for (int m = 0; m < 3; ++m)
                        {
                            file << cell2node[iCell].GetData(m) + 1 << "\n";
                        }
                        file << cell2node[iCell].GetData(2) + 1 << "\n";
                    }
                    else if (cell2node[iCell].size() == 4)
                    {
                        for (int m = 0; m < 4; ++m)
                        {
                            file << cell2node[iCell].GetData(m) + 1 << "\n";
                        }
                    }
                    else
                    {
                        std::cout << "Error: this function only support triangle and quad !\n";
                        std::exit(0);
                    }
                }
                delete[] cell2node;
                delete[] cell2face;


                file.close();
            }

            vector<RDouble *>().swap(phi);
        }
        else if (grid->GetDim() == 3)
        {
            Post_WriteTecplot generateVisualFile;
            generateVisualFile.Run();
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBOutput, PBOutputTecplot, OUTPUT_TECPLOT);
}