#include "Pre_GridBase.h"
#include "MixGrid.h"
#include "TK_Parse.h"
#include "IO_FileName.h"
#include "LinkStruct.h"
#include "Glb_Dimension.h"
#include "MultigridFactory.h"
#include "Geo_StructBC.h"
#include "Geo_UnstructBC.h"
#include "TK_Exit.h"
#include "GridType.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{

void Fantasy2Ustar()
{
    if (GetDim() == TWO_D)
    {
        Fantasy2Ustar2D();
    }
    else if (GetDim() == THREE_D)
    {
        Fantasy2Ustar3D();
    }
    else
    {
        TK_Exit::ExceptionExit("Error: Can't convert this dim to Ustar !\n");
    }

    PrintToWindow("HyperFlow grid has been converted to Ustar !\n");
}

void FantasyWriteToInp(const string &bnd_file, Grid **grids, int nBlocks)
{
    fstream GridgenBCfile;
    OpenFile(GridgenBCfile, bnd_file, ios_base::out);

    GridgenBCfile << 1 << endl;
    GridgenBCfile << nBlocks << endl;

    int wordwidth = 8;
    int nCoords = GetDim();

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        StructGrid *grid = StructGridCast(grids[iBlock]);

        int idim = grid->GetNI();
        int jdim = grid->GetNJ();
        int kdim = grid->GetNK();

        GridgenBCfile << setw(wordwidth) << idim;
        GridgenBCfile << setw(wordwidth) << jdim;
        if (nCoords == 3)
        {
            GridgenBCfile << setw(wordwidth) << kdim;
        }
        GridgenBCfile << endl;
        GridgenBCfile << iBlock + 1 << endl;

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int  numberOfBoundaryFaces = composite_bcregion->GetnBCRegion();
        GridgenBCfile << numberOfBoundaryFaces << endl;

        for (int ibcregion = 0; ibcregion < numberOfBoundaryFaces; ibcregion ++)
        {
            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();
            int bctype = structBC->GetBCType();
            
            for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
            {
                GridgenBCfile << setw(8) << s_st[iCoords];
                GridgenBCfile << setw(8) << s_ed[iCoords]; 
            }
            GridgenBCfile << setw(8) << bctype;
            GridgenBCfile << endl;

            if (bctype == -1)
            {
                int *t_st = structBC->GetTargetStart();
                int *t_ed = structBC->GetTargetEnd();
                int nbt = structBC->GetTargetRegionBlock();

                for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
                {
                    GridgenBCfile << setw(8) << t_st[iCoords];
                    GridgenBCfile << setw(8) << t_ed[iCoords];
                }
                GridgenBCfile << setw(8) << nbt + 1;
                GridgenBCfile << endl;
            }
        }
    }

    CloseFile(GridgenBCfile);
}

void FantasyWriteToGrd(const string &cel_file, Grid **grids, int nBlocks)
{
    fstream Gridgenfile;
    OpenFile(Gridgenfile, cel_file, ios_base::out|ios_base::binary|ios_base::trunc);
    Gridgenfile.write(reinterpret_cast<char *>(&nBlocks), sizeof(int));

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        StructGrid *grid = StructGridCast(grids[iBlock]);

        int idim = grid->GetNI();
        int jdim = grid->GetNJ();
        int kdim = grid->GetNK();

        Gridgenfile.write(reinterpret_cast<char *>(&idim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&jdim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&kdim), sizeof(int));
    }

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        StructGrid *grid = StructGridCast(grids[iBlock]);

        int ist, ied, jst, jed, kst, ked;
        grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

        grid->WriteXYZ(Gridgenfile, ist, ied, jst, jed, kst, ked);
    }

    CloseFile(Gridgenfile);
}

void ReProcessBCInfo(Grid **grids, int nBlocks)
{
    for (int iBlock = 0; iBlock < nBlocks; iBlock ++)
    {
        StructGrid *grid = StructGridCast(grids[iBlock]);
        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int numberOfBoundaryFaces = composite_bcregion->GetnBCRegion();
        for (int ibcregion = 0; ibcregion < numberOfBoundaryFaces; ibcregion ++)
        {
            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (structBC->GetFaceDirection() == m)
                {
                    if (s_st[m] != 1)
                    {
                        s_st[m] += 1;
                        s_ed[m] += 1;
                    }
                }
                else
                {
                    s_ed[m] += 1;
                }
            }
        }
    }
}

void ReadLnkFile(const string &fileName_in, Grid **grids, int nBlocks)
{
    string LnkFileName = ChangeExtensionOfFileName(fileName_in, "link");
    if (!SerialFileExist(LnkFileName))
    {
        TK_Exit::ExceptionExit("缺少link文件");
        return;
    }

    fstream file;
    OpenSerialFile(file, LnkFileName, ios_base::in|ios_base::binary);

    DataContainer *LnkData = new DataContainer();
    ReadFile(file, LnkData);

    VirtualFile *virtualFile = new VirtualFile(LnkData);
    virtualFile->BeginReadWork();

    for (int iBlock = 0; iBlock < nBlocks; iBlock ++)
    {
        StructGrid *grid = StructGridCast(grids[iBlock]);

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int nBCRegion = composite_bcregion->GetnBCRegion();
        for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
        {
            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

            int bctype = structBC->GetBCType();
            if (bctype >= 0) continue;

            int imin, imax, jmin, jmax, kmin, kmax, nbt;

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (GetDim() == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            structBC->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (GetDim() == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            PHRead(virtualFile, nbt);

            structBC->SetTargetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            structBC->SetTargetRegionBlock(nbt - 1);
            if (GetDim() == THREE_D)
            {
                structBC->ComputeRelativeParameters();
            }
        }
    }

    virtualFile->EndReadWork();

    delete virtualFile;
    delete LnkData;

    CloseFile(file);
}

void WriteLnkFile (const string &gridFileName, Grid **grids_in, int nBlocks)
{
    fstream file;

    Grid **grids = grids_in;

    string LnkFileName = ChangeExtensionOfFileName(gridFileName, "link");

    PHSPACE::OpenSerialFile(file, LnkFileName, ios_base::out|ios_base::binary);

    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    for (int izone = 0; izone < nBlocks; izone ++)
    {
        StructGrid *grid = StructGridCast(grids[izone]);

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int nBCRegion = composite_bcregion->GetnBCRegion();
        for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
        {
            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

            int bctype = structBC->GetBCType();
            if (bctype >= 0) continue;

            int dimension = GetDim();
            int *lnkInfo = structBC->GetInkInfo();

            PHWrite(virtualFile, lnkInfo[0]);
            PHWrite(virtualFile, lnkInfo[1]);
            PHWrite(virtualFile, lnkInfo[2]);
            PHWrite(virtualFile, lnkInfo[3]);

            if (dimension == THREE_D)
            {
                PHWrite(virtualFile, lnkInfo[4]);
                PHWrite(virtualFile, lnkInfo[5]);
            }

            PHWrite(virtualFile, lnkInfo[6]);
            PHWrite(virtualFile, lnkInfo[7]);
            PHWrite(virtualFile, lnkInfo[8]);
            PHWrite(virtualFile, lnkInfo[9]);

            if (dimension == THREE_D)
            {
                PHWrite(virtualFile, lnkInfo[10]);
                PHWrite(virtualFile, lnkInfo[11]);
            }

            int nbt = structBC->GetTargetRegionBlock() + 1;
            PHWrite(virtualFile, nbt);
        }
    }

    virtualFile->EndWriteWork();
    delete virtualFile;

    PHSPACE::CloseFile(file);
}

void CheckMeshMultigrid(Grid **grids_in, int nBlocks)
{
    Grid **grids = grids_in;

    int km, NN;
    int km_grid = 1024;

    for (int izone = 0; izone < nBlocks; izone ++)
    {
        StructGrid *grid = StructGridCast(grids[izone]);

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int  nBCRegion = composite_bcregion->GetnBCRegion();
        for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
        {
            km = 1;
            NN = 2;

            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);
            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();

            while ((((abs(s_st[0])-1)%NN)==0) && (((abs(s_ed[0])-1)%NN)==0) && (((abs(s_st[1])-1)%NN)==0) && (((abs(s_ed[1])-1)%NN)==0) && (((abs(s_st[2])-1)%NN)==0) && (((abs(s_ed[2])-1)%NN)==0))
            {
                NN = NN * 2;
                km = km + 1;
            }
            km_grid = min(km_grid, km);
        }
    }

    PrintToWindow("The most valid multi-grid level is: ", km_grid, "\n");
}

void WriteBcFileFromFTS (const string &gridFileName, Grid **grids_in, int nBlocks)
{
    Grid **grids = grids_in;

    fstream file;

    string BcFileName = ChangeExtensionOfFileName(gridFileName, "bc");
    PHSPACE::OpenSerialFile(file, BcFileName, ios_base::out|ios_base::binary);

    PHWrite(file, nBlocks);

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        VirtualFile *virtualFile = new VirtualFile(&file);
        virtualFile->BeginWriteWork();

        PHWrite(virtualFile, iZone);

        StructGrid *grid = StructGridCast(grids[iZone]);

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

        int nBCRegion = static_cast< int > (composite_bcregion->GetnBCRegion());

        PHWrite(virtualFile, nBCRegion);

        string *bcNameList = new string [nBCRegion];

        uint_long totalSize = 0;
        for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
        {
            StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

            string bcName = structBC->GetBCName();
            int bctype = structBC->GetBCType();
            if (bcName == "")
            {
                if (bctype == -1)
                {
                    bcName = "Connection";
                }
                else
                {
                    bcName = composite_bcregion->GetBCName(bctype);
                    if ((bctype >7) && (bctype <51))
                    {
                        bctype = 2;
                        structBC->SetBCType(bctype);
                    }
                }
                structBC->SetBCName(bcName);
            }
            
            bcNameList[ibcregion] = bcName;
            
            //totalSize += bcNameList[ibcregion].size();
            totalSize += static_cast< uint_long > (bcNameList[ibcregion].size());
            totalSize += 1;
        }

        char *bcNameChar = new char [totalSize];
        unsigned int count = 0;
        for (int boco = 0; boco < nBCRegion; ++ boco)
        {
            string &bcName = bcNameList[boco];
            streamsize nameSize = bcName.size();
            for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
            {
                bcNameChar[count ++] = bcName[iChar];
            }
            bcNameChar[count ++] = '\0';
        }

        PHWrite(virtualFile, totalSize);
        PHWrite(virtualFile, bcNameChar, static_cast<int>(totalSize));

        delete [] bcNameList;
        delete [] bcNameChar;

        virtualFile->EndWriteWork();
        delete virtualFile;
    }

    PHSPACE::CloseFile(file);
}

void Fantasy2Ustar2D()
{
    string from_gfile = "grid.fts";
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    string out_gfile = "grid.fts";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    fstream file;
    file.open(from_gfile.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(from_gfile);
    }

    Ustar2Fantasy u2fbc("support/fts2ustarbc.txt");

    using namespace PHSPACE;

    int nBlocks;
    file.read(reinterpret_cast<char*>(&nBlocks), sizeof(int));

    int *block_proc = new int [nBlocks];
    int *block_idx  = new int [nBlocks];
    int *block_type = new int [nBlocks];
    file.read(reinterpret_cast<char*>(block_proc), nBlocks * sizeof(int));
    file.read(reinterpret_cast<char*>(block_idx),  nBlocks * sizeof(int));
    file.read(reinterpret_cast<char*>(block_type), nBlocks * sizeof(int));
    delete [] block_proc;
    delete [] block_idx;
    delete [] block_type;

    streamsize nlen = 0;
    file.read(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    int nTotalNode, nTotalFace, nTotalCell, nBoundFace;

    cout << "reading the mesh files......\n";
    //! Read the total node number,face number and cell number.
    int sizeint = sizeof(int);
    int sizeRealGeom = sizeof(RDouble);

    file.read(reinterpret_cast<char*>(&nTotalNode), sizeint);
    file.read(reinterpret_cast<char*>(&nTotalFace), sizeint);
    file.read(reinterpret_cast<char*>(&nTotalCell), sizeint);

    cout << "the total numbers of nodes     ：" << nTotalNode;
    cout << "the total numbers of faces   ：" << nTotalFace;
    cout << "the total numbers of cells     ：" << nTotalCell << "\n";

    RDouble *x, *y, *z;
    RDouble xmin, ymin, xmax, ymax;
    xmin = ymin =   LARGE;
    xmax = ymax = - LARGE;
    x = new RDouble[nTotalNode];
    y = new RDouble[nTotalNode];
    z = new RDouble[nTotalNode];

    file.read(reinterpret_cast<char*>(x), nTotalNode * sizeRealGeom);
    file.read(reinterpret_cast<char*>(y), nTotalNode * sizeRealGeom);
    file.read(reinterpret_cast<char*>(z), nTotalNode * sizeRealGeom);

    cout << "mesh nodes reading end.\n";

    int *node_number_of_each_face = new int[nTotalFace];
    //! Assign the value to node_number_of_each_face.
    file.read(reinterpret_cast<char*>(node_number_of_each_face), nTotalFace * sizeint);

    int nsum = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        nsum += node_number_of_each_face[iface];
    }

    cout << "set the connection relationship of face to node......\n";
    //! Face to Node connection relationship.
    int *face2node = new int [nsum];

    file.read(reinterpret_cast<char*>(face2node), nsum * sizeint);

    cout << "set the connection relationship of face to cell......\n";

    int *left_cell_of_face  = new int [nTotalFace];
    int *right_cell_of_face = new int [nTotalFace];
    file.read(reinterpret_cast<char*>(left_cell_of_face), nTotalFace * sizeint);
    file.read(reinterpret_cast<char*>(right_cell_of_face), nTotalFace * sizeint);

    cout << "set the boundary conditions......\n";
    file.read(reinterpret_cast<char*>(&nBoundFace), sizeint);

    //! Set the boundary conditions.
    int *bctype = new int [nBoundFace];
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int tmp;
        file.read(reinterpret_cast<char*>(&tmp), sizeint);
        bctype[iface] = tmp;
    }

    file.close();
    file.clear();

    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int type = u2fbc.TranslateBC(bctype[iface]);
        bctype[iface] = type;
    }  

    file.open(out_gfile.c_str(), ios_base::out);
    file.setf(ios::scientific);
    file.precision(8);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(out_gfile);
    }

    int ngrids = 1;
    file << ngrids << "\n";
    file << 0 << "\n";
    file << "grid 1\n";
    file << nTotalNode << " " << nTotalFace << " " << nTotalCell << " " << nBoundFace << "\n";
 
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        file << x[inode] << " " << y[inode] << "\n";
    }

    int *nodePosi = new int [nTotalFace+1];

    nodePosi[0] = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        nodePosi[iface+1] = nodePosi[iface] + 2;
    }
    delete [] nodePosi;

    for (int i = 0; i < nsum; ++ i)
    {
        file << face2node[i] + 1;
        if (i % 2 == 0)
        {
            file << " ";
        }
        else
        {
            file << "\n";
        }
    }

    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        file << left_cell_of_face[ iface ] + 1 << " ";
        file << bctype[iface]   << "\n";
    }

    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        file << left_cell_of_face [ iface ] + 1 << " ";
        file << right_cell_of_face[ iface ] + 1 << "\n";
    }

    file.close();
    file.clear();

    delete [] bctype;
    delete [] node_number_of_each_face;
    delete [] face2node;
    delete [] left_cell_of_face;
    delete [] right_cell_of_face;
    delete [] x;
    delete [] y;
    delete [] z;
}

void Fantasy2Ustar3D()
{
    string from_gfile = "grid.fts";
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    string out_gfile = "mgrid.in";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    fstream file;
    file.open(from_gfile.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(from_gfile);
    }

    typedef Ustar2Fantasy FantasyToUstar;
    FantasyToUstar u2fbc("support/fts2ustarbc.txt");

    using namespace PHSPACE;

    int nBlocks;
    file.read(reinterpret_cast<char*>(&nBlocks), sizeof(int));

    int *block_proc = new int [nBlocks];
    int *block_idx  = new int [nBlocks];
    int *block_type = new int [nBlocks];
    file.read(reinterpret_cast<char*>(block_proc), nBlocks * sizeof(int));
    file.read(reinterpret_cast<char*>(block_idx),  nBlocks * sizeof(int));
    file.read(reinterpret_cast<char*>(block_type), nBlocks * sizeof(int));
    delete [] block_proc;
    delete [] block_idx;
    delete [] block_type;

    streamsize nlen = 0;
    file.read(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    int nTotalNode, nTotalFace, nTotalCell, nBoundFace;

    cout << "reading the mesh files......\n";
    //! Read the total node number, face number and cell number.
    int sizeint = sizeof(int);
    int sizeRealGeom = sizeof(RDouble);

    file.read(reinterpret_cast<char*>(&nTotalNode), sizeint);
    file.read(reinterpret_cast<char*>(&nTotalFace), sizeint);
    file.read(reinterpret_cast<char*>(&nTotalCell), sizeint);

    cout << "the total numbers of nodes     ：" << nTotalNode;
    cout << "the total numbers of faces   ：" << nTotalFace;
    cout << "the total numbers of cells     ：" << nTotalCell << "\n";

    RDouble *x, *y, *z;
    RDouble xmin, ymin, xmax, ymax;
    xmin = ymin =   LARGE;
    xmax = ymax = - LARGE;
    x = new RDouble[nTotalNode];
    y = new RDouble[nTotalNode];
    z = new RDouble[nTotalNode];

    file.read(reinterpret_cast<char*>(x), nTotalNode * sizeRealGeom);
    file.read(reinterpret_cast<char*>(y), nTotalNode * sizeRealGeom);
    file.read(reinterpret_cast<char*>(z), nTotalNode * sizeRealGeom);

    cout << "mesh nodes reading end.\n";

    int *node_number_of_each_face = new int[nTotalFace];
    //! Assign the value to node_number_of_each_face.
    file.read(reinterpret_cast<char*>(node_number_of_each_face), nTotalFace * sizeint);

    int nsum = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        nsum += node_number_of_each_face[iface];
    }

    cout << "set the connection relationship of face to node ......\n";
    //! Face to Node connection relationship.
    int *face2node = new int [nsum];

    file.read(reinterpret_cast<char*>(face2node), nsum * sizeint);

    cout << "set the connection relationship of face to cell ......\n";

    int *left_cell_of_face  = new int[nTotalFace];
    int *right_cell_of_face = new int[nTotalFace];
    file.read(reinterpret_cast<char*>(left_cell_of_face), nTotalFace * sizeint);
    file.read(reinterpret_cast<char*>(right_cell_of_face), nTotalFace * sizeint);

    cout << "set the boundary conditions......\n";
    file.read(reinterpret_cast<char*>(&nBoundFace), sizeint);

    //! Set the boundary conditions.
    int *bctype = new int [nBoundFace];
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int tmp;
        file.read(reinterpret_cast<char*>(&tmp), sizeint);
        bctype[iface] = tmp;
    }
    file.close();
    file.clear();

    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int type = u2fbc.TranslateBC(bctype[iface]);
        bctype[iface] = type;
    }  

    //! Output file in Ustar format.
    const int BINARY = 0;
    const int ASCII  = 1;

    int fileformat = 0;
    GlobalDataBase::GetData("fileformat", &fileformat, PHINT, 1);

    if (fileformat == ASCII)
    {
        file.open(out_gfile.c_str(), ios_base::out);
        file.setf(ios::scientific);
        file.precision(8);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(out_gfile);
        }

        file << nTotalNode << "  " << nTotalFace << "  " << nTotalCell << endl;

        for (int inode = 0; inode < nTotalNode; ++ inode)
        {
            file << x[inode] << "  " << y[inode] << "  " << z[inode] << endl;
        }

        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            file << node_number_of_each_face[iface] << endl;
        }

        int count = 0;
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                file << face2node[count++] + 1 << "  ";
            }
            file << endl;
        }

        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            file << left_cell_of_face[iface] + 1 << " ";
            file << bctype[iface] << "\n"; 
        }

        for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
        {
            file << left_cell_of_face [iface] + 1 << " ";
            file << right_cell_of_face[iface] + 1 << "\n";
        }

        file.close();
        file.clear();
    } 
    else
    {
        //! BINARY format.
        int count;
        int SizeInt    = sizeof(int);
        int SizeDouble = sizeof(RDouble);

        file.open(out_gfile.c_str(), ios_base::out|ios_base::binary);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(out_gfile);
        }

        file.write((char *)&nTotalNode, SizeInt);
        file.write((char *)&nTotalFace, SizeInt);
        file.write((char *)&nTotalCell, SizeInt);

        //! Output coordinates.
        count = 3 * nTotalNode;
        RDouble *coordinate = new RDouble [count];
        count = 0;
        for (int i = 0; i < nTotalNode; ++ i)
        {
            coordinate[count ++] = x[i];
            coordinate[count ++] = y[i];
            coordinate[count ++] = z[i];
        }

        count = 3 * nTotalNode;
        file.write((char *)coordinate, count * SizeDouble);
        delete [] coordinate;

        //! Output face to node relationship.
        file.write((char *)node_number_of_each_face, nTotalFace * SizeInt);

        count = 0;
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                ++ count;
            }
        }

        int *f2n_temp = new int [count];
        for (int j = 0; j < count; ++ j)
        {
            f2n_temp[j] = face2node[j] + 1;
        }
        
        file.write((char *)f2n_temp, count * SizeInt);
        delete [] f2n_temp;

        //! Output face to cell relationship.
        count = 2 * nTotalFace;
        int *f2c_temp = new int [count];
        
        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            f2c_temp[2 * iface] = left_cell_of_face[iface] + 1;
            f2c_temp[2 * iface + 1] = bctype[iface];
        }

        for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
        {
            f2c_temp[2 * iface] = left_cell_of_face[iface] + 1;
            f2c_temp[2 * iface + 1] = right_cell_of_face[iface] + 1;
        }

        file.write((char *)f2c_temp, count * SizeInt);
        delete [] f2c_temp;
    }    

    delete [] bctype;
    delete [] node_number_of_each_face;
    delete [] face2node;
    delete [] left_cell_of_face;
    delete [] right_cell_of_face;
    delete [] x;
    delete [] y;
    delete [] z;
}

Ustar2Fantasy::Ustar2Fantasy(string mapfile_in)
{
    mapfilename = mapfile_in;

    fstream file;
    string name = mapfilename;
    file.open(name.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(name);
    }

    int num_bc = 0;

    file >> num_bc;

    for (int i = 0; i < num_bc; ++ i)
    {
        int index, bcindex;
        file >> index  ;
        file >> bcindex;
        bcmap.insert(pair<int, int>(index, bcindex));
    }

    file.close();
    file.clear();
}

Ustar2Fantasy::~Ustar2Fantasy()
{
}

int Ustar2Fantasy::TranslateBC(int bctype)
{
    map <int, int>::iterator iter;
    iter = bcmap.find(bctype);

    if (iter != bcmap.end())
    {
        return iter->second;
    }

    return bctype;
}

//! Add binary ustar grid input.
void ReadUstarGrid(const string &gridfile, UnstructGrid *grid)
{
    const int BINARY = 0;
    const int ASCII  = 1;

    int fileformat = 0;
    GlobalDataBase::GetData("fileformat", &fileformat, PHINT, 1);

    if (fileformat == BINARY)
    {
        cout << "reading binary Ustar3D meshes ... \n";
        ReadUstarGrid3D_Binary(gridfile, grid);
        return;
    }
    else
    {
        cout << "reading ascii Ustar meshes ... \n";
    }

    string line, word;
    string separator = " =\t\r\n#$,;";

    fstream file;
    file.open(gridfile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridfile);
    }

    Ustar2Fantasy u2fbc("support/ustar2ftsbc.txt");

    if (GetDim() == TWO_D)
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        int ngrids = 1;
        from_string< int >(ngrids, word, std::dec);

        cout << "ngrids = " << ngrids << "\n";

        getline(file, line);
        getline(file, line);
    }

    int nTotalNode = 0, nTotalFace = 0, nTotalCell = 0, nBoundFace = 0;
    getline(file, line);

    line = FindNextWord(line, word, separator);
    from_string< int >(nTotalNode, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(nTotalFace, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(nTotalCell, word, std::dec);

    cout << "reading the unstructed data file......\n";
    //! Read the total node number,face number and cell number.

    if (GetDim() == TWO_D)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(nBoundFace, word, std::dec);
        cout << " the numbers of total nodes ："   << nTotalNode;
        cout << " the numbers of total faces ：" << nTotalFace;
        cout << " the numbers of total cells ："   << nTotalCell << "\n";
        cout << " the numbers of boundary faces ："   << nBoundFace << "\n";

    }
    else
    {
        cout << " the numbers of total nodes ："   << nTotalNode;
        cout << " the numbers of total faces ：" << nTotalFace;
        cout << " the numbers of total cells ："   << nTotalCell << "\n";
    }

    grid->SetNTotalNode(nTotalNode);
    grid->SetNTotalFace(nTotalFace);
    grid->SetNTotalCell(nTotalCell);

    RDouble *x, *y, *z;
    x = new RDouble [nTotalNode];
    y = new RDouble [nTotalNode];
    z = new RDouble [nTotalNode];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(x[inode], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(y[inode], word, std::dec);

        if (GetDim() == THREE_D)
        {
            line = FindNextWord(line, word, separator);
            from_string< RDouble >(z[inode], word, std::dec);
        }
        else
        {
            z[inode] = zero;
        }
    }

    int *node_number_of_each_face = new int[nTotalFace];
    //! Assign the value to node_number_of_each_face.
    grid->SetNodeNumberOfEachFace(node_number_of_each_face);

    if (GetDim() == TWO_D)
    {
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            node_number_of_each_face[iface] = 2;
        }
    }
    else
    {
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            file >> node_number_of_each_face[iface];
        }
    }

    int *nodePosi = new int[nTotalFace + 1];
    int nsum = 0;
    nodePosi[0] = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        nsum += node_number_of_each_face[iface];
        nodePosi[iface+1] = nsum;
    }

    cout << "set the connection relationship of face to node......\n";
    //! Face to Node connection relationship.
    int *face2node = new int [nsum];
    grid->SetFace2Node(face2node);

    for (int i = 0; i < nsum; ++ i)
    {
        file >> face2node[i];
        -- face2node[i];
    }

    int nBFaceReal = 0;
    int n_patch    = 0;

    cout << "set the connection relationship of face to cell......\n";
    int *left_cell_of_face  = new int[nTotalFace];
    int *right_cell_of_face = new int[nTotalFace];
    grid->SetLeftCellOfFace(left_cell_of_face);
    grid->SetRightCellOfFace(right_cell_of_face);

    if (GetDim() == TWO_D)
    {
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            file >> left_cell_of_face[iface];
            file >> right_cell_of_face[iface];

            if (left_cell_of_face[ iface ] <= 0)
            {
                -- left_cell_of_face [iface];
                -- right_cell_of_face[iface];
                // need to reverse the node ordering
                SWAP(face2node[2*iface], face2node[2*iface+1]);
                SWAP(left_cell_of_face[iface],right_cell_of_face[iface]);
                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else if (right_cell_of_face[iface] <= 0)
            {
                -- left_cell_of_face [iface];
                -- right_cell_of_face[iface];
                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else
            {
                -- left_cell_of_face [iface];
                -- right_cell_of_face[iface];
            }
        }
    }
    else
    {
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            file >> left_cell_of_face [iface];
            file >> right_cell_of_face[iface];
            if (left_cell_of_face[iface] < 0)
            {
                -- right_cell_of_face[iface];
                //! Need to reverse the node ordering.
                //! This problem is very hidden,mainly cause by the unclear definition of reverse.
                //! The bottom line is wrong.
                //std::reverse(face2node + nodePosi[iface],face2node + nodePosi[iface+1] - 1);
                std::reverse(face2node + nodePosi[iface], face2node + nodePosi[iface+1]);

                //! Now reverse le and re.
                SWAP(left_cell_of_face[iface], right_cell_of_face[iface]);

                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else if (right_cell_of_face[iface] < 0)
            {
                -- left_cell_of_face[iface];
                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else
            {
                -- left_cell_of_face [iface];
                -- right_cell_of_face[iface];
            }
        }
    }
    delete [] nodePosi;

    if (GetDim() == TWO_D)
    {
        if (nBFaceReal != nBoundFace)
        {
            TK_Exit::ExceptionExit("nBFaceReal != nBFace\n");
        }
    }
    nBoundFace = nBFaceReal;

    cout << "set the boundary conditions......\n";
    grid->SetNBoundFace(nBoundFace);

    cout << "nBoundFace = " << nBoundFace << "\n";

    //! Set the boundary conditions.
    CompositeUnstructBC **bcr = new CompositeUnstructBC * [nBoundFace];
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int bctype = u2fbc.TranslateBC(-right_cell_of_face[iface]);
        bcr[iface] = new CompositeUnstructBC();
        bcr[iface]->SetKey(bctype);
    }

    grid->SetBCRecord(bcr);

    InterfaceInfo *interfaceInfo = 0;
    grid->SetInterfaceInfo(interfaceInfo);

    grid->SpecifyRightCellofBC();

    file.close();
}

//! Read ustar 3d grid by binary.
void ReadUstarGrid3D_Binary(const string &gridfile, UnstructGrid *grid)
{
    if (GetDim() == TWO_D)
    {
        TK_Exit::ExceptionExit("暂时没考虑2D的二进制ustar文件读入 ... \n");
    }

    ifstream infile;
    infile.open(gridfile.c_str(), ios_base::in|ios_base::binary);
    if (!infile)
    {
        TK_Exit::FileOpenErrorExit(gridfile);
    }

    cout << "reading the unstructed data file......\n";

    Ustar2Fantasy u2fbc("support/ustar2ftsbc.txt");
        
    int SizeInt = sizeof(int);
    int SizeDouble = sizeof(RDouble);

    int nTotalNode, nTotalFace, nTotalCell, nBoundFace;

    infile.read((char *)&nTotalNode, SizeInt);
    infile.read((char *)&nTotalFace, SizeInt);
    infile.read((char *)&nTotalCell, SizeInt);

    if (GetDim() == TWO_D)
    {

    }
    else
    {
        cout << " the total numbers of nodes ："   << nTotalNode;
        cout << " the total numbers of faces ：" << nTotalFace;
        cout << " the total numbers of cells ："   << nTotalCell << "\n";
    }

    grid->SetNTotalNode(nTotalNode);
    grid->SetNTotalFace(nTotalFace);
    grid->SetNTotalCell(nTotalCell);

    RDouble *x, *y, *z;
    x = new RDouble[nTotalNode];
    y = new RDouble[nTotalNode];
    z = new RDouble[nTotalNode];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    //! Read the coordinatesof point.
    int count = 3 * nTotalNode;
    RDouble *coordinate = new RDouble [count];
    infile.read((char *)coordinate, count * SizeDouble);

    count = 0;
    for (int i = 0; i < nTotalNode; ++ i)
    {
        x[i] = coordinate[count ++];
        y[i] = coordinate[count ++];
        z[i] = coordinate[count ++];
    }

    delete [] coordinate;

    //! Read the face to node relationship.
    int *node_number_of_each_face = new int[nTotalFace];
    grid->SetNodeNumberOfEachFace(node_number_of_each_face);

    if (GetDim() == TWO_D)
    {
        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            node_number_of_each_face[iface] = 2;
        }
    }
    else
    {
        infile.read((char *)node_number_of_each_face, nTotalFace * SizeInt);
    }
    
    int *nodePosi = new int[nTotalFace + 1];
    int nsum = 0;
    nodePosi[0] = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        nsum += node_number_of_each_face[iface];
        nodePosi[iface+1] = nsum;
    }

    cout << "set the connection relationship of face to node......\n";

    int *face2node = new int [nsum];
    grid->SetFace2Node(face2node);

    infile.read((char *)face2node, nsum * SizeInt);
    for (int j = 0; j < nsum; ++ j)
    {
        face2node[j] = face2node[j] - 1;
    }

    //! Read the face to cell relationship.
    int nBFaceReal = 0;
    int n_patch    = 0;

    cout << "set the connection relationship of face to cells......\n";
    int *left_cell_of_face  = new int [nTotalFace];
    int *right_cell_of_face = new int [nTotalFace];
    grid->SetLeftCellOfFace(left_cell_of_face);
    grid->SetRightCellOfFace(right_cell_of_face);

    if (GetDim() == TWO_D)
    {

    }
    else
    {
        count = 2 * nTotalFace;
        int *f2c_temp = new int [count];
        infile.read((char *)f2c_temp, count * SizeInt);

        for (int iface = 0; iface < nTotalFace; ++ iface)
        {
            left_cell_of_face [iface] = f2c_temp[2 * iface];
            right_cell_of_face[iface] = f2c_temp[2 * iface + 1];

            if (left_cell_of_face[iface] < 0)
            {
                -- right_cell_of_face[iface];
                std::reverse(face2node + nodePosi[iface], face2node + nodePosi[iface+1]);

                //! Now reverse le and re.
                SWAP(left_cell_of_face[iface], right_cell_of_face[iface]);

                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else if (right_cell_of_face[iface] < 0)
            {
                -- left_cell_of_face[iface];
                n_patch = MAX(n_patch, -right_cell_of_face[iface]);
                ++ nBFaceReal;
            }
            else
            {
                -- left_cell_of_face [iface];
                -- right_cell_of_face[iface];
            }
        }

        delete [] f2c_temp;
    }
    delete [] nodePosi;

    nBoundFace = nBFaceReal;

    cout << "设置边界条件......\n";
    grid->SetNBoundFace(nBoundFace);

    cout << "nBoundFace = " << nBoundFace << "\n";

    //! Set the boundary conditions.
    CompositeUnstructBC **bcr = new CompositeUnstructBC * [nBoundFace];
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int bctype = u2fbc.TranslateBC(-right_cell_of_face[iface]);
        bcr[iface] = new CompositeUnstructBC();
        bcr[iface]->SetKey(bctype);
    }

    grid->SetBCRecord(bcr);

    InterfaceInfo *interfaceInfo = 0;
    grid->SetInterfaceInfo(interfaceInfo);

    grid->SpecifyRightCellofBC();

    infile.close();
    infile.clear();
}

void GetUnsIndex(StructGrid *grid, DataStruct_AdtTree<int, RDouble> *adtTree, int &count, RDouble tol, Int3D &unsindex, RDouble *xuns, RDouble *yuns, RDouble *zuns)
{
    typedef DataStruct_AdtTree<int, RDouble> AdtTree;
    typedef DataStruct_AdtNode<int, RDouble> AdtNode;
    AdtTree::AdtNodeList node_list;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &xs = * grid->GetStructX();
    RDouble3D &ys = * grid->GetStructY();
    RDouble3D &zs = * grid->GetStructZ();

    unsindex = 0;

    RDouble coor[3], minwin[3], maxwin[3];

    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                coor[0] = xs(i, j, k);
                coor[1] = ys(i, j, k);
                coor[2] = zs(i, j, k);

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                node_list.resize(0);
                adtTree->FindNodesInRegion(minwin, maxwin, node_list);

                if (node_list.size() == 0)
                {
                    //! This means that there is no point,can add them.
                    xuns[count] = coor[0];
                    yuns[count] = coor[1];
                    zuns[count] = coor[2];

                    AdtNode *node = new AdtNode(3, coor, count);
                    adtTree->AddNode(node);
                    unsindex(i, j, k) = count;
                    ++ count;
                }
                else
                {
                    if (node_list.size() > 1)
                    {
                        ostringstream oss;
                        oss << "FATAL ERROR : node_list.size() = " << node_list.size() << "\n";
                        TK_Exit::ExceptionExit(oss.str());
                    }
                    AdtNode *node =  node_list[0];
                    unsindex(i, j, k) = node->GetData();
                }
            }
        }
    }
}

void SetUnstructBC
(
    int ist,
    int ied,
    int jst,
    int jed,
    int kst,
    int ked,
    int bctype,
    cgsize_t *conn_list,
    int &start_pos,
    Int3D &unsindex
)
{
    cout << "ist,ied,jst,jed,kst,ked = " << ist << " " << ied << " " << jst << " " << jed << " " << kst << " " << ked << "\n";
    int nDim = GetDim();
    int numpt = 4;
    if (nDim == TWO_D) numpt = 2;

    if (ist == ied)
    {
        int i = ist;
        if (nDim == THREE_D)
        {
            for (int k = kst; k <= ked - 1; ++ k)
            {
                for (int j = jst; j <= jed - 1; ++ j)
                {
                    if (i == 1)
                    {
                        conn_list[start_pos + 0] = unsindex(i, j  , k ) + 1;
                        conn_list[start_pos + 1] = unsindex(i, j  , k+1) + 1;
                        conn_list[start_pos + 2] = unsindex(i, j+1, k+1) + 1;
                        conn_list[start_pos + 3] = unsindex(i, j+1, k ) + 1;
                    }
                    else
                    {
                        conn_list[start_pos + 0] = unsindex(i, j  , k ) + 1;
                        conn_list[start_pos + 1] = unsindex(i, j+1, k ) + 1;
                        conn_list[start_pos + 2] = unsindex(i, j+1, k+1) + 1;
                        conn_list[start_pos + 3] = unsindex(i, j  , k+1) + 1;
                    }
                    start_pos += numpt;
                }
            }
        }
        else
        {
            int k = kst;
            for (int j = jst; j <= jed - 1; ++ j)
            {
                if (i == 1)
                {
                    conn_list[start_pos + 0] = unsindex(i, j+1, k) + 1;
                    conn_list[start_pos + 1] = unsindex(i, j  , k) + 1;
                }
                else
                {
                    conn_list[start_pos + 0] = unsindex(i, j  , k) + 1;
                    conn_list[start_pos + 1] = unsindex(i, j+1, k) + 1;
                }
                start_pos += numpt;
            }
        }
        return;
    }

    if (jst == jed)
    {
        int j = jst;
        if (nDim == THREE_D)
        {
            for (int k = kst; k <= ked - 1; ++ k)
            {
                for (int i = ist; i <= ied - 1; ++ i)
                {
                    if (j == 1)
                    {
                        conn_list[start_pos + 0] = unsindex(i  , j, k ) + 1;
                        conn_list[start_pos + 1] = unsindex(i+1, j, k ) + 1;
                        conn_list[start_pos + 2] = unsindex(i+1, j, k+1) + 1;
                        conn_list[start_pos + 3] = unsindex(i  , j, k+1) + 1;
                    }
                    else
                    {
                        conn_list[start_pos + 0] = unsindex(i  , j, k ) + 1;
                        conn_list[start_pos + 1] = unsindex(i  , j, k+1) + 1;
                        conn_list[start_pos + 2] = unsindex(i+1, j, k+1) + 1;
                        conn_list[start_pos + 3] = unsindex(i+1, j, k ) + 1;
                    }
                    start_pos += numpt;
                }
            }
        }
        else
        {
            int k = kst;
            for (int i = ist; i <= ied - 1; ++ i)
            {
                if (j == 1)
                {
                    conn_list[start_pos + 0] = unsindex(i  , j, k) + 1;
                    conn_list[start_pos + 1] = unsindex(i+1, j, k) + 1;
                }   
                else
                {
                    conn_list[start_pos + 0] = unsindex(i+1, j, k) + 1;
                    conn_list[start_pos + 1] = unsindex(i  , j, k) + 1;
                }
                start_pos += numpt;
            }
        }
        return;
    }

    if (kst == ked)
    {
        int k = kst;
        for (int j = jst; j <= jed - 1; ++ j)
        {
            for (int i = ist; i <= ied - 1; ++ i)
            {
                if (k == 1)
                {
                    conn_list[start_pos + 0] = unsindex(i  , j  , k) + 1;
                    conn_list[start_pos + 1] = unsindex(i  , j+1, k) + 1;
                    conn_list[start_pos + 2] = unsindex(i+1, j+1, k) + 1;
                    conn_list[start_pos + 3] = unsindex(i+1, j  , k) + 1;
                }
                else
                {
                    conn_list[start_pos + 0] = unsindex(i  , j  , k) + 1;
                    conn_list[start_pos + 1] = unsindex(i+1, j  , k) + 1;
                    conn_list[start_pos + 2] = unsindex(i+1, j+1, k) + 1;
                    conn_list[start_pos + 3] = unsindex(i  , j+1, k) + 1;
                }
                start_pos += numpt;
            }
        }
        return;
    }

    TK_Exit::ExceptionExit(" error : ist != ied, jst != jed, kst != ked \n");
}

void FillSection(Grid **grids, int nBlocks, CGNSBase *base_cgns, BaseElement *base_elem, Int3D **unsindexlist)
{
    int nTotalCell, nBoundFace, nSections, nBCRegions;
    int istart[2], iend[2];
    ElementType_t etype[2];

    cgsize_t *conn_list;
    int start_pos;

    int nDim = GetDim();

    nTotalCell = 0;
    nBoundFace = 0;
    nBCRegions = 0;

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni-1);
        Range J(1, nj-1);
        Range K(1, nk-1);
        if (nk == 1) K.setRange(1, 1);
        nTotalCell += I.length() * J.length() * K.length();

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
        int nBCRegion = composite_bcregion->GetnBCRegion();

        for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
        {
            //! Read the source region imin,imax,jmin,jmax,kmin,kmax,bctype.
            int imin, imax, jmin, jmax, kmin, kmax, bctype;
            StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
            bcregion->GetNormalizeIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            bctype = bcregion->GetBCType();
            if (bctype != 7 || bctype / 10 != 7)
            {
                int di = MAX(imax - imin, 1);
                int dj = MAX(jmax - jmin, 1);
                int dk = MAX(kmax - kmin, 1);
                nBoundFace += di * dj * dk;
                nBCRegions ++;
            }
        }
    }

    base_cgns->SetNTotalCell(nTotalCell);
    base_cgns->SetNBoundFace(nBoundFace);

    nSections = 2;

    base_cgns->CreateElement(nSections);

    istart[0] = 1;
    iend[0] = nTotalCell;

    istart[1] = nTotalCell + 1;
    iend[1] = nTotalCell + 1 + nBoundFace;

    if (nDim == THREE_D)
    {
        etype[0] = HEXA_8;
        etype[1] = QUAD_4;
    }
    else
    {
        etype[0] = QUAD_4;
        etype[1] = BAR_2;
    }

    cgsize_t *base_cgns_istart = base_cgns->GetIndexOfStart();
    cgsize_t *base_cgns_iend   = base_cgns->GetIndexOfEnd();
    int *base_cgns_elem_type = base_cgns->GetElementType();

    cgsize_t **base_cgns_conn_list = base_cgns->GetElementConnectionList();

    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    for (int sec = 1; sec <= nSections; ++ sec)
    {
        base_cgns_istart[sec-1] = istart[sec-1];
        base_cgns_iend  [sec-1] = iend  [sec-1];
        base_cgns_elem_type[sec-1] = etype[sec-1];
    }

    for (int sec = 1; sec <= nSections; ++ sec)
    {
        int nelem = iend[sec-1] - istart[sec-1] + 1;
        int ptnum = base_elem_elem_pt[etype[sec-1]];
        base_cgns_conn_list[sec-1] = new cgsize_t[ptnum * nelem];
    }

    base_cgns->CreateBC(nBCRegions);

    conn_list = base_cgns_conn_list[0];
    start_pos = 0;
    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);

        Int3D &unsindex = *unsindexlist[izone];

        int ist, ied, jst, jed, kst, ked;
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

        int il1 = 1;
        int jl1 = 1;
        int kl1 = 1;

        if (nDim == TWO_D) kl1 = 0;

        int numpt = base_elem_elem_pt[etype[0]];
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    conn_list[start_pos + 0] = unsindex(i    , j    , k) + 1;
                    conn_list[start_pos + 1] = unsindex(i+il1, j    , k) + 1;
                    conn_list[start_pos + 2] = unsindex(i+il1, j+jl1, k) + 1;
                    conn_list[start_pos + 3] = unsindex(i    , j+jl1, k) + 1;
                    if (nDim == THREE_D)
                    {
                        conn_list[start_pos + 4] = unsindex(i    , j    , k+kl1) + 1;
                        conn_list[start_pos + 5] = unsindex(i+il1, j    , k+kl1) + 1;
                        conn_list[start_pos + 6] = unsindex(i+il1, j+jl1, k+kl1) + 1;
                        conn_list[start_pos + 7] = unsindex(i    , j+jl1, k+kl1) + 1;
                    }
                    start_pos += numpt;
                }
            }
        }
    }

    conn_list = base_cgns_conn_list[1];
    start_pos = 0;
    int boco  = 0;
    int ielem = nTotalCell;

    BaseBCType *bc_convert = new BaseBCType();

    int *base_cgns_nBCElem          = base_cgns->GetNumberOfBCElements();
    int *base_cgns_bc_elem_set_type = base_cgns->GetBoundaryElementType();
    int *base_cgns_bc_grid_location = base_cgns->GetBoundaryGridLocation();
    int *base_cgns_bc_type          = base_cgns->GetBCType();
    cgsize_t **base_cgns_bc_conn_list = base_cgns->GetBoundaryElementConnectionList();
    string *base_cgns_bc_name = base_cgns->GetBCName();

    cout << "    Warning: the cgsize_t type of conn_list has been convert to int type, when StrGridToUnsGrid()!!!" << endl;

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        Int3D &unsindex = *unsindexlist[izone];

        CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
        int nBCRegion = composite_bcregion->GetnBCRegion();
        for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
        {
            //! Read the source region imin,imax,jmin,jmax,kmin,kmax,bctype.
            int imin, imax, jmin, jmax, kmin, kmax, bctype;
            StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
             
            bcregion->GetNormalizeIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            bctype = bcregion->GetBCType();

            if (bctype != 7 || bctype / 10 != 7)
            {
                boco ++;
                int nBCElem = 2;

                int di = MAX(imax - imin, 1);
                int dj = MAX(jmax - jmin, 1);
                int dk = MAX(kmax - kmin, 1);

                base_cgns_nBCElem         [boco-1] = nBCElem;
                base_cgns_bc_conn_list    [boco-1] = new cgsize_t[nBCElem];
                base_cgns_bc_type         [boco-1] = bc_convert->Fantasy2CGNS(bctype);
                base_cgns_bc_elem_set_type[boco-1] = ElementRange;
                base_cgns_bc_grid_location[boco-1] = CellCenter;
                base_cgns_bc_name         [boco-1] = bcregion->GetBCName();
                base_cgns_bc_conn_list[boco-1][0]  = ielem + 1;
                base_cgns_bc_conn_list[boco-1][1]  = ielem + di * dj * dk;
                ielem += di * dj * dk;

                SetUnstructBC(imin, imax, jmin, jmax, kmin, kmax, bctype, conn_list, start_pos, unsindex);
            }
        }
    }

    delete bc_convert;
}

void StrGridToUnsGrid(Grid **grids, int nBlocks, UnstructGrid *uns_grid)
{
    using namespace PHSPACE;
    RDouble tol = 1.0e-8;

    RDouble mindis,maxdis;
    mindis =   LARGE;
    maxdis = - LARGE;

    int nTotalNode = 0;
    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble dismin, dismax;

        grid->GetMinMaxDS(dismin, dismax);

        mindis = MIN(mindis, dismin);
        maxdis = MAX(maxdis, dismax);
        nTotalNode += ni * nj * nk;
    }

    cout << " mindis = " << mindis << " maxdis = " << maxdis << "\n";

    tol = mindis / 4;

    //! Get the min max coordinates of the grid.
    RDouble pmin[3], pmax[3];

    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        RDouble *ptmin = grid->GetMinBox();
        RDouble *ptmax = grid->GetMaxBox();

        pmin[0] = MIN(pmin[0], ptmin[0]);
        pmin[1] = MIN(pmin[1], ptmin[1]);
        pmin[2] = MIN(pmin[2], ptmin[2]);

        pmax[0] = MAX(pmax[0], ptmax[0]);
        pmax[1] = MAX(pmax[1], ptmax[1]);
        pmax[2] = MAX(pmax[2], ptmax[2]);
    }

    ShiftMinMaxBox(pmin, pmax, two * tol);

    typedef DataStruct_AdtTree<int, RDouble> AdtTree;
    AdtTree::AdtNodeList nodeList;

    AdtTree *adtTree = new AdtTree(3, pmin, pmax);

    int count = 0;
    typedef DataStruct_AdtNode<int, RDouble> AdtNode;

    RDouble *xuns = new RDouble[nTotalNode];
    RDouble *yuns = new RDouble[nTotalNode];
    RDouble *zuns = new RDouble[nTotalNode];

    int nTotalCell = 0;
    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni-1);
        Range J(1, nj-1);
        Range K(1, nk-1);
        if (nk == 1) K.setRange(1, 1);
        nTotalCell += I.length() * J.length() * K.length();
    }

    Int3D **unsindexlist = new Int3D * [nBlocks];

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        StructGrid *grid = StructGridCast(grids[izone]);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        unsindexlist[izone] = new Int3D(Range(1, ni), Range(1, nj), Range(1, nk), fortranArray);
        Int3D &unsindex = *unsindexlist[izone];
        cout << " block = " << izone + 1 << "\n";
        GetUnsIndex(grid, adtTree, count, tol, unsindex, xuns, yuns, zuns);
    }

    delete adtTree;
    nTotalNode = count;

    cout << "first nTotalNode = " << nTotalNode << "\n";
    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        x[inode] = xuns[inode];
        y[inode] = yuns[inode];
        z[inode] = zuns[inode];
    }

    delete [] xuns;
    delete [] yuns;
    delete [] zuns;

    int nSingle_Blocks = 1;
    CGNSFactory *factory_cgns = new CGNSFactory(nSingle_Blocks);

    CGNSBase    *base_cgns = factory_cgns->GetCGNSBase(0);
    RawGrid     *rawgrid   = factory_cgns->GetRawGrid(0);
    BaseElement *base_elem = factory_cgns->GetBaseElement();

    rawgrid->SetNTotalNode(nTotalNode);
    rawgrid->SetX(x);
    rawgrid->SetY(y);
    rawgrid->SetZ(z);

    base_cgns->SetNTotalNode(nTotalNode);
    base_cgns->SetNTotalCell(nTotalCell);

    FillSection(grids, nBlocks, base_cgns, base_elem, unsindexlist);

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        delete unsindexlist[izone];
    }
    delete [] unsindexlist;

    factory_cgns->ConvertGrid2Fantasy();
    CGNS2UnsGrid(factory_cgns, 0, uns_grid);

    delete factory_cgns;
}

void StrGridToUnsGrid(StructGrid *grid, UnstructGrid *uns_grid)
{
    int nBlocks = 1;
    Grid **grids = new Grid *[nBlocks];
    grids[0] = grid;

    StrGridToUnsGrid(grids, nBlocks, uns_grid);

    delete [] grids;
}

void SetBcRegion(StructGrid *grid, VInt &imin, VInt &imax, VInt &jmin, VInt &jmax, VInt &kmin, VInt &kmax, VInt &bctype)
{
    //! The situation of multi block patched boundary is not considered temporary.
    int nBCRegion = static_cast<int>(imin.size());

    grid->CreateCompositeBCRegion(nBCRegion);
    CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();

    int nb = 0;
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = new StructBC(nb, ibcregion);
        composite_bcregion->SetBCRegion(ibcregion, bcregion);
        bcregion->SetIJKRegion(imin[ibcregion], imax[ibcregion], jmin[ibcregion], jmax[ibcregion], kmin[ibcregion], kmax[ibcregion]);
        bcregion->SetBCType(bctype[ibcregion]);
    }
}

}
