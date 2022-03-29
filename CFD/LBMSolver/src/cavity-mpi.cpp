
#include "CommonLibDll.h"
#include "D2Q9.h"
#include "D3Q15.h"
#include "TimeMarching.h"
#include "LatticeMPI.h"
#include "BoundaryFactory.h"
#include "GenBcMap.h"
#include "StlGeo.h"
#include <iostream> 

#include <GlobalDataBase.h>
#include <Region.h>
#include <Zone.h>
#include <Geometry.h>
#include <Geo_Grid.h>

#include <PHLBM.h>
#include <PHLBMGenerator.h>
#include <PHLBMGridGen.h>
#include <PHLBMRangeGen.h>
#include <PHLBMUnstructGrid.h>
#include "TK_Log.h"


#include <mpi.h> 

void test_linked_PdfCommunicator() {
  typedef nnwlbm::D2Q9 LC;
  constexpr const size_t dim = nnwlbm::Dim<2>::GetDim();
  //constexpr const size_t dim = 2;

  auto comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int commSize = 0;
  MPI_Comm_size(comm, &commSize);

  nnwlbm::Range<dim> r0(0,0, 2,3);
  nnwlbm::Range<dim> r1(2,0, 5,3);
  nnwlbm::Range<dim> r2(0,3, 3,5);
  nnwlbm::Range<dim> r3(3,3, 5,5);

  std::shared_ptr<RangeNeighbors<dim>> rangeNeighbors;
  switch(rank) {
  case 0:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r0);
    break;
  case 1:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r1);
    break;
  case 2:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r2);
    break;
  case 3:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r3);
    break;
  }

  rangeNeighbors->AddNeiborCandidate(r0);
  rangeNeighbors->AddNeiborCandidate(r1);
  rangeNeighbors->AddNeiborCandidate(r2);
  rangeNeighbors->AddNeiborCandidate(r3);

  std::map<nnwlbm::Range<dim>, int> rangeToRank;
  rangeToRank.emplace(r0, 0);
  rangeToRank.emplace(r1, 1);
  rangeToRank.emplace(r2, 2);
  rangeToRank.emplace(r3, 3);

  auto hostRange = rangeNeighbors->GetHostRange();

  hostRange.AddBuffer(1);
  nnwlbm::LinkedMaskRange<int, dim> mRange(hostRange, 0);
  nnwlbm::LinkedMaskData<double, dim> pdfData(nnwlbm::D2Q9::size(), mRange, nnwlbm::Color());

  int loc[dim];
  auto rgIndex = pdfData.GetNodeIndex();
  for(int carry = hostRange.IniCarry(loc); 
    carry == 0; carry = hostRange.AddOne(loc))
  {
    size_t nodeId = rgIndex.CalNodeID(loc);
    double* pdf = pdfData.GetVar(nodeId);
    for(size_t i = 0; i < pdfData.GetVarNum(); i++){
      pdf[i] = rank*100+double(i);
    }
  }

  PdfCommunicator<LC, nnwlbm::LinkedMaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);

  double origin[] = {0.0, 0.0};
  nnwlbm::UniformGrid<dim> grid(1.0, origin);
  std::string filename = "rank"+std::to_string(rank)+"before";
  nnwlbm::OutputCellCenteredDataHalfway(pdfData, grid, filename);

  pdfComm.SendAndRecv();
  pdfComm.wait();
  while(!pdfComm.update()) {
    std::cout << "update wrong " << rank << std::endl;
  }

  filename = "rank"+std::to_string(rank)+"after";
  nnwlbm::OutputCellCenteredDataHalfway(pdfData, grid, filename);
}

void test_PdfCommunicator() {
  typedef nnwlbm::D2Q9 LC;
  constexpr const size_t dim = nnwlbm::Dim<2>::GetDim();
  //constexpr const size_t dim = 2;

  auto comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int commSize = 0;
  MPI_Comm_size(comm, &commSize);
  /*if(commSize != 4) {
    if(rank == 0) {
     std::cout << "4 processes are required." << std::endl;
    }
    MPI_Abort(comm, 1);
  }*/

  nnwlbm::Range<dim> r0(0,0, 2,3);
  nnwlbm::Range<dim> r1(2,0, 5,3);
  nnwlbm::Range<dim> r2(0,3, 3,5);
  nnwlbm::Range<dim> r3(3,3, 5,5);

  std::shared_ptr<RangeNeighbors<dim>> rangeNeighbors;
  switch(rank) {
  case 0:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r0);
    break;
  case 1:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r1);
    break;
  case 2:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r2);
    break;
  case 3:
    rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(r3);
    break;
  }

  rangeNeighbors->AddNeiborCandidate(r0);
  rangeNeighbors->AddNeiborCandidate(r1);
  rangeNeighbors->AddNeiborCandidate(r2);
  rangeNeighbors->AddNeiborCandidate(r3);

  std::map<nnwlbm::Range<dim>, int> rangeToRank;
  rangeToRank.emplace(r0, 0);
  rangeToRank.emplace(r1, 1);
  rangeToRank.emplace(r2, 2);
  rangeToRank.emplace(r3, 3);

  auto hostRange = rangeNeighbors->GetHostRange();

  hostRange.AddBuffer(1);
  nnwlbm::MaskRange<int, dim> mRange(hostRange, 0);
  nnwlbm::MaskData<double, dim> pdfData(nnwlbm::D2Q9::size(), mRange, nnwlbm::Color());

  int loc[dim];
  auto rgIndex = pdfData.GetNodeIndex();
  for(int carry = hostRange.IniCarry(loc); 
    carry == 0; carry = hostRange.AddOne(loc))
  {
    size_t nodeId = rgIndex.CalNodeID(loc);
    double* pdf = pdfData.GetVar(nodeId);
    for(size_t i = 0; i < pdfData.GetVarNum(); i++){
      pdf[i] = rank*100+double(i);
    }
  }

  PdfCommunicator<LC, nnwlbm::MaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);

  double origin[] = {0.0, 0.0};
  nnwlbm::UniformGrid<dim> grid(1.0, origin);
  std::string filename = "rank"+std::to_string(rank)+"before";
  nnwlbm::OutputCellCenteredDataHalfway(pdfData, grid, filename);

  pdfComm.SendAndRecv();
  pdfComm.wait();
  pdfComm.update();

  filename = "rank"+std::to_string(rank)+"after";
  nnwlbm::OutputCellCenteredDataHalfway(pdfData, grid, filename);

}

std::pair<int, int> EqualFactoring(int in) {
  int max = (int)sqrt(in) + 1;
  if(max*max == in) return std::make_pair(max, max);
  for(int f1 = max-1; f1 > 0; f1--)
  {
    if(in%f1 == 0) return std::make_pair(f1, in/f1);
  }
  return std::make_pair(1, in);
}

//void cavity2d_mpi(){
//  // 物理参数设置
//  typedef nnwlbm::D2Q9 Ltc;  // 设置格子构型
//  Ltc lc(1.0);  // 生成格子，指定格子速度
//  const size_t dim = Ltc::GetDim();  // 获得维度信息
//  double wVel[] = {0.1, 0.0}; // 设置上板速度
//  double viscosity = 0.0001;  // 设置粘性
//  std::vector<double> iniMacroVar(dim+1); // 设置初条件
//  iniMacroVar[0] = 1.0; // 密度
//  iniMacroVar[1] = 0.0; // u速度
//  iniMacroVar[2] = 0.0; // v速度
//
//                        // 计算参数设置
//  const int num = 180; //  设置每条边的网格数
//  const double dx = 1.0/num;  // 计算网格大小
//  const double dt = dx/lc.velocity(); // 计算时间步长
//  double tau = nnwlbm::CalTau(lc, viscosity, dt, iniMacroVar[0]);
//
//  // BGK or MRT
//  typedef nnwlbm::Equilibrium<Ltc> EQ;
//  EQ equilibrium;
//  //typedef nnwlbm::BGK<Ltc, nnwlbm::Equilibrium<Ltc>> CLS;
//  //typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, nnwlbm::Equilibrium<Ltc>> CLS;
//  typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, EQ> CLS;
//  CLS kernel(lc, tau);
//
//
//  // 计算区域的几何形状
//  auto rec = std::make_shared<nnwlbm::Box>(dim);
//  (*rec)[0] = 0.0;
//  (*rec)[1] = 0.0;
//  (*rec)[2] = 1.0;
//  (*rec)[3] = 1.0;
//  // 外部形状
//  auto outsideBoundary = std::make_shared<nnwlbm::ComplementGeo>(rec);
//  // 上壁面形状
//  auto topWall = std::make_shared<nnwlbm::Box>(dim);
//  (*topWall)[0] = 0.0;
//  (*topWall)[1] = 1.0-1e-7;
//  (*topWall)[2] = 1.0;
//  (*topWall)[3] = 1.0+1e-7;
//
//  // 定义无滑移反弹边条件
//  auto bounceBack = std::make_shared<nnwlbm::BounceBack<Ltc>>();
//  // 定义几何体表示外部区域
//  std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<Ltc>*> bcBodyMap;
//  // 指定外部区域的边界为无滑移边条件
//  bcBodyMap.emplace(outsideBoundary.get(), bounceBack.get());
//
//  // 定义无滑移移动壁面反弹边条件
//  auto movingBounceBack = std::make_shared<nnwlbm::MovingBounceBack<Ltc>>(lc.velocity(), wVel);
//  // 定义几何体表示边界
//  std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<Ltc>*> bcSurfaceMap;
//  // 指定上壁面为无滑移边条件
//  bcSurfaceMap.emplace(topWall.get(), movingBounceBack.get());
//
//  // 初始化流场数据
//  nnwlbm::ConstFun initcond(std::make_shared<nnwlbm::UniversalGeo>(lc.GetDim()), iniMacroVar);
//
//  // 生成网格
//  double origin[2] = { 0.5*dx, 0.5*dx };
//  nnwlbm::UniformGrid<dim> grid(dx, origin);
//  origin[0] = origin[1] = 0.0; 
//  nnwlbm::UniformGrid<dim> gridC(dx, origin);
//
//  // 指定计算区域
//  std::vector<nnwlbm::Range<dim>> rangeVec;
//  auto comm = MPI_COMM_WORLD;
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  int commSize = 0;
//  MPI_Comm_size(comm, &commSize);
//  rangeVec.resize(commSize);
//
//  // partition
//  //int remainder = num%commSize;
//  //int quotient = num/commSize;
//  //for(int i = 0; i < remainder; i++) {
//  //  rangeVec[i] = nnwlbm::Range<dim>(i*(quotient+1), 0, (i+1)*(quotient+1), num);
//  //}
//  //for(int i = remainder; i < commSize; i++) {
//  //  rangeVec[i] = nnwlbm::Range<dim>(
//  //    remainder*(quotient+1)+(i-remainder)*quotient, 0, 
//  //    remainder*(quotient+1)+(i+1-remainder)*quotient, num
//  //    );
//  //}
//
//  // partition 2
//  if(commSize != 5) exit(1);
//  int nodeNumOfBlock = num*num/5;
//  int lengthOfCentralBlock = int(num/sqrt(5.0));
//  int widthOfOtherFourBlock = (num-lengthOfCentralBlock)/2;
//  int clb = widthOfOtherFourBlock;
//  int crt = widthOfOtherFourBlock+lengthOfCentralBlock;
//  
//  rangeVec[0] = nnwlbm::Range<dim>(0, 0, crt, clb);
//  rangeVec[1] = nnwlbm::Range<dim>(crt, 0, num, crt);
//  rangeVec[2] = nnwlbm::Range<dim>(clb, crt, num, num);
//  rangeVec[3] = nnwlbm::Range<dim>(0, clb, clb, num);
//  rangeVec[4] = nnwlbm::Range<dim>(clb, clb, crt, crt);
//
//  auto rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(rangeVec[rank]);
//  std::map<nnwlbm::Range<dim>, int> rangeToRank;
//
//  for(int i = 0; i < commSize; i++) {
//    rangeNeighbors->AddNeiborCandidate(rangeVec[i]);
//    rangeToRank.emplace(rangeVec[i], i);
//  }
//  auto hostRange = rangeNeighbors->GetHostRange();
//
//  hostRange.AddBuffer(1);
//  nnwlbm::LinkedMaskRange<int, dim> maskRange(hostRange, 0);
//  nnwlbm::LinkedMaskData<double, dim> maskData(3, maskRange);
//  nnwlbm::LinkedMaskData<double, dim> pdfData(nnwlbm::D2Q9::size(), maskRange, nnwlbm::Color());
//
//  PdfCommunicator<Ltc, nnwlbm::LinkedMaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);
//
//  //////////////////////////////////////////////////////////////////////////////////////
//  // 生成colorMap
//  nnwlbm::GenColorMap<dim, nnwlbm::SingleNodeBc<Ltc>>(maskRange, grid, bcBodyMap );
//
//  // 根据colorMap为分布函数和边界条件分配内存
//  //nnwlbm::LinkedMaskData<double, dim> pdfData(lc.size(), maskRange);
//  nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<Ltc>*, Ltc::GetDim(), nnwlbm::BcIdentifier > 
//    bcMapSingle(Ltc::size(), maskRange,  nnwlbm::BcIdentifier(nnwlbm::geo::boundary));
//
//  // 生成边界条件的映射
//  nnwlbm::gen_bcMap(bcMapSingle, grid, bcBodyMap, bcSurfaceMap);
//
//  //// 初始化流场
//  ini_field(maskData, initcond, grid);
//
//  // 初始化分布函数
//  auto loopRange = maskRange.GetRange();
//  loopRange.AddBuffer(-1);
//  nnwlbm::CalFieldPdf<Ltc>(pdfData, maskData, loopRange);
//
//  // 主循环 
//  CollisionAndStream<Ltc, CLS, EQ> mainLoop(lc, kernel, equilibrium, pdfData, maskData, bcMapSingle, gridC, pdfComm);
//  nnwlbm::Output<nnwlbm::LinkedMaskData<double, dim>, nnwlbm::UniformGrid<dim>> output(maskData, gridC);
//
//  time_marching_mpi(output, mainLoop, 100000, 10000, 10);
//}
//
//void cavity3d_mpi(){
//  // 物理参数设置
//  typedef nnwlbm::D3Q15 Ltc;  // 设置格子构型
//  Ltc lc(1.0);  // 生成格子，指定格子速度
//  const size_t dim = Ltc::GetDim();  // 获得维度信息
//  double wVel[] = {0.1, 0.0, 0.0}; // 设置上板速度
//  double viscosity = 0.0001;  // 设置粘性
//  std::vector<double> iniMacroVar(dim+1); // 设置初条件
//  iniMacroVar[0] = 1.0; // 密度
//  iniMacroVar[1] = 0.0; // u速度
//  iniMacroVar[2] = 0.0; // v速度
//  iniMacroVar[3] = 0.0; // w速度
//
//                        // 计算参数设置
//  const int num = 128; //  设置每条边的网格数
//  const double dx = 1.0/num;  // 计算网格大小
//  const double dt = dx/lc.velocity(); // 计算时间步长
//  double tau = nnwlbm::CalTau(lc, viscosity, dt, iniMacroVar[0]);
//
//  // BGK or MRT
//  typedef nnwlbm::Equilibrium<Ltc> EQ;
//  EQ equilibrium;
//  //typedef nnwlbm::BGK<Ltc, nnwlbm::Equilibrium<Ltc>> CLS;
//  typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, nnwlbm::Equilibrium<Ltc>> CLS;
//  CLS kernel(lc, tau);
//
//  // 计算区域的几何形状
//  auto rec = std::make_shared<nnwlbm::Box>(dim);
//  (*rec)[0] = (*rec)[1] = (*rec)[2] = 0.0;
//  (*rec)[3] = (*rec)[4] = (*rec)[5] = 1.0;
//  // 外部形状
//  auto outsideBoundary = std::make_shared<nnwlbm::ComplementGeo>(rec);
//  // 上壁面形状
//  auto topWall = std::make_shared<nnwlbm::Box>(dim);
//  (*topWall)[0] = 0.0-1e-5;
//  (*topWall)[1] = 0.0-1e-5;
//  (*topWall)[2] = 1.0-1e-5;
//  (*topWall)[3] = 1.0+1e-5;
//  (*topWall)[4] = 1.0+1e-5;
//  (*topWall)[5] = 1.0+1e-5;
//
//  // 定义无滑移反弹边条件
//  auto bounceBack = std::make_shared<nnwlbm::BounceBack<Ltc>>();
//  // 定义几何体表示外部区域
//  std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<Ltc>*> bcBodyMap;
//  // 指定外部区域的边界为无滑移边条件
//  bcBodyMap.emplace(outsideBoundary.get(), bounceBack.get());
//
//  // 定义无滑移移动壁面反弹边条件
//  auto movingBounceBack = std::make_shared<nnwlbm::MovingBounceBack<Ltc>>(lc.velocity(), wVel);
//  // 定义几何体表示边界
//  std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<Ltc>*> bcSurfaceMap;
//  // 指定上壁面为无滑移边条件
//  bcSurfaceMap.emplace(topWall.get(), movingBounceBack.get());
//
//  // 初始化流场数据
//  nnwlbm::ConstFun initcond(std::make_shared<nnwlbm::UniversalGeo>(lc.GetDim()), iniMacroVar);
//
//  // 生成网格
//  double origin[3] = { 0.5*dx, 0.5*dx, 0.5*dx };
//  nnwlbm::UniformGrid<dim> grid(dx, origin);
//  origin[0] = origin[1] = origin[2] = 0.0; 
//  nnwlbm::UniformGrid<dim> gridC(dx, origin);
//
//  // 指定计算区域
//  std::vector<nnwlbm::Range<dim>> rangeVec;
//  auto comm = MPI_COMM_WORLD;
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  int commSize = 0;
//  MPI_Comm_size(comm, &commSize);
//  rangeVec.resize(commSize);
//
//  // partition
//  int remainder = num%commSize;
//  int quotient = num/commSize;
//  for(int i = 0; i < remainder; i++) {
//    rangeVec[i] = nnwlbm::Range<dim>(i*(quotient+1), 0, 0, (i+1)*(quotient+1), num, num);
//  }
//  for(int i = remainder; i < commSize; i++) {
//    rangeVec[i] = nnwlbm::Range<dim>(
//      remainder*(quotient+1)+(i-remainder)*quotient, 0, 0,
//      remainder*(quotient+1)+(i+1-remainder)*quotient, num, num
//      );
//  }
//
//  auto rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(rangeVec[rank]);
//  std::map<nnwlbm::Range<dim>, int> rangeToRank;
//
//  for(int i = 0; i < commSize; i++) {
//    rangeNeighbors->AddNeiborCandidate(rangeVec[i]);
//    rangeToRank.emplace(rangeVec[i], i);
//  }
//  auto hostRange = rangeNeighbors->GetHostRange();
//
//  hostRange.AddBuffer(1);
//  nnwlbm::LinkedMaskRange<int, dim> maskRange(hostRange, 0);
//  nnwlbm::LinkedMaskData<double, dim> maskData(4, maskRange);
//  nnwlbm::LinkedMaskData<double, dim> pdfData(Ltc::size(), maskRange, nnwlbm::Color());
//
//  PdfCommunicator<Ltc, nnwlbm::LinkedMaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);
//
//  //////////////////////////////////////////////////////////////////////////////////////
//  // 生成colorMap
//  nnwlbm::GenColorMap<dim, nnwlbm::SingleNodeBc<Ltc>>(maskRange, grid, bcBodyMap );
//
//  // 根据colorMap为分布函数和边界条件分配内存
//  //nnwlbm::LinkedMaskData<double, dim> pdfData(lc.size(), maskRange);
//  nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<Ltc>*, Ltc::GetDim(), nnwlbm::BcIdentifier > 
//    bcMapSingle(Ltc::size(), maskRange,  nnwlbm::BcIdentifier(nnwlbm::geo::boundary));
//
//  // 生成边界条件的映射
//  nnwlbm::gen_bcMap(bcMapSingle, grid, bcBodyMap, bcSurfaceMap);
//
//  //// 初始化流场
//  ini_field(maskData, initcond, grid);
//
//  // 初始化分布函数
//  auto loopRange = maskRange.GetRange();
//  loopRange.AddBuffer(-1);
//  nnwlbm::CalFieldPdf<Ltc>(pdfData, maskData, loopRange);
//
//  // 主循环 
//  CollisionAndStream<Ltc, CLS, EQ> mainLoop(lc, kernel, equilibrium, pdfData, maskData, bcMapSingle, gridC, pdfComm);
//  nnwlbm::Output<nnwlbm::LinkedMaskData<double, dim>, nnwlbm::UniformGrid<dim>> output(maskData, gridC);
//
//  time_marching_mpi(output, mainLoop, 100000, 1000, 10);
//}
//
////void porescaleflow_mpi();
//
//void cavity2d_read_mpi() {
//  int whichType = PHSPACE::GlobalDataBase::GetIntParaFromDB("LBMSolver");
//
//  phlbm::BlockTextReader::ReadFromFile();
//  auto numBlock = phlbm::BlockTextReader::GetNumOfDataBlock();
//  auto data1 = phlbm::BlockTextReader::GetDataBlock(0);
//  auto dataBlock = phlbm::BlockTextReader::GetDataBlockVector();
//  auto dataBlockMap = phlbm::BlockTextReader::GetDataBlockMap();
//
//  // 物理参数设置
//  typedef nnwlbm::D2Q9 Ltc;  // 设置格子构型
//  Ltc lc(1.0);  // 生成格子，指定格子速度
//  const size_t dim = Ltc::GetDim();  // 获得维度信息
//  double viscosity = 0.0001;  // 设置粘性
//  std::vector<double> iniMacroVar(dim + 1); // 设置初条件
//  iniMacroVar[0] = 1.0; // 密度
//  iniMacroVar[1] = 0.0; // u速度
//  iniMacroVar[2] = 0.0; // v速度
//
//  const int num = 201; //  设置每条边的网格数
//  const double dx = 1.0 / num;  // 计算网格大小
//  const double dt = dx / lc.velocity(); // 计算时间步长
//  double tau = nnwlbm::CalTau(lc, viscosity, dt, iniMacroVar[0]);
//
//  // BGK or MRT
//  typedef nnwlbm::Equilibrium<Ltc> EQ;
//  EQ equilibrium;
//  typedef nnwlbm::BGK<Ltc, nnwlbm::Equilibrium<Ltc>> CLS;
//  //typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, nnwlbm::Equilibrium<Ltc>> CLS;
//  //typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, EQ> CLS;
//  CLS kernel(lc, tau);
//
//
//  // 计算区域的几何形状
//  phlbm::Generator<nnwlbm::Geometry> geoGen;
//  auto rec = geoGen(dataBlockMap.at("box"), dataBlockMap);
//  // 外部形状
//  auto outsideBoundary = std::make_shared<nnwlbm::ComplementGeo>(rec);
//  // 上壁面形状
//  auto topWall = geoGen(dataBlockMap.at("topWall"), dataBlockMap);
//
//  // 定义无滑移反弹边条件
//  phlbm::Generator<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>> factoryGen;
//  auto bounceBackFactory = factoryGen(dataBlockMap.at("bounceback"), dataBlockMap);
//
//  std::map<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<Ltc>*> bcBodyMap;
//  // 指定外部区域的边界为无滑移边条件
//  bcBodyMap.emplace(outsideBoundary.get(), bounceBackFactory.get());
//
//  // 定义无滑移移动壁面反弹边条件
//  auto movingBounceBackFactory = factoryGen(dataBlockMap.at("movingbounceback"), dataBlockMap);
//  // 定义几何体表示边界
//  std::map<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<Ltc>*> bcSurfaceMap;
//  // 指定上壁面为无滑移边条件
//  bcSurfaceMap.emplace(topWall.get(), movingBounceBackFactory.get());
//
//  // 初始化流场数据
//  nnwlbm::ConstFun initcond(std::make_shared<nnwlbm::UniversalGeo>(lc.GetDim()), iniMacroVar);
//
//  // 生成网格
//  double origin[2] = { 0.5*dx, 0.5*dx };
//  nnwlbm::UniformGrid<dim> grid(dx, origin);
//  auto gridC = phlbm::gen_grid_t<2>(dataBlockMap.at("grid"));
//  std::cout << "dx=" << gridC->get_dx() << std::endl;
//  // 指定计算区域
//  std::vector<nnwlbm::Range<dim>> rangeVec;
//  auto comm = MPI_COMM_WORLD;
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  int commSize = 0;
//  MPI_Comm_size(comm, &commSize);
//  rangeVec.resize(commSize);
//
//  /* partition
//  int remainder = num%commSize;
//  int quotient = num/commSize;
//  for(int i = 0; i < remainder; i++) {
//    rangeVec[i] = nnwlbm::Range<dim>(i*(quotient+1), 0, (i+1)*(quotient+1), num);
//  }
//  for(int i = remainder; i < commSize; i++) {
//    rangeVec[i] = nnwlbm::Range<dim>(
//      remainder*(quotient+1)+(i-remainder)*quotient, 0,
//      remainder*(quotient+1)+(i+1-remainder)*quotient, num
//      );
//  }*/
//
//  //partition 2
//  if (commSize != 5) exit(1);
//  int nodeNumOfBlock = num * num / 5;
//  int lengthOfCentralBlock = int(num / sqrt(5.0));
//  int widthOfOtherFourBlock = (num - lengthOfCentralBlock) / 2;
//  int clb = widthOfOtherFourBlock;
//  int crt = widthOfOtherFourBlock + lengthOfCentralBlock;
//
//  rangeVec[0] = nnwlbm::Range<dim>(0, 0, crt, clb);
//  rangeVec[1] = nnwlbm::Range<dim>(crt, 0, num, crt);
//  rangeVec[2] = nnwlbm::Range<dim>(clb, crt, num, num);
//  rangeVec[3] = nnwlbm::Range<dim>(0, clb, clb, num);
//  rangeVec[4] = nnwlbm::Range<dim>(clb, clb, crt, crt);
//
//  auto rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(rangeVec[rank]);
//  std::map<nnwlbm::Range<dim>, int> rangeToRank;
//
//  for (int i = 0; i < commSize; i++) {
//    rangeNeighbors->AddNeiborCandidate(rangeVec[i]);
//    rangeToRank.emplace(rangeVec[i], i);
//  }
//  auto hostRange = rangeNeighbors->GetHostRange();
//
//  hostRange.AddBuffer(1);
//  nnwlbm::LinkedMaskRange<int, dim> maskRange(hostRange, 0);
//  nnwlbm::LinkedMaskData<double, dim> maskData(3, maskRange);
//  nnwlbm::LinkedMaskData<double, dim> pdfData(nnwlbm::D2Q9::size(), maskRange, nnwlbm::Color());
//
//  PdfCommunicator<Ltc, nnwlbm::LinkedMaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);
//
//  //////////////////////////////////////////////////////////////////////////////////////
//  // 生成colorMap
//  nnwlbm::GenColorMap<dim, nnwlbm::BoundaryFactory<Ltc>>(maskRange, *gridC, bcBodyMap);
//
//  // 根据colorMap为分布函数和边界条件分配内存
//
//  // 初始化分布函数
//  ini_field(maskData, initcond, *gridC);
//
//  auto loopRange = maskRange.GetRange();
//  loopRange.AddBuffer(-1);
//  nnwlbm::CalFieldPdf<Ltc>(pdfData, maskData, loopRange);
//
//  nnwlbm::Output<nnwlbm::LinkedMaskData<double, dim>, nnwlbm::UniformGrid<dim>> output1(pdfData, *gridC);
//  output1.SetFileName(std::to_string(30000) + "pdf" + std::to_string(rank));
//  output1();
//  //nnwlbm::LinkedMaskData<double, dim> pdfData(lc.size(), maskRange);
//  nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SimpleBC>, Ltc::GetDim(), nnwlbm::BcIdentifier>bcMapSingle(Ltc::size(), maskRange, nnwlbm::BcIdentifier(nnwlbm::geo::boundary));
//
//  // 生成边界条件的映射
//  nnwlbm::gen_bcMap<Ltc>(bcMapSingle, grid, bcBodyMap, bcSurfaceMap, pdfData);
//
//  // 主循环 
//  CollisionAndStream1<Ltc, CLS, EQ> mainLoop(lc, kernel, equilibrium, pdfData, maskData, bcMapSingle, *gridC, pdfComm);
//  nnwlbm::Output<nnwlbm::LinkedMaskData<double, dim>, nnwlbm::UniformGrid<dim>> output(maskData, *gridC);
//
//  time_marching_mpi(output, mainLoop, 100000, 20000, 10);
//}

void stl_mpi() {
  //int whichType = PHSPACE::GlobalDataBase::GetIntParaFromDB("LBMSolver");

  phlbm::BlockTextReader::ReadFromFile();
  auto dataBlockMap = phlbm::BlockTextReader::GetDataBlockMap();

  //! calPara-calculatePara:计算物理参数设置
  typedef nnwlbm::D3Q15 Ltc;
  double latticeVel;
  dataBlockMap.at("calPara").GetData("latticeVel", &latticeVel, 1);
  Ltc lc(latticeVel);
  const size_t dim = Ltc::GetDim();
  double viscosity;
  dataBlockMap.at("calPara").GetData("viscosity", &viscosity, 1);
  std::vector<double> iniMacroVar(dim + 1);
  dataBlockMap.at("calPara").GetData("iniMacroVar", iniMacroVar.data(), dim + 1);
  int stepEnd, stepOutput;
  dataBlockMap.at("calPara").GetData("stepEnd", &stepEnd, 1);
  dataBlockMap.at("calPara").GetData("stepOutput", &stepOutput, 1);


  //! 读入初始化的格点以及网格数设置
  double calArea[2 * Ltc::GetDim()];
  dataBlockMap.at("calBox").GetData("box", calArea, 2 * Ltc::GetDim());
  auto range = phlbm::GenRangeT<3>(dataBlockMap.at("range"));
  int left[Ltc::GetDim()], right[Ltc::GetDim()];
  range->GetLeft(left); range->GetRight(right);
  //！计算离散步长
  const double dx = (calArea[Ltc::GetDim()] - calArea[0]) / (right[0] - left[0]);
  //！计算时间步长
  const double dt = dx / lc.velocity();
  //！计算松弛时间
  double tau = nnwlbm::CalTau(lc, viscosity, dt, iniMacroVar[0]);

  //! 读入实际的物理区间（时间步长和网格初始化坐标）
  double origin[Ltc::GetDim()];
  for (int counter = 0; counter < Ltc::GetDim(); counter++)
    origin[counter] = calArea[counter];
  nnwlbm::UniformGrid<dim> gridC(dx, origin);
  auto grid = gridC.GetStaggeredGrid();

  //! BGK和MRT选择
  typedef nnwlbm::Equilibrium<Ltc> EQ;
  EQ equilibrium;
  //typedef nnwlbm::BGK<Ltc, nnwlbm::Equilibrium<Ltc>> CLS;
  typedef nnwlbm::MRT<Ltc, nnwlbm::TransMatrix<Ltc::size()>, EQ> CLS;
  CLS kernel(lc, tau);


  //! 计算区域的几何形状
  phlbm::Generator<nnwlbm::Geometry> geoGen;
  //! calBox存储参数文件计算区域box
  auto calBox = geoGen(dataBlockMap.at("calBox"), dataBlockMap);
  //! outsideBoundary作为同样存储box
  auto outsideBoundary = std::make_shared<nnwlbm::ComplementGeo>(calBox);

  //! 入口区域范围读取
  auto inLet = geoGen(dataBlockMap.at("inLet"), dataBlockMap);
  std::make_shared<nnwlbm::Box>(dim);
  //! 出口区域范围读取
  auto outLet = geoGen(dataBlockMap.at("outLet"), dataBlockMap); 
  std::make_shared<nnwlbm::Box>(dim);

  //! 读入无滑移边界条件，设置bcBodyMap
  phlbm::Generator<nnwlbm::BoundaryFactory<Ltc>> factoryGen;
  auto bounceBackFactory = factoryGen(dataBlockMap.at("bounceback"), dataBlockMap);
  std::map<int,std::pair<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<Ltc>*>> bcBodyMap;
  bcBodyMap.emplace(0,std::make_pair(outsideBoundary.get(), bounceBackFactory.get()));

  //! 读入压力非平衡外推边界条件，包括入口与出口，设置bcSurfaceMap
  auto pressureFactoryIn = factoryGen(dataBlockMap.at("nonEquilibriumOfPressureIn"), dataBlockMap);
  auto pressureFactoryOut = factoryGen(dataBlockMap.at("nonEquilibriumOfPressureOut"), dataBlockMap);
  std::map<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<Ltc>*> bcSurfaceMap;
  bcSurfaceMap.emplace(inLet.get(), pressureFactoryIn.get());
  bcSurfaceMap.emplace(outLet.get(), pressureFactoryOut.get());

  //! 设置宏观初始量，
  nnwlbm::ConstFun initcond(std::make_shared<nnwlbm::UniversalGeo>(lc.GetDim()), iniMacroVar);


  //! MPI并行
  std::vector<nnwlbm::Range<dim>> rangeVec;
  auto comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int commSize = 0;
  MPI_Comm_size(comm, &commSize);
  rangeVec.resize(commSize);

  //! 区域划分设置
  int remainder = (right[0] - left[0]) % commSize;
  int quotient = (right[0] - left[0]) / commSize;
  for (int i = 0; i < remainder; i++) {
    rangeVec[i] = nnwlbm::Range<dim>(i*(quotient + 1), left[1], left[2], (i + 1) * (quotient + 1), right[1], right[2]);
  }
  for (int i = remainder; i < commSize; i++) {
    rangeVec[i] = nnwlbm::Range<dim>(
      remainder * (quotient + 1) + (i - remainder) * quotient, left[1], left[2],
      remainder * (quotient + 1) + (i + 1 - remainder) * quotient, right[1], right[2]
      );
  }

  auto rangeNeighbors = std::make_shared<RangeNeighbors<dim>>(rangeVec[rank]);
  std::map<nnwlbm::Range<dim>, int> rangeToRank;
  for (int i = 0; i < commSize; i++) {
    rangeNeighbors->AddNeiborCandidate(rangeVec[i]);
    rangeToRank.emplace(rangeVec[i], i);
  }
  auto hostRange = rangeNeighbors->GetHostRange();

  hostRange.AddBuffer(1);
  nnwlbm::LinkedMaskRange<int, dim> maskRange(hostRange, 0);
  nnwlbm::LinkedMaskData<double, dim> maskData(Ltc::GetDim() + 1, maskRange);
  nnwlbm::LinkedMaskData<double, dim> pdfData(Ltc::size(), maskRange, nnwlbm::Color());

  //! 读入对应的模型
  auto stlGeo = geoGen(dataBlockMap.at("StlGeo"), dataBlockMap);
  dynamic_cast<nnwlbm::StlGeo*>(stlGeo.get())->GenBoundingBox(dx);
  //dynamic_cast<nnwlbm::StlGeo*>(stlGeo.get())->sharpen_boundary();
  bcBodyMap.emplace(1,std::make_pair(stlGeo.get(), bounceBackFactory.get()));

  //! 更新colorMap 
  nnwlbm::GenColorMap<dim, nnwlbm::BoundaryFactory<Ltc>>(maskRange, grid, bcBodyMap);
  //! 初始化流场
  ini_field(maskData, initcond, grid);
  //! 初始化分布函数
  auto loopRange = maskRange.GetRange();
  loopRange.AddBuffer(-1);
  nnwlbm::CalFieldPdf<Ltc>(pdfData, maskData, loopRange);

  PdfCommunicator<Ltc, nnwlbm::LinkedMaskData> pdfComm(&pdfData, *rangeNeighbors, rangeToRank, comm);

  //! 初始化边界条件bcMapSingle
  nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SimpleBC>, Ltc::GetDim(), nnwlbm::BcIdentifier>bcMapSingle(Ltc::size(), maskRange, nnwlbm::BcIdentifier(nnwlbm::geo::boundary));
  nnwlbm::gen_bcMap<Ltc>(bcMapSingle, grid, bcBodyMap, bcSurfaceMap, pdfData);

  //! 主循环 
  CollisionAndStream1<Ltc, CLS, EQ> mainLoop(lc, kernel, equilibrium, pdfData, maskData, bcMapSingle, gridC, pdfComm);
  nnwlbm::Output<nnwlbm::LinkedMaskData<double, dim>, nnwlbm::UniformGrid<dim>> output(maskData, gridC);

  PHSPACE::PHMakeDir("results");
  time_marching_mpi(output, mainLoop, stepEnd, stepOutput, 10);
}

//void porescaleflow_mpi();

namespace phlbm {

  //std::vector<double> gen_triangleVec_from() {
  std::vector<double> gen_triangleVec_from() {
    PHSPACE::Region region;
    region.ReadGrid();
    PHSPACE::Zone* zPtr = region.GetZone(0);
    PHSPACE::Grid* phGrid = zPtr->GetGeometry()->GetGrid();
    PHSPACE::UnstructGrid *unsPhGrid = PHSPACE::UnstructGridCast(phGrid);

    int nBoundFace = unsPhGrid->GetNBoundFace();
    int *face2node = unsPhGrid->GetFace2Node();
    int *nodeNumberOfEachFace  = unsPhGrid->GetNodeNumberOfEachFace();

    PHSPACE::RDouble *x = unsPhGrid->GetX();
    PHSPACE::RDouble *y = unsPhGrid->GetY();
    PHSPACE::RDouble *z = unsPhGrid->GetZ();

    //! Boundary faces.
    int nodePosition = 0;
    int triNum = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
      triNum += nodeNumberOfEachFace[iFace] - 2;
      nodePosition += nodeNumberOfEachFace[iFace];
    }
    std::vector<double> ans(12 * triNum);
    nodePosition = 0;
    triNum = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
      for (int jNode = 2; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
      {
        double *tri = ans.data() + 12*triNum;
        int id = face2node[nodePosition];
        tri[0 + 3 * 1] = x[id];
        tri[1 + 3 * 1] = y[id];
        tri[2 + 3 * 1] = z[id];
      
        id = face2node[nodePosition + jNode - 1];
        tri[0 + 3 * 2] = x[id];
        tri[1 + 3 * 2] = y[id];
        tri[2 + 3 * 2] = z[id];
        id = face2node[nodePosition + jNode];
        tri[0 + 3 * 3] = x[id];
        tri[1 + 3 * 3] = y[id];
        tri[2 + 3 * 3] = z[id];

        double v1[3], v2[3];
        for(unsigned iDim = 0; iDim < 3; iDim++) {
          v1[iDim] = tri[iDim + 3 * 2] - tri[iDim + 3 * 1];
          v2[iDim] = tri[iDim + 3 * 3] - tri[iDim + 3 * 2];
        }
        Geo3d::crossproduct(tri, v1, v2);
        Geo3d::normalize(tri);
        triNum++;
      }
      nodePosition += nodeNumberOfEachFace[iFace];
    }
    return std::move(ans);
  }
}


void RunLBMSolver()
{
  //auto stlGeo = std::make_shared<nnwlbm::StlGeo>();
  //stlGeo->ReadFile("Plane_t.stl");
  //stlGeo->GenBoundingBox(0.5);
  //cavity2d_read_mpi();
#ifdef _DEBUG
  std::cout << "running DEBUG version" << std::endl;
#else
  std::cout << "running RELEASE version" << std::endl;
#endif
  stl_mpi();
  //cavity3d_mpi();
  //porescaleflow_mpi();

  //MPI_Finalize();
  return;
}

