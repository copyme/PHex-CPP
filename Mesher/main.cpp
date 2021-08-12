#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/edge_topology.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/dcel.h>
#include <hedra/polygonal_write_OFF.h>
#include <igl/serialize.h>

#include "Mesh.h"
#include "../common/setup_parameterization.h"
#include "../common/CLI11.hpp"

Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VMeshCut;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd paramFuncsd, wholeCornerParamFuncsN;   //#F x N parameterization functions
Eigen::MatrixXd genV;
Eigen::MatrixXi genF, genT;
Eigen::VectorXi genD;

void read_compact_func(const std::string& filename, Eigen::MatrixXd & paramFuncsd) {
    std::ifstream file;
    int lineCounter = 0;
    std::string line;
    std::string token;

    file.open(filename);
    if (!file.is_open())
    {
        std::cout << "The compact function file could not be read!" << std::endl;
        exit(0);
    }

    while(std::getline(file, line))
    {
        std::istringstream iss(line);
        lineCounter++;
        std::vector<double> vals;
        while (std::getline(iss, token, ',')) {
            vals.push_back(std::stod(token));
        }
        paramFuncsd.conservativeResize(lineCounter, vals.size());
        for(int i =0; i < vals.size(); i++){
            paramFuncsd(lineCounter - 1, i) = vals[i];
        }
    }
    file.close();
}

void read_per_corner_fuc(const std::string& filename, Eigen::MatrixXd & wholeCornerParamFuncsN) {
    std::ifstream file;
    int lineCounter = 0;
    std::string line;
    std::string token;

    file.open(filename);
    if (!file.is_open())
    {
        std::cout << "The per-corner function file could not be read!" << std::endl;
        exit(0);
    }

    file.flags(std::ios::scientific);
    file.precision(std::numeric_limits<double>::digits10 + 1);

    while(std::getline(file, line))
    {
        std::istringstream iss(line);
        lineCounter++;
        std::vector<double> vals;
        while (std::getline(iss, token, ',')) {
            vals.push_back(std::stod(token));
        }
        wholeCornerParamFuncsN.conservativeResize(lineCounter, vals.size());
        for(int i =0; i < vals.size(); i++)
            wholeCornerParamFuncsN(lineCounter - 1, i) = vals[i];
    }
    file.close();
}

double min_triangle_angle(const Eigen::MatrixXd & V, const Eigen::RowVectorXi & F)
{
  double min = 100000;
  for(int i = 0; i < 3; i++)
  {
    Eigen::RowVector3d v0 = V.row(F(0, i));
    Eigen::RowVector3d vm1 = V.row(F(0, (3 + (i - 1)) % 3));
    Eigen::RowVector3d vp1 = V.row(F(0, (i + 1) % 3));

    Eigen::RowVector3d vv1 = (vm1 - v0).normalized();
    Eigen::RowVector3d vv2 = (vp1 - v0).normalized();
    double tcos = vv1.dot(vv2);
    if(tcos > 1.0)
      tcos = 1.0;
    if(tcos < -1.0)
      tcos = -1.0;
    double angle = std::acos(tcos);
    if(min > angle)
      min = angle;
  }
  return min;
}

void triangulate(const Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & D)
{
  int OLD_SIZE = F.rows();
  for (int i = 0; i < OLD_SIZE; i++)
  {
    if(D(i) != 4)
      continue;

    Eigen::RowVectorXi tmp = F.row(i);
    Eigen::RowVectorXi f0, f1, f02, f12;

    //first set to pick from
    f0.resize(1, 3);
    f1.resize(1, 3);
    f0(0) = tmp(0);
    f0(1) = tmp(1);
    f0(2) = tmp(3);

    f1(0) = tmp(1);
    f1(1) = tmp(2);
    f1(2) = tmp(3);

    //second set to pick from
    f02.resize(1, 3);
    f12.resize(1, 3);
    f02(0) = tmp(0);
    f02(1) = tmp(1);
    f02(2) = tmp(2);

    f12(0) = tmp(0);
    f12(1) = tmp(2);
    f12(2) = tmp(3);

    //find min angles
    double minf0 = min_triangle_angle(V, f0);
    double minf1 = min_triangle_angle(V, f1);
    double minmax0 = minf1;
    if(minf0 > minf1)
      minmax0 = minf0;

    double minf02 = min_triangle_angle(V, f02);
    double minf12 = min_triangle_angle(V, f12);

    double minmax2 = minf12;
    if(minf02 > minf12)
      minmax2 = minf02;

    D(i) = 3;
    D.conservativeResize(D.rows() + 1, 1);
    F.conservativeResize(F.rows() + 1, F.cols());
    D(D.rows() - 1) = 3;

    if(minmax0 > minmax2) {
      F.block<1, 3>(i, 0) = f0;
      F.block<1, 3>(F.rows() - 1, 0) = f1;
    }
    else
    {
      F.block<1, 3>(i, 0) = f02;
      F.block<1, 3>(F.rows() - 1, 0) = f12;
    }
  }
}

int main(int argc, char* argv[])
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputMeshFileName;
  std::string inputCutMeshFileName;
  std::string inputSerializedPDFileName;
  std::string inputPerCornerFFileName;
  std::string inputCompactFFFileName;

  std::string outputNHMeshFileName;
  std::string outputSHMeshFileName;
  std::string outputDebugLogFileName;

  app.description("This program re-meshes the input mesh using the parametrization.");

  //add all the options supported.
  app.add_option("-i,--input, 1", inputMeshFileName, "The path to the OFF file containing the input mesh.")->required()->check(CLI::ExistingFile);
  app.add_option("-c,--input_cutmesh, 2", inputCutMeshFileName, "The path to the OFF file containing the input cut mesh.")->required()->check(CLI::ExistingFile);
  app.add_option("-u,--input_serialized_pd, 3", inputSerializedPDFileName, "The path to the serialized parametrization data.")->required()->check(CLI::ExistingFile);
  app.add_option("-p,--input_per_cor_func, 4", inputPerCornerFFileName, "Input per corner functions.")->required()->check(CLI::ExistingFile);
  app.add_option("-C,--input_comp_func, 5", inputCompactFFFileName, "Input per corner functions.")->required()->check(CLI::ExistingFile);

  app.add_option("-n,--output_nsmesh, 6", outputNHMeshFileName, "The path to the output non-simplified output mesh.");
  app.add_option("-o,--output_smesh, 7", outputSHMeshFileName, "The path to the output simplified mesh.");
  app.add_option("-d,--output_debug, 8", outputDebugLogFileName, "The path to the output debug log file.");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  if (!igl::readOFF(inputMeshFileName, VMeshWhole, FMeshWhole))
  {
      std::cout << "The mesh file could not be read!" << std::endl;
      exit(0);
  }
  if (!igl::readOFF(inputCutMeshFileName, VMeshCut, FMeshCut))
  {
      std::cout << "The cut mesh file could not be read!" << std::endl;
      exit(0);
  }
  
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

  directional::ParameterizationData pd;
  if (!igl::deserialize(pd, inputSerializedPDFileName))
  {
      std::cout << "The data could not be deserialized!" << std::endl;
      exit(0);
  }

  read_per_corner_fuc(inputPerCornerFFileName, paramFuncsd);
  read_compact_func(inputCompactFFFileName, wholeCornerParamFuncsN);

  std::ofstream DebugLog;
  std::ofstream dummy;
  dummy.setstate(std::ios::failbit) ;

  if(!outputDebugLogFileName.empty())
  {
      DebugLog.open(outputDebugLogFileName);
      if (!DebugLog.is_open())
      {
          std::cout << "Debug log file could not be open in the write mode. Logging will be omitted." << std::endl;
          DebugLog.setstate(std::ios::failbit); //somehow redundant
      }
  }
  Mesh TMesh(dummy), HMesh(DebugLog);

    //creating TMesh
  Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  Eigen::MatrixXd FEsPoly;
  hedra::polygonal_edge_topology(Eigen::VectorXi::Constant(FMeshWhole.rows(),3), FMeshWhole,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  hedra::dcel(Eigen::VectorXi::Constant(FMeshWhole.rows(),3),FMeshWhole,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,
              HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
  //creating TMesh
  TMesh.fromHedraDCEL(Eigen::VectorXi::Constant(FMeshWhole.rows(),3),VMeshWhole, FMeshWhole, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly, VHPoly,
                      EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, VMeshCut, FMeshCut, paramFuncsd, pd.d,
                      pd.vertexTrans2CutMatInteger*pd.symmMatInteger,  pd.constraintMatInteger*pd.symmMatInteger, wholeCornerParamFuncsN, pd.integerVars);
 
  TMesh.GenerateMesh(pd.N, HMesh);
  
  Eigen::MatrixXi genFfuncNum,gePrescribedAnglesInt;
  Eigen::MatrixXd cornerAngles;
  HMesh.toHedra(genV,genD,  genF, genFfuncNum, gePrescribedAnglesInt,cornerAngles);

  if(!outputNHMeshFileName.empty())
  {
      if (!hedra::polygonal_write_OFF(std::string(argv[3]), genV, genD, genF, true))
      {
          std::cout << "The output for non-simplified mesh file could not be written!" << std::endl;
          exit(0);
      }
  }

  std::cout << "Simplifying Mesh." << std::endl;
  bool success = HMesh.SimplifyHexMesh(pd.N);
  if (!success){
    std::cout << "Mesh simplification failed at some point! " << std::endl;
    exit(0);
  } else {
    std::cout << "Mesh simplification seems to have succeeded!" << std::endl;
  }
  
  Eigen::MatrixXi genFfuncNum2, gePrescribedAnglesInt2;
  Eigen::MatrixXd cornerAngles2;
  HMesh.toHedra(genV,genD,  genF, genFfuncNum2, gePrescribedAnglesInt2, cornerAngles2);

  if(pd.N == 6)
    triangulate(genV, genF, genD);


  if(!outputSHMeshFileName.empty())
  {
      if (!hedra::polygonal_write_OFF(outputSHMeshFileName, genV, genD, genF))
      {
          std::cout << "The output mesh file could not be written!" << std::endl;
          exit(0);
      }
  }
  std::cout << "Done!" << std::endl;
 
  return 0; 
}



