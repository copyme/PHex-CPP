#include <iostream>
#include <Eigen/Core>
#include <igl/edge_topology.h>

#include <igl/serialize.h>
#include <igl/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/write_singularities.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/combing.h>
#include <directional/read_matching.h>
#include <directional/write_matching.h>
#include <directional/cut_mesh_with_singularities.h>
#include <hedra/polygonal_write_OFF.h>

#include "parameterize.h"
#include "../common/setup_parameterization.h"
#include "../common/CLI11.hpp"

Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VMeshCut;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd paramFuncsd, paramFuncsN, wholeCornerParamFuncsN;   //#F x N parameterization functions

int N;

void effort_from_matching(const Eigen::MatrixXd& V,
                          const Eigen::MatrixXi& F,
                          const Eigen::MatrixXi& EV,
                          const Eigen::MatrixXi& EF,
                          const Eigen::MatrixXd& rawField,
                          const Eigen::VectorXi& matching,
                          Eigen::VectorXd& effort)
{
  typedef std::complex<double> Complex;

  Eigen::MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);

  int N = rawField.cols() / 3;

  Eigen::VectorXcd edgeTransport(EF.rows()); //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
  Eigen::MatrixXd edgeVectors(EF.rows(), 3);
  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1)
      continue;
    edgeVectors.row(i) = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
    Complex ef(edgeVectors.row(i).dot(B1.row(EF(i, 0))), edgeVectors.row(i).dot(B2.row(EF(i, 0))));
    Complex eg(edgeVectors.row(i).dot(B1.row(EF(i, 1))), edgeVectors.row(i).dot(B2.row(EF(i, 1))));
    edgeTransport(i) = eg / ef;
  }

  effort = Eigen::VectorXd::Zero(EF.rows());
  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1)
      continue;
//computing the full effort for 0->indexMinFromZero, and readjusting the matching to fit principal effort
    double currEffort=0;
    for (int j = 0; j < N; j++) {
      Eigen::RowVector3d vecjf = rawField.block(EF(i, 0), 3*j, 1, 3);
      Complex vecjfc = Complex(vecjf.dot(B1.row(EF(i, 0))), vecjf.dot(B2.row(EF(i, 0))));
      Eigen::RowVector3d vecjg = rawField.block(EF(i, 1), 3 * ((matching(i)+j+N)%N), 1, 3);
      Complex vecjgc = Complex(vecjg.dot(B1.row(EF(i, 1))), vecjg.dot(B2.row(EF(i, 1))));
      Complex transvecjfc = vecjfc*edgeTransport(i);
      currEffort+= arg(vecjgc / transvecjfc);
    }

    effort(i) = currEffort;
  }
}

directional::ParameterizationData &
setFieldInformation(const std::string &inputMatchingFileName,
                    const std::string &inputSerializedFileName,
                    const std::string &inputSingFileName,
                    bool & isFieldCombed,
                    directional::ParameterizationData &pd,
                    Eigen::VectorXd &curl) {
  if(inputMatchingFileName.empty()) {
    if(isFieldCombed)
    {
      std::cout << "No matching file provided but the field marked as combed: isFieldCombed set to false." << std::endl;
      isFieldCombed = false;
    }
    igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
    directional::curl_matching(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, effort, curl);
  }
  else
  {
    if(isFieldCombed)
      directional::read_matching(inputMatchingFileName, combedMatching, EF, EV, FE, N);
    else {
      directional::read_matching(inputMatchingFileName, matching, EF, EV, FE, N);
      effort_from_matching(VMeshWhole, FMeshWhole, EV, EF, rawField, matching, effort);
    }
  }
  if(isFieldCombed) {
    combedField = rawField;
    if(inputSingFileName.empty())
    {
      directional::effort_to_indices(VMeshWhole, FMeshWhole, EV, EF, effort, matching, N, singVertices, singIndices);
    }
    if(inputSerializedFileName.empty())
      directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  }
  else
  {
    if(!inputSerializedFileName.empty())
    {
      std::cout << "Input serialized PD data provided but the input field is not combed: PD data IGNORED!" << std::endl;
    }

    if(!inputSingFileName.empty())
    {
      std::cout << "Input singularity data provided but the input field is not combed: singularity data IGNORED!" << std::endl;
    }

    pd = {}; // set to zero;
    directional::effort_to_indices(VMeshWhole, FMeshWhole, EV, EF, effort, matching, N, singVertices, singIndices);
    directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
    directional::combing(VMeshWhole, FMeshWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField,combedMatching);
  }
  return pd;
}

int main(int argc, char* argv[])
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputMeshFileName;
  std::string inputFieldFileName;
  std::string inputMatchingFileName;
  std::string inputCutMeshFileName;
  std::string inputSerializedFileName;
  std::string inputSingFileName;
  std::string outputFieldFileName;
  std::string outputMatchingFileName;
  std::string outputSingularFileName;
  std::string outputCutMeshFileName;
  std::string outputSerliazedFileName;
  double lengthRatio = 1;
  bool normalizeField = false;
  bool isInteger = false;
  bool isFieldCombed = false;
  std::string matname("poisson.mat");

  app.description("This program combs a directional field defined on a triangular mesh.");

  //add all the options supported.
  app.add_option("-i,--input, 1", inputMeshFileName, "The path to the OFF file containing the input mesh.")->required()->check(CLI::ExistingFile);
  app.add_option("-f,--input_field, 2", inputFieldFileName, "The path to the rawfield file containing the vector field.")->required()->check(CLI::ExistingFile);
  app.add_option("-U,--input_serialized_pd, 3", inputSerializedFileName, "Input PD do be deserialized.")->check(CLI::ExistingFile);
  app.add_option("-C,--input_cut_mesh, 4", inputCutMeshFileName, "Input cut mesh")->check(CLI::ExistingFile);
  app.add_option("-s,--input_sing, 5", inputSingFileName, "The path to the sing file (singularities).")->check(CLI::ExistingFile);
  app.add_option("-R,--length_ratio, 6", lengthRatio , "Length ratio.")->required()->check(CLI::NonNegativeNumber);
  app.add_option("-M,--input_matching, 8", inputMatchingFileName, "File with predefined matching.")->check(CLI::ExistingFile);
  app.add_option("-o,--output, 7", outputSerliazedFileName, "The path to the serialized PD.");
  app.add_option("-L,--output_matlab, 9", matname, "File name for the matlab output.")->default_val("poisson.mat");
  app.add_option("-m,--output_matching, 10", outputMatchingFileName, "Output combed matching.");
  app.add_option("-F,--output_combed_filed, 11", outputFieldFileName, "Output combed filed.");
  app.add_option("-S,--output_sing_filed, 12", outputSingularFileName, "Path to singularities information of the combed field. (sign file format).");
  app.add_option("-c,--output_cut_mesh, 13", outputCutMeshFileName, "Output cut mesh.");
  app.add_flag("-N,--normalize_field", normalizeField, "Keep the boundary point at the boundary of the original triangulation.");
  app.add_flag("-I,--use_integers", isInteger, "Rounding on/off.");
  app.add_flag("-A,--is_combed", isFieldCombed, "Assume the input field and matching is already combed.");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  igl::readOFF(inputMeshFileName, VMeshWhole, FMeshWhole);
  if(!inputCutMeshFileName.empty())
    igl::readOFF(inputCutMeshFileName, VMeshCut, FMeshCut);
  directional::read_raw_field(inputFieldFileName, N, rawField);

  if(normalizeField){
    double totalLengthSum=0.0;
    for (int i = 0; i < FMeshWhole.rows(); i++)
        for (int j = 0; j < N; j++)
            totalLengthSum += rawField.block(i, 3 * j, 1, 3).norm();
    rawField.array() *= (double)(N * FMeshWhole.rows()) / totalLengthSum;
  }

  //combing and cutting
  directional::ParameterizationData pd;
  Eigen::VectorXd curl;

  if(!inputSerializedFileName.empty())
    igl::deserialize(pd, inputSerializedFileName);

  if(!inputSingFileName.empty())
    directional::read_singularities(inputSingFileName, N, singVertices, singIndices);

  pd = setFieldInformation(inputMatchingFileName, inputSerializedFileName, inputSingFileName, isFieldCombed, pd, curl);

  if(!outputFieldFileName.empty())
    directional::write_raw_field(outputFieldFileName, combedField, true);

  if(!outputMatchingFileName.empty())
    directional::write_matching(outputMatchingFileName, combedMatching, EF, EV, FE, N);

  if(!outputSingularFileName.empty())
    directional::write_singularities(outputSingularFileName, N, singVertices, singIndices);

  Eigen::MatrixXi symmFunc;
  if (N!=6){
  //generic sign symmetry
  symmFunc.resize(N,N/2);
  symmFunc.block(0,0,N/2,N/2)=Eigen::MatrixXi::Identity(N/2,N/2);
  symmFunc.block(N/2,0,N/2,N/2)=-Eigen::MatrixXi::Identity(N/2,N/2);
  } else {
  symmFunc.resize(6,2);
  symmFunc.row(0)<<1,0;
  symmFunc.row(1)<<0,1;
  symmFunc.row(2)<<-1,1;
  symmFunc.block(3,0,3,2)=-symmFunc.block(0,0,3,2);
  }

  Eigen::MatrixXi intLatticeFunc = Eigen::MatrixXi::Identity(N/2., N/2.);

  if(!isFieldCombed || inputCutMeshFileName.empty() || inputSerializedFileName.empty())
    directional::setup_parameterization(symmFunc, intLatticeFunc, VMeshWhole, FMeshWhole, EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);

  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, matname, paramFuncsd, paramFuncsN, wholeCornerParamFuncsN);

  if(!outputSerliazedFileName.empty())
    igl::serialize(pd, outputSerliazedFileName);

  if(!outputCutMeshFileName.empty()) {
    Eigen::VectorXi DMeshCut = Eigen::VectorXi::Constant(FMeshCut.rows(), 3);
    hedra::polygonal_write_OFF(outputCutMeshFileName, VMeshCut, DMeshCut, FMeshCut, true);
  }
  return 0;
}


