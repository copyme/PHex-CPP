#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/edge_topology.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/write_matching.h>
#include <directional/read_matching.h>
#include <directional/read_singularities.h>
#include <directional/write_singularities.h>
#include <directional/principal_matching.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/combing.h>
#include <directional/cut_mesh_with_singularities.h>
#include <hedra/polygonal_write_OFF.h>
#include <map>

#include "../common/setup_parameterization.h"
#include "../common/CLI11.hpp"

Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VMeshCut;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;


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
  using namespace Eigen;
  using namespace std;

  MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);

  int N = rawField.cols() / 3;

  VectorXcd edgeTransport(EF.rows()); //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
  MatrixXd edgeVectors(EF.rows(), 3);
  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1)
      continue;
    edgeVectors.row(i) = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
    Complex ef(edgeVectors.row(i).dot(B1.row(EF(i, 0))), edgeVectors.row(i).dot(B2.row(EF(i, 0))));
    Complex eg(edgeVectors.row(i).dot(B1.row(EF(i, 1))), edgeVectors.row(i).dot(B2.row(EF(i, 1))));
    edgeTransport(i) = eg / ef;
  }

  effort = VectorXd::Zero(EF.rows());
  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1)
      continue;
    //computing the full effort for 0->indexMinFromZero, and readjusting the matching to fit principal effort
    double currEffort=0;
    for (int j = 0; j < N; j++) {
      RowVector3d vecjf = rawField.block(EF(i, 0), 3*j, 1, 3);
      Complex vecjfc = Complex(vecjf.dot(B1.row(EF(i, 0))), vecjf.dot(B2.row(EF(i, 0))));
      RowVector3d vecjg = rawField.block(EF(i, 1), 3 * ((matching(i)+j+N)%N), 1, 3);
      Complex vecjgc = Complex(vecjg.dot(B1.row(EF(i, 1))), vecjg.dot(B2.row(EF(i, 1))));
      Complex transvecjfc = vecjfc*edgeTransport(i);
      currEffort+= arg(vecjgc / transvecjfc);
    }
    effort(i) = currEffort;
  }
}

int main(int argc, char* argv[])
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputMeshFileName;
  std::string inputFieldFileName;
  std::string outputCutMeshFileName;
  std::string inputMatchingFileName;
  std::string inputSingFileName;

  std::string outputFieldFileName;
  std::string outputMatchingFileName;
  std::string outputCombedMatchingFileName;
  std::string outputSingularFileName;
  std::string outputSerializedFileName;

  app.description("This program combs a directional field defined on a triangular mesh.");

  //add all the options supported.
  app.add_option("-i,--input, 1", inputMeshFileName, "The path to the OFF file containing the input mesh.")->required()->check(CLI::ExistingFile);

  app.add_option("-f,--input_field, 2", inputFieldFileName, "The path to the rawfield file containing the vector field.")->required()->check(CLI::ExistingFile);

  app.add_option("-M,--input_matching, 3", inputMatchingFileName, "File with predefined matching.")->check(CLI::ExistingFile);

  app.add_option("-s,--input_sing, 4", inputSingFileName, "The path to the sing file (singularities).")->check(CLI::ExistingFile);

  app.add_option("-o,--output, 5", outputCutMeshFileName, "The path to the output OFF file.");

  app.add_option("-U,--output_serialized_pd, 6", outputSerializedFileName, "The output serialized pd");

  app.add_option("-m,--output_combed_matching, 7", outputCombedMatchingFileName, "Output combed matching.");

  app.add_option("-a,--output_matching, 8", outputMatchingFileName, "Output uncombed matching.");

  app.add_option("-F,--output_combed_filed, 9", outputFieldFileName, "Output combed filed.");

  app.add_option("-S,--output_sing_filed, 10", outputSingularFileName, "Path to singularities information. (sign file format).");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  igl::readOFF(inputMeshFileName, VMeshWhole, FMeshWhole);
  directional::read_raw_field(inputFieldFileName, N, rawField);
  
  //combing and cutting
  directional::ParameterizationData pd;
  Eigen::VectorXd curl;

  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  if(inputMatchingFileName.empty()) {
    directional::curl_matching(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, effort, curl);
  }
  else
  {
    directional::read_matching(inputMatchingFileName, matching, EF, EV, FE, N);
    effort_from_matching(VMeshWhole, FMeshWhole, EV, EF, rawField, matching, effort);
  }

  if(inputSingFileName.empty())
    directional::effort_to_indices(VMeshWhole, FMeshWhole, EV, EF, effort, matching, N, singVertices, singIndices);
  else
    directional::read_singularities(inputSingFileName, N, singVertices, singIndices);
  directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  directional::combing(VMeshWhole, FMeshWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField,combedMatching);

  if (!outputMatchingFileName.empty())
      directional::write_matching(outputMatchingFileName, matching, EF, EV, FE, N);

  if(!outputCombedMatchingFileName.empty())
    directional::write_matching(outputCombedMatchingFileName, combedMatching, EF, EV, FE, N);

  if(!outputFieldFileName.empty())
    directional::write_raw_field(outputFieldFileName, combedField, true);

  if(!outputSingularFileName.empty())
    directional::write_singularities(outputSingularFileName, N, singVertices, singIndices);

  Eigen::MatrixXi symmFunc;
  if (N!=6){
  //generic sign symmetry
  symmFunc.resize(N,N/2);
  symmFunc.block(0,0,N/2,N/2) = Eigen::MatrixXi::Identity(N/2,N/2);
  symmFunc.block(N/2,0,N/2,N/2) = -Eigen::MatrixXi::Identity(N/2,N/2);
} else {
  symmFunc.resize(6,2);
  symmFunc.row(0)<<1,0;
  symmFunc.row(1)<<0,1;
  symmFunc.row(2)<<-1,1;
  symmFunc.block(3,0,3,2)=-symmFunc.block(0,0,3,2);
}
  Eigen::MatrixXi intLatticeFunc = Eigen::MatrixXi::Identity(N/2., N/2.);
  directional::setup_parameterization(symmFunc, intLatticeFunc, VMeshWhole, FMeshWhole, EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);

  if(!outputSerializedFileName.empty()) {
    igl::serialize(pd, outputSerializedFileName);
  }

  if(!outputCutMeshFileName.empty()) {
    Eigen::VectorXi DMeshCut = Eigen::VectorXi::Constant(FMeshCut.rows(), 3);
    hedra::polygonal_write_OFF(outputCutMeshFileName, VMeshCut, DMeshCut, FMeshCut, true);
  }
  return 0;
}

