#include <fstream>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <hedra/dual_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>

#include "../common/CLI11.hpp"

Eigen::MatrixXi F;
Eigen::MatrixXd V;
Eigen::VectorXi D;

Eigen::MatrixXi dualF;
Eigen::MatrixXd dualV;
Eigen::VectorXi dualD;


template <typename T>
void removeRow(T& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

void remove_border(const Eigen::MatrixXi& EF, Eigen::MatrixXi& F, Eigen::VectorXi& D, Eigen::MatrixXd& V)
{
    std::vector<int> face2remove;

    for (int i = 0; i < EF.rows(); i++)
    {
        if (EF(i, 0) == -1)
            face2remove.push_back(EF(i, 1));

        if (EF(i, 1) == -1)
            face2remove.push_back(EF(i, 0));
    }

    std::sort(face2remove.begin(), face2remove.end(), std::greater<>());
    face2remove.erase(std::unique(face2remove.begin(), face2remove.end()), face2remove.end());

    for (unsigned int fid : face2remove)
    {
        removeRow<Eigen::MatrixXi>(F, fid);
        removeRow<Eigen::VectorXi>(D, fid);
    }

    std::vector<bool> markV(V.rows(), false);
    for (int i = 0; i < F.rows(); i++)
    {
        for (int d = 0; d < D(i); d++)
        {
            markV[F(i, d)] = true;
        }
    }
    std::vector<int> vertex2remove;
    for (size_t i = 0; i < markV.size(); i++)
    {
        if (!markV[i])
            vertex2remove.push_back(i);
    }
    std::sort(vertex2remove.begin(), vertex2remove.end());
    vertex2remove.erase(std::unique(vertex2remove.begin(), vertex2remove.end()), vertex2remove.end());

    std::sort(vertex2remove.begin(), vertex2remove.end(), std::greater<>());

    for (unsigned int vid : vertex2remove)
    {
        removeRow<Eigen::MatrixXd>(V, vid);
    }

    for (int i = 0; i < F.rows(); i++)
    {
        for (int d = 0; d < D(i); d++)
        {
            for (int vid : vertex2remove)
            {
                if (F(i, d) > vid) {
                    F(i, d)--;
                }
            }
        }
    }
}

int main(int argc, char* argv[])
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputMeshFileName;
  std::string outputMeshFileName;
  std::string boundaryMode;

  app.description("Program for computing a dual mesh.");

  //add all the options supported.
  app.add_option("-i,--input, 1", inputMeshFileName, "The path to the OFF file containing the input mesh.")->required()->check(CLI::ExistingFile);

  app.add_option("-o,--output, 2", outputMeshFileName, "The path to the output OFF file.")->required();

  app.add_option("-b,--boundary_mode, 3", boundaryMode, "How to treat the boundary.")->required()->check(CLI::IsMember({"NOTHING", "CLIPPED", "CLIPPED_TRIG", "CLIPPED_TRIG_DUAL"}));

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  if(! hedra::polygonal_read_OFF(inputMeshFileName, V, D, F) )
    throw std::runtime_error("The input file could not be read!");

  Eigen::MatrixXi EV, FE, EF, EFi;
  Eigen::MatrixXd FEs;
  Eigen::VectorXi InnerEdges;
  hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, InnerEdges);

  if (boundaryMode == "CLIPPED")
  {
      hedra::dual_mesh(V, D, F, hedra::LINEAR_SUBDIVISION, dualV, dualD, dualF, true);
      hedra::polygonal_write_OFF(outputMeshFileName, dualV, dualD, dualF, true);
  }
  else if (boundaryMode == "CLIPPED_TRIG")
  {
      remove_border(EF, F, D, V);
      hedra::dual_mesh(V, D, F, hedra::LINEAR_SUBDIVISION, dualV, dualD, dualF, false);
      hedra::polygonal_write_OFF(outputMeshFileName, dualV, dualD, dualF, true);
  }
  else if (boundaryMode == "CLIPPED_TRIG_DUAL")
  {
      remove_border(EF, F, D, V);
      hedra::dual_mesh(V, D, F, hedra::LINEAR_SUBDIVISION, dualV, dualD, dualF, true);
      hedra::polygonal_write_OFF(outputMeshFileName, dualV, dualD, dualF, true);
  }
  else if(boundaryMode == "NOTHING")
  {
      hedra::dual_mesh(V, D, F, hedra::LINEAR_SUBDIVISION, dualV, dualD, dualF, false);
      hedra::polygonal_write_OFF(outputMeshFileName, dualV, dualD, dualF, true);
  }

  return 0;
}





