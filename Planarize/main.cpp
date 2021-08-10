#include <iostream>
#include <fstream>
#include "common.h"
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/planarity.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/dcel.h>
#include <igl/per_face_normals.h>
#include <igl/unique_edge_map.h>
#include <map>
#include <set>
#include <igl/copyleft/cgal/submesh_aabb_tree.h>
#include "CeresPlanarize.h"
#include <igl/local_basis.h>

#include "../common/CLI11.hpp"

//#define EXTRA_OUTPUT

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


Eigen::MatrixXd V, VO, N, NTrig, centers, VOrg;
Eigen::MatrixXi F, T, FOrg;
Eigen::VectorXi D,  DOrg, hexCenter2tri;
Eigen::MatrixXd FEs;
Eigen::MatrixXi EF, FE, EV, EV_diagonal, EFi, TEdges;
Eigen::VectorXi innerEdges;
Eigen::VectorXd planarity;
Eigen::VectorXi innerFacesFlags;

std::vector<Triangle> triangles;
std::map<int, Eigen::Vector3d> hexVer2Closest;
std::map<int, std::vector<int> > vert2star;
std::map<int, Eigen::Matrix3d> trig2S;
Tree tree;

double LC, LP;
int iter;

void findClosestElements() {
  hexCenter2tri.resize(centers.rows());
  for (int i = 0; i < centers.rows(); i++)
  {
    Point_3 query(centers(i,0), centers(i,1), centers(i,2));
    auto projection = tree.closest_point_and_primitive(query);
    size_t fid = projection.second - triangles.begin();
    hexCenter2tri(i) = fid;
  }

  hexVer2Closest.clear();
  for (int i = 0; i < VO.rows(); i++)
  {
    Point_3 query(VO(i,0), VO(i,1), VO(i,2));
    auto projection = tree.closest_point_and_primitive(query);
    Point_3 closest_point = projection.first;
    hexVer2Closest.insert({i, Eigen::Vector3d(CGAL::to_double(closest_point.x()), CGAL::to_double(closest_point.y()), CGAL::to_double(closest_point.z()))});
  }
}

void computeEdgeLength(Eigen::VectorXd & lengths)
{
  lengths.resize(EV.rows(), 1);
  for(int curr = 0; curr < EV.rows(); curr++)
  {
    lengths(curr) = (VO.row(EV(curr, 0)) - VO.row(EV(curr, 1))).norm();
  }
}

void findDiagonals(const Eigen::MatrixXi & F, const Eigen::VectorXi & innerFacesFlags, const Eigen::VectorXi & D, Eigen::MatrixXi & EV_diagonal)
{
  EV_diagonal.resize(0, 2);
  for(int i = 0; i < F.rows(); i++) {
    if (((D(i) % 2) != 0) || (! innerFacesFlags(i)))
      continue;
    for(int j = 0; j < D(i); j++) {
      EV_diagonal.conservativeResize(EV_diagonal.rows() + 1, 2);
      EV_diagonal.row(EV_diagonal.rows() - 1) << Eigen::RowVector2i(F(i, j), F(i, (j + D(i) / 2) % D(i)));
    }
  }
}

void getPlaneEnergyValues(double currSolution[], std::vector<double> & planeEnergyValues)
{
  PlaneEnergy pe;
  double residual = 0;

  for(int i = 0; i < F.rows(); i++)
  {
    for (int j = 0; j < D(i); j++)
    {
     pe(currSolution + 3 * F(i, j), currSolution + 3 * F(i, (j + 1) % D(i)), currSolution + 3 * (V.rows() + i), &residual);
     planeEnergyValues.push_back(residual);
    }
  }
}


void findInnerFaces(const Eigen::MatrixXi & F, const Eigen::MatrixXi & EF, Eigen::VectorXi & innerFacesFlags)
{
  innerFacesFlags.resize(F.rows(), 1);
  for(int i = 0; i < EF.rows(); i++)
  {
   int fid0 = EF(i, 0);
   int fid1 = EF(i, 1);
   if(fid0 != -1 && fid1 != -1)
   {
     innerFacesFlags(fid0) = 1;
     innerFacesFlags(fid1) = 1;
   }
   else if (fid0 != -1)
     innerFacesFlags(fid0) = 0;
   else
     innerFacesFlags(fid1) = 0;
  }
}


void find_true_boundary_edges(const Eigen::MatrixXd & V,
                              const Eigen::MatrixXi & F,
                              const Eigen::VectorXi & D,
                              const Eigen::VectorXi & innerEdges,
                              const Eigen::MatrixXi & EF,
                              const Eigen::MatrixXi & EV,
                              const Eigen::MatrixXi & EFi,
                              std::vector<std::pair<int, int> > & true_boundary_edges,
                              Eigen::VectorXi & vdeg)
{
  Eigen::VectorXi VH, HE, HF, nextH, prevH, twinH, HV;
  Eigen::MatrixXi EH, FH;
  hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH,  HV,  HE, HF, nextH, prevH, twinH);

  vdeg = Eigen::VectorXi::Zero(V.rows());

  for (int i = 0; i < EV.rows(); i++)
  {
    vdeg(EV(i, 0))++;
    vdeg(EV(i, 1))++;
  }
  for (int vid = 0; vid < VH.rows(); vid++)
  {
    int hid = VH(vid);
    if (twinH(hid) != -1 || vdeg(vid) < 3)
      continue;

    int vid0 = vid;
    int vid1 = -1;

    int next = nextH(hid);
    do
    {
      if(vdeg(HV(next)) > 2 && twinH(prevH(next)) == -1)
      {
        vid1 = HV(next);
        break;
      }

      next = nextH(next);
    } while(next != hid);

    if(vid1 != -1)
      true_boundary_edges.emplace_back(vid0, vid1);
  }
}

void find_closest_point_on_boundary(const Eigen::MatrixXd & V, const Eigen::MatrixXd & VH, const std::vector<std::pair<int, int> > & org_boundry_edges, const Eigen::VectorXi & innerVertices, std::map<int, Eigen::Vector3d> & hexVer2Closest)
{
  for (int i = 0; i < VH.rows(); i++)
  {
    double mind = 100000;
    Eigen::Vector3d closest;
    Eigen::Vector3d P = VH.row(i);
    if(innerVertices(i))
      continue;
    for(auto e : org_boundry_edges)
    {
      Eigen::Vector3d A = V.row(e.first);
      Eigen::Vector3d B = V.row(e.second);
      auto AP = P - A;
      auto AB = B - A;
      auto tmp = A + AP.dot(AB) / AB.dot(AB) * AB;
      double dist = (P - tmp).norm();
      if(dist < mind)
      {
        mind = dist;
        closest = tmp;
      }
    }
    hexVer2Closest[i] = closest;
  }
}

void findVertexStar(const Eigen::MatrixXi & F, const Eigen::VectorXi & D, const Eigen::VectorXi & innerVertices, std::map<int, std::vector<int> > & v2star)
{
  // re-do the internal vertices to have a correct ordering
  Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  Eigen::MatrixXd FEsPoly;
  hedra::polygonal_edge_topology(D, F,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  hedra::dcel(D,F, EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);

  std::map<int, std::set<int> > tmp_v2star;
  for(int i = 0; i < EVPoly.rows(); i++)
  {
    int v0 = EV(i, 0);
    int v1 = EV(i, 1);

    tmp_v2star[v0].insert(v1);
    tmp_v2star[v1].insert(v0);
  }
  for(const auto& v2s : tmp_v2star)
  {
    for (auto val : v2s.second)
      v2star[v2s.first].push_back(val);
  }

  for (int i = 0; i < innerVertices.rows(); i++)
  {
    if(!innerVertices(i))
      continue;
    int count = 0;
    auto e = VHPoly(i);
    auto next = e;
    do
    {
      v2star[i][count++] = HVPoly(twinHPoly(next));
      next = nextHPoly(twinHPoly(next));
    } while(next != e);
  }

}

int main(int argc, char ** argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string inputTrigFileName;
  std::string outputFileName;
  std::string outputNormalsFileName;
  std::string outputPlanarityFileName;
  bool snapBoundary = false;
  bool quadMode = false;
  bool singSym = false;

  double planarityThreshold;
  int maxInter = 100;

  app.description("This program planarizes a polyhedral mesh where the faces can be of any degreee higher than three.");

  //add all the options supported.
  app.add_option("-i,--input, 1", inputFileName, "The path to the OFF file containing the input mesh.")->required()->check(CLI::ExistingFile);

  app.add_option("-t,--input_trig, 2", inputTrigFileName, "The path to the OFF file containing a triangulation of the mesh to be planarized.")->required()->check(CLI::ExistingFile);

  app.add_option("-o,--output, 3", outputFileName, "The path to the output OFF file.")->required();

  app.add_option("-T,--threshold, 4", planarityThreshold , "Max. planarity error threshold (in percent) for stopping the optimization.")->required()->check(CLI::NonNegativeNumber);

  app.add_option("-n,--output_normals, 5", outputNormalsFileName, "The path to the output CSV file containing normals.");

  app.add_option("-P,--output_planarity, 6", outputPlanarityFileName, "The path to a CSV containing planarity per face for each iteration.");

  app.add_option("-m,--max_iter, 7", maxInter , "Max. number of iterations.", true)->default_val(100);

  app.add_flag("-S,--snap_boundary", snapBoundary , "Keep the boundary points at the boundary of the original triangulation.")->default_val(false);
  app.add_flag("-Q,--quad_mode", quadMode , "Add extra energy part that refines the quad mesh.")->default_val(false);
  app.add_flag("-s,--sing_symm_mode", singSym, "Add extra energy part that enforces symmetry of odd degree faces.")->default_val(false);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  iter = 0;

  // Load a mesh in OFF format
  if(!hedra::polygonal_read_OFF(inputFileName, V, D, F))
    throw std::runtime_error("File could not be read!");

  Eigen::MatrixXd planarityPerFaceIter;
  planarityPerFaceIter.resize(F.rows(), 1);

  if(!hedra::polygonal_read_OFF(inputTrigFileName, VOrg, DOrg, FOrg))
    throw std::runtime_error("File could not be read!");

  hedra::polygonal_edge_topology(DOrg, FOrg, EV,FE, EF, EFi, FEs, innerEdges);

  std::vector<std::pair<int, int> > org_boundry_edges;
  for(int i = 0; i < EF.rows(); i++)
  {
    int fid0 = EF(i, 0);
    int fid1 = EF(i, 1);
    if (fid0 != -1 && fid1 != -1)
      continue;
    org_boundry_edges.emplace_back(EV(i, 0), EV(i, 1));
  }

  Eigen::MatrixXd B1, B2, B3;
  igl::local_basis(VOrg, FOrg, B1, B2, B3);
  
  LC = 0.01;
  LP = 0.01;
  VO = V;
  hedra::polygonal_edge_topology(D, F, EV,FE, EF, EFi, FEs, innerEdges);

  std::vector<std::pair<int, int> > true_boundary_edges;
  Eigen::VectorXi vdeg;
  find_true_boundary_edges(V, F, D, innerEdges, EF, EV, EFi, true_boundary_edges, vdeg);

  Eigen::VectorXi innerVertices = Eigen::VectorXi::Constant(V.rows(), 1);
  findInnerFaces(F,EF, innerFacesFlags);

  for(int i = 0; i < EF.rows(); i++)
  {
    if(EF(i, 0) != -1 && EF(i, 1) != -1)
      continue;
    innerVertices(EV(i, 0)) = 0;
    innerVertices(EV(i, 1)) = 0;
  }
  find_closest_point_on_boundary(VOrg, V, org_boundry_edges, innerVertices, hexVer2Closest);
  if(quadMode)
    findVertexStar(F, D, innerVertices, vert2star);

  hedra::planarity(V, D, F, planarity);

  if(!outputPlanarityFileName.empty())
  {
    for(int i = 0; i < F.rows(); i++)
    {
      planarityPerFaceIter(i, 0) = planarity(i);
    }
  }

  std::cout << "INITIAL MAX PLANARITY: " << planarity.maxCoeff() << std::endl;

  hedra::polygonal_face_centers(VO, D, F, centers);
  Eigen::MatrixXi EOrg, uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<size_t> > uE2E;
  igl::unique_edge_map(FOrg, EOrg, uE, EMAP, uE2E);
  Eigen::VectorXi closest_facet_orientations;

  std::vector<bool> in_I;
  Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi>(FOrg.rows(), 0, FOrg.rows()-1);
  igl::copyleft::cgal::submesh_aabb_tree(VOrg, FOrg, I , tree,triangles,in_I);

  findClosestElements();

  igl::per_face_normals(VOrg, FOrg, NTrig);


  double planarize = planarity.maxCoeff();

  Eigen::VectorXd edge_lengths;
  computeEdgeLength(edge_lengths);

  std::cout << "Min. initial edge length: " << edge_lengths.minCoeff() << ", Max. initial edge length: " << edge_lengths.maxCoeff() << std::endl;
  double min_len = edge_lengths.minCoeff();


  findDiagonals(F, innerFacesFlags,  D, EV_diagonal);

  CeresPlanarize cf;
  while(planarize > planarityThreshold && iter <= maxInter)
  {
    hedra::polygonal_face_centers(VO, D, F, centers);
    polyhedral_face_normals(VO, D, F, N);
    findClosestElements();
    if(snapBoundary)
      find_closest_point_on_boundary(VOrg, VO, org_boundry_edges, innerVertices, hexVer2Closest);

    cf.init(F, VO, V, N, D, centers, hexVer2Closest, hexCenter2tri, trig2S, EF, FE, innerEdges, EV, EV_diagonal, LP, LC, innerVertices, innerFacesFlags, true_boundary_edges, vdeg, vert2star, singSym);
    LP *= 2.0;
    cf.solve(false);

    // write the function values for plotting
    std::vector<double> planeEnergyValues;
    getPlaneEnergyValues(cf.currSolution, planeEnergyValues);

    iter++;
    if(iter % 5 == 0)
      LC *= 2.0;

    VO.resize(V.rows(), V.cols());
    for (int i = 0; i < V.rows(); i++) {
        Eigen::Vector3d tmp(cf.currSolution[3 * i], cf.currSolution[3 * i + 1], cf.currSolution[3 * i + 2]);
        VO.row(i) = tmp;
    }

    hedra::planarity(VO, D, F, planarity);
    planarize = planarity.maxCoeff();

    if(!outputPlanarityFileName.empty())
    {
      planarityPerFaceIter.conservativeResize(F.rows(), planarityPerFaceIter.cols()+1);
      int idx = planarityPerFaceIter.cols() - 1;
      for(int i = 0; i < F.rows(); i++)
      {
        planarityPerFaceIter(i, idx) = planarity(i);
      }
    }

    std::cout << "Iter : " << iter << ", max planarity " << planarize << std::endl;
    computeEdgeLength(edge_lengths);
    std::cout << "Iter : " << iter << ", min. edge length: " << edge_lengths.minCoeff() << ", max. edge length: " << edge_lengths.maxCoeff() << std::endl;
  }

  // export the normals
  if (!outputNormalsFileName.empty())
  {
    std::ofstream file;
    file.flags (std::ios::scientific);
    file.precision (std::numeric_limits<double>::digits10 + 1);
    file.open(outputNormalsFileName);
    for (int i = 0; i < F.rows(); i++) {
      Eigen::Vector3d tmp(cf.currSolution[3 * (V.rows() + i)], cf.currSolution[3 * (V.rows() + i) + 1], cf.currSolution[3 * (V.rows() + i) + 2]);
      tmp.normalize();
      file << tmp(0) << ", " << tmp(1) << ", " << tmp(2) << std::endl;
    }
    file.close();
  }
  hedra::polygonal_write_OFF(outputFileName, VO, D, F, true);

  if(!outputPlanarityFileName.empty()) {
    std::ofstream file;
    file.flags (std::ios::scientific);
    file.precision (std::numeric_limits<double>::digits10 + 1);
    file.open(outputPlanarityFileName);
    for (int i = 0; i < F.rows(); i++) {
      for(int j = 0; j < planarityPerFaceIter.cols() - 1; j++)
      {
        file << planarityPerFaceIter(i, j) << ", ";
      }
      file << planarityPerFaceIter(i, planarityPerFaceIter.cols() - 1) << std::endl;
    }
    file.close();
  }
    return 0;
}