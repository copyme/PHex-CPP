
//
// Created by kacper on 12/16/19.
//

#ifndef PLANARIZE_CERESPLANARIZE_H
#define PLANARIZE_CERESPLANARIZE_H

#include <cstring>
#include <cmath>
#include <ceres/ceres.h>
#include <ceres/loss_function.h>
#include <ceres/dynamic_autodiff_cost_function.h>
#include <glog/logging.h>
#include <memory>
#include <utility>
#include <valarray>

#define MAX_LOG_LEVEL -100

struct PlaneEnergy
{
  PlaneEnergy() = default;

  template <typename T>
  bool operator()(const T * const  v0, const T * const  v1, const T * const  n, T * residuals) const
  {
    T v0mv1[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    T lenSQR = ceres::sqrt(v0mv1[0] * v0mv1[0] + v0mv1[1] * v0mv1[1] + v0mv1[2] * v0mv1[2]);
    T INNER = (v0mv1[0] * n[0] + v0mv1[1] * n[1] + v0mv1[2] * n[2]);
    residuals[0] = INNER/* / lenSQR*/;
    return true;
  }
};

struct ClosenessEnergy
{
  ClosenessEnergy(const Eigen::Vector3d & n, const Eigen::Vector3d & org) : no {n(0), n(1), n(2)}, vo {org(0), org(1), org(2)} {}

  template <typename T>
  bool operator()(const T * const  vi, T * residuals) const
  {
    T v0mv1[3] = {vi[0] - vo[0], vi[1] - vo[1], vi[2] - vo[2]};
    T lenSQR = v0mv1[0] * v0mv1[0] + v0mv1[1] * v0mv1[1] + v0mv1[2] * v0mv1[2];
    T INNER = v0mv1[0] * no[0] + v0mv1[1] * no[1] + v0mv1[2] * no[2];
    residuals[0] = INNER;
    return true;
  }

private:
  double no[3]{};
  double vo[3]{};
};


struct VertexSimilarityOddEnergy
{
  typedef ceres::DynamicAutoDiffCostFunction<VertexSimilarityOddEnergy, 5> VertexSimilarityOddEnergyDynamicCostFunction;
  explicit VertexSimilarityOddEnergy(int pDeg) : deg(pDeg) {};

  template <typename T>
  bool operator()(T const* const* points, T* residuals) const {
    //we want the vertices to have a projection on the middle of the opposite edge and minimize the distance from the projection to the edge center
    for(int i = 0; i < deg; i++) {
      int idxVE0 = (int)(i + std::floor(deg / 2.)) % deg;
      int idxVE1 = (idxVE0 + 1) % deg;

      T ap[3] = {points[i][0] - points[idxVE0][0], points[i][1] - points[idxVE0][1], points[i][2] - points[idxVE0][2]};
      T ab[3] = {points[idxVE1][0] - points[idxVE0][0], points[idxVE1][1] - points[idxVE0][1], points[idxVE1][2] - points[idxVE0][2]};
      T ABdotAB = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
      T APdotAB = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2];
      T projPontoAB[3] = {points[idxVE0][0] + APdotAB / ABdotAB  * ab[0], points[idxVE0][1] + APdotAB / ABdotAB  * ab[1], points[idxVE0][2] + APdotAB / ABdotAB  * ab[2]};
      T edgeC[3] = {(points[idxVE1][0] + points[idxVE0][0])/ T(2.), (points[idxVE1][1] + points[idxVE0][1]) / T(2.), (points[idxVE1][2] + points[idxVE0][2]) / T(2.)};
      residuals[3 * i] = projPontoAB[0] - edgeC[0];
      residuals[3 * i + 1] = projPontoAB[1] - edgeC[1];
      residuals[3 * i + 2] = projPontoAB[2] - edgeC[2];
    }
    return true;
  }
private:
  int deg;
public:
  // Factory method to create a CostFunction from a VertexSimilarityOddEnergyDynamicCostFunction to
  // conveniently add to a ceres problem.
  static VertexSimilarityOddEnergyDynamicCostFunction* Create(const int deg,
                                             const std::vector<int> & indexes,
                                             double* system,
                                             std::vector<double*>* parameter_blocks) {
    auto* constraint = new VertexSimilarityOddEnergy(deg);
    auto* cost_function = new VertexSimilarityOddEnergyDynamicCostFunction(constraint);
    // Add all the parameter blocks that affect this constraint.
    parameter_blocks->clear();
    for (int i = 0; i < deg; i++) {
      parameter_blocks->push_back(system + 3 * indexes[i]);
      cost_function->AddParameterBlock(3);
    }
    cost_function->SetNumResiduals(3 * deg);
    return (cost_function);
  }
};


struct VertexSimilarityEvenEnergy
{
  explicit VertexSimilarityEvenEnergy(const Eigen::Vector3d & center) : c {center(0), center(1), center(2)} {};

  template <typename T>
  bool operator()(const T * const  v0, const T * const  v1, T * residuals) const
  {
    T vf[3] = {c[0] - v0[0], c[1] - v0[1], c[2] - v0[2]};
    T vg[3] = {v1[0] - c[0], v1[1] - c[1], v1[2] - c[2]};
    residuals[0] = vf[0] - vg[0];
    residuals[1] = vf[1] - vg[1];
    residuals[2] = vf[2] - vg[2];
    return true;
  }

private:
  double c[3]{};
};

struct DistanceToSurface {
  explicit DistanceToSurface(const Eigen::Vector3d & p) : sp {p(0), p(1), p(2)} {}

  template <typename T>
  bool operator()(const T * const  v0, T * residuals) const
  {
    residuals[0] = v0[0] - sp[0];
    residuals[1] = v0[1] - sp[1];
    residuals[2] = v0[2] - sp[2];
    return true;
  }
private:
  double sp[3]{};
};


struct LengthEnergyBarrierGeo
{
  explicit LengthEnergyBarrierGeo(double SL) : S(SL * SL) {};

  template <typename T>
  bool operator()(const T * const  v0, const T * const  v1, T * residuals) const
  {
    T v0mv1[3] = {v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]};
    T LEN = ceres::sqrt(v0mv1[0] * v0mv1[0] + v0mv1[1] * v0mv1[1] + v0mv1[2] * v0mv1[2]);
    T ST = T(S);

    if(LEN <= T(10e-20))
      residuals[0] = T(10000000);
    else if(LEN >= ST)
      residuals[0] = T(0.);
    else
      residuals[0] = (T(1.) / (ceres::pow(LEN / ST, T(3.)) - ((T(3.) * LEN * LEN) / (ST * ST)) + (T(3.) * LEN / ST))) - T(1.);
    return true;
  }

private:
  double S;
};


struct FairnessDynamic
{
  typedef ceres::DynamicAutoDiffCostFunction<FairnessDynamic, 5> FarinessDynamicCostFunction;
  explicit FairnessDynamic(int pN) : N(pN) {}

  template <typename T>
    bool operator()(T const* const* points, T* residuals) const {
      T S[3] = {T(0.)};
      for (int i = 1; i < N; i++) {
        for (int j = 0; j < 3; j++)
          S[j] += points[i][j];
      }

    residuals[0] = points[0][0] - S[0] / T(N - 1);
    residuals[1] = points[0][1] - S[1] / T(N - 1);
    residuals[2] = points[0][2] - S[2] / T(N - 1);
    return true;
  }
private:
  int N;
public:
  // Factory method to create a CostFunction from a FarinessDynamicCostFunction to
  // conveniently add to a ceres problem.
  static FarinessDynamicCostFunction* Create(const int N,
                                             const std::vector<int> & indexes,
                                             double* system,
                                             std::vector<double*>* parameter_blocks) {
    auto* constraint = new FairnessDynamic(N);
    auto* cost_function = new FarinessDynamicCostFunction(constraint);
    // Add all the parameter blocks that affect this constraint.
    parameter_blocks->clear();
    for (int i = 0; i < N; i++) {
      parameter_blocks->push_back(system + 3 * indexes[i]);
      cost_function->AddParameterBlock(3);
    }
    cost_function->SetNumResiduals(3);
    return (cost_function);
  }
};

class CeresPlanarize {
public:
  CeresPlanarize() : currSolution(nullptr), problem(nullptr) {
  }

  ~CeresPlanarize() {
    delete[] currSolution;
    delete problem;
  }

private:
  ceres::Problem *problem;

public:
  double *currSolution;

public:
  void init(const Eigen::MatrixXi & F, const Eigen::MatrixXd & V,
            const Eigen::MatrixXd & VINIT,
            const Eigen::MatrixXd & N, const Eigen::MatrixXi & D,
            const Eigen::MatrixXd & centers, const std::map<int, Eigen::Vector3d> & hexVer2Closest,
            const Eigen::VectorXi & hexCenter2tr, const std::map<int, Eigen::Matrix3d> & trig2S,
            const Eigen::MatrixXi & EF,
            const Eigen::MatrixXi & FE,
            const Eigen::VectorXi & innerEdges,
            const Eigen::MatrixXi & EV,
            const Eigen::MatrixXi & EV_diagonal,
            double LP, double LC,
            const Eigen::VectorXi & innerVertices,
            const Eigen::VectorXi & innerFacesFlags,
            const std::vector<std::pair<int, int> > & true_boundary_edges,
            const Eigen::VectorXi & vdeg,
            const std::map<int, std::vector<int>>& v2star,
            bool singSym
            )
  {
    delete[] currSolution;
    delete problem;
    problem = new ceres::Problem;
    currSolution = new double[V.size() + N.size()]();
    for(int i = 0; i < V.rows(); i++)
    {
      currSolution[i * 3] = V(i, 0);
      currSolution[i * 3 + 1] = V(i, 1);
      currSolution[i * 3 + 2] = V(i, 2);
    }
    for(int i = 0; i < N.rows(); i++)
    {
      if(N.row(i).norm() < 10e-6 || std::isnan(N.row(i).norm()))
      {
        std::cout << "Corrupted! " << std::endl;
        continue;
      }
      currSolution[3 * (V.rows() + i)] = N(i, 0);
      currSolution[3 * (V.rows() + i) + 1] = N(i, 1);
      currSolution[3 * (V.rows() + i) + 2] = N(i, 2);
    }

    for(int i = 0; i < F.rows(); i++) {

      problem->AddParameterBlock(currSolution + 3 * (V.rows() + i), 3, new ceres::HomogeneousVectorParameterization(3));
      for(int j = 0; j < D(i); j++) {
        ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<PlaneEnergy, 1, 3, 3, 3>(new PlaneEnergy);
        auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), LP, ceres::TAKE_OWNERSHIP);
        problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * F(i, j), currSolution + 3 * F(i, (j + 1) % D(i)), currSolution + 3 * (V.rows() + i));
      }
    }


    for (int i = 0; i < V.rows(); i++)
    {
      double w = 0.5;
      if(!innerVertices(i))
        w =  10;
      ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<DistanceToSurface, 3, 3>(new DistanceToSurface(hexVer2Closest.at(i)));
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), w, ceres::TAKE_OWNERSHIP);
      problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * i);
    }


    for(int i = 0; i < F.rows(); i++) {
      if (((D(i) % 2) != 0) || (! innerFacesFlags(i)))
        continue;
      for(int j = 0; j < D(i) / 2.; j++) {
        ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<VertexSimilarityEvenEnergy, 3, 3, 3>(new VertexSimilarityEvenEnergy(centers.row(i)));
        auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), 1, ceres::TAKE_OWNERSHIP);
        problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * F(i, j),currSolution + 3 * F(i, (j + D(i) / 2) % D(i)));
      }
    }
    // for singularities
    for(int i = 0; (i < F.rows()) && singSym; i++) {
      std::vector<int> indexes;
      if (((D(i) % 2) == 0) || (! innerFacesFlags(i)))
        continue;
      for(int j = 0; j < D(i); j++) {
        indexes.push_back(F(i, j));
      }
      std::vector<double*> parameter_blocks;
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), 0.25, ceres::TAKE_OWNERSHIP);
      VertexSimilarityOddEnergy::VertexSimilarityOddEnergyDynamicCostFunction * range_cost_function = VertexSimilarityOddEnergy::Create(D(i), indexes, currSolution, &parameter_blocks);
      problem->AddResidualBlock(range_cost_function, loss_function, parameter_blocks);
    }

    Eigen::MatrixXd VN;
    average_onto_vertices(V, F,D, N, VN);
    for(int i = 0; i < V.rows(); i++) {

      Eigen::Vector3d n = VN.row(i).normalized();
      ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<ClosenessEnergy, 1, 3>(new ClosenessEnergy(n, V.row(i)));
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), 0.01, ceres::TAKE_OWNERSHIP);
      problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * i);
    }

    for (int curr = 0; curr < EV.rows(); curr++)
    {
      double len = (VINIT.row(EV(curr, 0)) - VINIT.row(EV(curr, 1))).norm();

      ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<LengthEnergyBarrierGeo, 1, 3, 3>(new LengthEnergyBarrierGeo(len));
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), LC/10, ceres::TAKE_OWNERSHIP);
      problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * EV(curr, 0), currSolution + 3 * EV(curr, 1));
    }

    for (int curr = 0; curr < EV_diagonal.rows(); curr++)
    {
      double len = (VINIT.row(EV_diagonal(curr, 0)) - VINIT.row(EV_diagonal(curr, 1))).norm();
      ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<LengthEnergyBarrierGeo, 1, 3, 3>(new LengthEnergyBarrierGeo(len));
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), LC/10, ceres::TAKE_OWNERSHIP);
      problem->AddResidualBlock(cost_function, loss_function, currSolution + 3 * EV_diagonal(curr, 0), currSolution + 3 * EV_diagonal(curr, 1));
    }

    for (auto v2s : v2star)
    {
      if (!innerVertices(v2s.first))
        continue;
      double w = 5.0;
      std::vector<double*> parameter_blocks;
      std::vector<int> indexes;
      indexes.push_back(v2s.first);
      indexes.insert(indexes.end(), v2s.second.begin(), v2s.second.end());
      auto loss_function = new ceres::ScaledLoss(new ceres::TrivialLoss(), w, ceres::TAKE_OWNERSHIP);
      FairnessDynamic::FarinessDynamicCostFunction * range_cost_function = FairnessDynamic::Create(indexes.size(), indexes, currSolution, &parameter_blocks);
      problem->AddResidualBlock(range_cost_function, loss_function, parameter_blocks);
    }
  }

  void solve(const bool outputProgress)
  {
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = outputProgress;
    options.max_num_iterations = 250;
    options.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.num_threads = 4;
    options.logging_type = ceres::SILENT;

    ceres::Solver::Summary summary;
    ceres::Solve(options, problem, &summary);
    if (outputProgress)
      std::cout << summary.FullReport() << "\n";
  }
};

#endif //PLANARIZE_CERESPLANARIZE_H
