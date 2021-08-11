// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PARAMETERIZE_H
#define DIRECTIONAL_PARAMETERIZE_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/matlab_format.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/matlab/MatlabWorkspace.h>

#include "branched_gradient.h"
#include "../common/setup_parameterization.h"

namespace directional
{

  // Creates a parameterization of (currently supported) (u,v, -u,-v) functions from a directional field by solving the Poisson equation, with custom edge weights
  // Input:
  //  wholeV:          #V x 3 vertex coordinates of the original mesh.
  //  wholeF:          #F x 3 face vertex indices of the original mesh.
  //  FE:              #F x 3 faces to edges indices.
  //  rawField:        #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  //  lengthRatio      #edgeLength/bounding_box_diagonal of quad mesh (scaling the gradient).
  //  pd:              Parameterization data obtained from directional::setup_parameterization.
  //  cutV:            #cV x 3 vertices of the cut mesh.
  //  cutF:            #F x 3 faces of the cut mesh.
  //  roundIntegers:   which variables (from #V+#T) are rounded iteratively to double integers. for each "x" entry that means that the [4*x,4*x+4] entries of vt will be double integer.
  // matname:          name of the output MATLAB file
  // Output:
  //  paramFuncsd             #cV x d parameterization functions per cut vertex (the compact version)
  //  paramFuncsN             #cV x N parameterization functions per cut vertex (full version with all symmetries unpacked)
  // wholeCornerParamFuncsN   (3*N) x #F parameterization functions per corner of whole mesh
  IGL_INLINE void parameterize(const Eigen::MatrixXd& wholeV,
                               const Eigen::MatrixXi& wholeF,
                               const Eigen::MatrixXi& FE,
                               const Eigen::MatrixXd rawField,
                               const double lengthRatio,
                               const ParameterizationData& pd,
                               const Eigen::MatrixXd& cutV,
                               const Eigen::MatrixXi& cutF,
                               const bool roundIntegers,
                               const std::string matname,
                               Eigen::MatrixXd& paramFuncsd,
                               Eigen::MatrixXd& paramFuncsN,
                               Eigen::MatrixXd& wholeCornerParamFuncsN)
  
  
  {
    Eigen::VectorXd edgeWeights = Eigen::VectorXd::Constant(FE.maxCoeff() + 1, 1.0);
    double length = igl::bounding_box_diagonal(wholeV) * lengthRatio;

    //TODO: in vertex space, not corner...
    int N = pd.N; //rawField.cols() / 3;
    long numVars = pd.symmMat.cols();
    //constructing face differentials
    std::vector<Eigen::Triplet<double> >  d0Triplets;
    std::vector<Eigen::Triplet<double> > M1Triplets;
    Eigen::VectorXd gamma(3 * N * wholeF.rows());
    for(int i = 0; i < cutF.rows(); i++)
    {
      for(int j = 0; j < 3; j++)
      {
        for(int k = 0; k < N; k++)
        {
          d0Triplets.emplace_back(3 * N * i + N * j + k, N * cutF(i, j) + k, -1.0);
          d0Triplets.emplace_back(3 * N * i + N * j + k, N * cutF(i, (j + 1) % 3) + k, 1.0);
          Eigen::Vector3d edgeVector = (cutV.row(cutF(i, (j + 1) % 3)) - cutV.row(cutF(i, j))).transpose();
          gamma(3 * N * i + N * j + k) = (rawField.block(i, 3 * k, 1, 3) * edgeVector)(0, 0) / length;
          M1Triplets.emplace_back(3 * N * i + N * j + k, 3 * N * i + N * j + k, edgeWeights(FE(i, j)));
        }
      }
    }
    Eigen::SparseMatrix<double> d0(3 * N * wholeF.rows(), N * cutV.rows());
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
    Eigen::SparseMatrix<double> M1(3 * N * wholeF.rows(), 3 * N * wholeF.rows());
    M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
    Eigen::SparseMatrix<double> d0T = d0.transpose();


    //The variables that should be fixed in the end
    Eigen::VectorXi fixedMask(numVars);
    fixedMask.setZero();

    for (int i=0;i<pd.fixedIndices.size();i++)
      fixedMask(pd.fixedIndices(i)) = 1;

    if(roundIntegers)
      for(int i = 0; i < pd.integerVars.size(); i++)
        for (int j=0;j<pd.d;j++)
          fixedMask(pd.d*pd.integerVars(i)+j) = 1;

    //the variables that were already fixed to begin with
    Eigen::VectorXi alreadyFixed(numVars);
    alreadyFixed.setZero();

    for (int i=0;i<pd.fixedIndices.size();i++)
      alreadyFixed(pd.fixedIndices(i)) = 1;

    //the values for the fixed variables (size is as all variables)
    Eigen::VectorXd fixedValues(numVars);
    fixedValues.setZero();  //for everything but the originally fixed values
    for (int i=0;i<pd.fixedValues.size();i++)
      fixedValues(pd.fixedIndices(i))=pd.fixedValues(i);

    Eigen::SparseMatrix<double> Efull = d0 * pd.vertexTrans2CutMat * pd.symmMat * pd.intSpanMat;
    Eigen::VectorXd x, xprev;

    // until then all the N depedencies should be resolved?

    //reducing constraintMat
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > qrsolver;
    Eigen::SparseMatrix<double> Cfull = pd.constraintMat * pd.symmMat * pd.intSpanMat;

    long CRank = 0;
    if (Cfull.rows() > 0)
    {
        qrsolver.compute(Cfull.transpose());
        CRank = qrsolver.rank();
    }
    //creating sliced permutation matrix
    Eigen::VectorXi PIndices;
    if(Cfull.rows() > 0)
        PIndices = qrsolver.colsPermutation().indices();

    std::vector<Eigen::Triplet<double> > CTriplets;
    for(int k = 0; k < Cfull.outerSize(); ++k)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Cfull, k); it; ++it)
      {
        for(int j = 0; j < CRank; j++)
          if(it.row() == PIndices(j))
            CTriplets.emplace_back(j, it.col(), it.value());
      }
    }

    Cfull.resize(CRank, Cfull.cols());
    Cfull.setFromTriplets(CTriplets.begin(), CTriplets.end());
    Eigen::SparseMatrix<double> var2AllMat;
    Eigen::VectorXd fullx(numVars); fullx.setZero();
    for(int intIter = 0; intIter < fixedMask.sum(); intIter++)
    {
      //the non-fixed variables to all variables
      var2AllMat.resize(numVars, numVars - alreadyFixed.sum());
      int varCounter = 0;
      std::vector<Eigen::Triplet<double> > var2AllTriplets;
      for(int i = 0; i < numVars; i++)
      {
        if (!alreadyFixed(i)){
          var2AllTriplets.emplace_back(i, varCounter++, 1.0);
        }

      }
      var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());

      Eigen::SparseMatrix<double> Epart = Efull * var2AllMat;
      Eigen::VectorXd torhs = -Efull * fixedValues;
      Eigen::SparseMatrix<double> EtE = Epart.transpose() * M1 * Epart;
      Eigen::SparseMatrix<double> Cpart = Cfull * var2AllMat;

      //reducing rank on Cpart
      long CpartRank = 0;
      if (Cfull.rows() > 0)
      {
          qrsolver.compute(Cpart.transpose());
          CpartRank = qrsolver.rank();
      }
      //creating sliced permutation matrix
      if(Cfull.rows() > 0)
          PIndices = qrsolver.colsPermutation().indices();

      std::vector<Eigen::Triplet<double> > CPartTriplets;
      Eigen::VectorXd bpart(CpartRank);
      for(int k = 0; k < Cpart.outerSize(); ++k)
      {
       for (Eigen::SparseMatrix<double>::InnerIterator it(Cpart, k); it; ++it)
        {
          for (int j = 0; j < CpartRank; j++)
            if (it.row() == PIndices(j))
              CPartTriplets.emplace_back(j, it.col(), it.value());
        }
      }

      Cpart.resize(CpartRank, Cpart.cols());
      Cpart.setFromTriplets(CPartTriplets.begin(), CPartTriplets.end());
      Eigen::SparseMatrix<double> A(EtE.rows()+ Cpart.rows(), EtE.rows() + Cpart.rows());

      std::vector<Eigen::Triplet<double>> ATriplets;
      for(int k = 0; k < EtE.outerSize(); ++k)
      {
          for (Eigen::SparseMatrix<double>::InnerIterator it(EtE, k); it; ++it)
          ATriplets.emplace_back(it.row(), it.col(), it.value());
      }

      for(int k = 0; k < Cpart.outerSize(); ++k)
      {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Cpart, k); it; ++it)
        {
          ATriplets.emplace_back(it.row() + EtE.rows(), it.col(), it.value());
          ATriplets.emplace_back(it.col(), it.row() + EtE.rows(), it.value());
        }
      }

      A.setFromTriplets(ATriplets.begin(), ATriplets.end());

      //Right-hand side with fixed values
      Eigen::VectorXd b = Eigen::VectorXd::Zero(EtE.rows() + Cpart.rows());
      b.segment(0, EtE.rows())= Epart.transpose() * M1 * (gamma + torhs);
      Eigen::VectorXd bfull = -Cfull * fixedValues;
      for(int k = 0; k < CpartRank; k++)
        bpart(k)=bfull(PIndices(k));
      b.segment(EtE.rows(), Cpart.rows()) = bpart;

      Eigen::SparseLU<Eigen::SparseMatrix<double> > lusolver;
      lusolver.compute(A);
      if(lusolver.info() != Eigen::Success)
        throw std::runtime_error("LU decomposition failed!");
      x = lusolver.solve(b);

      fullx = var2AllMat * x.head(numVars - alreadyFixed.sum()) + fixedValues;

//#ifndef NDEBUG
      if (Cfull.rows() != 0)
      {
          std::cout << "(Cfull * fullx).lpNorm<Infinity>(): " << (Cfull * fullx).lpNorm<Eigen::Infinity>() << std::endl;
      }
      std::cout << "Poisson error: " << (Efull * fullx - gamma).lpNorm<Eigen::Infinity>() << std::endl;
//#endif

      if((alreadyFixed - fixedMask).sum() == 0)
        break;
      
      double minIntDiff = std::numeric_limits<double>::max();
      int minIntDiffIndex = -1;
      for (int i = 0; i < numVars; i++)
      {
        if ((fixedMask(i)) && (!alreadyFixed(i)))
        {
          double currIntDiff =0;
          double func = fullx(i); 
            currIntDiff += std::fabs(func - std::round(func));
          if (currIntDiff < minIntDiff)
          {
            minIntDiff = currIntDiff;
            minIntDiffIndex = i;
          }
        }
      }

//#ifndef NDEBUG
      std::cout << "variable index: " << minIntDiffIndex << std::endl;
      std::cout << "variable Integer error: " << minIntDiff << std::endl;
//#endif
      
      if (minIntDiffIndex != -1)
      {
        alreadyFixed(minIntDiffIndex) = 1;
        double func = fullx(minIntDiffIndex) ;
        double funcInteger=std::round(func);
        fixedValues(minIntDiffIndex) = funcInteger;
      }
      
      xprev.resize(x.rows() - 1);
      varCounter = 0;
      for(int i = 0; i < numVars; i++)
        if (!alreadyFixed(i))
          xprev(varCounter++) = fullx(i);
      
      xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());
    }

    //the results are packets of N functions for each vertex, and need to be allocated for corners
    Eigen::VectorXd paramFuncsVec = pd.vertexTrans2CutMat * pd.symmMat * pd.intSpanMat * fullx;
    paramFuncsN.conservativeResize(cutV.rows(), pd.N);
    for(int i = 0; i < paramFuncsN.rows(); i++)
      paramFuncsN.row(i) << paramFuncsVec.segment(pd.N * i, N).transpose();
    
    paramFuncsd = fullx;
    
    //allocating per corner
    wholeCornerParamFuncsN.conservativeResize(wholeF.rows(), N*3);
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<3;j++)
        wholeCornerParamFuncsN.block(i, N*j, 1, N) = paramFuncsN.row(cutF(i,j));
    
    //exporting things to MATLAB
    
    Eigen::SparseMatrix<double> G;
    Eigen::MatrixXd FN;
    igl::per_face_normals(cutV, cutF, FN);
    branched_gradient(cutV,cutF, pd.N, G);

    Eigen::SparseMatrix<double> Gd=G*pd.vertexTrans2CutMat * pd.symmMat * pd.intSpanMat;
    Eigen::SparseMatrix<double> x2CornerMat=pd.vertexTrans2CutMat * pd.symmMat * pd.intSpanMat;
    
    igl::matlab::MatlabWorkspace mw;
    Eigen::VectorXi integerIndices(pd.integerVars.size()*pd.d);
    for(int i = 0; i < pd.integerVars.size(); i++)
      for (int j=0;j<pd.d;j++)
       integerIndices(pd.d*i+j) = pd.d*pd.integerVars(i)+j;
    
    mw.save(Efull, "A");
    mw.save(fullx,"x0");
    mw.save(rawField, "rawField");
    mw.save_index(pd.fixedIndices, "fixedIndices");
    mw.save(pd.fixedValues, "fixedValues");
    mw.save_index(pd.singularIndices, "singularIndices");
    mw.save(gamma, "b");
    mw.save(Cfull, "C");
    mw.save(Gd, "G");
    mw.save(FN, "FN");
    mw.save(pd.N, "N");
    mw.save(lengthRatio, "lengthRatio");
    mw.save(pd.d, "n");
    mw.save(pd.intSpanMat, "intSpanMat");
    mw.save(cutV, "V");
    mw.save_index(cutF, "F");
    mw.save(x2CornerMat, "x2CornerMat");
    mw.save_index(integerIndices, "integerIndices");
    mw.save(pd.vertexTrans2CutMat, "vertexTrans2CutMat");
    mw.save(pd.symmMat, "symmMat");
    mw.write(matname);
  }

}
  
#endif

  
