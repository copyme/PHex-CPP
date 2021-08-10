//
// Created by kacper on 5/3/19.
//

#ifndef INC_COMMON_H
#define INC_COMMON_H
#include <Eigen/Core>
#include <Eigen/Sparse>


bool polyhedral_face_normals(const Eigen::MatrixXd& V,
                             const Eigen::VectorXi& D,
                             const Eigen::MatrixXi& F,
                            Eigen::MatrixXd& faceNormals)
{
  using namespace Eigen;
  faceNormals.resize(D.rows(),3);
  for (int i=0;i<D.rows();i++){
    RowVector3d faceNormal; faceNormal<<0.0,0.0,0.0;
    for (int j=0;j<D(i);j++){
      RowVector3d vn=V.row(F(i,(j+D(i)-1)%D(i)));
      RowVector3d v0=V.row(F(i,j));
      RowVector3d v1=V.row(F(i,(j+1)%D(i)));
      faceNormal += (v1-v0).cross(vn-v0);
    }
    faceNormals.row(i)=faceNormal.normalized();
  }

  return true;
}



template<typename DerivedV,typename DerivedF,typename DerivedD,typename DerivedS>
void average_onto_vertices(const Eigen::MatrixBase<DerivedV> &V,
                                      const Eigen::MatrixBase<DerivedF> &F,
                                      const Eigen::MatrixBase<DerivedD> &D,
                                      const Eigen::MatrixBase<DerivedS> &S,
                                      Eigen::MatrixBase<DerivedS> &SV)
{
  SV = DerivedS::Zero(V.rows(),S.cols());
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> COUNT(V.rows());
  COUNT.setZero();
  for (int i = 0; i <F.rows(); ++i)
  {
    for (int j = 0; j<D(i); ++j)
    {
      SV.row(F(i,j)) += S.row(i);
      COUNT[F(i,j)] ++;
    }
  }
  for (int i = 0; i <V.rows(); ++i)
    SV.row(i) /= COUNT[i];
};


#endif //INC_COMMON_H
