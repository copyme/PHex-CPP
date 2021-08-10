#ifndef TRIANGLE_MESH_HEADER_FILE
#define TRIANGLE_MESH_HEADER_FILE


#include <iosfwd>
#include <vector>
#include <deque>
#include <set>
#include <cmath>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>

#include "Definitions.h"

class Vertex{
public:
  int ID;
  Point3D Coordinates;
  EPoint3D ECoordinates;
  int AdjHalfedge;
  
  bool isHex;
  bool Valid;
  
  Vertex():ID(-1), AdjHalfedge(-1), isHex(false), Valid(true){}
  ~Vertex()= default;
};

class Halfedge{
public:
  int ID;
  int Origin;
  int Next;
  int Prev;
  int Twin;
  int AdjFace;
  Eigen::VectorXd paramFuncs;
  std::vector<ENumber> exactParamFuncs;
  bool isHex;
  bool Valid;
  
  //Parametric function values
  int OrigHalfedge;
  int OrigParamFunc;  //the original parameteric function assooicated with this edge
  int prescribedAngleDiff;
  double prescribedAngle;  //the actual prescribed angle
  
  
  Halfedge():ID(-1), Origin(-1), Next(-1), Prev(-1), Twin(-1), AdjFace(-1), isHex(false), Valid(true), OrigHalfedge(-1), OrigParamFunc(-1), prescribedAngleDiff(-1), prescribedAngle(-1.0){}
  ~Halfedge()= default;
};


class Face{
public:
  int ID;
  int AdjHalfedge;
  bool Valid;
  
  Face():ID(-1), AdjHalfedge(-1), Valid(true){}
  ~Face()= default;
};

class Mesh{
private:
    std::vector<Vertex> Vertices;
    std::vector<Halfedge> Halfedges;
    std::vector<Face> Faces;
    std::vector<int> TransVertices;
    std::ostream & DebugLog;

public:
  explicit Mesh(std::ostream & log) : DebugLog (log) {}
  void UnifyEdges(int heindex);
  bool CheckMesh(bool checkHalfedgeRepetition, bool CheckTwinGaps, bool checkPureBoundary);
  void CleanMesh();
  void WalkBoundary(int &CurrEdge);
  void GenerateMesh(int numParamFuncs, Mesh& HexMesh);
  bool SimplifyHexMesh(int N);

  void Allocate(int NumofVertices, int NumofFaces, int NumofHEdges)
  {
    Vertices.resize(NumofVertices);
    Faces.resize(NumofFaces);
    Halfedges.resize(NumofHEdges);
  }
  
  //produces y = M*x
  static void exactSparseMult(const Eigen::SparseMatrix<int> & M, const std::vector<ENumber>& x,std::vector<ENumber>& y){
    y.resize(M.rows());
    
    for (auto & i : y)
      i=ENumber(0);
    
    for (int k=0; k<M.outerSize(); ++k)
      for (Eigen::SparseMatrix<int>::InnerIterator it(M,k); it; ++it)
        y[it.row()]+=ENumber((long)it.value())*x[it.col()];
  }
  
  void fromHedraDCEL(const Eigen::VectorXi& D,
                     const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const Eigen::MatrixXi& EV,
                     const Eigen::MatrixXi& FE,
                     const Eigen::MatrixXi& EF,
                     const Eigen::MatrixXi& EFi,
                     const Eigen::MatrixXd& FEs,
                     const Eigen::VectorXi& innerEdges,
                     const Eigen::VectorXi& VH,
                     const Eigen::MatrixXi& EH,
                     const Eigen::MatrixXi& FH,
                     const Eigen::VectorXi& HV,
                     const Eigen::VectorXi& HE,
                     const Eigen::VectorXi& HF,
                     const Eigen::VectorXi& nextH,
                     const Eigen::VectorXi& prevH,
                     const Eigen::VectorXi& twinH,
                     const Eigen::MatrixXd& cutV,
                     const Eigen::MatrixXi& cutF,
                     const Eigen::MatrixXd& paramFuncsd,
                     const int d,
                     const Eigen::SparseMatrix<int>& d2NMat,
                     const Eigen::SparseMatrix<int>& constraintMatInteger,
                     const Eigen::MatrixXd& cornerParamFuncs,
                     const Eigen::VectorXi& integerVars){
    
    using namespace std;
    using namespace CGAL;
    Vertices.resize(V.rows());
    Halfedges.resize(HE.rows());
    Faces.resize(F.rows());
    
    int N = cornerParamFuncs.cols()/3;
    
    for (int i=0;i<V.rows();i++){
      Vertices[i].Coordinates=Point3D(V(i,0), V(i,1), V(i,2));
      Vertices[i].AdjHalfedge=VH(i);
      Vertices[i].ID=i;
    }
    
    for (int i=0;i<HE.rows();i++){
      Halfedges[i].ID=i;
      Halfedges[i].Origin=HV(i);
      Halfedges[i].Next=nextH(i);
      Halfedges[i].Prev=prevH(i);
      Halfedges[i].Twin=twinH(i);
      Halfedges[i].AdjFace=HF(i);
    }
    
    for (int i=0;i<FH.rows();i++)
      for (int j=0;j<FH.cols();j++)
        Halfedges[FH(i,j)].paramFuncs = cornerParamFuncs.block(i, N*j, 1, N).transpose();
    
    for (int i=0;i<F.rows();i++){
      Faces[i].ID=i;
      Faces[i].AdjHalfedge=FH(i);
    }
    
    //computing exact rational corner values by quantizing the free variables d and then manually performing the sparse matrix multiplication
    
    //resolution is set to 10e-10 of bounding box of mesh
    vector<Point3D> coordList;
    for (auto & Vertice : Vertices)
      coordList.push_back(Vertice.Coordinates);
    
    Bbox_3 boundBox = CGAL::bbox_3  ( coordList.begin(), coordList.end());
    
    double minRange = 3276700.0;
    for (int i=0;i<2;i++)
      minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));
    
    double Resolution=1e8;

    vector<ENumber> exactParamFuncsd(paramFuncsd.size());
    for (int i=0;i<paramFuncsd.size();i++){
      exactParamFuncsd[i]=ENumber(Gmpz(round(paramFuncsd(i)* Resolution)), Gmpz(Resolution));
      cout<<"rounding diff of var "<<i<<" is "<<exactParamFuncsd[i].to_double()-paramFuncsd(i)<<endl;;
    }
    
    for (int i=0;i<integerVars.size();i++){
      for (int j=0;j<d;j++){
        exactParamFuncsd[d*integerVars(i)+j]=ENumber(Gmpz(round(paramFuncsd(d*integerVars(i)+j))));
        cout<<"rounding diff of integer var "<<d*integerVars(i)+j<<" is "<<exactParamFuncsd[d*integerVars(i)+j].to_double()-paramFuncsd(d*integerVars(i)+j)<<endl;
      }
    }
    
    vector<ENumber> exactParamFuncsVec;
    exactSparseMult(d2NMat, exactParamFuncsd,exactParamFuncsVec);
    
    vector<ENumber> constraintError;
    exactSparseMult(constraintMatInteger, exactParamFuncsd,constraintError);
    
    ENumber MaxError(0);
    
    for (auto & i : constraintError)
      if (i>MaxError)
        i=MaxError;
    
    cout<<"constraintMatInteger*exactParamFuncsd MaxError: "<<MaxError<<endl;
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    vector<vector<ENumber> > exactParamFuncsN(cutV.rows());
    for(int i = 0; i < cutV.rows(); i++){
      exactParamFuncsN[i].resize(N);
      for (int j = 0; j < N; j++) {
          exactParamFuncsN[i][j] = exactParamFuncsVec[N * i + j];
      }
    }
    
    //allocating per corner
    vector<vector<ENumber> > exactWholeCornerParamFuncsN(F.rows());
    
    for (int i=0;i<F.rows();i++){
      exactWholeCornerParamFuncsN[i].resize(N*3);
      for (int j=0;j<3;j++)
        for (int k=0;k<N;k++)
          exactWholeCornerParamFuncsN[i][N*j+k] = exactParamFuncsN[cutF(i,j)][k];
    }
    
    
    for (int i=0;i<FH.rows();i++)
      for (int j=0;j<FH.cols();j++){
        Halfedges[FH(i,j)].exactParamFuncs.resize(N);
        for (int k=0;k<N;k++)
          Halfedges[FH(i,j)].exactParamFuncs[k] = exactWholeCornerParamFuncsN[i][N * j + k];
      }
    
    //sanity check
    double maxError = -32767000.0;
    for (auto & Halfedge : Halfedges){
      for (int j=0;j<N;j++){
        double fromExact = Halfedge.exactParamFuncs[j].to_double();
        if (std::fabs(fromExact-Halfedge.paramFuncs[j])>maxError)
          maxError =std::fabs(fromExact-Halfedge.paramFuncs[j]);
      }
      
    }
    cout<<"maxError: "<<maxError<<endl;
  }
  
  
  //corner angles is per vertex in each F
  void toHedra(Eigen::MatrixXd& generatedV, Eigen::VectorXi& generatedD, Eigen::MatrixXi& generatedF, Eigen::MatrixXi& generatedFfuncNum, Eigen::MatrixXi& prescribedAnglesInt, Eigen::MatrixXd& cornerAngles){
    generatedV.resize(Vertices.size(),3);
    
    generatedD.resize(Faces.size());
    
    for (int i=0;i<Vertices.size();i++)
      generatedV.row(i)<<Vertices[i].Coordinates.x(), Vertices[i].Coordinates.y(),Vertices[i].Coordinates.z();
    
    for (int i=0;i<Faces.size();i++){
      int hebegin = Faces[i].AdjHalfedge;
      //reseting to first vertex
      int vCount=0;
      int heiterate=hebegin;
      do{
        vCount++;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      generatedD(i)=vCount;
    }
    
    generatedF.resize(Faces.size(),generatedD.maxCoeff());
    for (int i=0;i<Faces.size();i++){
       int hebegin = Faces[i].AdjHalfedge;
      int vCount=0;
      int heiterate=hebegin;
      do{
        generatedF(i,vCount++)=Halfedges[heiterate].Origin;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      
    }
    
    generatedFfuncNum.resize(Faces.size(),generatedD.maxCoeff());
    cornerAngles.resize(Faces.size(),generatedD.maxCoeff());
    prescribedAnglesInt.resize(Faces.size(),generatedD.maxCoeff());
    for (int i=0;i<Faces.size();i++){
      int hebegin = Faces[i].AdjHalfedge;
      int vCount=0;
      int heiterate=hebegin;
      do{
        generatedFfuncNum(i,vCount)=Halfedges[heiterate].OrigParamFunc;
        cornerAngles(i,vCount)=Halfedges[heiterate].prescribedAngle;
        prescribedAnglesInt(i,vCount++)=Halfedges[heiterate].prescribedAngleDiff;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
    }
    
    
  }
  
  Mesh()= default;
  ~Mesh()= default;
  
};


struct EdgeData{
  int ID;
  bool isHex;
  int OrigHalfedge;
  bool isBoundary;
  int funcNum;  //in case of function segment
  
  EdgeData():ID(-1), isHex(false), OrigHalfedge(-1), isBoundary(false), funcNum(-1){}
  ~EdgeData()= default;
};

namespace CGAL {

template <class ArrangementA, class ArrangementB, class ArrangementR>
class Arr_hex_overlay_traits
{
public:
  
  typedef typename ArrangementA::Face_const_handle    Face_handle_A;
  typedef typename ArrangementB::Face_const_handle    Face_handle_B;
  typedef typename ArrangementR::Face_handle          Face_handle_R;
  
  typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
  typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
  typedef typename ArrangementR::Vertex_handle          Vertex_handle_R;
  
  typedef typename ArrangementA::Halfedge_const_handle    Halfedge_handle_A;
  typedef typename ArrangementB::Halfedge_const_handle    Halfedge_handle_B;
  typedef typename ArrangementR::Halfedge_handle          Halfedge_handle_R;
  
public:
  
  virtual void create_face (Face_handle_A f1,
                            Face_handle_B f2,
                            Face_handle_R f) const
  {
    // Overlay the data objects associated with f1 and f2 and store the result
    // with f.
    f->set_data (f1->data());
 }
  
  //-1 - triangle vertex (non-hex vertex)
  //-2 - hex vertex
  virtual void	create_vertex ( Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  
  virtual void create_vertex ( Vertex_handle_A v1, Halfedge_handle_B e2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  
  virtual void create_vertex ( Vertex_handle_A v1, Face_handle_B f2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  virtual void create_vertex ( Halfedge_handle_A e1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  virtual void create_vertex ( Face_handle_A f1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  virtual void create_vertex ( Halfedge_handle_A e1, Halfedge_handle_B e2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  
  virtual void create_edge ( Halfedge_handle_A e1, Halfedge_handle_B e2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-2;
    data.OrigHalfedge=e1->data().OrigHalfedge;
    data.funcNum = e2->data().funcNum;
    e->set_data(data);
    e->twin()->set_data(data);
  }
  
  
  virtual void create_edge ( Halfedge_handle_A e1, Face_handle_B f2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-1;
    data.OrigHalfedge=e1->data().OrigHalfedge;
    data.funcNum = -1;  //triangle edge
    e->set_data(data);
    e->twin()->set_data(data);
  }
  
  virtual void create_edge ( Face_handle_A f1, Halfedge_handle_B e2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-2;
    data.funcNum = e2->data().funcNum;
    e->set_data(data);
    e->twin()->set_data(data);
  }
};

} //namespace CGAL



#endif
