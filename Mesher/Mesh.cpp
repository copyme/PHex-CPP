#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <igl/PI.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>
#include <cmath>


#include "Definitions.h"
#include "Mesh.h"

using namespace std;

typedef CGAL::Arr_linear_traits_2<EKernel>                    Traits_2;
typedef Traits_2::Point_2                                     Point_2;
typedef Traits_2::Segment_2                                   Segment_2;
typedef Traits_2::Ray_2                                       Ray_2;
typedef Traits_2::Line_2                                      Line_2;
typedef Traits_2::X_monotone_curve_2                          X_monotone_curve_2;
typedef CGAL::Arr_extended_dcel<Traits_2, int,EdgeData,int>   Dcel;
typedef CGAL::Arrangement_with_history_2<Traits_2, Dcel>      Arr_2;
typedef Arr_2::Face_iterator							      Face_iterator;
typedef Arr_2::Face_handle								      Face_handle;
typedef Arr_2::Edge_iterator							      Edge_iterator;
typedef Arr_2::Halfedge_iterator						      Halfedge_iterator;
typedef Arr_2::Vertex_iterator							      Vertex_iterator;
typedef Arr_2::Vertex_handle							      Vertex_handle;
typedef Arr_2::Halfedge_handle							      Halfedge_handle;
typedef Arr_2::Ccb_halfedge_circulator					      Ccb_halfedge_circulator;
typedef CGAL::Arr_hex_overlay_traits <Arr_2,Arr_2,Arr_2>      Overlay_traits;


void Mesh::GenerateMesh(const int numParamFuncs, Mesh& HexMesh)
{
  HexMesh.Vertices.clear();
  HexMesh.Halfedges.clear();
  HexMesh.Faces.clear();
  
  //resolution is set to 10e-6 of bounding box of mesh
  vector<Point3D> coordList;
  for (auto & Vertice : Vertices)
    coordList.push_back(Vertice.Coordinates);
  
  CGAL::Bbox_3 boundBox = CGAL::bbox_3  ( coordList.begin(), coordList.end());

  double Resolution=1e8;

  std::cout << "Face loop" << std::endl;
  for (auto & findex : Faces){
    
    //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face
    
    int ebegin=findex.AdjHalfedge;
    int eiterate=ebegin;

    vector<vector<ENumber> > funcValues(3);
    
    vector<ENumber> minFuncs(numParamFuncs);
    vector<ENumber> maxFuncs(numParamFuncs);
    for (int k=0;k<numParamFuncs;k++){
      minFuncs[k]=ENumber(327600);
      maxFuncs[k]=ENumber(-327600);
    }
    
    Arr_2 ParamArr,TriangleArr, FullArr;
    ebegin=findex.AdjHalfedge;
    eiterate=ebegin;
    int currVertex=0;
    do{
      for(int i=0;i<numParamFuncs;i++){
          if (Halfedges[eiterate].exactParamFuncs[i]>maxFuncs[i]) maxFuncs[i]=Halfedges[eiterate].exactParamFuncs[i];
          if (Halfedges[eiterate].exactParamFuncs[i]<minFuncs[i]) minFuncs[i]=Halfedges[eiterate].exactParamFuncs[i];
      }
      funcValues[currVertex++]=Halfedges[eiterate].exactParamFuncs;
      eiterate=Halfedges[eiterate].Next;
    }while (eiterate!=ebegin);
    
    //building the one-triangle arrangement
    ebegin=findex.AdjHalfedge;
    eiterate=ebegin;
    vector<EPoint2D> ETriPoints2D;
    vector<EPoint3D> ETriPoints3D;
    vector<EdgeData> EdgeDatas;
    ETriPoints2D.emplace_back(0,0);
    ETriPoints2D.emplace_back(1,0);
    ETriPoints2D.emplace_back(0,1);

    do{
      Point3D Position=Vertices[Halfedges[eiterate].Origin].Coordinates;
      ENumber x=ENumber(CGAL::Gmpz(std::round(Position.x()*Resolution)), CGAL::Gmpz(Resolution));
      ENumber y=ENumber(CGAL::Gmpz(std::round(Position.y()*Resolution)), CGAL::Gmpz(Resolution));
      ENumber z=ENumber(CGAL::Gmpz(std::round(Position.z()*Resolution)), CGAL::Gmpz(Resolution));
      ETriPoints3D.emplace_back(x,y,z);
      int DomEdge;
      
      if ((Halfedges[eiterate].Twin<0)||(Halfedges[eiterate].Twin>eiterate))
        DomEdge=eiterate;
      else
        DomEdge=Halfedges[eiterate].Twin;
      EdgeData ed; ed.OrigHalfedge=DomEdge;
      ed.isBoundary=(Halfedges[eiterate].Twin<0);
      EdgeDatas.push_back(ed);
      eiterate=Halfedges[eiterate].Next;
    }while(ebegin!=eiterate);
    
    for (int i=0;i<3;i++){
      X_monotone_curve_2 c =ESegment2D(ETriPoints2D[i],ETriPoints2D[(i+1)%3]);
      Halfedge_handle he=CGAL::insert_non_intersecting_curve(TriangleArr,c);
      he->set_data(EdgeDatas[i]);
      if (EdgeDatas[i].isBoundary)
        he->source()->data()=he->target()->data()=0;
      else
        he->source()->data()=he->target()->data()=1;
      
      he->twin()->set_data(EdgeDatas[i]);
    }
    
    for (Face_iterator fi= TriangleArr.faces_begin(); fi != TriangleArr.faces_end(); fi++){
      if (fi->is_unbounded())
        fi->data()=0;
      else
        fi->data()=1;
    }
    
    //creating the primal arrangement of lines
    vector<ELine2D> paramLines;
    vector<EDirection2D> isoDirections(numParamFuncs);
    int jumps = (numParamFuncs%2==0 ? 2 : 1);
    for (int funcIter=0;funcIter<numParamFuncs/jumps;funcIter++){
      
      vector<EInt> isoValues;
      EInt q,r;
      CGAL::div_mod(minFuncs[funcIter].numerator(), minFuncs[funcIter].denominator(), q, r);
      EInt minIsoValue = q + (r<0 ? -1 : 0);
      CGAL::div_mod(maxFuncs[funcIter].numerator(), maxFuncs[funcIter].denominator(), q, r);
      EInt maxIsoValue = q + (r<0  ? 0 : -1);

      for (EInt isoValue=minIsoValue-2;isoValue <=maxIsoValue+2;isoValue++){
        isoValues.push_back(isoValue);
      }
      
      //computing gradient of function in plane
      EVector2D e01 =ETriPoints2D[1] - ETriPoints2D[0];
      EVector2D e12 =ETriPoints2D[2] - ETriPoints2D[1];
      EVector2D e20 =ETriPoints2D[0] - ETriPoints2D[2];
      
      //a and b values of lines
      EVector2D gradVector = funcValues[2][funcIter]*EVector2D(-e01.y(), e01.x())+
      funcValues[0][funcIter]*EVector2D(-e12.y(), e12.x())+
      funcValues[1][funcIter]*EVector2D(-e20.y(), e20.x());
      
      isoDirections[funcIter]=EDirection2D(gradVector);
      
      //TODO: find c = z1*u+z2 of ax+by+c(u) ad then use it to generate all values between floor and ceil.
      
      //pinv of [a 1;b 1;c 1] is [           2*a - b - c,           2*b - a - c,           2*c - b - a]
      //[ b^2 - a*b + c^2 - a*c, a^2 - b*a + c^2 - b*c, a^2 - c*a + b^2 - c*b]/(2*a^2 - 2*a*b - 2*a*c + 2*b^2 - 2*b*c + 2*c^2)
      
      ENumber a=funcValues[0][funcIter];
      ENumber b=funcValues[1][funcIter];
      ENumber c=funcValues[2][funcIter];
      if ((a==b)&&(b==c))
        continue;  //that means a degenerate function on the triangle
      
      ENumber rhs[3];
      rhs[0]=-gradVector[0]*ETriPoints2D[0].x()-gradVector[1]*ETriPoints2D[0].y();
      rhs[1]=-gradVector[0]*ETriPoints2D[1].x()-gradVector[1]*ETriPoints2D[1].y();
      rhs[2]=-gradVector[0]*ETriPoints2D[2].x()-gradVector[1]*ETriPoints2D[2].y();
      
      ENumber invM[2][3];
      invM[0][0]= 2*a-b-c;
      invM[0][1]= 2*b-a-c;
      invM[0][2]= 2*c-b-a;
      invM[1][0]=b*b - a*b + c*c - a*c;
      invM[1][1]=a*a - b*a + c*c - b*c;
      invM[1][2]=a*a - c*a + b*b - c*b;
      for (auto & row : invM)
        for (auto & col : row)
          col/=(ENumber(2)*(a*a - a*b - a*c + b*b- b*c + c*c));
      
      ENumber x[2];
      x[0] = invM[0][0]*rhs[0]+invM[0][1]*rhs[1]+invM[0][2]*rhs[2];
      x[1] = invM[1][0]*rhs[0]+invM[1][1]*rhs[1]+invM[1][2]*rhs[2];
      
      //sanity check
      ENumber error[3];
      error[0]=x[0]*a+x[1]-rhs[0];
      error[1]=x[0]*b+x[1]-rhs[1];
      error[2]=x[0]*c+x[1]-rhs[2];

      
      //generating all lines
      for (auto & isoValue : isoValues){
        ENumber currc = isoValue*x[0]+x[1];
        paramLines.emplace_back(gradVector[0],gradVector[1],currc);
      }
    }

    CGAL::insert(ParamArr, paramLines.begin(), paramLines.end());
    
    //giving edge data to curve arrangement
    Arr_2::Edge_iterator                  eit;
    Arr_2::Originating_curve_iterator     ocit;
    for (eit = ParamArr.edges_begin(); eit != ParamArr.edges_end(); ++eit) {
      for (ocit = ParamArr.originating_curves_begin(eit);
           ocit != ParamArr.originating_curves_end(eit); ++ocit){
        EDirection2D thisDirection =  EDirection2D(ocit->supporting_line().a(), ocit->supporting_line().b());
        for (int paramIter = 0;paramIter<numParamFuncs/jumps;paramIter++){
          if ((thisDirection==isoDirections[paramIter])||(thisDirection==-isoDirections[paramIter])){
            eit->data().funcNum=paramIter;
            eit->twin()->data().funcNum=paramIter;
          }
        }
      }
    }

    //
    //creating the overlay

    Overlay_traits ot;
    overlay (TriangleArr, ParamArr, FullArr, ot);

    for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
      if (!fi->data())
        continue;  //not participating
      
      Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
      Ccb_halfedge_circulator heiterate=hebegin;
      do{
        
        if (heiterate->source()->data()<0){  //new vertex
          Vertex NewVertex;
          NewVertex.ID=HexMesh.Vertices.size();
          NewVertex.isHex=(heiterate->source()->data()==-2);
          HexMesh.Vertices.push_back(NewVertex);
          heiterate->source()->data()=NewVertex.ID;
        }
        
        if (heiterate->data().ID<0){  //new halfedge
          Halfedge NewHalfedge;
          NewHalfedge.ID=HexMesh.Halfedges.size();
          NewHalfedge.isHex=(heiterate->data().ID==-2);
          NewHalfedge.Origin=heiterate->source()->data();
          NewHalfedge.OrigHalfedge=heiterate->data().OrigHalfedge;
          NewHalfedge.OrigParamFunc=heiterate->data().funcNum;
          //cout<<"NewHalfedge.OrigParamFunc :"<<NewHalfedge.OrigParamFunc<<endl;
          HexMesh.Vertices[heiterate->source()->data()].AdjHalfedge=NewHalfedge.ID;
          HexMesh.Halfedges.push_back(NewHalfedge);
          heiterate->data().ID=NewHalfedge.ID;
        }
        heiterate++;
      }while(heiterate!=hebegin);
      

      //now assigning nexts and prevs
      do{
        HexMesh.Halfedges[heiterate->data().ID].Next=heiterate->next()->data().ID;
        HexMesh.Halfedges[heiterate->data().ID].Prev=heiterate->prev()->data().ID;
        HexMesh.Halfedges[heiterate->data().ID].Twin=heiterate->twin()->data().ID;
        if (heiterate->twin()->data().ID>=0)
          HexMesh.Halfedges[heiterate->twin()->data().ID].Twin=heiterate->data().ID;
        
        heiterate++;
      }while (heiterate!=hebegin);
    }
    
    //constructing the actual vertices
    for (Vertex_iterator vi=FullArr.vertices_begin();vi!=FullArr.vertices_end();vi++){
      if (vi->data()<0)
        continue;
      
      //finding out barycentric coordinates
      ENumber BaryValues[3] = { 0 };
      ENumber Sum=0;
      for (int i=0;i<3;i++){
        ETriangle2D t(vi->point(), ETriPoints2D[(i+1)%3], ETriPoints2D[(i+2)%3]);
        BaryValues[i]=t.area();
        Sum+=BaryValues[i];
      }
      for (auto & BaryValue : BaryValues)
        BaryValue/=Sum;

      EPoint3D ENewPosition(0,0,0);
      for (int i=0;i<3;i++)
        ENewPosition=ENewPosition+(ETriPoints3D[i]-CGAL::ORIGIN)*BaryValues[i];
      
      Point3D NewPosition(to_double(ENewPosition.x()), to_double(ENewPosition.y()), to_double(ENewPosition.z()));
      
      HexMesh.Vertices[vi->data()].Coordinates=NewPosition;
      HexMesh.Vertices[vi->data()].ECoordinates=ENewPosition;
      
    }
    
    for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
      if (!fi->data())
        continue;
      
      int FaceSize=0;
      Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
      Ccb_halfedge_circulator heiterate=hebegin;
      do{ FaceSize++;  heiterate++; }while(heiterate!=hebegin);
      int CurrPlace=0;
      
      Face NewFace;
      NewFace.ID=HexMesh.Faces.size();
      //NewFace.NumVertices=FaceSize;
      NewFace.AdjHalfedge=hebegin->data().ID;
      
      do{
        HexMesh.Halfedges[heiterate->data().ID].AdjFace=NewFace.ID;
        heiterate++;
      }while(heiterate!=hebegin);
      HexMesh.Faces.push_back(NewFace);
    }
  }
  
  //devising angles from differences in functions
  int ratio = (numParamFuncs%2==0 ? 1 : 2);
  for (int hi=0;hi<HexMesh.Halfedges.size();hi++){
    if ((HexMesh.Halfedges[hi].OrigParamFunc==-1)||(HexMesh.Halfedges[HexMesh.Halfedges[hi].Prev].OrigParamFunc==-1))
      HexMesh.Halfedges[hi].prescribedAngleDiff=-1;  //one of the edges is a triangle edge, and it will be devised later.
    else{
      int func1 =(ratio*(HexMesh.Halfedges[hi].OrigParamFunc)) % (numParamFuncs/(3-ratio));
      int func2 =(ratio*(HexMesh.Halfedges[HexMesh.Halfedges[hi].Prev].OrigParamFunc)) % (numParamFuncs/(3-ratio));
      HexMesh.Halfedges[hi].prescribedAngleDiff=(func2-func1+numParamFuncs/(3-ratio)) % (numParamFuncs/(3-ratio));
    }
    
  }
}


struct PointPair{
  int Index1, Index2;
  ENumber Distance;
  
  PointPair(int i1, int i2, ENumber d):Index1(i1), Index2(i2), Distance(d){}
  ~PointPair()= default;
  
  bool operator<(const PointPair& pp) const {
    if (Distance>pp.Distance) return false;
    if (Distance<pp.Distance) return true;
    
    if (Index1>pp.Index1) return false;
    if (Index1<pp.Index1) return true;
    
    if (Index2>pp.Index2) return false;
    if (Index2<pp.Index2) return true;
    
    return false;
    
  }
};

vector<pair<int,int>> FindVertexMatch(ostream& DebugLog, vector<EPoint3D>& Set1, vector<EPoint3D>& Set2)
{
  set<PointPair> PairSet;
  for (int i=0;i<Set1.size();i++)
    for (int j=0;j<Set2.size();j++)
      PairSet.insert(PointPair(i,j,squared_distance(Set1[i],Set2[j])));
  
  DebugLog<<"Matching set ";
  for (auto & i : Set1)
    DebugLog<<i.x().to_double()<<" "<<i.y().to_double()<<" "<<i.z().to_double()<<", "<<endl;
  
  DebugLog<<"to set" <<endl;
  for (auto & i : Set2)
    DebugLog<<i.x().to_double()<<" "<<i.y().to_double()<<" "<<i.z().to_double()<<", "<<endl;
  
  if (Set1.size()!=Set2.size())  //should not happen anymore
    DebugLog<<"The two sets are of different sizes!! "<<endl;
  
  //adding greedily legal connections until graph is full
  vector<bool> Set1Connect(Set1.size());
  vector<bool> Set2Connect(Set2.size());
  
  vector<pair<int, int> > Result;
  
  for (int i=0;i<Set1.size();i++)
    Set1Connect[i]=false;
  
  for (int i=0;i<Set2.size();i++)
    Set2Connect[i]=false;

  int NumConnected=0;
  
  //categorically match both ends
  Result.emplace_back(0,0);
  Result.emplace_back(Set1.size()-1, Set2.size()-1);
  for (const auto& CurrPair : PairSet)
  {
    //checking legality - if any of one's former are connected to ones latters or vice versa
    bool FoundConflict=false;
    for (auto & i : Result){
      if (((i.first>CurrPair.Index1)&&(i.second<CurrPair.Index2))||
          ((i.first<CurrPair.Index1)&&(i.second>CurrPair.Index2))){
        FoundConflict=true;
        break;
      }
    }
    
    if (FoundConflict)
      continue;
    
    //if both are already matched, this matching is redundant
    if ((Set1Connect[CurrPair.Index1])&&(Set2Connect[CurrPair.Index2]))
      continue;  //there is no reason for this matching
    
    //otherwise this edge is legal, so add it
    Result.emplace_back(CurrPair.Index1, CurrPair.Index2);
    if (!Set1Connect[CurrPair.Index1]) NumConnected++;
    if (!Set2Connect[CurrPair.Index2]) NumConnected++;
    Set1Connect[CurrPair.Index1]=Set2Connect[CurrPair.Index2]=true;
  }
  
  for (int i=0;i<Set1.size();i++)
    if (!Set1Connect[i])
      DebugLog<<"Relative Vertex "<<i<<" in Set1 is unmatched!"<<endl;
  
  for (int i=0;i<Set2.size();i++)
    if (!Set2Connect[i])
      DebugLog<<"Relative Vertex "<<i<<" in Set2 is unmatched!"<<endl;

  DebugLog<<"matching points are ";
  for (auto & i : Result){
    DebugLog<<"("<<i.first<<","<<i.second<<") with dist "<<squared_distance(Set1[i.first],Set2[i.second])<<endl;
    if (squared_distance(Set1[i.first],Set2[i.second])>0)
      DebugLog<<"Distance is abnormally not zero!"<<endl;
  }

  return Result;
  
}

struct TwinFinder{
  int index;
  int v1,v2;
  
  TwinFinder(int i, int vv1, int vv2):index(i), v1(vv1), v2(vv2){}
  ~TwinFinder()= default;
  
  bool operator<(const TwinFinder& tf) const
  {
    if (v1<tf.v1) return false;
    if (v1>tf.v1) return true;
    
    if (v2<tf.v2) return false;
    if (v2>tf.v2) return true;
    
    return false;
  }
};

//edge has to be sourcing a 2-valence vertex!!
void Mesh::UnifyEdges(int heindex)
{	
  if (Halfedges[heindex].Twin<0){
    DebugLog<<"Unifying edge "<<heindex<<" with source "<<Halfedges[heindex].Origin<<"\n";
    DebugLog<<"new source "<<Halfedges[Halfedges[heindex].Prev].Origin<<"\n";
    
    DebugLog<<"old positions: ("<<Vertices[Halfedges[Halfedges[heindex].Prev].Origin].Coordinates<<")->("<<Vertices[Halfedges[heindex].Origin].Coordinates<<")->("<<Vertices[Halfedges[Halfedges[heindex].Next].Origin].Coordinates<<")\n";
  }
  //adjusting source
  Vertices[Halfedges[heindex].Origin].Valid=false;
  Halfedges[heindex].Origin=Halfedges[Halfedges[heindex].Prev].Origin;
  if (Halfedges[heindex].prescribedAngleDiff==-1)
    Halfedges[heindex].prescribedAngleDiff=Halfedges[Halfedges[heindex].Prev].prescribedAngleDiff;
  Vertices[Halfedges[heindex].Origin].AdjHalfedge=heindex;
  
  Faces[Halfedges[heindex].AdjFace].AdjHalfedge=Halfedges[heindex].Next;
  //Faces[Halfedges[heindex].AdjFace].NumVertices--;
  
  if (Halfedges[heindex].Twin<0)
    DebugLog<<"Removing edge "<<Halfedges[heindex].Prev<<" and connecting "<<Halfedges[Halfedges[heindex].Prev].Prev<<"->"<<heindex<<"\n";
  
  //adjusting halfedges
  Halfedges[Halfedges[heindex].Prev].Valid=false;
  Halfedges[heindex].Prev=Halfedges[Halfedges[heindex].Prev].Prev;
  Halfedges[Halfedges[heindex].Prev].Next=heindex;
  
  if (Halfedges[heindex].Twin<0)
    DebugLog<<"new positions: ("<<Vertices[Halfedges[heindex].Origin].Coordinates<<")->("<<Vertices[Halfedges[Halfedges[heindex].Next].Origin].Coordinates<<")\n";
  
  
  //adjusting twin, if exists
  if (Halfedges[heindex].Twin>=0){
    if (Halfedges[Halfedges[heindex].Twin].prescribedAngleDiff==-1)
      Halfedges[Halfedges[heindex].Twin].prescribedAngleDiff=Halfedges[Halfedges[Halfedges[heindex].Twin].Next].prescribedAngleDiff;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Valid=false;
    Halfedges[Halfedges[heindex].Twin].Next=Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Next;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[heindex].Twin;
    Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
  }
}

void Mesh::CleanMesh()
{
  DebugLog<<"Cleaning mesh (removing invalid components"<<endl;
  //removing nonvalid vertices
  vector<int> TransVertices(Vertices.size());
  vector<Vertex> NewVertices;
  for (int i=0;i<Vertices.size();i++){
    if (!Vertices[i].Valid)
      continue;
    
    Vertex NewVertex=Vertices[i];
    NewVertex.ID=NewVertices.size();
    NewVertices.push_back(NewVertex);
    TransVertices[i]=NewVertex.ID;
  }
  
   DebugLog<<"ok here 1"<<endl;
  
  Vertices=NewVertices;
  for (auto & Halfedge : Halfedges)
    Halfedge.Origin=TransVertices[Halfedge.Origin];
  
  DebugLog<<"ok here 2"<<endl;

  //removing nonvalid faces
  vector<Face> NewFaces;
  vector<int> TransFaces(Faces.size());
  for (int i=0;i<Faces.size();i++){
    if (!Faces[i].Valid)
      continue;
    
    Face NewFace=Faces[i];
    NewFace.ID=NewFaces.size();
    NewFaces.push_back(NewFace);
    TransFaces[i]=NewFace.ID;
  }
  Faces=NewFaces;
  for (auto & Halfedge : Halfedges)
    Halfedge.AdjFace=TransFaces[Halfedge.AdjFace];
  
   DebugLog<<"ok here 4"<<endl;
  
  //removing nonvalid halfedges
  vector<Halfedge> NewHalfedges;
  vector<int> TransHalfedges(Halfedges.size());
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    
    Halfedge NewHalfedge=Halfedges[i];
    NewHalfedge.ID=NewHalfedges.size();
    NewHalfedges.push_back(NewHalfedge);
    TransHalfedges[i]=NewHalfedge.ID;
  }
  
   DebugLog<<"ok here 5"<<endl;
  
  Halfedges=NewHalfedges;
  for (auto & Face : Faces)
    Face.AdjHalfedge=TransHalfedges[Face.AdjHalfedge];
  
   DebugLog<<"ok here 6"<<endl;
  
  for (auto & Vertice : Vertices)
    Vertice.AdjHalfedge=TransHalfedges[Vertice.AdjHalfedge];
  
   DebugLog<<"ok here 7"<<endl;
  
  for (auto & Halfedge : Halfedges){
    if (Halfedge.Twin!=-1)
      Halfedge.Twin=TransHalfedges[Halfedge.Twin];
    Halfedge.Next=TransHalfedges[Halfedge.Next];
    Halfedge.Prev=TransHalfedges[Halfedge.Prev];
  }
  
}

bool Mesh::CheckMesh(bool checkHalfedgeRepetition, bool CheckTwinGaps, bool checkPureBoundary)
{
  for (int i=0;i<Vertices.size();i++){
    if (!Vertices[i].Valid)
      continue;
    
    if (Vertices[i].AdjHalfedge==-1){
      DebugLog<<"Valid Vertex "<<i<<" points to non-valid value -1 "<<endl;
      return false;
    }
    
    if (!Halfedges[Vertices[i].AdjHalfedge].Valid){
      DebugLog<<"Valid Vertex "<<i<<" points to non-valid AdjHalfedge "<<Vertices[i].AdjHalfedge<<endl;
      return false;
    }
    
    
    if (Halfedges[Vertices[i].AdjHalfedge].Origin!=i){
      DebugLog<<"Adjacent Halfedge "<<Vertices[i].AdjHalfedge<<" of vertex "<<i<<"does not point back"<<endl;
      return false;
    }
    
  }
  
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    
    
    if (Halfedges[i].Next==-1){
      DebugLog<<"Valid Halfedge "<<i<<"points to Next non-valid value -1"<<endl;
      return false;
    }
    
    if (Halfedges[i].Prev==-1){
      DebugLog<<"Valid Halfedge "<<i<<"points to Prev non-valid value -1"<<endl;
      return false;
    }
    
    
    if (Halfedges[i].Origin==-1){
      DebugLog<<"Valid Halfedge "<<i<<"points to Origin non-valid value -1"<<endl;
      return false;
    }
    
    if (Halfedges[i].AdjFace==-1){
      DebugLog<<"Valid Halfedge "<<i<<"points to AdjFace non-valid value -1"<<endl;
      return false;
    }
    
    if (Halfedges[Halfedges[i].Next].Prev!=i){
      DebugLog<<"Halfedge "<<i<<"Next is "<<Halfedges[i].Next<<" which doesn't point back as Prev"<<endl;
      return false;
    }
    
    
    if (Halfedges[Halfedges[i].Prev].Next!=i){
      DebugLog<<"Halfedge "<<i<<"Prev is "<<Halfedges[i].Prev<<" which doesn't point back as Next"<<endl;
      return false;
    }
    
    if (!Vertices[Halfedges[i].Origin].Valid){
      DebugLog<<"The Origin of halfedges "<<i<<", vertex "<<Halfedges[i].Origin<<" is not valid"<<endl;
      return false;
    }
    
    if (!Faces[Halfedges[i].AdjFace].Valid){
      DebugLog<<"The face of halfedges "<<i<<", face "<<Halfedges[i].AdjFace<<" is not valid"<<endl;
      return false;
    }
    
    if (Halfedges[Halfedges[i].Next].Origin==Halfedges[i].Origin){  //a degenerate edge{
      DebugLog<<"Halfedge "<<i<<" with twin"<<Halfedges[i].Twin<<" is degenerate with vertex "<<Halfedges[i].Origin<<endl;
      return false;
    }
    
    if (Halfedges[i].Twin>=0){
      if (Halfedges[Halfedges[i].Twin].Twin!=i){
        DebugLog<<"Halfedge "<<i<<"twin is "<<Halfedges[i].Twin<<" which doesn't point back"<<endl;
        return false;
      }
      
      if (!Halfedges[Halfedges[i].Twin].Valid){
        DebugLog<<"halfedge "<<i<<" is twin with invalid halfedge"<<Halfedges[i].Twin<<endl;
        return false;
      }
    }
    
    if (!Halfedges[Halfedges[i].Next].Valid){
      DebugLog<<"halfedge "<<i<<" has next invalid halfedge"<<Halfedges[i].Next<<endl;
      return false;
    }
    
    if (!Halfedges[Halfedges[i].Prev].Valid){
      DebugLog<<"halfedge "<<i<<" has prev invalid halfedge"<<Halfedges[i].Prev<<endl;
      return false;
    }
    
    if (Halfedges[i].isHex){  //checking that it is not left alone
      if (Halfedges[i].Prev==Halfedges[i].Twin){
        DebugLog<<"Hex halfedge "<<i<<" has Halfedge "<<Halfedges[i].Prev<<" and both prev and twin"<<endl;
        return false;
      }
      
      
      if (Halfedges[i].Next==Halfedges[i].Twin){
        DebugLog<<"Hex halfedge "<<i<<" has Halfedge "<<Halfedges[i].Next<<" and both next and twin"<<endl;
        return false;
      }
    }
  }
  
  vector<set<int>> halfedgesinFace(Faces.size());
  vector<set<int>> verticesinFace(Faces.size());
  for (int i=0;i<Faces.size();i++){
    if (!Faces[i].Valid)
      continue;
    int hebegin=Faces[i].AdjHalfedge;
    int heiterate=hebegin;
    int NumEdges=0;
    int actualNumVertices=0;
    
    do{
      if (verticesinFace[i].find(Halfedges[heiterate].Origin)!=verticesinFace[i].end())
        DebugLog<<"Warning: Vertex "<<Halfedges[heiterate].Origin<<" appears more than once in face "<<i<<endl;
      
      verticesinFace[i].insert(Halfedges[heiterate].Origin);
      halfedgesinFace[i].insert(heiterate);
      actualNumVertices++;
      if (!Halfedges[heiterate].Valid)
        return false;
      
      if (Halfedges[heiterate].AdjFace!=i){
        DebugLog<<"Face "<<i<<" has halfedge "<<heiterate<<" that does not point back"<<endl;
        return false;
      }
    
      heiterate=Halfedges[heiterate].Next;
      NumEdges++;
      if (NumEdges>Halfedges.size()){
        DebugLog<<"Infinity loop!"<<endl;
        return false;
      }
       
      
    }while (heiterate!=hebegin);
  }
  
  //checking if all halfedges that relate to a face are part of its recognized chain (so no floaters)
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    int currFace = Halfedges[i].AdjFace;
    if (halfedgesinFace[currFace].find(i)==halfedgesinFace[currFace].end()){
      DebugLog<<"Halfedge "<<i<<" is floating in face "<<currFace<<endl;
      return false;
    }
  }
  
  //check if mesh is a manifold: every halfedge appears only once
  if (checkHalfedgeRepetition){
    std::set<TwinFinder> HESet;
    for (int i=0;i<Halfedges.size();i++){
      if (!Halfedges[i].Valid)
        continue;
      auto HESetIterator=HESet.find(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
      if (HESetIterator!=HESet.end()){
        DebugLog<<"Warning: the halfedge ("<<Halfedges[i].Origin<<","<<Halfedges[Halfedges[i].Next].Origin<<") appears at least twice in the mesh"<<endl;
        DebugLog<<"for instance halfedges "<<i<<" and "<<HESetIterator->index<<endl;
        return false;
        //return false;
      }else{
        HESet.insert(TwinFinder(i,Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
      }
    }
  }
  
  if (CheckTwinGaps){
    std::set<TwinFinder> HESet;
    //checking if there is a gap: two halfedges that share the same opposite vertices but do not have twins
    for (int i=0;i<Halfedges.size();i++){
      if (!Halfedges[i].Valid)
        continue;
      
      auto HESetIterator=HESet.find(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
      if (HESetIterator==HESet.end()){
        HESet.insert(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
        continue;
      }
      
      HESetIterator=HESet.find(TwinFinder(i, Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
      if (HESetIterator!=HESet.end()){
        
        if (Halfedges[i].Twin==-1){
          DebugLog<<"Halfedge "<<i<<"has no twin although halfedge "<<HESetIterator->index<<" can be a twin"<<endl;
          return false;
        }
        if (Halfedges[HESetIterator->index].Twin==-1){
          DebugLog<<"Halfedge "<<HESetIterator->index<<"has no twin although halfedge "<<i<<" can be a twin"<<endl;
          return false;
        }
      }
    }
  }
  
  //checking if there are pure boundary faces (there shouldn't be)
  if (checkPureBoundary){
   for (int i=0;i<Halfedges.size();i++){
     if (!Halfedges[i].Valid)
       continue;
     if ((Halfedges[i].Twin<0)&&(Halfedges[i].isHex))
       DebugLog<<"WARNING: Halfedge "<<i<<" is a hex edge without twin!"<<endl;
     
     if (Halfedges[i].Twin>0)
       continue;
     
     bool pureBoundary=true;
     int hebegin=i;
     int heiterate=hebegin;
     do{
       if (Halfedges[heiterate].Twin>0){
         pureBoundary=false;
         break;
       }
       heiterate = Halfedges[heiterate].Next;
     }while (heiterate!=hebegin);
     if (pureBoundary){
       DebugLog<<"Face "<<Halfedges[i].AdjFace<<"is a pure boundary face!"<<endl;
       return false;
     }
   }
     
    //checking for latent valence 2 faces
    vector<int> Valences(Vertices.size());
    for (int i=0;i<Vertices.size();i++)
      Valences[i]=0;
    
    for (auto & Halfedge : Halfedges){
      if (Halfedge.Valid){
        Valences[Halfedge.Origin]++;
        if (Halfedge.Twin<0)  //should account for the target as well
          Valences[Halfedges[Halfedge.Next].Origin]++;
      }
    }
     
    int countThree;
    for (int i=0;i<Faces.size();i++){
      if (!Faces[i].Valid)
        continue;
      countThree=0;
      int hebegin = Faces[i].AdjHalfedge;
      int heiterate=hebegin;
      int numEdges=0;
      do{
        if (Valences[Halfedges[heiterate].Origin]>2)
          countThree++;
        heiterate=Halfedges[heiterate].Next;
        numEdges++;
        if (numEdges>Halfedges.size()){
          DebugLog<<"Infinity loop in face "<<i<<"!"<<endl;
          return false;
        }
      }while (heiterate!=hebegin);
      if (countThree<3){
        DebugLog<<"Face "<<i<<" is a latent valence 2 face!"<<endl;
        DebugLog<<"Its vertices are "<<endl;
        do{
          DebugLog<<"Vertex "<<Halfedges[heiterate].Origin<<" halfedge "<<heiterate<<" valence "<<Valences[Halfedges[heiterate].Origin]<<endl;
          
          if (Valences[Halfedges[heiterate].Origin]>2)
            countThree++;
          heiterate=Halfedges[heiterate].Next;
          numEdges++;
          if (numEdges>Halfedges.size()){
            DebugLog<<"Infinity loop in face "<<i<<"!"<<endl;
            return false;
          }
        }while (heiterate!=hebegin);
      }
    }
  }
  
  DebugLog<<"Mesh is clear according to given checks"<<endl;
  return true;  //most likely the mesh is solid
}

void Mesh::WalkBoundary(int &CurrEdge)
{
  do{
    CurrEdge=Halfedges[CurrEdge].Next;
    if (Halfedges[CurrEdge].Twin<0)
      break;  //next boundary over a 2-valence vertex
    CurrEdge=Halfedges[CurrEdge].Twin;
  }while (Halfedges[CurrEdge].Twin>=0);
  
}

typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;


bool Mesh::SimplifyHexMesh(int N)
{
  //unifying vertices which are similar

  if (!CheckMesh(false, false, false))
    return false;
  
  int MaxOrigHE=-3276700.0;
  for (auto & Halfedge : Halfedges)
    MaxOrigHE=std::max(MaxOrigHE, Halfedge.OrigHalfedge);
  
  vector<bool> visitedOrig(MaxOrigHE+1);
  for (int i=0;i<MaxOrigHE+1;i++) visitedOrig[i]=false;
  for (int i=0;i<Halfedges.size();i++){
    if (Halfedges[i].OrigHalfedge<0)
      continue;
    if (visitedOrig[Halfedges[i].OrigHalfedge])
      continue;
    
    int hebegin = i;
    int heiterate = hebegin;
    DebugLog<<"Walking original triangle boundary"<<endl;
    do{
      visitedOrig[Halfedges[heiterate].OrigHalfedge]=true;
      DebugLog<<"Walking boundary halfedge "<<heiterate<<" with vertex "<<Halfedges[heiterate].Origin<<" on original halfedge "<<Halfedges[heiterate].OrigHalfedge<<endl;
      WalkBoundary(heiterate);
    }while (heiterate!=hebegin);
    
  }
  
  vector< vector<int> > BoundEdgeCollect1(MaxOrigHE+1);
  vector< vector<int> > BoundEdgeCollect2(MaxOrigHE+1);
  vector<bool> Marked(Halfedges.size());
  for (int i=0;i<Halfedges.size();i++) Marked[i]=false;
  //finding out vertex correspondence along twin edges of the original mesh by walking on boundaries
  for (int i=0;i<Halfedges.size();i++){
    if ((Halfedges[i].OrigHalfedge<0)||(Marked[i]))
      continue;
    
    //find the next beginning of a boundary
    int PrevOrig;
    int CurrEdge=i;
    do{
      PrevOrig=Halfedges[CurrEdge].OrigHalfedge;
      WalkBoundary(CurrEdge);
    }while(PrevOrig==Halfedges[CurrEdge].OrigHalfedge);
    
    //filling out strips of boundary with the respective attached original halfedges
    int BeginEdge=CurrEdge;
    vector<pair<int,int> > CurrEdgeCollect;
    do{
      CurrEdgeCollect.emplace_back(Halfedges[CurrEdge].OrigHalfedge, CurrEdge);
      Marked[CurrEdge]=true;
      WalkBoundary(CurrEdge);
    }while (CurrEdge!=BeginEdge);
    
    PrevOrig=-1000;
    bool In1;
    for (auto & j : CurrEdgeCollect){
      if (j.first!=PrevOrig)
        In1=BoundEdgeCollect1[j.first].empty();
      
      if (In1)
        BoundEdgeCollect1[j.first].push_back(j.second);
      else
        BoundEdgeCollect2[j.first].push_back(j.second);
      PrevOrig=j.first;
    }
  }
  
  //editing the edges into two vector lists per associated original edge
  vector< vector<int> > VertexSets1(MaxOrigHE+1), VertexSets2(MaxOrigHE+1);
  for (int i=0;i<MaxOrigHE+1;i++){
    for (int j=0;j<BoundEdgeCollect1[i].size();j++)
      VertexSets1[i].push_back(Halfedges[BoundEdgeCollect1[i][j]].Origin);
    
    if (!BoundEdgeCollect1[i].empty())
      VertexSets1[i].push_back(Halfedges[Halfedges[BoundEdgeCollect1[i][BoundEdgeCollect1[i].size()-1]].Next].Origin);
    
    for (int j=0;j<BoundEdgeCollect2[i].size();j++)
      VertexSets2[i].push_back(Halfedges[BoundEdgeCollect2[i][j]].Origin);
    
    if (!BoundEdgeCollect2[i].empty())
      VertexSets2[i].push_back(Halfedges[Halfedges[BoundEdgeCollect2[i][BoundEdgeCollect2[i].size()-1]].Next].Origin);
    
    std::reverse(VertexSets2[i].begin(),VertexSets2[i].end());
  }
  
  //finding out vertex matches
  vector<pair<int, int> > VertexMatches;
  for (int i=0;i<MaxOrigHE+1;i++){
    vector<EPoint3D> PointSet1(VertexSets1[i].size());
    vector<EPoint3D> PointSet2(VertexSets2[i].size());
    for (int j=0;j<PointSet1.size();j++)
      PointSet1[j]=Vertices[VertexSets1[i][j]].ECoordinates;
    
    for (int j=0;j<PointSet2.size();j++)
      PointSet2[j]=Vertices[VertexSets2[i][j]].ECoordinates;
    
    vector<pair<int, int> > CurrMatches;
    if ((!PointSet1.empty())&&(!PointSet2.empty()))
      CurrMatches=FindVertexMatch(DebugLog, PointSet1, PointSet2);
    
    for (auto & CurrMatche : CurrMatches){
      CurrMatche.first =VertexSets1[i][CurrMatche.first];
      CurrMatche.second=VertexSets2[i][CurrMatche.second];
      DebugLog<<"Vertex "<<CurrMatche.first<<" is matched with vertex "<<CurrMatche.second<<endl;
    }
    
    VertexMatches.insert(VertexMatches.end(), CurrMatches.begin(), CurrMatches.end() );
  }
  
  //finding connected components, and uniting every component into a random single vertex in it (it comes out the last mentioned)
  Graph MatchGraph;
  for (int i=0;i<Vertices.size();i++)
    add_vertex(MatchGraph);
  for (auto & VertexMatche : VertexMatches)
    add_edge(VertexMatche.first, VertexMatche.second, MatchGraph);
  
  double MaxDist=-327670000.0;
  for (auto & VertexMatche : VertexMatches)
    MaxDist=std::max(MaxDist, (Vertices[VertexMatche.first].Coordinates-Vertices[VertexMatche.second].Coordinates).squared_length());
  
  cout<<"Max matching distance: "<<MaxDist<<endl;
  
  TransVertices.resize(Vertices.size());
  int NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);
  for (int i=0;i<NumNewVertices;i++){
    DebugLog<<"TransVertices group "<<i<<": ";
    for (int j=0;j<TransVertices.size();j++)
      if (TransVertices[j]==i)
        DebugLog<<j<<", ";
    DebugLog<<endl;
  }
  
  for (int i=0;i<Faces.size();i++){
    int hebegin = Faces[i].AdjHalfedge;
    int heiterate = hebegin;
    DebugLog<<"Face "<<i<<" is initially of ";
    do{
      DebugLog<<heiterate<<"["<<Halfedges[heiterate].Origin<<"], ";
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
    DebugLog<<endl;
  }
  
  for (int i=0;i<Halfedges.size();i++){
    DebugLog<<"Halfedge "<<i<<" is initially a ";
    if (Halfedges[i].isHex)
      DebugLog<<"Hex edge with ";
    else
      DebugLog<<"Triangle edge with ";
    
    DebugLog<<"Origin: "<<Halfedges[i].Origin<<"("<<TransVertices[Halfedges[i].Origin]<<")\n";
    DebugLog<<"Prev: "<<Halfedges[i].Prev<<"\n";
    DebugLog<<"Next: "<<Halfedges[i].Next<<"\n";
    DebugLog<<"Twin: "<<Halfedges[i].Twin<<"\n";
    DebugLog<<"Face: "<<Halfedges[i].AdjFace<<"\n";
  }

  if (!CheckMesh(false, false, false))
    return false;
  
  vector<bool> transClaimed(NumNewVertices);
  for (int i=0;i<NumNewVertices;i++)
    transClaimed[i]=false;
  //unifying all vertices into the TransVertices
  vector<Vertex> NewVertices(NumNewVertices);
  for (int i=0;i<Vertices.size();i++){  //redundant, but not terrible
    if (!Vertices[i].Valid)
      continue;
    Vertex NewVertex=Vertices[i];
    NewVertex.ID=TransVertices[i];
    transClaimed[TransVertices[i]]=true;
    NewVertices[TransVertices[i]]=NewVertex;
  }
  
  for (int i=0;i<NumNewVertices;i++)
    if (!transClaimed[i]){
      DebugLog<<"TransVertex "<<i<<" not claimed!"<<endl;
      DebugLog<<"Group of vertices: ";
      for (int j=0;j<TransVertices.size();j++)
        if (TransVertices[j]==i)
          DebugLog<<j<<", ";
      DebugLog<<endl;
      NewVertices[i].Valid=false;  //this vertex is dead to begin with
    }
  
  Vertices=NewVertices;

  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
    Vertices[Halfedges[i].Origin].AdjHalfedge=i;
  }

  if (!CheckMesh(true, false, false))
    return false;
  
  //twinning up edges
  set<TwinFinder> Twinning;
  for (int i=0;i<Halfedges.size();i++){
    if ((Halfedges[i].Twin>=0)||(!Halfedges[i].Valid))
      continue;
    
    auto Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
    if (Twinit!=Twinning.end()){
      DebugLog<<"Twinning halfedge "<<i<<" and halfedge "<<Twinit->index<<"\n";
      if (Halfedges[Twinit->index].Twin!=-1)
        DebugLog<<"warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<endl;
      if (Halfedges[i].Twin!=-1)
        DebugLog<<"warning: halfedge "<<i<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<endl;
      Halfedges[Twinit->index].Twin=i;
      Halfedges[i].Twin=Twinit->index;
      
      if (Halfedges[i].isHex){
        DebugLog<<"Halfedge "<<i<<" is hex, infecting the other\n";
        Halfedges[Twinit->index].isHex = true;
      } else if (Halfedges[Twinit->index].isHex){
        DebugLog<<"Halfedge "<<Twinit->index<<" is hex, infecting the other\n";
        Halfedges[i].isHex = true;
      }
      Twinning.erase(*Twinit);
    } else {
      Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
    }
  }
  
  //check if there are any non-twinned edge which shouldn't be in a closed mesh
  for (int i=0;i<Halfedges.size();i++){
    if (Halfedges[i].Twin==-1)
      DebugLog<<"Halfedge "<<i<<" does not have a twin!"<<endl;
  }
  

  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    if (Halfedges[i].isHex)
      DebugLog<<"Hex edge "<<i<<"\n";
    else
      DebugLog<<"Triangle edge "<<i<<"\n";
    
    DebugLog<<"Origin: "<<Halfedges[i].Origin<<"\n";
    DebugLog<<"Prev: "<<Halfedges[i].Prev<<"\n";
    DebugLog<<"Next: "<<Halfedges[i].Next<<"\n";
    DebugLog<<"Twin: "<<Halfedges[i].Twin<<"\n";
    DebugLog<<"Face: "<<Halfedges[i].AdjFace<<"\n";
  }
  
  if (!CheckMesh(true, true, true))
    return false;
  
  //removing triangle components
  
  //starting with pure triangle vertices
  std::vector<bool> isPureTriangle(Vertices.size());
  std::vector<bool> isBoundary(Vertices.size());
  for (int i=0;i<Vertices.size();i++){
    isPureTriangle[i]=true;
    isBoundary[i]=false;
  }
  for (auto & Halfedge : Halfedges){
    if ((Halfedge.isHex)&&(Halfedge.Valid)){
      isPureTriangle[Halfedge.Origin]=isPureTriangle[Halfedges[Halfedge.Next].Origin]=false;  //adjacent to at least one hex edge
    }
    if (Halfedge.Twin==-1){
      isBoundary[Halfedge.Origin]=true;
      isPureTriangle[Halfedge.Origin]=false;  //this shouldn't be removed
    }
  }
  
  std::vector<bool> isEar(Vertices.size());
  for (int i=0;i<Vertices.size();i++){
    isEar[i] = (Halfedges[Vertices[i].AdjHalfedge].Twin==-1)&&(Halfedges[Halfedges[Vertices[i].AdjHalfedge].Prev].Twin==-1);
    if (isEar[i]) isPureTriangle[i]=false;
  }
  
  //realigning halfedges in hex vertices to only follow other hex edges
  DebugLog<<"Realigning edges around vertices"<<endl;
  for (int i=0;i<Vertices.size();i++){
    if ((isPureTriangle[i])||(!Vertices[i].Valid))
      continue;
    
    vector<int> hexHEorder;
    int hebegin = Vertices[i].AdjHalfedge;
    if (isBoundary[i]){
      //finding the first hex halfedge
      while (Halfedges[Halfedges[hebegin].Prev].Twin!=-1)
        hebegin =Halfedges[Halfedges[hebegin].Prev].Twin;
    }
    
    int heiterate=hebegin;
    do{
      if ((Halfedges[heiterate].isHex)||(Halfedges[heiterate].Twin==-1))
        hexHEorder.push_back(heiterate);
      if (Halfedges[heiterate].Twin==-1)
        break;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    
    
    for (int j=0;j<hexHEorder.size();j++){
      if ((isBoundary[i])&&(j==hexHEorder.size()-1))
        continue;
      Halfedges[hexHEorder[(j+1)%hexHEorder.size()]].Prev =Halfedges[hexHEorder[j]].Twin;
      Halfedges[Halfedges[hexHEorder[j]].Twin].Next =hexHEorder[(j+1)%hexHEorder.size()];
      Vertices[Halfedges[hexHEorder[j]].Origin].AdjHalfedge=hexHEorder[j];
    }
    
    if (isBoundary[i]){ //connect first to the prev
      Halfedges[hexHEorder[0]].Prev = Halfedges[hebegin].Prev;
      Halfedges[Halfedges[hebegin].Prev].Next =hexHEorder[0];
      Vertices[Halfedges[hexHEorder[0]].Origin].AdjHalfedge=hexHEorder[0];
    }
  }

  //invalidating all triangle vertices and edges
  DebugLog<<"Invalidating triangle vertices and edges"<<endl;
  for (int i=0;i<Vertices.size();i++)
    if (isPureTriangle[i])
      Vertices[i].Valid=false;
  
  for (auto & Halfedge : Halfedges)
    if ((!Halfedge.isHex)&&(Halfedge.Twin!=-1))
      Halfedge.Valid=false;
  
  //realigning faces
  DebugLog<<"Realigning faces"<<endl;
  Eigen::VectorXi visitedHE=Eigen::VectorXi::Zero(Halfedges.size());
  Eigen::VectorXi usedFace=Eigen::VectorXi::Zero(Faces.size());
  for (int i=0;i<Halfedges.size();i++){
    if ((!Halfedges[i].Valid)||(visitedHE[i]!=0))
      continue;
    
    //following the loop and reassigning face
    int currFace=Halfedges[i].AdjFace;
    Faces[currFace].AdjHalfedge=i;
    usedFace[currFace]=1;
    int hebegin=i;
    int heiterate=hebegin;
    int infinityCounter=0;
    do{
      infinityCounter++;
      if (infinityCounter>Halfedges.size()){
        DebugLog<<"Infinity loop in realigning faces on halfedge "<<i<<endl;
        return false;
      }
      Halfedges[heiterate].AdjFace=currFace;
      heiterate=Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
  }
  
  int countThree=0;
  DebugLog<<"Invalidating remainder faces"<<endl;
  for (int i=0;i<Faces.size();i++)
    if (!usedFace[i])
      Faces[i].Valid=false;
  
  
  //killing perfect ear faces (not doing corners atm)
  //counting valences
  vector<int> Valences(Vertices.size());
  for (int i=0;i<Vertices.size();i++)
    Valences[i]=0;
  
  for (auto & Halfedge : Halfedges){
    if (Halfedge.Valid){
      Valences[Halfedge.Origin]++;
      if (Halfedge.Twin<0)  //should account for the target as well
        Valences[Halfedges[Halfedge.Next].Origin]++;
    }
  }
     
  DebugLog<<"Invalidating ear (latent valence 2) faces"<<endl;
  for (int i=0;i<Faces.size();i++){
    if (!Faces[i].Valid)
      continue;
    countThree=0;
    int hebegin = Faces[i].AdjHalfedge;
    int heiterate=hebegin;
    do{
      if (Valences[Halfedges[heiterate].Origin]>2)
        countThree++;
      heiterate=Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
    if (countThree<3){
      DebugLog<<"Invalidating ear Face "<<i<<endl;
      DebugLog<<"Its vertices are "<<endl;
      do{
        Halfedges[heiterate].Valid=false;
        if (Halfedges[heiterate].Twin!=-1)
          Halfedges[Halfedges[heiterate].Twin].Twin=-1;
        if ((Halfedges[heiterate].Twin==-1)&&(Halfedges[Halfedges[heiterate].Prev].Twin==-1))  //origin is a boundary vertex
          Vertices[Halfedges[heiterate].Origin].Valid=false;
        
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      Faces[i].Valid=false;
    }
  }
  
  //need to realign all vertices pointing
  for (int i=0;i<Halfedges.size();i++)
    if (Halfedges[i].Valid)
      Vertices[Halfedges[i].Origin].AdjHalfedge=i;
    
  if (!CheckMesh(true, true, true))
    return false;

  
  for (int i=0;i<Valences.size();i++)
    if ((Vertices[i].Valid)&&(Valences[i]<2))
      Vertices[i].Valid=false;
  
  DebugLog<<"Starting unifying edges" <<endl;
  for (int i=0;i<Vertices.size();i++){
    
    if ((Vertices[i].Valid)&&(Valences[i]<=2)&&(!isEar[i]))
      UnifyEdges(Vertices[i].AdjHalfedge);
  }

  if (!CheckMesh(true, true, true))
    return false;
  
  //remove non-valid components
  CleanMesh();
  
  //checking if mesh is valid
  if (!CheckMesh(true, true, true))
    return false;

  //completing angle diffs - currently not working with a boundary
  for (auto & Vertice : Vertices){
    int hebegin = Vertice.AdjHalfedge;
    //weeding out boundaries - WHAT TO DO WITH THEM?
    int heiterate = hebegin;
    bool isBoundary=false;
    do{
      if (Halfedges[heiterate].Twin==-1){
        isBoundary=true;
        break;
      }
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    
    if (isBoundary)
      continue;
    
    heiterate = hebegin;
    int sumPrescribedDiff=0;
    int missingPrescribed=0;
    do{
      if (Halfedges[heiterate].prescribedAngleDiff!=-1)
        sumPrescribedDiff+=Halfedges[heiterate].prescribedAngleDiff;
      else
        missingPrescribed++;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);

    double completionAngle = 2*igl::PI*(1-sumPrescribedDiff/N);
    do{
      if (Halfedges[heiterate].prescribedAngleDiff!=-1)
        Halfedges[heiterate].prescribedAngle=2*igl::PI*(double)Halfedges[heiterate].prescribedAngleDiff/N;
      else
        Halfedges[heiterate].prescribedAngle = completionAngle;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);

  }
  
  return true;
}


