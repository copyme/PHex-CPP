#ifndef DEFINITIONS_HEADER_FILE
#define DEFINITIONS_HEADER_FILE

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Cartesian.h> 
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <list>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpq>  EKernel;
typedef CGAL::Gmpz  EInt;

typedef EKernel::RT ENumber;
typedef Kernel::Point_3 Point3D;
typedef EKernel::Line_2 ELine2D;
typedef EKernel::Segment_2 ESegment2D;
typedef EKernel::Point_2 EPoint2D;
typedef EKernel::Vector_2 EVector2D;
typedef EKernel::Point_3 EPoint3D;
typedef EKernel::Triangle_2 ETriangle2D;
typedef EKernel::Direction_2 EDirection2D;

#endif
