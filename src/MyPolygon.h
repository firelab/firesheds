#pragma once

#include <vector>

using std::vector;

/***************************************************************************
// Data Structures
***************************************************************************/

//data structure for points
struct MyPoint
{
    double x;
    double y;
};

//data structure for a polygon part of a multiPolygon
struct MyPolygonPart
{
    vector<MyPoint> pointString;
};

//data structure for a multiPolygon
struct MyMultiPolygon
{
    vector<MyPolygonPart> polygonPart;
};

//data structure to hold bounding box
struct SBoundingBox
{
    double maxX;
    double maxY;
    double minX;
    double minY;
};

class MyPolygonUtility
{
public:
	static bool IsOverlapping(const MyPoint& gridPoint, vector<MyPoint>& CurrentFirePerimeter);
private:
    static double Direction(const MyPoint& gridPoint, const MyPoint& FirePerimeterPoint);
};
