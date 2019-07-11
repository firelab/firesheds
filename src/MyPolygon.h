#pragma once

#include <vector>

using std::vector;

/***************************************************************************
// Data Structures
***************************************************************************/

//data structure for points
typedef struct MyPoint2D
{
    double X;
    double Y;
}MyPoint2D;

//data structure for a multipoint feature
typedef struct MultipointFeature
{
    vector<MyPoint2D>PointsOfFeature;
}MultipointFeature;

//data structure for lines
typedef struct MyLine2D
{
    vector<MyPoint2D> LineString;
}MyLine2D;

//data structure for a line feature
typedef struct LineFeature
{
    vector<MyLine2D>LinesOfFeature;
}LineFeature;

//data structure for rings
typedef struct MyRing2D
{
    vector<MyPoint2D> RingString;
    bool IsClockwised;
}MyRing2D;

//data structure for polygons
typedef struct MyPolygon2D
{
    vector<MyRing2D>Polygon;
}MyPolygon2D;

//data structure for a polygon feature
typedef struct PolygonFeature
{
    vector<MyPolygon2D>PolygonsOfFeature;
}PolygonFeature;

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
	static bool IsOverlapping(const MyPoint2D& gridPoint, vector<MyPoint2D>& CurrentFirePerimeter);
private:
    static double Direction(const MyPoint2D& gridPoint, const MyPoint2D& FirePerimeterPoint);
};
