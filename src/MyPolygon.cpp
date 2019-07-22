#include "MyPolygon.h"

#include <math.h>

static const double PI = 3.14159265358979323846;

double MyPolygonUtility::Direction(const MyPoint& gridPoint, const MyPoint& FirePerimeterPoint)
{
    // calculates sweep direction for angle determination
    double zangle = 999.9, xDiff, ydiff;

    xDiff = FirePerimeterPoint.x - gridPoint.x;
    ydiff = FirePerimeterPoint.y - gridPoint.y;

    if (fabs(xDiff) < 1e-9)
    {
        xDiff = 0.0;
    }
    if (fabs(ydiff) < 1e-9)
    {
        ydiff = 0.0;
    }
    if (xDiff != 0.0)
    {
        zangle = atan(ydiff / xDiff);
        if (xDiff > 0.0)
        {
            zangle = (PI / 2.0) - zangle;
        }
        else
        {
            zangle = (3.0*PI / 2.0) - zangle;
        }
    }
    else
    {
        if (ydiff >= 0.0)
        {
            zangle = 0;
        }
        else if (ydiff < 0.0)
        {
            zangle = PI;
        }
    }

    return zangle;
}

bool MyPolygonUtility::IsOverlapping(const MyPoint& gridPoint, vector<MyPoint>& CurrentFirePerimeter)
{
    // determines if point is inside or outside a fire polygon (CurrentFirePerimeter)
    long NumVertex = CurrentFirePerimeter.size();
    long count = 0;
    long count1 = 0;
    long count2 = 0;
    bool inside = false;
    double angleA = 0.0, angleB;
    double cumulativeAngle = 0.0, angleDifference;

    MyPoint firePerimeterPoint;
    while (count < NumVertex) // make sure that startx, starty != x[0]y[0]
    {
        firePerimeterPoint.x = CurrentFirePerimeter[count].x;
        firePerimeterPoint.y = CurrentFirePerimeter[count].y;
        angleA = Direction(gridPoint, firePerimeterPoint);
        count++;
        if (angleA != 999.9)
        {
            break;
        }
    }

    for (count1 = count; count1 <= NumVertex; count1++)
    {
        if (count1 == NumVertex)
        {
            count2 = count - 1;
        }
        else
        {
            count2 = count1;
        }

        firePerimeterPoint.x = CurrentFirePerimeter[count2].x;
        firePerimeterPoint.y = CurrentFirePerimeter[count2].y;

        angleB = Direction(gridPoint, firePerimeterPoint);
        if (angleB != 999.9)
        {
            angleDifference = angleB - angleA;
            if (angleDifference > PI)
            {
                angleDifference = -(2.0*PI - angleDifference);
            }
            else if (angleDifference < -PI)
            {
                angleDifference = (2.0*PI + angleDifference);
            }
            cumulativeAngle -= angleDifference;
            angleA = angleB;
        }
    }
    if (fabs(cumulativeAngle) >= PI) // if absolute value of cumulative angle is >= PI
    {
        inside = true; // then point is inside polygon
    }

    return inside;
}
