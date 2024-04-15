#include <iostream>
#include "Vector.hpp"
#include "Utils.cpp"


int main()
{
    {
        Point<double> pt = Point<double>(1, 2, 3);
        std::cout << pt.X() << "\n";
        pt.X()++;
        std::cout << pt.X() << "\n";
        Point<double> pt2 = pt;
        std::cout << pt2.Z();
        Vector<double> vec = Vector<double>(1, 2, 3);
        Vector<double> vec2 = Vector<double>(5, 5, 5);
        normalize(vec);
        std::cout << vec.X() << std::endl;
        vec = vec2;
        vec += vec2;
        Line<double> lin = Line<double>(vec, pt2);
        Line<double> lin2 = Line<double>(pt, pt2);
        std::cout << vec.X() << std::endl;
        vec *= 1000.00;
        std::cout << vec.Z();
    }
    {
        Point<double> p1 = Point<double>(5, 10, 0);
        Point<double> p2 = Point<double>(5, 1, 0);
        Point<double> p3 = Point<double>(1, 5, 0);
        Point<double> p4 = Point<double>(10, 5, 0);
        Segment<double> seg1 = Segment<double>(p2, p1);
        Segment<double> seg2 = Segment<double>(p4, p3);
        std::cout << isPointsCreatePlane(seg1.getPoints().first, seg1.getPoints().second, seg2.getPoints().first, seg2.getPoints().second) << std::endl;
        std::cout << "Segment crossed " << isSegmentCrossed(seg1, seg2) << std::endl;
        auto res = intersectPointOfSegments(seg1, seg2);
        if (res.has_value()) {
            auto pt = res.value();
            std::cout << pt.X() << " " << pt.Y() << " " << pt.Z() << std::endl;
        }
        else {
            std::cout << "nullOpt\n";
        }
    }
    {
        Line<double> line = Line<double>(Point<double>(-1.55, -1.38, 0), Point<double>(-2, 1, 1.21));
        Plane<double> plane = Plane<double>(Point<double>(0.33, 0.49, 0), Point<double>(0.88, 1.69, 0), Point<double>(-2.33, 2.5, 1));
        std::cout << plane.getD() << std::endl;
        std::cout << getAnglePlaneLine(plane, line);
    }
    {
        Plane<double> plane = Plane<double>(Vector<double>(2, 3, 4), 8.0);
        Point<double> pt1 = Point<double>(2, 2, 3);
        Point<double> pt2 = Point<double>(4, 3, 6);
        Line<double> line = Line<double>(pt1, pt2);

        auto pt = PointOfLineandPlaneIntersect(line, plane);
        if (pt.has_value()) {
            auto res = pt.value();
            std::cout << res.X() << " " << res.Y() << " " << res.Z() << std::endl;
        }
        else {
            std::cout << "No intersect" << std::endl;
        }
    }
    return 0;
}
