#pragma once 
#include "Line.hpp"
#include "Plane.hpp"
#include "Point.hpp"
#include "Segment.hpp"
#include "Vector.hpp"

#include <optional>
#include <math.h>
#include <vector>
/**
  * @brief Find if vector`s is colinear. Pass 2 vectors. x1/x2 = y1/y2/ = z1/z2 - statement for colinear
  */
template<typename T>
bool isVectorCollinear(const Vector<T>& vec1, const Vector<T>& vec2) {
    auto x = vec1.X() / vec2.X();
    auto y = vec1.Y() / vec2.Y();
    auto z = vec1.Z() / vec2.Z();
    if (fabs(x - y) < TOLERANCE && fabs(y - z) < TOLERANCE) {
        return true;
    }
    else {
        return false;
    }
}
/**
  * @brief Find a magnitude of a vector. Pass vector
  */
template<typename T>
double magnitude(const Vector<T>& vec) {
    double res = pow(vec.X(), 2) + pow(vec.Y(), 2) + pow(vec.Z(), 2);
    return sqrt(res);
}
/**
  * @brief Normalize a vector. Pass a Vector&. Function edited Vector and modify fields of Vector
  */
template<typename T>
void normalize(Vector<T>& vec) {
    double mag = magnitude(vec);
    vec.X() /= mag;
    vec.Y() /= mag;
    vec.Z() /= mag;
}
/**
  * @brief Find a dot product of 2 vector: Pass (Vector, Vector). return a double - result of a dot production
  */
template<typename T>
double dotProduct(const Vector<T>& vec1, const Vector<T>& vec2) {
    double res = vec1.X() * vec2.X() + vec1.Y() * vec2.Y() + vec1.Z() * vec2.Z();
    return res;
}
/**
  * @brief Find a dot product of 2 vector and a angle between them: use magnitude and cos(angle) to find an dot prduct. Pass an angle in radian
  */
template<typename T>
double dotProduct(const Vector<T>& vec1, const Vector<T>& vec2, const double angle) {
    double res = magnitude(vec1) * magnitude(vec2) * cos(angle);
    return res;
}
/**
  * @brief Find a croos product of a 2 vectors. Pass (Vector, Vector). Return an Normal Vector to a plane, in which lies 2 given Vector
  */
template<typename T>
Vector<T> crossProduct(const Vector<T>& vec1, const Vector<T>& vec2) {
    double x = vec1.Y() * vec2.Z() - vec2.Y() * vec1.Z();
    double y = -(vec1.X() * vec2.Z() - vec2.X() * vec1.Z());
    double z = vec1.X() * vec2.Y() - vec2.X() * vec1.Y();
    return Vector<T>(x, y, z);
}
/**
  * @brief Find a croos product of a 2 vectors. Pass (Vector, Vector, angle(in radian)). Return an double, which provides a possible if 2 vector is right side or left side, and we can deduce direction of a normal
  */
template<typename T>
double crossProduct(const Vector<T>& vec1, const Vector<T>& vec2, const double angle) {
    double res = magnitude(vec1) * magnitude(vec2) * sin(angle);
    return res;
}
/**
  * @brief Find an angle(in radian) betweeen 2 Vectors. Pass (Vector, Vector). Using dotProduct finds an angle: returns an angle in radian
  */
template<typename T>
double findAngle(const Vector<T>& v1, const Vector<T>& v2) {
    double dtPro = dotProduct(v1, v2);
    double dtProS = magnitude(v1) * magnitude(v2);
    double angle = acos(dtPro / dtProS);
    
    const double pi = acos(-1);

    if (angle >= pi) {
        return 2 * pi - angle;
    }
    else {
        return angle;
    }
}


/**
 * @brief enum show possible states for 2 segments: crossed, or no touch
 */
enum POSITION_OFSEGMENTS {
    CROSSED, NOTOUCH
};
/**
 * @brief Return if points created plane. Pass (Point, Point, Point, Point). Return bool true - if creates, false - not
 */
template<typename T>
bool isPointsCreatePlane(const Point<T>& p1, const Point<T>& p2, const Point<T>& p3, const Point<T>& p4) {
    Vector<T> v1 = Vector<T>(p1, p2);
    Vector<T> v2 = Vector<T>(p1, p3);

    Vector<T> n = crossProduct(v1, v2);

    if (fabs(n.X() * (p4.X() - p1.X()) + n.Y() * (p4.Y() - p1.Y()) + n.Z() * (p4.Z() - p1.Z())) <= (0 + TOLERANCE)) {
        return true;
    }
    return false;
}
/**
 * @brief Return if 2 segment crossed. Pass (Segment, Segment). Return bool true - if crossed, false - not.
 * @brief do it finding angles between vector created by points and finding a sum of it. IF sum equals, return true. Calculate it for all points
 */
template<typename T>
bool isSegmentCrossed(const Segment<T>& seg1, const Segment<T>& seg2) {
    Point<T> p1 = seg1.getPoints().first;
    Point<T> p2 = seg1.getPoints().second;

    Point<T> p3 = seg2.getPoints().first;
    Point<T> p4 = seg2.getPoints().second;

    if (!isPointsCreatePlane(p1, p2, p3, p4))
        return false;


    double angleP3P2P1 = findAngle(Vector<T>(p2, p3), Vector<T>(p2, p1));
    double angleP1P2P4 = findAngle(Vector<T>(p2, p1), Vector<T>(p2, p4));
    double angleP3P2P4 = findAngle(Vector<T>(p2, p3), Vector<T>(p2, p4));
    if (!(fabs(angleP3P2P4 - (angleP3P2P1 + angleP1P2P4)) < TOLERANCE))
        return false;

    double angleP3P4P1 = findAngle(Vector<T>(p4, p3), Vector<T>(p4, p1));
    double angleP3P4P2 = findAngle(Vector<T>(p4, p3), Vector<T>(p4, p2));
    double angleP1P4P2 = findAngle(Vector<T>(p4, p1), Vector<T>(p4, p2));
    if (!(fabs(angleP1P4P2 - (angleP3P4P1 + angleP3P4P2)) < TOLERANCE))
        return false;

    double angleP4P3P1 = findAngle(Vector<T>(p3, p4), Vector<T>(p3, p1));
    double angleP4P3P2 = findAngle(Vector<T>(p3, p4), Vector<T>(p3, p2));
    double angleP1P3P2 = findAngle(Vector<T>(p3, p1), Vector<T>(p3, p2));
    if (!(fabs(angleP1P3P2 - (angleP4P3P1 + angleP4P3P2)) < TOLERANCE))
        return false;

    double angleP3P1P2 = findAngle(Vector<T>(p1, p3), Vector<T>(p1, p2));
    double angleP2P1P4 = findAngle(Vector<T>(p1, p2), Vector<T>(p1, p4));
    double angleP3P1P4 = findAngle(Vector<T>(p1, p3), Vector<T>(p1, p4));
    if (!(fabs(angleP3P1P4 - (angleP3P1P2 + angleP2P1P4)) < TOLERANCE))
        return false;

    return true;
}
/**
 * @brief gaussElimination code to solve matrix. Generated by chat gpt. Pass (std::vectot). Return std::vector with solution
 */
template<typename T>
std::vector<T> gaussElimination(std::vector<std::vector<T>>& matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size() - 1; // Exclude the last column (the constants column)

    // Forward elimination
    for (int i = 0; i < rows; ++i) {
        // Find the pivot row
        int maxRowIndex = i;
        for (int k = i + 1; k < rows; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRowIndex][i])) {
                maxRowIndex = k;
            }
        }

        // Swap current row with the pivot row
        std::swap(matrix[i], matrix[maxRowIndex]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < rows; ++k) {
            double factor = -matrix[k][i] / matrix[i][i];
            for (int j = i; j <= cols; ++j) {
                matrix[k][j] += matrix[i][j] * factor;
            }
        }
    }

    // Back substitution
    std::vector<T> solution(cols, 0);
    for (ptrdiff_t i = rows - 1; i >= 0; --i) {
        solution[i] = matrix[i][cols];
        for (int j = i + 1; j < cols; ++j) {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }
    return solution;
}
/**
 * @brief Finds a point of intersection segments. Pass (Segment, Segment). Returns std::nullopt if point does not exist, or Point, if point exist
 */
template<typename T>
std::optional<Point<T>> intersectPointOfSegments(const Segment<T>& seg1, const Segment<T>& seg2) {
    if (!isSegmentCrossed(seg1, seg2))
        return std::nullopt;

    Point<T> pt = seg1.getPoints().first;
    Point<T> pt2 = seg1.getPoints().second;
    Vector<T> ln = Vector<T>(pt, pt2);

    Point<T> pt_ = seg2.getPoints().first;
    Point<T> pt_2 = seg2.getPoints().second;
    Vector<T> ln_ = Vector<T>(pt_, pt_2);

    double det = -ln.X() * ln_.Y() + ln.Y() * ln_.X();

    if (det == 0)
        return std::nullopt;
    double matrix[2][3];

    matrix[0][0] = ln.X();
    matrix[0][1] = -ln_.X();
    matrix[0][2] = pt_.X() - pt.X();

    matrix[1][0] = ln.Y();
    matrix[1][1] = -ln_.Y();
    matrix[1][2] = pt_.Y() - pt.Y();

    std::vector<std::vector<T>> vec = {
        {matrix[0][0], matrix[0][1], matrix[0][2] },
        {matrix[1][0], matrix[1][1], matrix[1][2]}
    };
    std::vector<T> res = gaussElimination(vec);

    double Xcheck = fabs(pt.X() + ln.X() * res[0] - (pt_.X() + ln_.X() * res[1]));

    double Ycheck = fabs(pt.Y() + ln.Y() * res[0] - (pt_.Y() + ln_.Y() * res[1]));

    double Zcheck = fabs(pt.Z() + ln.Z() * res[0] - (pt_.Z() + ln_.Z() * res[1]));

    if (Xcheck < TOLERANCE && Ycheck < TOLERANCE && Zcheck < TOLERANCE)
        return Point<T>((pt.X() + ln.X() * res[0]), (pt.Y() + ln.Y() * res[0]), (pt.Z() + ln.Z() * res[0]));
    else
        return std::nullopt;
}
/**
 * @brief Function which returns a state of a 2 segment - crossed, or no touch
 */
template<typename T>
POSITION_OFSEGMENTS getPositionOfSegments(const Segment<T>& seg1, const Segment<T>& seg2) {
    if (isSegmentCrossed(seg1, seg2))
        return CROSSED;
    return NOTOUCH;
}

/**
 * @brief returns a distance between 2 points. Pass (Point, Point). return double distance
 */
template<typename T>
double distBetweenPoints(const Point<T>& p1, const Point<T>& p2) {
    double dist = sqrt(pow(p1.X() - p2.X(), 2) + pow(p1.Y() - p2.Y(), 2) + pow(p1.Z() - p2.Z(), 2));
    return dist;
}
/**
 * @brief returns if a point lies on a line of a segment. Pass(Point, Segment) Returns a true, if point lies on a segment. False - does not
 */
template <typename T>
bool isPointOnSegment(const Point<T>& point, const Segment<T>& seg) {
    std::pair<Point<T>, Point<T>> pair = seg.getPoints();
    double dist_P1_Pgiven = distBetweenPoints(pair.first, point);
    double dist_P2_Pgiven = distBetweenPoints(pair.second, point);
    double dist_P2_P1 = distBetweenPoints(pair.first, pair.second);
    if (fabs(dist_P2_P1 - (dist_P1_Pgiven + dist_P2_Pgiven)) < TOLERANCE)
        return true;
    return false;
}
/**
 * @brief returns if a point lies on a line of a line. Pass(Point, Line) Returns a true, if point lies on a segment. False - does not
 */
template <typename T>
bool isPointOnLine(const Point<T>& pt, const Line<T>& line) {
    Point<T> ptL = line.getPoint();
    Vector<T> vec = line.getVector();

    double xDif = (pt.X() - ptL.X()) / vec.X();
    double yDif = (pt.Y() - ptL.Y()) / vec.Y();
    double zDif = (pt.Z() - ptL.Z()) / vec.Z();

    if (fabs(xDif - yDif) < TOLERANCE && fabs(yDif - zDif) < TOLERANCE && fabs(xDif - zDif) < TOLERANCE)
        return true;
    return false;
}
/**
 * @brief Returns true if planes are parallel: else false. Pass (Plane, Plane).
 */
template <typename T>
bool ifPlanesparallel(const Plane<T>& plane1, const Plane<T>& plane2) {
    if (findAngle(plane1.getVector(), plane2.getVector()) == 0)
        return true;
    return false;
}

/**
 * @brief returns a line - intersection of 2 planes. Pass (Plane, Plane). Returns an std::nullopt if line does not exist(planes are parallel). Returns a line if exist.
 */
template <typename T>
Line<T> getLineOfInterPlanes(const Plane<double>& plane1, const Plane<double>& plane2) {
    if (ifPlanesparallel(plane1, plane2))
        return std::nullopt;

    Vector<T> dir = crossProduct(plane1.getVector(), plane2.getVector());

    double a = (plane2.getD() * dotProduct(plane2.getVector(), plane1.getVector()) - plane1.getD() * pow(magnitude(plane2.getVector()), 2)) / (pow(dotProduct(plane1.getVector(), plane2.getVector()), 2) - pow(magnitude(plane2.getVector()), 2) * pow(magnitude(plane1.getVector()), 2));
    double b = (plane1.getD() * dotProduct(plane2.getVector(), plane1.getVector()) - plane2.getD() * pow(magnitude(plane1.getVector()), 2)) / (pow(dotProduct(plane1.getVector(), plane2.getVector()), 2) - pow(magnitude(plane2.getVector()), 2) * pow(magnitude(plane1.getVector()), 2));
    Point<T> resP = plane1.getVector() * a + plane2.getVector() * b;
    return Line<double>(dir, resP);
}
/**
 * @brief returns true if line is parallel to a plane. Pass (Line, Plane) returns true - if parallel. False - NOT parallel.
 */
template <typename T>
bool isLineParallelToPlane(const Line<T>& line, const Plane<double>& plane) {
    if (dotProduct(line.getVector(), plane.getVector()) == 0) {
        return true;
    }
    else {
        return false;
    }
}
/**
 * @brief Returns a Point of Intersection of line and plane. Pass(Line, Plane). Returns std::nullopt if line is parallel to a plane. Else return a Point of intersection.
 */
template <typename T>
std::optional<Point<T>> PointOfLineandPlaneIntersect(const Line<T>& line, const Plane<T>& plane) {
    if (isLineParallelToPlane(line, plane))
        return std::nullopt;
    double x = line.getPoint().X() + line.getVector().X() * ((-dotProduct(plane.getVector(), Vector<double>(line.getPoint())) + plane.getD()) / dotProduct(plane.getVector(), line.getVector()));
    double y = line.getPoint().Y() + line.getVector().Y() * ((-dotProduct(plane.getVector(), Vector<double>(line.getPoint())) + plane.getD()) / dotProduct(plane.getVector(), line.getVector()));
    double z = line.getPoint().Z() + line.getVector().Z() * ((-dotProduct(plane.getVector(), Vector<double>(line.getPoint())) + plane.getD()) / dotProduct(plane.getVector(), line.getVector()));
    return Point<double>(x, y, z);
}
/**
 * @brief Finds a volume of tetrahedron(призми), which is built on this 3 vectors. Pass(Vector, Vector, Vector). REturns double - volume.
 */
template <typename T>
double volumeOfTetrahedron(const Vector<T>& vec1, const Vector<T>& vec2, const Vector<T>& vec3) {
    double volume = 1 / 6 * dotProduct(vec1, crossProduct(vec2, vec3));
    return volume;
}
/**
 * @brief Check if 3 vectors are complanar. IF 3 vectors lies in same plane - return true; else false. Uses a volumeOfTetrahedron to find a volume. Pass(Vector, Vector, Vector)
 */
template <typename T>
bool isComplanarVector(const Vector<T>& vec1, const Vector<T>& vec2, const Vector<T>& vec3) {
    if (volumeOfTetrahedron(vec1, vec2, vec3) < TOLERANCE) {
        return true;
    }
    else {
        return false;
    }
}
/**
 * @brief Check if 4 points are complanar. Creates a 3 vectors. IF 3 vectors lies in same plane - return true; else false. Uses a volumeOfTetrahedron to find a volume. Pass(Point, Point, Point, Point)
 */
template <typename T>
bool isComplanarPoints(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3, const Point<T>& pt4) {
    if (volumeOfTetrahedron(Vector<double>(pt1, pt2), Vector<double>(pt1, pt3), Vector<double>(pt1, pt4)) < TOLERANCE) {
        return true;
    }
    else {
        return false;
    }
}


/**
 * @brief Enum show a possible state of 2 lines - parallel, intersected, Passing.
 */
enum POSITION_OF_LINES {
    PARALLEL, INTERSECTED, PASSING
};
/**
 * @brief Finds if 2 lines line in one plane. Pass(Line, Line). Returns true - if lying; false - if not lying
 */
template <typename T>
bool isLinesOnOnePLane(const Line<T>& line1, const Line<T>& line2) {
    Point<T> p1 = line1.getPoint();
    Point<T> p2 = line2.getPoint();

    Vector<T> vec1 = Vector<T>(p2.X() - p1.X(), p2.Y() - p1.Y(), p2.Z() - p1.Z());
    Vector<T> vec2 = line1.getVector();
    Vector<T> vec3 = line2.getVector();

    double res = vec1.X() * vec2.Y() * vec3.Z() + vec1.Y() * vec2.Z() * vec3.X() + vec1.z() * vec2.X() * vec3.Y();
    res -= vec3.X() * vec2.Y() * vec1.Z() + vec1.X() * vec2.Z() * vec3.Y() + vec1.Y() * vec2.X() * vec3.Z();
    if (res == 0)
        return true;
    return false;
}
/**
 * @brief Checks if 2 lines are parallel. Pass(Line, Line). Returns true if parallel. False - not
 */
template<typename T>
bool isLinesParallel(const Line<T>& line1, const Line<T>& line2) {
    Vector<T> vec1 = line1.getVector();
    Vector<T> vec2 = line2.getVector();

    if (fabs(vec1.X() / vec2.X() - vec1.Y() / vec2.Y()) < TOLERANCE && fabs(vec1.Y() / vec2.Y() - vec1.Z() / vec2.Z()) < TOLERANCE)
        return true;
    return false;
}
/**
 * @brief REturns a Point of intersection of line. IF lines are in same plane and not parallel - returns Point. Else - return std::nullopt. Pass(Line, Line)
 */
template <typename T>
std::optional<Point<T>> intersectPointOfLine(const Line<T>& line1, const Line<T>& line2) {
    if (isLinesParallel(line1, line2))
        return std::nullopt;
    if (!isLinesOnOnePLane(line1, line2))
        return std::nullopt;



    Point<T> pt = line1.getPoint();
    Vector<T> ln = line1.getPoint();

    Point<T> pt_ = line2.getPoint();
    Vector<T> ln_ = line2.getVector();

    double det = -ln.X() * ln_.Y() + ln.Y() * ln_.X();

    if (det == 0)
        return std::nullopt;
    double matrix[2][3];

    matrix[0][0] = ln.X();
    matrix[0][1] = -ln_.X();
    matrix[0][2] = pt_.X() - pt.X();

    matrix[1][0] = ln.Y();
    matrix[1][1] = -ln_.Y();
    matrix[1][2] = pt_.Y() - pt.Y();

    std::vector<std::vector<T>> vec = {
        {matrix[0][0], matrix[0][1], matrix[0][2] },
        {matrix[1][0], matrix[1][1], matrix[1][2]}
    };
    std::vector<T> res = gaussElimination(vec);

    double Xcheck = fabs(pt.X() + ln.X() * res[0] - (pt_.X() + ln_.X() * res[1]));

    double Ycheck = fabs(pt.Y() + ln.Y() * res[0] - (pt_.Y() + ln_.Y() * res[1]));

    double Zcheck = fabs(pt.Z() + ln.Z() * res[0] - (pt_.Z() + ln_.Z() * res[1]));

    if (Xcheck < TOLERANCE && Ycheck < TOLERANCE && Zcheck < TOLERANCE)
        return Point<T>((pt.X() + ln.X() * res[0]), (pt.Y() + ln.Y() * res[0]), (pt.Z() + ln.Z() * res[0]));
    else
        return std::nullopt;
}
/**
 * @brief Function to get the state of 2 line: parallel, intersected, or passing
 */
template<typename T>
POSITION_OF_LINES getPositionOfLines(const Line<T>& line1, const Line<T>& line2) {
    if (isLinesParallel(line1, line2))
        return PARALLEL;
    auto res = intersectPointOfLine(line1, line2);
    if (res.has_value())
        return INTERSECTED;
    return PASSING;
}

/**
 * @brief Finds a distance between Point and plane. Returns a double - distance. Pass (Point, Plane)
 */
template <typename T>
double distBetweenPointAndPlane(const Point<T>& pt, const Plane<T>& plane) {
    double dist = abs(plane.getVector().X() * pt.X() + plane.getVector().Y() * pt.Y() + plane.getVector().Z() * pt.Z() + plane.getD()) / sqrt(pow(plane.getVector().X(), 2) + pow(plane.getVector().Y(), 2) + pow(plane.getVector().Z(), 2));
    return dist;
}
/**
 * @brief Finds a distance between Point and line. Returns a double - distance. Pass (Point, Line)
 */
template <typename T>
double distBetweenPointAndLine(const Point<T>& pt, const Line<T>& line) {
    auto AC = pt - line.getPoint();
    auto t = dotProduct(line.getVector(), AC);

    auto xt = line.getPoint() + line.getVector() * t;

    auto dist_vec = xt - pt;
    return dist_vec.magnitude();
}

/**
 * @brief Returns an angle between plane and line. Pass(line, plane)
 */
template<typename T>
double getAnglePlaneLine(const Plane<T>& plane, const Line<T>& line) {
    double angle = findAngle(plane.getVector(), line.getVector());
    const double pi = acos(-1);
    if (angle >= (pi / 2)) {
        return angle - pi / 2;
    }
    else {
        return angle;
    }
}
/**
 * @brief Returns an angle between plane and plane. Pass(plane, plane)
 */
template <typename T>
double getAnglePlanes(const Plane<double>& plane1, const Plane<double>& plane2) {
    double angle = findAngle(plane1.getVector(), plane2.getVector());
    const double pi = acos(-1);
    if (angle >= (pi / 2)) {
        return angle - pi / 2;
    }
    else {
        return angle;
    }
}