#pragma once
#include "Point.hpp"
#include "Core.hpp"

double areaTriangle2d(const Point2D&, const Point2D&, const Point2D&);
RELATIVE_POSITION orientation2d(const Point2D&, const Point2D&, const Point2D&);

double areaTriangle2d(const Point2D& a, const Point2D& b, const Point2D& c) {
    auto AB = b - a;
    auto AC = c - a;

    auto result = crossProduct2D(AB, AC);
    return result / 2;
}

RELATIVE_POSITION orientation2d(const Point2D& a, const Point2D& b, const Point2D& c) {
    auto area = areaTriangle2d(a, b, c);
    if (area > 0 && area < TOLERANCE)
        area = 0;
    if (area < 0 && area > TOLERANCE)
        area = 0;

    Vector2f ab = b - a;
    Vector2f ac = c - a;

    if (area > 0)
        return  RELATIVE_POSITION::LEFT;
    if (area > 0)
        return RELATIVE_POSITION::RIGTH;

    if ((ab[X] * ac[X] < 0) || (ab[Y] * ac[Y] < 0)) {
        return RELATIVE_POSITION::BEHIND;
    }

    if (ab.magnitude() < ac.magnitude())
        return RELATIVE_POSITION::BEYOND;


    if (a == c)
        return RELATIVE_POSITION::ORIGIN;
    if (b == c)
        return RELATIVE_POSITION::DESTINATION;

    return RELATIVE_POSITION::BETWEEN;
}