#pragma once
#include  "Vector.hpp"
#include "Point.hpp"
/**
  * @brief Line class
  * @brief (x - x1) / a = (y - y1)/ b = (z - z1) /c = t - parametric formula of a line in 3D
  * @brief (a,b,c) - directional vector | напрямний вектор
  * @brief (x1,x2,x3) - point of a line
  */
template<typename T>
class Line {
private:
    /**
      * @brief directional vector of a line
      */
    Vector<T> _vector;
    /**
      * @brief Point
      */
    Point<T> _point;
public:
    /**
     @brief Point p1, Point p2. Use p1 as point for line. Vector calculated as p2 - p1;
     @return Created Line
    */
    Line(const Point<T>& p1, const Point<T>& p2);

    /**
     @brief Vector v, Point p. Use p as point for line. Vector v used as directional;
     @return Created Line
    */
    Line(const Vector<T>& dir, const Point<T>& p) : _vector(dir), _point(p) {};

    /**
     @brief Copy constructor
    */
    Line(const Line<T>& other);
    Line& operator=(const Line<T>& other);

    /**
     @return return point for a line
    */
    Point<T> getPoint() const;
    /**
     @return return vector of a line
    */
    Vector<T> getVector() const;
};
template<typename T>
Point<T> Line<T>::getPoint() const {
    return this->_point;
}
template<typename T>
Vector<T> Line<T>::getVector() const {
    return this->_vector;
}
template<typename T>
Line<T>::Line(const Point<T>& p1, const Point<T>& p2) {
    Vector<T> newVec = Vector<T>(p1, p2);
    _vector = newVec;
    _point = p1;
}
template <typename T>
Line<T>::Line(const Line<T>& other)
{
    this->_vector = other._vector;
    this->_point = other._point;
}
template <typename T>
Line<T>& Line<T>::operator=(const Line<T>& other)
{
    if (this != &other) {
        this->_vector = other._vector;
        this->_point = other._point;
    }
    return *this;
}