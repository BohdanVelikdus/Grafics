#pragma once
#include "Vector.hpp"
#include "Point.hpp"
/**
  * @brief Plane Class
  * @brief AX+BY+CZ = D, (A,B,C) - normal vector
  * @brief A(x-x1) + B(y - y1) + C(z - z1) = 0
  * @brief D is calculated a A*x1 + B*y1 + C*z1
  */
template<typename T>
class Plane {
private:
    /**
      * @brief Noraml vector of a plane
      */
    Vector<T> _NVec;
    /**
      * @brief D argument of a plane
      */
    T _point;
public:
    /**
      * @brief Construct of a plane using Normal Vector and Point
      */
    Plane(const Vector<T>& vec, const Point<T>& pt);
    /**
      * @brief Construct of a plane using Normal Vector and D
      */
    Plane(const Vector<T>& vec, T pt) : _NVec(vec), _point(pt) {}
    /**
      * @brief Construct a plane using 2 Vector`s in plane. Find an crossProduct of a vectors. D find as dotproduct of normal vector and given vector
      */
    Plane(const Vector<T>& vec1, const Vector<T>& vec2);
    /**
      * @brief Construct a plane using 3 points. Find a 2 vector to use it in constructor with 2 vectors
      */
    Plane(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3);
    /**
      * @brief Copy constructor for a plane
      */
    Plane(const Plane<T>& other);
    /**
      * @brief get normal vector of a plane (const)
      */
    Vector<T> getVector() const;
    /**
      * @brief get D (const)
      */
    T getD() const;
};
template <typename T>
Plane<T>::Plane(const Vector<T>& vec, const Point<T>& pt) {
    _NVec = vec;
    _point = vec.X() * pt.X() + vec.Y() * pt.Y() + vec.Z() * pt.Z();
}
template<typename T>
Vector<T> Plane<T>::getVector() const {
    return _NVec;
}
template<typename T>
T Plane<T>::getD() const {
    return _point;
}
template <typename T>
Plane<T>::Plane(const Vector<T>& vec1, const Vector<T>& vec2)
{
    _NVec = crossProduct(vec1, vec2);
    _point = dotProduct(_NVec, vec1);
}
template <typename T>
Plane<T>::Plane(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3)
{
    Vector<T> vec1 = Vector<T>(pt1, pt2);
    Vector<T> vec2 = Vector<T>(pt3, pt2);
    _NVec = crossProduct(vec1, vec2);
    _point = dotProduct(_NVec, vec1);
}