#pragma once
/**
  * @brief Point Class. X,Y,Z=0 - coordinates of point
  */
template<typename T>
class Point {
private:
    /**
      * @brief X,Y,Z coordinates
      */
    T _x = 0;
    T _y = 0;
    T _z = 0;
public:
    Point() {}
    /**
      * @brief X,Y,Z to create a point
      */
    Point(T x, T y, T z) : _x(x), _y(y), _z(z) {}
    /**
      * @brief X,Y - to create a point. Z = 0
      */
    Point(T x, T y) : _x(x), _y(y) {}
    /**
      * @brief return Point.X by reference(can be edited)
      */
    T& X() {
        return _x;
    }
    /**
      * @brief return Point.Y by reference(can be edited)
      */
    T& Y() {
        return _y;
    }
    /**
      * @brief return Point.Z by reference(can be edited)
      */
    T& Z() {
        return _z;
    }
    /**
      * @brief return Point.X by value(can NOT be edited)
      */
    T X() const {
        return _x;
    }
    /**
      * @brief return Point.Y by value(can NOT be edited)
      */
    T Y() const {
        return _y;
    }
    /**
      * @brief return Point.Z by value(can NOT be edited)
      */
    T Z() const {
        return _z;
    }
    /**
      * @brief Copy constructor
      */
    Point(const Point& other) {
        this->_x = other._x;
        this->_y = other._y;
        this->_z = other._z;
    }
    /**
      * @brief Assignment operator
      */
    Point& operator=(const Point& other);
    /**
      * @brief -= operator for Point -= Point
      */
    Point& operator -=(const Point<T>& other);
    /**
      * @brief - operator: Point = Point - Point(using -=)
      */
    Point operator-(const Point<T>& other);
};
template<typename T>
Point<T>& Point<T>::operator=(const Point<T>& other) {
    if (this != &other) {
        this->_x = other._x;
        this->_y = other._y;
        this->_z = other._z;
    }
    return *this;
}
template <typename T>
Point<T>& Point<T>::operator-=(const Point<T>& other)
{
    this->_x -= other._x;
    this->_y -= other._y;
    this->_z -= other._z;
    return *this;
}
template <typename T>
Point<T> Point<T>::operator-(const Point<T>& other)
{
    Point<T> nPt = *this;
    return nPt -= other;
}