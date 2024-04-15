#pragma once
#include <math.h>
#include "Point.hpp"
#include <type_traits>
/**
  * @brief Tolerance for measuring floats
  */
#define TOLERANCE 0.00001
  /**
   * @brief Vector class
   * @brief (x,y,z) - coords of a vector
   */
template<typename T>
class Vector {
private:
    /**
      * @brief X,Y,Z coords
      */
    T _x = 0;
    T _y = 0;
    T _z = 0;
public:
    Vector() {}
    /**
      * @brief Point start, Point end. Vector coords = end - start
      * @return create a Vector
      */
    Vector(const Point<T>& start, const Point<T>& end);
    /**
      * @brief X, Y, Z as coordinates for Vector
      * @return create a Vector
      */
    Vector(T x, T y, T z) : _x(x), _y(y), _z(z) {}
    /**
      * @brief X, Y, Z=0 as coordinates for Vector
      * @return create a Vector
      */
    Vector(T x, T y) : _x(x), _y(y) {}
    /**
      * @brief Copy constructor
      */
    Vector(const Point<T>& other);

    /**
      * @return X coordinate of Vector by reference(can be edited)
      */
    T& X() {
        return _x;
    }
    /**
      * @return Y coordinate of Vector by reference(can be edited)
      */
    T& Y() {
        return _y;
    }
    /**
      * @return Z coordinate of Vector by reference(can be edited)
      */
    T& Z() {
        return _z;
    }
    /**
      * @return X const coordinates of Vector by value
      */
    T X() const {
        return _x;
    }
    /**
      * @return Y const coordinates of Vector by value
      */
    T Y() const {
        return _y;
    }
    /**
      * @return Z const coordinates of Vector by value
      */
    T Z() const {
        return _z;
    }
    /**
      * @brief Copy constructor
      */
    Vector(const Vector& other);
    /**
      * @brief Copy assignment
      */
    Vector& operator=(const Vector& other);
    /**
      * @brief += operator
      */
    Vector& operator+=(const Vector& other);
    /**
      * @brief + operator(using +=)
      */
    Vector operator+(const Vector& other);
    /**
      * @brief -= operator
      */
    Vector& operator-=(const Vector& other);
    /**
      * @brief - operator(using -=)
      */
    Vector operator-(const Vector& other);
    /**
      * @brief Multiplying Vector to numeric value
      * @brief *= operator
      */
    template<typename U>
    Vector& operator*=(const U other);
    /**
      * @brief Multiplying Vector to numeric value
      * @brief * operator(using *=)
      */
    template<typename U>
    Vector operator*(const U other);
    /**
      * @brief Equality of Vector
      */
    bool operator==(const Point<T>& other) const;
    /**
      * @brief !Equality of Vector
      */
    bool operator!=(const Point<T>& other) const;
};
template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
    if (this != &other) {
        this->_x = other._x;
        this->_y = other._y;
        this->_z = other._z;
    }
    return *this;
}
template<typename T>
Vector<T>::Vector(const Vector<T>& other) {
    this->_x = other._x;
    this->_y = other._y;
    this->_z = other._z;
}
template<typename T>
Vector<T>::Vector(const Point<T>& start, const Point<T>& end) {
    this->_x = end.X() - start.X();
    this->_y = end.Y() - start.Y();
    this->_z = end.Z() - start.Z();
}
template <typename T>
Vector<T>::Vector(const Point<T>& other)
{
    this->_x = other.X();
    this->_y = other.Y();
    this->_z = other.Z();
}
template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& other) {
    this->_x += other._x;
    this->_y += other._y;
    this->_z += other._z;
    return *this;
}
template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& other) {
    Vector<T> vec = *this;
    return vec += other;
}
template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& other) {
    this->_x -= other._x;
    this->_y -= other._y;
    this->_z -= other._z;
}
template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& other) {
    Vector<T> vec = *this;
    return vec -= other;
}
template <typename T>
bool Vector<T>::operator==(const Point<T>& other) const
{
    if (fabs(this->_x - other._x) < TOLERANCE && fabs(this->_y - other._y) < TOLERANCE && fabs(this->_z - other._z) < TOLERANCE)
        return true;
    return false;
}
template <typename T>
bool Vector<T>::operator!=(const Point<T>& other) const
{
    return !(*this == other);
}
template<typename T>
template<typename U>
Vector<T>& Vector<T>::operator*=(const U other) {
    this->_x *= static_cast<T>(other);
    this->_y *= static_cast<T>(other);
    this->_z *= static_cast<T>(other);
    return *this;
}

template<typename T>
template<typename U>
Vector<T> Vector<T>::operator*(const U other) {
    Vector<T> vec = *this;
    return vec *= other;
}