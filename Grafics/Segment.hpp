#pragma once
#include "Point.hpp"
#include <utility>

/**
  * @brief Segment class. Defined by 2 points Point, Point
  */
template<typename T>
class Segment {
private:
    Point<T> p1;
    Point<T> p2;
public:
    /**
      * @brief Create a Segment using Point p1, Point p2.
      */
    Segment(const Point<T>& p1, const Point<T>& p2) : p1(p1), p2(p2) {}
    /**
      * @brief Copy constructor for Segment
      */
    Segment(const Segment<T>& other);
    /**
      * @brief Assignment operator
      */
    Segment& operator=(const Segment<T>& other);
    /**
      * @brief Return a pair(std::pair<Point, Point>) of points of given segments
      */ 
    std::pair<Point<T>, Point<T>> getPoints() const {
        return std::make_pair(p1, p2);
    }
};

template <typename T>
Segment<T>::Segment(const Segment<T>& other)
{
    this->p1 = other.p1;
    this->p2 = other.p2;
}

template <typename T>
Segment<T>& Segment<T>::operator=(const Segment<T>& other)
{
    if (this != &other) {
        this->p1 = other.p1;
        this->p2 = other.p2;
    }
    return *this;
}