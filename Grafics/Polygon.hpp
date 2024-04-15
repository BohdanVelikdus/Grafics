#pragma once
#include "Point.hpp"
#include "Vector.hpp"
#include <vector>
#include <list>

template <typename T>
struct Vertex
{
	Vector<T> point;
	Vertex* next = nullptr;
	Vertex* prev = nullptr;

	Vertex(const Vector<T>& _point, const Vertex<T>& _next, const Vertex<T>& _prev) :
		point(_point), next(_next), prev(_prev) {
	};

	Vertex(const Vector<T>& _point) :
		point(_point) {
	};
};


template <typename T>
class Polygon
{
private:
	std::vector<Vertex<T>*> vertex_list;
public:
	Polygon(std::list<Vector<T>>& points) {
		const int size = points.size();

		if (size < 3) {
			std::count << "Not enough points\n";
			return;
		}

		for (auto _point : points) {
			vertex_list.push_back(new Vertex(_point));
		}

		for (size_t i = 0; i < size; i++) {
			vertex_list[i].next = &vertex_list[(i + 1) % size];
			
			if (i != 0) {
				vertex_list.prev = &vertex_list[i - 1];
			}
			else {
				vertex_list.prev = &vertex_list[size - 1];
			}
		}
	}

	std::vector<Vertex<T>> getVertices() {
		return vertex_list;
	}
};
