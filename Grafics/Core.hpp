#pragma once
#include<math.h>
#define TOLERANCE 0.0000001

static bool IsEqualD(double a, double b) {
	return fabs(a - b) < TOLERANCE;
};

enum class RELATIVE_POSITION {
	LEFT = 0, RIGTH, BEHIND, BEYOND, BETWEEN, ORIGIN, DESTINATION
};

static bool _xor(bool x, bool y) {
	return x ^ y;
}