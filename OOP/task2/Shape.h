#ifndef SHAPE_H
#define SHAPE_H

#include <iostream>

#include "Point.h"

struct Shape {
	virtual ~Shape() = default;

	virtual double getArea() const = 0;
	virtual Point getCenter() const = 0;
	virtual std::string getName() const = 0;

	virtual void scale(double k) = 0;
};

#endif