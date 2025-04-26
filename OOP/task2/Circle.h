#ifndef CIRCLE_H
#define CIRCLE_H

#include "Shape.h"

class Circle : public Shape {
	Point center_;
	double r_;
public:
	Circle(Point, double);

	double getArea() const override;
	Point getCenter() const override;
	void scale(double k) override;
	std::string getName() const override;
};

#endif