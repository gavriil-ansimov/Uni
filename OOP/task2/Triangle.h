#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Shape.h"

struct Triangle : public Shape {
	Triangle(Point, Point, Point);

	double getArea() const override;
	Point getCenter() const override;
	std::string getName() const override;

	void scale(double k) override;
private:
	Point p1_;
	Point p2_;
	Point p3_;
};

#endif