#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "Shape.h"

struct Rectangle : public Shape {
	Rectangle(Point, Point);

	double getArea() const override;
	Point getCenter() const override;
	std::string getName() const override;

	void scale(double k) override;
private:
	Point leftBottom_;
	Point rightTop_;
};	
#endif
