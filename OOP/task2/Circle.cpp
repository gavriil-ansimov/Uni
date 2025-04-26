#include "Circle.h"

double PI = atan(1);

Circle::Circle(Point center, double r) : center_(center), r_(r)
{
	if (r < 0)
		throw std::exception("Incorrect values!");
}

double Circle::getArea() const {
	return PI * r_ * r_;
}

Point Circle::getCenter() const {
	return center_;
}

void Circle::scale(double k) {
	if (k < 0) {
		center_.x *= k;
		center_.y *= k;
	}
	r_ *= abs(k);
}

std::string Circle::getName() const {
	return "CIRCLE";
}