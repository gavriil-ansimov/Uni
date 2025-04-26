#include "Triangle.h"

Triangle::Triangle(Point p1, Point p2, Point p3) {
	if ((p1.x == p2.x && p1.y == p2.y) || (p1.x == p3.x && p1.y == p3.y)
		|| (p2.x == p3.x && p2.y == p3.y))
		throw std::exception("Wrong points!");
	p1_ = p1;
	p2_ = p2;
	p3_ = p3;
}

double Triangle::getArea() const {
	return 0.5 * ((p1_.x - p3_.x) * (p2_.y - p3_.y) -
		(p2_.x - p3_.x) * (p1_.y - p3_.y));
}

Point Triangle::getCenter() const {
	Point center = { (p1_.x + p2_.x + p3_.x) / 3 , (p1_.y + p2_.y + p3_.y) / 3 };
	return center;
}

std::string Triangle::getName() const {
	return "TRIANGLE";
}

void Triangle::scale(double k) {
	Point c = this->getCenter();
	p1_.x *= k;
	p1_.x -= c.x;
	p1_.y *= k;
	p1_.y -= c.y;
	p2_.x *= k;
	p2_.x -= c.x;
	p2_.y *= k;
	p2_.y -= c.y;
	p3_.x *= k;
	p3_.x -= c.x;
	p3_.y *= k;
	p3_.y -= c.y;
}