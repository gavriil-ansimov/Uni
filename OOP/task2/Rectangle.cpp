#include "Rectangle.h"

Rectangle::Rectangle(Point lb, Point rt)
{
	if (lb.x >= rt.x || lb.y >= rt.y)
		throw std::exception("Wrong points!");
	leftBottom_ = lb;
	rightTop_ = rt;
}

double Rectangle::getArea() const {
	return (rightTop_.x - leftBottom_.x) * (rightTop_.y - leftBottom_.y);
}

Point Rectangle::getCenter() const {
	Point center = { (leftBottom_.x + rightTop_.x) / 2, (leftBottom_.y + rightTop_.y) / 2 };
	return center;
}

std::string Rectangle::getName() const {
	return "RECTANGLE";
}

void Rectangle::scale(double k) {
	if (k >= 0) 
	{
		Point c = this->getCenter();
		leftBottom_.x *= k;
		leftBottom_.x -= c.x;
		leftBottom_.y *= k;
		leftBottom_.y -= c.y;
		rightTop_.x *= k;
		rightTop_.x -= c.x;
		rightTop_.y *= k;
		rightTop_.y -= c.y;
	}
}

