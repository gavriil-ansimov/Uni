#include "Fraction.h"

Fraction::Fraction(): numenator_(0), denomenator_(1)
{}
Fraction::Fraction(int num, int den)
{
	if (den == 0)
		throw std::exception("Division by zero!");
	if ((num < 0 && den < 0) || (num > 0 && den > 0)) {
		numenator_ = abs(num);
		denomenator_ = abs(den);
	}
	else {
		numenator_ = -abs(num);
		denomenator_ = abs(den);
	}
}
Fraction::Fraction(int num) : numenator_(num), denomenator_(1)
{}

int Fraction::getNum()const {
	return numenator_;
}
int Fraction::getDen()const {
	return denomenator_;
}

void Fraction::reduce() {
	int gcd = 1;
	int num = numenator_;
	int den = denomenator_;
	if (num < den) {
		std::swap(num, den);
	}
	while (num % den != 0) {
		num = num % den;
		std::swap(num, den);
	}
	gcd = den;
	numenator_ /= gcd;
	denomenator_ /= gcd;
}

std::ostream& operator<<(std::ostream& os, const Fraction& fraction) {
	os << fraction.getNum() << '/' << fraction.getDen();
	return os;
}

Fraction operator+(const Fraction& l, const Fraction& r) {
	int den = r.getDen() * l.getDen();
	int num = r.getNum() * l.getDen() + l.getNum() * r.getDen();
	Fraction res(num, den);
	res.reduce();
	return res;
}
Fraction operator-(const Fraction& l, const Fraction& r) {
	int den = r.getDen() * l.getDen();
	int num = r.getNum() * l.getDen() - l.getNum() * r.getDen();
	Fraction res(num, den);
	res.reduce();
	return res;
}
Fraction operator*(const Fraction& l, const Fraction& r) {
	int den = r.getDen() * l.getDen();
	int num = r.getNum() * l.getNum();
	Fraction res(num, den);
	res.reduce();
	return res;
}
Fraction operator/(const Fraction& l, const Fraction& r) {
	int den = r.getNum() * l.getDen();
	int num = r.getDen() * l.getNum();
	Fraction res(num, den);
	res.reduce();
	return res;
}

bool operator<(const Fraction& l, const Fraction& r) {
	int num1 = l.getNum() * r.getDen();
	int num2 = r.getNum() * l.getDen();
	if (num1 < num2)
		return 1;
	return 0;
}
bool operator>(const Fraction& l, const Fraction& r) {
	return r < l;
}
bool operator<=(const Fraction& l, const Fraction& r) {
	return !(l > r);
}
bool operator>=(const Fraction & l, const Fraction & r){
	return !(l < r);
}
bool operator==(const Fraction& l, const Fraction& r) {
	int num1 = l.getNum() * r.getDen();
	int num2 = r.getNum() * l.getDen();
	if (num1 == num2)
		return 1;
	return 0;
}
bool operator!=(const Fraction& l, const Fraction& r) {
	return !(r == l);
}