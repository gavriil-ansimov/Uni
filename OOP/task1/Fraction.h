#ifndef FRACTION_H
#define FRACTION_H

#include <iostream>
class Fraction {
	int numenator_;
	int denomenator_;
public:
	Fraction();
	Fraction(int, int);
	Fraction(int);

	int getNum()const;
	int getDen()const;

	void reduce();
};

std::ostream& operator<<(std::ostream&, const Fraction&);

Fraction operator+(const Fraction&, const Fraction&);
Fraction operator-(const Fraction&, const Fraction&);
Fraction operator*(const Fraction&, const Fraction&);
Fraction operator/(const Fraction&, const Fraction&);

bool operator<(const Fraction&, const Fraction&);
bool operator>(const Fraction&, const Fraction&);
bool operator<=(const Fraction&, const Fraction&);
bool operator>=(const Fraction&, const Fraction&);
bool operator==(const Fraction&, const Fraction&);
bool operator!=(const Fraction&, const Fraction&);

#endif