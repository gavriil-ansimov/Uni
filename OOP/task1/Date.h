#pragma once
#include <iostream>

struct Date {
	
private:
	int day_;
	int month_;
	int year_;
public:
	Date();
	Date(int day, int month, int year);

	Date& operator+=(int d);
	Date& operator-=(int d);

	int maxday()const;
	int get_day()const;
	int get_month()const;
	int get_year()const;
	bool vys(int year)const;
	int days_count() const;
	int difference(const Date&)const;
	
	
	friend bool operator<(const Date& l, const Date& r);
	friend bool operator>(const Date& l, const Date& r);
	friend bool operator<=(const Date& l, const Date& r);
	friend bool operator>=(const Date& l, const Date& r);
	friend bool operator==(const Date& l, const Date& r);
	friend bool operator!=(const Date& l, const Date& r);
	
};

std::ostream& operator<<(std::ostream&, const Date&);
