#include "Date.h"
#include <iostream>


Date::Date() : day_(0), month_(0), year_(0)
    {}
Date::Date(int day, int month, int year) {
    if (month < 1 || month > 12)
        throw std::exception("Incorrect values");
    if (year < 0)
        throw std::exception("Incorrect values");
    int days[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };
    if (vys(year))days[1] = 29;
    if (day < 0 || day > days[month])
        throw std::exception("Incorrect values");
    day_ = day;
    month_ = month;
    year_ = year;
}

Date& Date::operator+=(int d) {
    Date tmp(this->day_, this->month_, this->year_);
    int dd = d;
    if (day_ + d > maxday()) {
        d -= maxday() - day_;
        if (++month_ > 12) {
            month_ = 1;
            year_++;
        }
        while (d / maxday()) {
            if (++month_ > 12) {
                month_ = 1;
                year_++;
            }
            d -= maxday();
        }
        day_ = d;
    }
    else day_ += d;
    while (abs(tmp.difference(*this)) > dd) {
        *this -= 1;
    }
    return *this;
}
Date& Date::operator-=(int d) {
    Date tmp(this->day_, this->month_, this->year_);
    int dd = d;
    if (day_ - d < 1) {
        d -= day_;
        if (--month_ == 0) {
            month_ = 12;
            year_--;
        }
        while (d / maxday()) {
            d -= maxday();
            if (--month_ == 0) {
                month_ = 12;
                year_--;
            }
        }
        day_ = maxday() - d;
    }
    else day_ -= d;
    while (abs(this->difference(tmp)) > dd) {
        *this += 1;
    }
    return *this;
}


    int Date::maxday()const {
        int days[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };
        if (vys(year_))days[1] = 29;
        return days[month_ - 1];
    }

    int Date::get_day() const {
        return this->day_;
    }
    int Date::get_month() const {
        return this->month_;
    };
    int Date::get_year() const {
        return this->year_;
    };

    bool Date::vys(int year) const
    {
        bool res = false;
        if (year % 4 == 0)
            res = true;
        if (year % 100 == 0)
            res = false;
        if (year % 400 == 0)
            res = true;
        return res;
    }

    int Date::days_count() const
    {
        int d = this->day_, m = this->month_, y = this->year_;
        int k = d;
        switch (m - 1)
        {
        case 12: k += 31;
        case 11: k += 30;
        case 10: k += 31;
        case  9: k += 30;
        case  8: k += 31;
        case  7: k += 31;
        case  6: k += 30;
        case  5: k += 31;
        case  4: k += 30;
        case  3: k += 31;
        case  2: if (vys(y)) k += 29; else k += 28;
        case  1: k += 31;
        }
        k += (y - 1) * 365 + ((y - 1) / 4);
        return k;
    }
    int Date::difference(const Date& other) const {

        int k = this->days_count() - other.days_count();
        return k;
    };

    bool operator<(const Date& l, const Date& r) {
        if (l.difference(r) < 0) {
            return 1;
        }
        else if (l.difference(r) > 0) {
            return 0;
        }
    }
    bool operator>(const Date& l, const Date& r) {
        return r < l;
    }
    bool operator<=(const Date& l, const Date& r) {
        return !(l > r);
    }
    bool operator>=(const Date& l, const Date& r) {
        return !(l < r);
    }
    bool operator==(const Date& l, const Date& r) {
        if (l.difference(r) == 0) {
            return 1;
        }
        else {
            return 0;
        }
    }
    bool operator!=(const Date& l, const Date& r) {
        return !(l == r);
    }

std::ostream& operator << (std::ostream& os, const Date& date)
{
    os << date.get_day() << '.' << date.get_month() << '.' << date.get_year();
    return os;
}
