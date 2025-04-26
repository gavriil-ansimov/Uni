#include <iostream>

#include "Rectangle.h"
#include "Triangle.h"

bool operator<(const Shape& l, const Shape& r) {
	return (l.getArea() < r.getArea());
}

void insertionSort(Shape* arrayPtr[], int length)
{
	Shape* temp = nullptr;
	int	item = 0;
	for (int counter = 1; counter < length; counter++)
	{
		temp = arrayPtr[counter];
		item = counter - 1;
		while (item >= 0 && *temp < *arrayPtr[item])
		{
			arrayPtr[item + 1] = arrayPtr[item];
			arrayPtr[item] = temp;
			item--;
		}
	}
}

int main() {
	Rectangle rec1({ 10, 10 }, { 23.5, 51.2 });
	Rectangle rec2({ 46.8, 85.9 }, { 48.6, 123.4 });
	Rectangle rec3({ -64.2, -1.5 }, { 44.5, 59.1 });
	Triangle tri1({ 0, 0 }, { 10, 10 }, { 100, 100 });
	Triangle tri2({ -15.3, 18.2 }, { 38.5, -29.1 }, { 50, 68.3 });
	Shape* arr[5] = { &rec1, &rec2 , &rec3, &tri1 , &tri2 };
	insertionSort(arr, 5);
	for (int i = 0; i < 5; ++i) {
		Point center = (*arr[i]).getCenter();
		std::cout << (*arr[i]).getName() << ", Center = (" << center.x << ',' << center.y << "), " << (*arr[i]).getArea() << std::endl;
		(*arr[i]).scale(2);
	}
	std::cout << std::endl << "k = -1.5" << std::endl;
	for (int i = 0; i < 5; ++i) {
		Point center = (*arr[i]).getCenter();
		std::cout << (*arr[i]).getName() << ", Center = (" << center.x << ',' << center.y << "), " << (*arr[i]).getArea() << std::endl;
	}
}