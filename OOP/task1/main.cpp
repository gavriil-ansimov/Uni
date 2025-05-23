﻿//Разработать указанные классы.
//• Можно создать более одного класса, если вы считаете, что это логично для вашей задачи.
//• Реализовать необходимые конструкторы.
//• Конструкторы копирования / перемещения, операторы присваивания и деструкторы в этом
//задании не требуются, т.к.классы не должны напрямую владеть ресурсами(в т.ч.
//	динамической памятью)
//	• Реализовать методы доступа к полям(если нужно).
//	• Реализовать перегруженный оператор << для вывода информации в поток.
//	• Классы должны контролировать свои данные и не допускать создания объектов с
//	некорректным состоянием(можно выбрасывать исключения с помощью throw).
//	• Некоторые методы можно реализовать в виде перегруженных операторов(если подходят
//		по смыслу).
//	• Объявление класса должно быть в заголовочном файле(.h, не забывайте про header guard),
//	а определения методов в файле реализации(.cpp).
//	• В main продемонстрировать работу всех(!) написанных методов и конструкторов.
//Разработать класс для работы с обыкновенными дробями(числитель / знаменатель).Методы
//для сложения, умножения, вычитания, деления, а также для сравнения дробей.Метод для
//сокращения дроби.
#include <iostream>
#include "Fraction.h"
#include "Date.h"
int main()
{
	Fraction f1;
	std::cout << "f1 = " << f1 << std::endl;
	f1 = f1 + static_cast<Fraction>(5);
	std::cout << "f1 + 5 = " << f1 << std::endl;
	Fraction f2(-8, -15);
	std::cout << "f2 = " << f2 << std::endl;
	Fraction f3;
	f3 = f1 - f2;
	std::cout << "f3 = f1 - f2 = " << f3 << std::endl;
	Fraction f4;
	f4 = f3 * f2;
	std::cout << "f4 = f3 * f2 = " << f4 << std::endl;
	f4 = f4 / f2;
	std::cout << "f4 / f2 = " << f4 << std::endl;
	Fraction f5(-20, 30);
	std::cout << "f5 = " << f5 << std::endl;
	f5.reduce();
	std::cout << "Reduced f5 = " << f5 << std::endl;
	std::cout << "f1 < f2: " << ((f1 < f2) ? "True": "False") << std::endl;
	std::cout << "f3 > f4: " << ((f3 > f4) ? "True" : "False") << std::endl;
	std::cout << "f1 <= f3: " << ((f1 <= f3) ? "True" : "False") << std::endl;
	std::cout << "f1 >= f3: " << ((f1 >= f3) ? "True" : "False") << std::endl;
	std::cout << "f2 == f5: " << ((f2 == f5) ? "True" : "False") << std::endl;
	std::cout << "f2 != f5: " << ((f2 != f5) ? "True" : "False") << std::endl;


}