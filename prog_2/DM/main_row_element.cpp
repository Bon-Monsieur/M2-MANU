#include <iostream>
#include "class_row.hpp"
#include "class_row_element.hpp"
#include "header/class_dynamic_vector.hpp"

int main()
{
    std::cout << "=============================\n";
    std::cout << " TEST DE row_element\n";
    std::cout << "=============================\n";

    row_element<double> e1(2.0, 1);
    row_element<double> e2(3.0, 1);

    std::cout << "e1 = " << e1 << std::endl;
    std::cout << "e2 = " << e2 << std::endl;

    e1 += 5.0;
    std::cout << "e1 += 5 -> " << e1 << std::endl;

    e1 -= 2.0;
    std::cout << "e1 -= 2 -> " << e1 << std::endl;

    e1 *= 3.0;
    std::cout << "e1 *= 3 -> " << e1 << std::endl;

    e1 /= 3.0;
    std::cout << "e1 /= 3 -> " << e1 << std::endl;

    e1 += e2;
    std::cout << "e1 += e2 -> " << e1 << std::endl;

    row_element<double> e3 = e1 * 2.0;
    std::cout << "e3 = e1 * 2 -> " << e3 << std::endl;

    std::cout << "Comparaisons de colonnes:\n";
    std::cout << "(e1 == e2) : " << (e1 == e2) << std::endl;
    std::cout << "(e1 <  e2) : " << (e1 <  e2) << std::endl;
    std::cout << "(e1 >  e2) : " << (e1 >  e2) << std::endl;

    return 0;
}