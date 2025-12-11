#include <iostream>
#include "class_row_element.hpp"  

int main() {
    std::cout << "=== Tests de row_element ===\n\n";

    // --- Constructeurs ---
    row_element<double> a(3.5, 2);
    row_element<double> b(1.5, 5);
    row_element<double> c;  // valeur=0, colonne=-1

    std::cout << "Constructeurs :\n";
    print(a);
    print(b);
    print(c);
    std::cout << "\n";

    // --- Opérateur d'affectation ---
    std::cout << "Affectation c = a\n";
    c = a;
    print(c);
    std::cout << "\n";

    // --- Getters / Setters ---
    std::cout << "Setters sur c : valeur=10, colonne=7\n";
    c.set_value(10);
    c.set_column(7);
    print(c);
    std::cout << "\n";

    // --- Opérateurs arithmétiques internes ---
    std::cout << "Tests des opérateurs internes :\n";
    row_element<double> x(4.0, 1);

    std::cout << "x += 2 : ";
    x += 2;
    print(x);

    std::cout << "x -= 1.5 : ";
    x -= 1.5;
    print(x);

    std::cout << "x *= 3 : ";
    x *= 3;
    print(x);

    std::cout << "x /= 2 : ";
    x /= 2;
    print(x);

    std::cout << "x += a : ";
    x += a;
    print(x);

    std::cout << "x -= b : ";
    x -= b;
    print(x);

    std::cout << "x *= b : ";
    x *= b;
    print(x);

    std::cout << "x /= a : ";
    x /= a;
    print(x);

    std::cout << "\n";

    // --- Opérateurs arithmétiques externes ---
    std::cout << "Opérateurs externes :\n";
    print(a + 2.);
    print(2. + a);
    print(a - 5.0);
    print(5.0 - a);
    print(a * 3.);
    print(3. * a);
    print(a / 2.);
    std::cout << "\n";

    // --- Comparateurs ---
    std::cout << "Comparateurs :\n";
    std::cout << "a < b ? " << (a < b) << "\n";
    std::cout << "a > b ? " << (a > b) << "\n";
    std::cout << "a == b ? " << (a == b) << "\n";

    row_element<double> d(99, a.get_column());
    std::cout << "a == d ? (mêmes colonnes) " << (a == d) << "\n";

    std::cout << "\n=== Fin des tests ===\n";

    return 0;
}
