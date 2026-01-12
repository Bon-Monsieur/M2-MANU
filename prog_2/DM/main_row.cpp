#include <iostream>
#include "class_row.hpp"
#include "class_row_element.hpp"
#include "header/class_dynamic_vector.hpp"


int main(){

    // Création d'une ligne avec un premier élément
    row<double> r(1.0, 0);   // (1, col 0)

    r.append(2.0, 2);
    r.append(3.0, 4);

    std::cout << "Row r : " << std::endl;
    print(r);
    std::cout << std::endl;

    // Accès par colonne
    std::cout << "r[0] = " << r[0] << std::endl;
    std::cout << "r[1] = " << r[1] << std::endl;
    std::cout << "r[2] = " << r[2] << std::endl;
    std::cout << "r[4] = " << r[4] << std::endl;

    // Somme de la ligne
    std::cout << "Somme de la ligne = " << r.row_sum() << std::endl;
    std::cout << std::endl;

    // Multiplication scalaire
    r *= 2.0;
    std::cout << "r *= 2 : " << std::endl;
    print(r);
    std::cout << std::endl;

    // Division scalaire
    r /= 2.0;
    std::cout << "r /= 2 : " << std::endl;
    print(r);
    std::cout << std::endl;

    // Produit scalaire avec un vecteur
    dynamic_vector<double> v(5, 1.0); // vecteur [1,1,1,1,1]
    double dot = r * v;
    std::cout << "Produit scalaire r.v = " << dot << std::endl;

    // Renumérotation des colonnes
    dynamic_vector<int> renum(5);
    renum(0) = 4;
    renum(2) = 1;
    renum(4) = 3;

    r.renumber_columns(renum);
    std::cout << "Apres renumber_columns : " << std::endl;
    print(r);
    std::cout << std::endl;

    // Drop items
    dynamic_vector<int> mask(5, 0);
    mask(1) = 1;
    mask(3) = 1;

    r.drop_items(mask);
    std::cout << "Apres drop_items : " << std::endl;
    print(r);
    std::cout << std::endl;

    std::cout << "\n=== FIN DES TESTS ===\n";

    return 0;
}