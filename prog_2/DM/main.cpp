#include <iostream>
#include <fstream>
#include "class_row.hpp"
#include "header/class_dynamic_vector.hpp"
#include "class_sparse_matrix.hpp"

int main()
{
    try
    {
        std::cout << "===== TEST 1 : Constructeur par défaut =====" << std::endl;
        sparse_matrix<double> M0(3);
        std::cout << "Nombre de lignes : " << M0.row_number() << std::endl;
        std::cout << "Nombre de colonnes : " << M0.column_number() << std::endl;
        std::cout << "Ordre : " << M0.order() << std::endl;
        print(M0);

        std::cout << "\n===== TEST 2 : Constructeur diagonale =====" << std::endl;
        sparse_matrix<double> M1(4, 2.0);
        std::cout << "Nombre de lignes : " << M1.row_number() << std::endl;
        std::cout << "Nombre de colonnes : " << M1.column_number() << std::endl;
        std::cout << "Ordre : " << M1.order() << std::endl;
        print(M1);

        std::cout << "\n===== TEST 3 : operator()(i,j) =====" << std::endl;
        std::cout << "M1(0,0) = " << M1(0,0) << std::endl;
        std::cout << "M1(1,1) = " << M1(1,1) << std::endl;

        try {
            std::cout << "M1(0,1) = " << M1(0,1) << std::endl;
        }
        catch (std::exception const& e) {
            std::cout << "Exception attendue : " << e.what() << std::endl;
        }

    }
    catch (std::exception const& e)
    {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }

    return 0;
}
