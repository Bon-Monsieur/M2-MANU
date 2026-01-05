#include <iostream>

#include "class_row.hpp"
#include "header/class_dynamic_vector.hpp"

int main()
{
    std::cout << "===== TEST drop_items(mask) =====\n";

    // 1) Construction de la row
    row<double> r(1.0, 0);
    r.append(2.0, 1);
    r.append(3.0, 2);
    r.append(4.0, 3);
    r.append(5.0, 4);

    std::cout << "\nRow AVANT drop_items :\n";
    print(r);
    std::cout << "\n";

    // 2) Création du masque
    // mask[i] = 1 -> on garde la colonne i
    // mask[i] = 0 -> on supprime la colonne i
    dynamic_vector<int> mask(5);
    mask(0) = 1;
    mask(1) = 0;
    mask(2) = 1;
    mask(3) = 0;
    mask(4) = 1;
    // 3) Application du masking
    r.drop_items(mask);

    // 4) Résultat
    std::cout << "\nRow APRES drop_items :\n";
    print(r);
    std::cout << "\n";

    std::cout << "===== FIN DU TEST =====\n";
    return 0;
}
