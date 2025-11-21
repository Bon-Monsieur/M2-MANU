#include "class_linked_list.hpp"


template<typename T>
void print(linked_list<T> const& head) {
    linked_list<T> const* current = &head;

    while (current != nullptr) {
        std::cout << current->item() << " ";
        current = current->p_next();
    }

    std::cout << std::endl;
}




int main() {
    // Création dynamique de chaque nœud
    linked_list<double>* node1 = new linked_list<double>(3.21);
    linked_list<double>* node2 = new linked_list<double>(4.56, node1);
    linked_list<double>* node3 = new linked_list<double>(7.89, node2);

    print(*node3);
    std::cout << "Last item: " << node3->last().item() << std::endl;

    
    delete node2; // Supprime toute la liste après le node2 (non compris)

    std::cout << "After deletion of node1:" << std::endl;
    print(*node3); // Affiche la liste après la suppression
    
    return 0;
}