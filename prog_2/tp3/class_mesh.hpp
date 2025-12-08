#pragma once
#include "class_element.hpp"
#include "../tp2/class_linked_list.hpp"
#include "class_node.hpp"
#include "class_point.hpp"
#include <vector>
#include <string>

template<class T>
class mesh : public linked_list<T>{
    protected:
        // Rien pour l'instant
    public:
        mesh(T& e):linked_list<T>(e){}   // constructor

        int indexing();   // indexing the nodes in the mesh
}; 

typedef mesh<triangle> triangulation;



template<class T>
int mesh<T>::indexing(){
    for(mesh<T>*p_scan = this ; p_scan; p_scan=(mesh<T>*)p_scan->p_next())  // on reset tous les indices
        p_scan->item().reset_indices();
    int count=1;
    for(mesh<T>*p_scan = this ; p_scan; p_scan=(mesh<T>*)p_scan->p_next())  // on indexe les indices dans l'ordre ensuite
        p_scan->item().indexing(count);
    return count;
} 

