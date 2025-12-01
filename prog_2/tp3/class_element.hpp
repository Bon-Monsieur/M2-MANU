#pragma once
#include "class_node.hpp"
#include "class_point.hpp"
// On utilise des adresses de point car les points peuvent être partagés entre plusieurs éléments

template<class T, int NVERTICES>
class element{
    protected: 
        node<T>* vertices_[NVERTICES];  // array of pointers to the nodes of the element
    public :
        element(node<T>& v1, node<T>& v2, node<T>& v3); // constructeur avec 3 adresses de points
        element(const element<T, NVERTICES>& e); // constructeur de copie 
        element<T, NVERTICES>& operator=(const element<T, NVERTICES>& e); // opérateur d'affectation
        ~element(){
            for (int i = 0; i < NVERTICES; ++i) {
                vertices_[i]->less_sharing_elements();  // quand on détruit un élément , il faut diminuer le nombre d'éléments partageant chaque noeud (car ils ne savent pas à quels éléments ils appartiennent)
            }
        }

};


template<class T, int NVERTICES>
element<T, NVERTICES>::element(node<T>& v1, node<T>& v2, node<T>& v3){
    vertices_[0] = &v1;
    vertices_[1] = &v2;
    vertices_[2] = &v3;
    v1.more_sharing_elements(); // augmente le nombre d'éléments partageant ce noeud
    v2.more_sharing_elements(); // augmente le nombre d'éléments partageant ce noeud
    v3.more_sharing_elements(); // augmente le nombre d'éléments partageant ce noeud
}

template<class T, int NVERTICES>
element<T, NVERTICES>::element(const element<T, NVERTICES>& e){
    for (int i = 0; i < NVERTICES; ++i) {
        vertices_[i] = e.vertices_[i];
        vertices_[i]->more_sharing_elements(); // augmente le nombre d'éléments partageant ce noeud
    }
}

template<class T, int NVERTICES>
element<T, NVERTICES>& element<T, NVERTICES>::operator=(const element<T, NVERTICES>& e){
    if(this != &e){
        // Diminuer le nombre d'éléments partageant les anciens noeuds
        for (int i = 0; i < NVERTICES; ++i) {
            vertices_[i]->less_sharing_elements();
        }
        // Copier les nouveaux noeuds et augmenter leur compteur
        for (int i = 0; i < NVERTICES; ++i) {
            vertices_[i] = e.vertices_[i];
            vertices_[i]->more_sharing_elements();
        }
    }
    return *this;
}