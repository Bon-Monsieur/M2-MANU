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
        ~element(){
            for (int i = 0; i < NVERTICES; ++i) {
                vertices_[i]->less_sharing_elements();  // quand on détruit un élément , il faut diminuer le nombre d'éléments partageant chaque noeud (car ils ne savent pas à quels éléments ils appartiennent)
            }
        }
        element<T, NVERTICES>& operator=(const element<T, NVERTICES>& e); // opérateur d'affectation
        node<T>& operator()(size_t i) { return *vertices_[i];} // renvoie la i-ème node de l'élément et peut la modifier
        node<T> const& operator[](size_t i) const { return *vertices_[i];} // read-only

        node<T>* vertex(size_t i) { return vertices_[i];}   // renvoie le pointeur vers le i-ème noeud de l'élément
        node<T> const* vertex(size_t i) const { return vertices_[i]; } 

        void reset_indices(); // met les indices des noeuds à -1 (leur valeur initiale)
        void indexing(int& count); // assigne des indices aux noeuds de l'élément en partant de la valeur count et incrémente count en conséquence

        
        friend std::ostream& operator<<(std::ostream& os, const element<T,NVERTICES>& e) {
            for (int i = 0; i < NVERTICES; ++i) {
                os << *(e.vertices_[i])<< " ";
            }
            return os;
        }
        
};

typedef element<point2d,3> triangle;


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


template<class T, int NVERTICES>
void element<T, NVERTICES>::reset_indices(){
    for (int i = 0; i < NVERTICES; ++i) {
        vertices_[i]->index() = -1; // réinitialise l'indice du noeud à -1
    }
}

template<class T, int NVERTICES>
void element<T, NVERTICES>::indexing(int& count){
    for (int i = 0; i < NVERTICES; ++i) {
        if (vertices_[i]->index() < 0) { // si le noeud n'a pas encore d'indice
            vertices_[i]->index() = count++; // on lui assigne l'indice courant et augmente count
        }
    }
}



template<class T, int NVERTICES>
size_t operator<(node<T> const& n, element<T,NVERTICES> const& e){
    for (int i = 0; i < NVERTICES; ++i) {
        if (&n == e.vertex(i)) {
            return true;
        }
    }
    return false;
}

