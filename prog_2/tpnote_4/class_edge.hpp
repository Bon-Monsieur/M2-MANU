#pragma once

#include "class_node.hpp"
#include "class_face.hpp"
#include "class_element.hpp"

template <typename T>
class edge
{
protected:
    node<T>* n1_; //sommet
    node<T>* n2_; //sommet
    int index_; //index des arêtes
    std::array<int,2> neighbors_; //élément voisins de l'arêtes. on prend pour convetion -1 s'il n'y a pas de voisins
public:
    //edge(node<T>* n1, node<T>* n2): n1_(n1), n2_(n2) {} //définition arête, 2 pointeurs vers des nodes, le constructeur

    edge(node<T>*  n1, node<T>* n2): n1_(n1), n2_(n2), index_(-1) {neighbors_[0]=-1, neighbors_[1]=-1} //nouveau constructeur qui prend en compte index et neighbors

    node<T>* first() const { return n1_; } //accesseurs pour lire les sommets
    node<T>* second() const { return n2_; }

    bool is_same(const edge<T>& other) const; //accesseur pour savoir si 2 arêtes sont égales

    double length() const; //longueur arête


    template <typename S> //operateur d'affichage
    friend std::ostream& operator<<(std::ostream& os, const edge<S>& e);

    int index() const; //accesseur pour l'indice global de l'arête
    void set_index(int i); //acceseur pour indexer les arêtes

    int neighbor(std::size_t i) const; //accesseur pour lire un voisin
    void set_neighbor(std::size_t i, int elem_index); //acceseur pour modifier un voisin

    point2d midpoint() const; //accesseur pour le milieu d'une arête 

    point2d normal() const; //accesseur pour la normale unitaire à l'arête

    bool contains(const node<T>* n) const; //acceseur test d'appartenance à un sommet
};

template <typename T> //savoir si 2 arêtes sont égales
bool edge<T>::is_same(const edge<T>& other) const
{
    return ( (n1_ == other.n1_ && n2_ == other.n2_) ||
             (n1_ == other.n2_ && n2_ == other.n1_) );
}


template <typename T> //longueur de l'arête
double edge<T>::length() const
{
    double dx = n1_->x() - n2_->x();
    double dy = n1_->y() - n2_->y();
    return std::sqrt(dx*dx + dy*dy);
}


template <typename T> //opérateur affichage
std::ostream& operator<<(std::ostream& os, const edge<T>& e)
{
    os << "Edge("
       << e.first()->index()
       << ", "
       << e.second()->index()
       << ")";
    return os;
}


template <typename T> //indice arête
int edge<T>::index() const
{
   return index_;
}

template <typename T> //indexer arêtes
void edge<T>::set_index(int i)
{
    index_=i;
}

template <typename T> //lire voisin
int edge<T>::neighbor(std::size_t i) const
{
    return neighbors_[i];
}

template <typename T> //modifier voisin
void edge<T>::set_neighbor(std::size_t i, int elem_index)
{
    neighbors_[i]=elem_index;
}

template <typename T> //calculer le milieu de l'arête
point2d edge<T>::midpoint() const
{
    double xm=0.5*(n1_->x()+n2_->x());
    double ym=0.5*(n1_->y()+n2_->y());
    return point2d(xm,ym);

}

template <typename T> //calculer la normale sortante à une arête
point2d edge<T>::normal() const
{
    double dx=n2_->x()-n1_->x();
    double dy=n2_->y()-n1_->y();

    //rotation à 90°

    double nx=dx;
    double ny=-dy;

    double L=std::sqrt(nx*nx+ny*ny);

    if (L<1e-14) //zéro machine
       return point2d(0.,0.);
    
    return point2d(nx / L, ny / L);

}

template <typename T> //test d'appartenance d'un sommet à une arête
bool edge<T>::contains(const node<T>* n) const
{
    return (n == n1_ || n == n2_);
}
