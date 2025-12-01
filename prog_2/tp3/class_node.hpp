#pragma once 
#include "class_point.hpp"

template<typename POINT_T>
class node{
    protected: 
        POINT_T location_;        // pour créer des node<point2d> ou node<point3d>
        int index_;     // index among the nodes of the mesh
        size_t sharing_elements_;  // number of elements sharing this node
        size_t sharing_faces_; // number of faces sharing this node
    public :
        node(POINT_T const& loc=0., int ind=-1, size_t sharing_elements=0, size_t sharing_faces=0);
        node(const node<POINT_T>& n);
        node<POINT_T>& operator=(const node<POINT_T>& n);

        inline int& index() {return index_;}
        inline int index() const {return index_;}
        inline size_t& sharing_elements() {return sharing_elements_;}
        inline size_t sharing_elements() const {return sharing_elements_;}
        inline size_t& sharing_faces() {return sharing_faces_;}
        inline size_t sharing_faces() const {return sharing_faces_;}
        inline void more_sharing_elements() { ++sharing_elements_;}
        inline void less_sharing_elements() { --sharing_elements_;}
        
        friend std::ostream& operator<<(std::ostream& os, const node<POINT_T>& n) {
            os << n.location_;
            return os;
        }

        inline void print() const { std::cout << "index: " << index() << std::endl; }
            
        
};


template<typename POINT_T>
node<POINT_T>::node(POINT_T const& loc, int ind, size_t sharing_elements, size_t sharing_faces) : location_(loc), index_(ind), sharing_elements_(sharing_elements), sharing_faces_(sharing_faces){}


template<typename POINT_T>
node<POINT_T>::node(const node<POINT_T>& n) : location_(n.location_), index_(n.index_), sharing_elements_(n.sharing_elements_), sharing_faces_(n.sharing_faces_){}

template<typename POINT_T>
node<POINT_T>& node<POINT_T>::operator=(const node<POINT_T>& n){
    if(this != &n){
        this->location_ = n.location_;
        this->index_ = n.index_;
        this->sharing_elements_ = n.sharing_elements_;
        this->sharing_faces_ = n.sharing_faces_;
    }
    return *this;
}