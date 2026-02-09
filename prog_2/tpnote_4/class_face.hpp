#pragma once

#include "class_element.hpp"
#include "class_node.hpp"
#include "class_point.hpp"
#include <cmath>

#include <iostream>

template <typename T, int NVERTICES>
class face
{
protected:
    node<T>* vertices_[NVERTICES];
    int index_;
    std::array<int, NVERTICES> neighbors_;         // indices des triangles adjacents (-1 si bord)
    std::array<triangle*, NVERTICES> p_neighbors_; // adresses des triangles adjacents (nullptr si bord)
    point2d centroid_;
    double length_;
    point2d normal_;

public:
    face();                                                                                      // default constructor
    face(node<T>&, node<T>&, int ind = -1, int tri_1 = -1, int tri_2 = -1, triangle* = nullptr); // constructor for segments
    face(face<T, NVERTICES> const&);                                                             // copy constructor
    const face<T, NVERTICES>& operator=(face<T, NVERTICES> const&);
    ~face();
    node<T>& operator()(size_t i) { return *(this->vertices_[i]); }             //  read/write ith vertex
    const node<T>& operator[](size_t i) const { return *(this->vertices_[i]); } //  read only ith vertex

    inline int index() const { return index_; };
    inline int& index() { return index_; };

    inline int neighbor(std::size_t i) const { return this->neighbors_[i]; };
    inline int& neighbor(std::size_t i) { return this->neighbors_[i]; };

    inline triangle* p_neighbor(std::size_t i) const { return this->p_neighbors_[i]; };
    inline void set_p_neighbor(std::size_t i, triangle* pt) { p_neighbors_[i] = pt; };

    // TP4
    inline point2d centroid() const { return centroid_; };
    inline point2d& centroid() { return centroid_; };

    inline double length() const { return length_; };
    inline double& length() { return length_; };

    point2d compute_centroid();
    double compute_length();

    inline point2d normal() const { return normal_; };
    inline point2d& normal() { return normal_; };

    point2d compute_unit_normal();
    
};

// face for d=2
typedef face<point2d, 2> edge;


// méthode rajoutees:

template <typename T, int NVERTICES>
point2d face<T, NVERTICES>::compute_unit_normal(){
    point2d temp;

        double dx=vertices_[1]->x()-vertices_[0]->x();
        double dy=vertices_[1]->y()-vertices_[0]->y();

    //rotation à 90°
    double nx=-dy;
    double ny=dx;
    
    normal_ = point2d(nx / length_, ny / length_);


    return normal_;
}


template <typename T, int NVERTICES>
double face<T, NVERTICES>::compute_length(){
    double dx = (this->vertices_[0])->x() - (this->vertices_[1])->x();
    double dy = (this->vertices_[0])->y() - (this->vertices_[1])->y();

    double l = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    length_ = l;
    return l;
}


template <typename T, int NVERTICES>
point2d face<T, NVERTICES>::compute_centroid(){
    double xm=0.;
    double ym=0.;

    for (auto ii=0; ii<NVERTICES; ii++){
        xm+= this->vertices_[ii]->x();
        ym+= this->vertices_[ii]->y();
    }

    xm/=NVERTICES; 
    ym/=NVERTICES;

    centroid_= point2d(xm, ym);
    return centroid_;
}



template <typename T, int NVERTICES>
face<T, NVERTICES>::face()
{
    for (auto i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i] = nullptr;
    }
} //  default constructor

template <typename T, int NVERTICES>
face<T, NVERTICES>::face(node<T>& a, node<T>& b, int ind, int tri_1, int tri_2, triangle* pt) : index_(ind), neighbors_{tri_1, tri_2}, p_neighbors_{pt, nullptr}
{
    this->vertices_[0] = &a;
    this->vertices_[1] = &b;

    for (auto i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i]->more_sharing_faces();
    }

} //  constructor

template <typename T, int NVERTICES>
face<T, NVERTICES>::face(face<T, NVERTICES> const& e) : index_(e.index_), neighbors_(e.neighbors_), p_neighbors_(e.p_neighbors_)
{
    for (auto i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i] = e.vertices_[i];
        this->vertices_[i]->more_sharing_faces();
        // std::cout << "moreSharingFaces called in copy constructor" << std::endl;
    }
} //  copy constructor

template <typename T, int NVERTICES>
const face<T, NVERTICES>& face<T, NVERTICES>::operator=(face<T, NVERTICES> const& e)
{
    if (this != &e)
    {
        for (auto i = 0; i < NVERTICES; i++)
        {
            this->vertices_[i] = e.vertices_[i];
            this->vertices_[i]->more_sharing_faces();
        }

        this->index_       = e.index();
        this->neighbors_   = e.neighbors_;
        this->p_neighbors_ = e.p_neighbors_;
    }
    return *this;
} //  assignment operator

template <typename T, int NVERTICES>
face<T, NVERTICES>::~face()
{
    for (auto i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i]->less_sharing_faces();
    }
} //   destructor
