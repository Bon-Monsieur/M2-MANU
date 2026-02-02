#pragma once

#include "class_node.hpp"
#include "class_point.hpp"
#include <cmath>

#include <iostream>

template <typename T, int NVERTICES>
class element
{
protected:
    node<T>* vertices_[NVERTICES];
    int index_;
    // geometry (element dependent)
    double area_; // value of the triangle area
    point2d centroid_; // coordinates of the centroid
    std::array<point2d, NVERTICES> normals_; // coordinates of the outgoing unit normals to faces
    std::array<double, NVERTICES> faces_length_; // faces lengths

    std::array<int, NVERTICES> neighbors_; // indices globaux parmis les {triangles} des adjacents
    std::array<int, NVERTICES> faces_; // indices globaux des edges de l'élément 
    
public:
    element();                             //  default constructor
    element(node<T>&, node<T>&, node<T>&); // for triangles
    element(element<T, NVERTICES> const&);
    
    node<T>& operator()(size_t i) { return *(this->vertices_[i]); }             //  read/write ith vertex
    const node<T>& operator[](size_t i) const { return *(this->vertices_[i]); } //  read only ith vertex

    node<T>* vertex(std::size_t i) const { return this->vertices_[i]; }; //  get ith vertex address (pointer)
    node<T>*& vertex(std::size_t i) { return this->vertices_[i]; };      //  get ith vertex address (pointer)

    inline double area() const { return area_; };
    inline void set_area(double a) { area_ = a; };
    inline point2d centroid() {return centroid_;};
    inline std::array<double, NVERTICES> faces_length() { return faces_length_;};
    inline point2d normal(size_t ii) const { return normals_[ii];};
    inline int neighbor(std::size_t i) const {return this->neighbors_[i];};
    inline void set_neighbor_index(std::size_t ii, int jj) {this->neightbors_[ii] = jj;}
    inline int face(std::size_t i) const {return this->faces_[i];};
    inline int& face(std::size_t i){return this->faces_[i];};

    const element<T, NVERTICES>& operator=(element<T, NVERTICES>&);
    ~element();

    void reset_indices();               //  reset vertices indices to -1
    void indexing_vertices(int& count); //  indexing the vertices

    inline int index() const { return index_;};
    inline void set_index(int ind){ index_ = ind;};

    template <typename S, int MVERTICES>
    friend std::ostream& operator<<(std::ostream& os, element<S, MVERTICES> const& e);
    void print_vertices_indices();


    void compute_area();
    void compute_centroid();
    void compute_normals();
    void compute_faces_length();
};

typedef element<point2d, 3> triangle;

template <typename T, int NVERTICES>
void element<T, NVERTICES>::compute_faces_length(){
    double dx1 = vertices_[1]->x() - vertices_[0]->x();
    double dy1 = vertices_[1]->y() - vertices_[0]->y();
    double dx2 = vertices_[2]->x() - vertices_[1]->x();
    double dy2 = vertices_[2]->y() - vertices_[1]->y();
    double dx3 = vertices_[0]->x() - vertices_[2]->x();
    double dy3 = vertices_[0]->y() - vertices_[2]->y();

    faces_length_[0] = std::sqrt(dx1*dx1 + dy1*dy1);
    faces_length_[1] = std::sqrt(dx2*dx2 + dy2*dy2);
    faces_length_[2] = std::sqrt(dx3*dx3 + dy3*dy3);
} // compute_faces_length

template <typename T, int NVERTICES>
void element<T, NVERTICES>::compute_normals(){
    for (auto i = 0; i <NVERTICES; i++){
        //scan vers les sommets de l'arête
        node<T>* pA = vertices_[i];
        node<T>* pB = vertices_[(i+1)%3];

        // 1. Vecteur directeur de l'arête
        double dx = pB->x() - pA->x();
        double dy = pB->y() - pA->y();

        //2. Normale sortante (orthogonale à l'arête)
        // Rotation de 90° dans le sens horaire
        double nx = dy;
        double ny = -dx;

        // 3. Calcul de la longueur (mesure de la face)
        double length = std::sqrt(nx*nx + ny*ny);
        faces_length_[i] = length;

        // 4. Normalisation unitaire
        if (length > 1e-15)
        {
            double invL = 1.0 / length;
            normals_[i] = point2d(nx * invL, ny * invL);
        }
    }
} // compute_normals

template <typename T, int NVERTICES>
void element<T, NVERTICES>::compute_centroid()
{
    double x_sum = 0.;
    double y_sum = 0.;

    for (auto ii = 0; ii < NVERTICES; ii++)
    {
        x_sum += this->vertices_[ii]->x();
        y_sum += this->vertices_[ii]->y();
    }

    centroid_.x() = x_sum / NVERTICES;
    centroid_.y() = y_sum / NVERTICES;
} // compute_centroid

template <typename T, int NVERTICES>
void element<T, NVERTICES>::compute_area() // area for triangles
{
    double x0 = this->vertices_[0]->x();
    double y0 = this->vertices_[0]->y();
    double x1 = this->vertices_[1]->x();
    double y1 = this->vertices_[1]->y();
    double x2 = this->vertices_[2]->x();
    double y2 = this->vertices_[2]->y();

    area_ = 0.5 * std::abs(x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
} // compute_area


template <typename T, int NVERTICES>
std::ostream& operator<<(std::ostream& os, element<T, NVERTICES> const& e)
{
    for (auto ii = 0; ii < NVERTICES; ii++)
        os << e[ii];
    return os;
} //

template <typename T, int NVERTICES>
void element<T, NVERTICES>::print_vertices_indices()
{
    for (auto ii = 0; ii < NVERTICES; ii++)
        std::cout << "this->vertices_[i]->index()= " << this->vertices_[ii]->index() << std::endl;
}

template <typename T, int NVERTICES>
element<T, NVERTICES>::element(node<T>& a, node<T>& b, node<T>& c)
{
    this->vertices_[0] = &a;
    this->vertices_[1] = &b;
    this->vertices_[2] = &c;

    for (auto i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i]->more_sharing_elements();
    }
} //  constructor

template <typename T, int NVERTICES>
element<T, NVERTICES>::element(element<T, NVERTICES> const& e)
{
    for (int i = 0; i < NVERTICES; i++)
    {
        this->vertices_[i] = e.vertices_[i];
        this->vertices_[i]->more_sharing_elements();
    }
} //  copy constructor

template <typename T, int NVERTICES>
const element<T, NVERTICES>& element<T, NVERTICES>::operator=(element<T, NVERTICES>& e)
{
    if (this != &e)
    {
        for (int i = 0; i < NVERTICES; i++)
            this->vertices_[i]->less_sharing_elements();
        // delete this->vertices_[i];

        for (int i = 0; i < NVERTICES; i++)
        {
            this->vertices_[i] = e.vertices_[i];
            this->vertices_[i]->more_sharing_elements();
        }
        this->set_index(e.index());
    }
    return *this;
} //  assignment operator

template <typename T, int NVERTICES>
element<T, NVERTICES>::~element()
{
    for (int i = 0; i < NVERTICES; i++)
        this->vertices_[i]->less_sharing_elements();
} //   destructor

template <typename T, int NVERTICES>
void element<T, NVERTICES>::reset_indices()
{
    for (int i = 0; i < NVERTICES; i++)
        this->vertices_[i]->index() = -1;
} //  reset indices to -1

template <typename T, int NVERTICES>
void element<T, NVERTICES>::indexing_vertices(int& count)
{
    for (int i = 0; i < NVERTICES; i++)
    {
        if (this->vertices_[i]->index() < 0)
        {
            this->vertices_[i]->index() = count++;
        }
    }
} //  indexing the vertices

template <typename T, int NVERTICES>
int operator<(const node<T>& n, const element<T, NVERTICES>& e)
{
    for (int i = 0; i < NVERTICES; i++)
        if (&n == &(e[i]))
            return i + 1;

    return 0;
} //  check whether a node n is in a finite element e
