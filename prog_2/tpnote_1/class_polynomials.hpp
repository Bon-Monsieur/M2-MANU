#pragma once
#include <iostream>
#include "class_list.hpp"
#include <cmath>

// nupmber = nb coef = degree + 1

template<typename T>
class polynomial : public list<T>{

public:
    // EXO 1
    polynomial(size_t degree);
    polynomial( size_t degree, T value );

    inline size_t degree() const {return this->number()-1;} 

    T operator()(T x);

    // Surcharge de l’opérateur <<
    friend std::ostream& operator<<(std::ostream& os, const polynomial<T>& p) {
        for (int ii = p.degree(); ii >= 1; ii--) {
            os << p[ii]<< "*x^"<< ii << " + ";
        }
        os << p[0] << std::endl;
        return os;
    }

    // EXO2 
    friend const polynomial<T> operator*(T lambda,  polynomial<T> const& p){
        polynomial<T> temp = p;
        for (size_t ii = 0; ii <= p.degree(); ii++) {
            *(temp.item(ii)) *= lambda;
        }
        return temp;
    }

    friend const polynomial<T> operator*(polynomial<T> const& p,T lambda){
        return lambda * p;
    }
    
    friend const polynomial<T> operator+(polynomial<T> const& p1, polynomial<T> const& p2){
        size_t max_degree = std::max(p1.degree(), p2.degree()); // la somme prend le plus grand degre
        polynomial<T> temp(max_degree);
        for (size_t ii = 0; ii <= max_degree; ii++) {   // on doit check que les coefs existent avant de les sommer
            T coef1 = (ii <= p1.degree()) ? *(p1.item(ii)) : T(0);
            T coef2 = (ii <= p2.degree()) ? *(p2.item(ii)) : T(0);
            *(temp.item(ii)) = coef1 + coef2;
        }
        return temp;
    }


    // EXO 3
    friend const T horner_polynomial(polynomial<T> const& p, T x){
        T res = 0;
        for (int ii = p.degree(); ii >= 0; ii--) {
            res = res * x + *(p.item(ii));
        }
        return res;
    }

};

template<typename T>
polynomial<T>::polynomial(size_t degree) : list<T>(degree + 1){
    for (size_t i = 0; i <= degree; ++i) {
        this->item(i) = new T(0);
    }
}

template<typename T>
polynomial<T>::polynomial(size_t degree, T value) : list<T>(degree + 1){
    for (size_t i = 0; i <= degree; ++i) {
        this->item(i) = new T(value);
    }
}


template<typename T>
T polynomial<T>::operator()(T x) {
    T res = 0;
    for (size_t ii = 0; ii <= degree(); ii++) {
        res += *this->item(ii) * std::pow(x, ii);
    }
    return res;
}
