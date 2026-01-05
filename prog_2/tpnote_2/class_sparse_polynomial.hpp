#pragma once
#include "class_monomial.hpp"
#include "class_linked_list_modified.hpp"
#include <cmath>

// sparse_polynomial = linked list de monomes 
class sparse_polynomial :  linked_list_modified<monomial<double>>{
    public:
        sparse_polynomial();
        sparse_polynomial(monomial<double> const& monom , linked_list_modified<monomial<double>>*pt);
        ~sparse_polynomial();
        size_t degree() const;
        void append(monomial<double> const& monom){linked_list_modified<monomial<double>>::append(monom);};
        void insert_first_item(monomial<double> const& monom){linked_list_modified<monomial<double>>::insert_first_item(monom);};
        friend void print(const sparse_polynomial& l);
        friend ostream& operator<<(ostream& os, const sparse_polynomial& sp);
        // La fonction + utilise la fonction += qui ne compile pas, j'ai donc commenté l'appel à cette fonction dans le corps de +
        friend const sparse_polynomial operator+(const sparse_polynomial& p, const sparse_polynomial& q);
        
        // la fonction += ne compile pas, j'ai modifié la fonction dans class_linked_list_modified.hpp mais je n'ai donc pas pu la tester
        //sparse_polynomial& operator+=(sparse_polynomial& L){return linked_list_modified<monomial<double>>::operator+=(L);};

        double operator()(double const& x);
};


sparse_polynomial::sparse_polynomial() : linked_list_modified<monomial<double>>() {};

sparse_polynomial::sparse_polynomial(monomial<double> const& monom , linked_list_modified<monomial<double>>*pt = nullptr){
    this->item() = monom;
    this->p_next() = pt;
}

sparse_polynomial::~sparse_polynomial(){};

size_t sparse_polynomial::degree() const {
    return this->item().power();        // par convention le terme de plus haut degre est au début
}

void print(const sparse_polynomial& sp){
    std::cout << sp.item() << "+";
    if (sp.p_next())
    {
        print(*sp.p_next());
    }
    
}

ostream& operator<<(ostream& os, const sparse_polynomial& sp)
{
    print(sp);
    
    return os;
}


const sparse_polynomial operator+(const sparse_polynomial& p, const sparse_polynomial& q){
    sparse_polynomial temp(p);
    //temp+=q; devrai compiler si l'opérateur += est ok
    return temp;
}


double sparse_polynomial::operator()(double const& x){
    double sum = 0.;
    linked_list_modified<monomial<double>>* scan = this;
    for (; scan; scan = scan->p_next()){
        sum += scan->item().value()*std::pow(x, scan->item().power());
    }
    return sum;
};