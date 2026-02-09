#pragma once
#include <cmath>
#include <iostream>

#include "class_list.hpp"

using namespace std;

template<typename T>
class polynomial : public list<T>
{
public:
    //  constructor
    polynomial(size_t degree) : list<T>(degree + 1) {}

    //  constructor with a double argument
    polynomial(size_t degree, T value) : list<T>(degree + 1, value) {}

    ~polynomial() {} //  destructor

    size_t degree() const { return this->number_ - 1; } //  degree of polynomial

    template <typename S>
    friend S calculate_polynomial(polynomial<S> const& p, S x);

    template <typename S>
    friend S horner_polynomial(polynomial<S> const& p, S x);

    T operator()(T xx)
    {
        return calculate_polynomial(*this, xx);
    };

    template <typename S>
    friend ostream& operator<<(ostream& inout, polynomial<S> const& polyn);
};

template <typename S>
ostream& operator<<(ostream& inout, polynomial<S> const& polyn)
{
    for (auto ii = polyn.number()-1; ii > 0; --ii)
    {

        inout << *polyn.item_[ii] << "x^" << ii << " + ";
    }
    inout << *polyn.item_[0] << endl;
    return inout;
}

//  add two polynomials
template <typename S>
const polynomial<S> operator+(const polynomial<S>& p, const polynomial<S>& q)
{
    if (p.number() >= q.number())
    {
        polynomial<S> sump(p);
        for (auto ii = 0; ii < q.number(); ++ii)
        {
            sump.list::operator()(ii) += q[ii];
        }
        return sump;
    }
    else
    {
        polynomial<S> sump(q);
        for (auto ii = 0; ii < p.number(); ++ii)
        {
            sump.list::operator()(ii) += p[ii];
        }
        return sump;
    }
}

//  scalar times polynomial
template <typename S>
const polynomial<S> operator*(S a, polynomial<S> const& p)
{
    polynomial<S> tmp = p;
    for (auto ii = 0; ii < p.number(); ii++)
    {
        tmp.list::operator()(ii) *= a;
    }
    return tmp;
}

//  scalar times polynomial
template <typename S>
const polynomial<S> operator*(polynomial<S> const& p, S a)
{
    polynomial<S> tmp = p;
    for (auto ii = 0; ii < p.number(); ii++)
    {
        tmp.list::operator()(ii) *= a;
    }
    return tmp;
}

//  calculate a polynomial
template <typename S>
S calculate_polynomial(polynomial<S> const& p, S x)
{
    S power_of_X = 1;
    S sum      = 0.;

    for (auto i = 0; i < p.number(); i++)
    {
        sum += p[i] * power_of_X;
        power_of_X *= x;
    }
    return sum;
}

template <typename S>
S horner_polynomial(polynomial<S> const& p, S x)
{
    S result = p[p.degree()];
    for(auto i=p.degree(); i>0; i--)
    {
      result *= x;
      result += p[i-1];
    }
    return result;
}  //  Horner algorithm to calculate a polynomial
