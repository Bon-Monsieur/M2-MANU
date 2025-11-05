#pragma once
#include "class_dynamic_vector.hpp"
#include "class_list.hpp"


template<typename T>
class binomial : public list<DynamicVector<T>> {
    public: 
        binomial(int M);
        T operator()(int n, int k) {
            if (k < 0 || k > n) {
                throw std::out_of_range("Index k out of range");
            }
            return (*(this->item_(n)))[k];
        }
};


template<typename T>
binomial<T>::binomial(int M) : list<DynamicVector<T>>(M+1) {
    this->item(0) = new DynamicVector<T>(1,T(1));
    this->item(1) = new DynamicVector<T>(2,T(1));
    for (int mm = 2; mm <= M; mm++){
        DynamicVector<T> temp(mm+1,T(1));
        for (int kk=1; kk < mm; kk++){
            temp[kk] = (*(this->item(mm-1)))[kk-1] + (*(this->item(mm-1)))[kk];
        }
        this->item(mm) = new DynamicVector<T>(temp);
    }
}