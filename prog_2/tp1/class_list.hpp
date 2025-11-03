#include <iostream>
#include <cassert>
#include <cmath>
#pragma once


template<typename T>
class list {
protected:
    int number_;
    T** item_;

public:
    list(int number);
    list(int number, T tab);
    ~list() {
        for (int ii = 0; ii < number_; ii++)
            delete item_[ii];
        delete[] item_;
    }
    list ( list <T> const&); // copy constructor

    list<T>& operator=(const list<T>& l1);

    const T& operator[](int i) const;   // read-only 
    T& operator()(int i);               // read-write

    int number() const { return number_; }
    int& number() { return number_; }

    T* item(int i) const { return item_[i]; }   // read-only
    T*& item(int i) { return item_[i]; }        // read-write


    // Surcharge de l’opérateur <<
    friend std::ostream& operator<<(std::ostream& os, const list<T>& list) {
        os << "List with " << list.number_ << " items:\n";
        for (int ii = 0; ii < list.number_; ii++) {
            os << "  item_[" << ii << "] address = " << list[ii] << "\n";
        }
        return os;
    }

};

template<typename T>
list<T>::list(int number){
    if (number>0){
        number_ = number;
        item_ = new T*[number_];

    }
    else{
        number_ = 0;
        item_ = nullptr;
    }
}

template<typename T>
list<T>::list(int number, T tab){
    number_ = number;
    item_ = new T*[number_];
    for (int ii = 0; ii < number_; ii++){
        item_[ii] = new T(tab);   // Créer une copie de tab pour chaque élément afin de ne pas partager 
                                //la même référence mémoire 
    }
}

template<typename T>
list<T>::list ( list <T> const& tab) : number_(tab.number_){
    item_ = new T*[number_];
    for (int ii = 0; ii < number_; ii++){
        item_[ii] = new T(*(tab.item_[ii])); // Créer une copie de chaque élément
    }
}


template<typename T>
list<T>& list<T>::operator=(const list<T>& l1) {
    if (this == &l1)  // auto-assignment protection
        return *this;

    // Libérer l'ancienne mémoire
    for (int ii = 0; ii < number_; ii++)
        delete item_[ii];
    delete[] item_;

    // Copier les données
    number_ = l1.number_;
    item_ = new T*[number_];
    for (int ii = 0; ii < number_; ii++)
        item_[ii] = new T(*(l1.item_[ii]));  // copie profonde

    return *this;  
}

template<typename T>
const T& list<T>::operator[](int i) const {
    assert(i >= 0 && i < number_ && "Index out of bounds");
    return *(item_[i]);
}

template<typename T>
T& list<T>::operator()(int i) {
    assert(i >= 0 && i < number_ && "Index out of bounds");
    return *(item_[i]);
}