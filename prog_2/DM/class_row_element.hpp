#pragma once 


template<typename T>
class row_element{
    protected: 
        T value_;
        int column_;

    public:
        // constructors
        row_element(T value=0, int column=-1): value_(value), column_(column){}; 
        row_element<T>& operator=(const row_element<T>& other);

        // getters and setters
        inline T const& get_value() const { return value_; }
        inline int const& get_column() const { return column_;}
        inline void set_value(T value){ value_ = value;}
        inline void set_column(int column){ column_ = column;}

        // basic operators
        row_element<T>& operator+=(T value);
        row_element<T>& operator+=(const row_element<T>& other);
        row_element<T>& operator-=(T value);
        row_element<T>& operator-=(const row_element<T>& other);
        row_element<T>& operator*=(T value);
        row_element<T>& operator*=(const row_element<T>& other);
        row_element<T>& operator/=(T value);
        row_element<T>& operator/=(const row_element<T>& other);

};


// copy constructor 
template<typename T>
row_element<T>& row_element<T>::operator=(const row_element<T>& other){
    if(this != &other){
        this->value_ = other.get_value();
        this->column_ = other.get_column();
    }
    return *this;
}



// basic operators
template<typename T>
row_element<T>& row_element<T>::operator+=(T value){
    this->set_value(this->get_value() + value);
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator+=(const row_element<T>& other){
    this->set_value(this->get_value() + other.get_value());
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator-=(T value){
    this->set_value(this->get_value() - value);
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator-=(const row_element<T>& other){
    this->set_value(this->get_value() - other.get_value());
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator*=(T value){
    this->set_value(this->get_value() * value);
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator*=(const row_element<T>& other){
    this->set_value(this->get_value() * other.get_value());
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator/=(T value){
    if (value == 0){
        throw std::runtime_error("Error: division by zero");
    }
    this->set_value(this->get_value() / value);
    return *this;
}

template<typename T>
row_element<T>& row_element<T>::operator/=(const row_element<T>& other){
    if (other.get_value() == 0){
        throw std::runtime_error("Error: division by zero");
    }
    this->set_value(this->get_value() / other.get_value());
    return *this;
}



// column comparator
template<typename T>
int operator<(const row_element<T>& re1, const row_element<T>& re2){
    return re1.get_column() < re2.get_column();
}

template<typename T>
int operator>(const row_element<T>& re1, const row_element<T>& re2){
    return re1.get_column() > re2.get_column();
}

template<typename T>
int operator==(const row_element<T>& re1, const row_element<T>& re2){
    return re1.get_column() == re2.get_column();
}




// arithmetic operators

template<typename T>
const row_element<T> operator+(row_element<T> const& e, T const& t){
    row_element<T> res = e;
    res += t;
    return res;
}

template<typename T>
const row_element<T> operator+(T const& t , row_element<T> const& e){
    row_element<T> res = e;
    res += t;
    return res;
}

template<typename T>
const row_element<T> operator-(row_element<T> const& e, T const& t){
    row_element<T> res = e;
    res -= t;
    return res;
}

template<typename T>
const row_element<T> operator-(T const& t , row_element<T> const& e){
    row_element<T> res = e;
    res.set_value(t - res.get_value());
    return res;
}

template<typename T>
const row_element<T> operator*(row_element<T> const& e, T const& t){
    row_element<T> res = e;
    res *= t;
    return res;
}

template<typename T>
const row_element<T> operator*(T const& t , row_element<T> const& e){
    row_element<T> res = e;
    res *= t;
    return res;
}


template<typename T>
const row_element<T> operator/(row_element<T> const& e, T const& t){
    row_element<T> res = e;
    if (t == 0){
        throw std::runtime_error("Error: division by zero");
    }
    res /= t;
    return res;
}


template<typename T>
void print(const row_element<T>& re){
    std::cout << "(val: " << re.get_value() << ", col: " << re.get_column() << ") ";
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const row_element<T>& re)
{
    os << "(val: " << re.get_value() << ", col: " << re.get_column() << ") ";
    return os;
}