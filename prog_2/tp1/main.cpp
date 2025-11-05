#include "class_list.hpp"
#include "class_dynamic_vector.hpp"
#include "class_binomial.hpp"

template<typename T>
T inverse(T x){
    return T(1)/x;
}

template<typename T>
const list<T> derive_inverse(T const& a, int n){
    list<T> temp(n+1,0);
    temp(0) = inverse(a);
    for (int k = 1; k < n+1; k++){
        temp(k) = -k/a*temp(k-1);
    }
    return temp;
}


template<class T>
const T Taylor ( list<T> const& f , T const& h ){
    T sum = f[0];
    T h_temp = 1;
    T k_fact = 1;
    for (int k=1; k < f.number(); k++){
        sum += f[k]*(h/k)*(h_temp/k_fact);
        h_temp *= h;
        k_fact *= k;
    }
    return sum;
}



int main(){
    
    /*
    // Essaie liste des dérivées de 1/x en 2 jusqu'à l'ordre 5
    double a = 2.0;
    int n = 5;
    list<double> my_list(derive_inverse(a,n));
    std::cout << my_list;

    // Essaie de l'approximation de 1/(2.1) par le polynôme de Taylor d'ordre 5
    double h = 0.1;
    double approx = Taylor(my_list,h);
    std::cout << "Approximation de 1/" << a+h << " : " << approx << std::endl;
    */

    // Essaie de la classe binomiale
    int M = 10;
    binomial<double> my_binom(M);
    std::cout << my_binom;
    std::cout << "Binomial coefficient C(10,3) = " << my_binom(10,3) << std::endl;

    return 0;
}