#include "point.h"

#include "operations.h"

#include <fstream>

template<typename type>
point<type>::point(type x, type y) : x(x), y(y) {}

template<typename type>
point<type> point<type>::operator+(const point<type>& p){
    return point(x + p.x, y + p.y);
}

template<typename type>
point<type> point<type>::operator*(const point<type>& b){
    return point(x * b.x, y * b.y);
}

template struct point<float>;
template struct point<double>;

template<typename out_type, typename type>
out_type& operator<<(out_type& out, const point<type>& p){
    out << '(' << p.x << ',' << p.y << ')';
    return out;
}

template<typename type>
type dot(const point<type>& a, const point<type>& b){
    return a.x * b.x + a.y * b.y;
}

template std::ofstream& operator<<(std::ofstream& out, const point<float>& p);
template std::ofstream& operator<<(std::ofstream& out, const point<double>& p);

template float dot(const point<float>& a, const point<float>& b);
template double dot(const point<double>& a, const point<double>& b);

template<typename type>
points<type>::points(size_t n_x, size_t n_y) : n_x(n_x), n_y(n_y) {}

template<typename type>
void points<type>::generate(type x_0, type x_n, type y_0, type y_n){
    point<type> p_0 = point<type>(x_0,y_0);
    type dx = (x_n - x_0) / (n_x - 1);
    type dy = (y_n - y_0) / (n_y - 1);

    data = ::generate<point<type>>(n_x, n_y);
    for(size_t i = 0; i < n_x; i++){
        for(size_t j = 0; j < n_y; j++){
            data[i][j] = p_0 + point<type>(dx * i, dy * j);
        }
    }
}

template<typename type>
point<type>* points<type>::operator[](const size_t i){
    return data[i];
}

template<typename type>
const point<type>* points<type>::operator[](const size_t i) const {
    return data[i];
}

template<typename type>
points<type>::~points(){
    destroy(data, n_x);
}

template struct points<float>;
template struct points<double>;
