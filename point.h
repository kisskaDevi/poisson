#ifndef POINT_H
#define POINT_H

#include <ostream>

#include "operations.h"

template<typename type>
struct point{
    type x{0};
    type y{0};

    point() = default;
    point(type x, type y) : x(x), y(y) {}

    point operator+(const point& p){
        return point(x + p.x, y + p.y);
    }

    point operator*(const point<type>& b){
        return point(x * b.x, y * b.y);
    }
};

template<typename type>
std::ostream& operator<<(std::ostream& out, const point<type>& p){
    out << '(' << p.x << ',' << p.y << ')';
    return out;
}

template<typename type>
type dot(const point<type>& a, const point<type>& b){
    return a.x * b.x + a.y * b.y;
}

template<typename type>
struct points{
    point<type>** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    points(size_t n_x, size_t n_y) : n_x(n_x), n_y(n_y) {}

    void generate(type x_0, type x_n, type y_0, type y_n){
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

    point<type>* operator[](const size_t i){
        return data[i];
    }

    const point<type>* operator[](const size_t i) const {
        return data[i];
    }

    ~points(){
        destroy(data, n_x);
    }
};

#endif
