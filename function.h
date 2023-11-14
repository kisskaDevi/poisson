#ifndef FUNCTION_H
#define FUNCTION_H

#include "point.h"

template<typename type>
struct function
{
    type** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    ~function(){
        destroy(data, n_x);
    }

    void generate(const points<type>& grid, type (*f)(type x, type y)){
        n_x = grid.n_x;
        n_y = grid.n_y;
        data = ::generate<type>(n_x, n_y);
        for(size_t i = 0; i < n_x; i++){
            for(size_t j = 0; j < n_y; j++){
                data[i][j] = f(grid.data[i][j].x, grid.data[i][j].y);
            }
        }
    }

    type* operator[](const size_t i){
        return data[i];
    }

    const type* operator[](const size_t i) const {
        return data[i];
    }
};

#endif // FUNCTION_H
