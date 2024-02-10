#include "function.h"

#include "operations.h"

template<typename type>
function<type>::~function(){
    destroy(data, n_x);
}

template<typename type>
void function<type>::generate(const points<type>& grid, type (*f)(type x, type y)){
    n_x = grid.n_x;
    n_y = grid.n_y;
    data = ::generate<type>(n_x, n_y);
    for(size_t i = 0; i < n_x; i++){
        for(size_t j = 0; j < n_y; j++){
            data[i][j] = f(grid.data[i][j].x, grid.data[i][j].y);
        }
    }
}

template<typename type>
type* function<type>::operator[](const size_t i){
    return data[i];
}

template<typename type>
const type* function<type>::operator[](const size_t i) const {
    return data[i];
}

template struct function<float>;
template struct function<double>;
