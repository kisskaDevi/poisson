#ifndef CUDA_GAUSS_SEIDEL_H
#define CUDA_GAUSS_SEIDEL_H

#include "../field.h"
#include "../function.h"

namespace cuda {

template<typename type>
struct gauss_seidel_info{
    type eps{0};
    size_t max_it{0};
};

template<typename type>
void gauss_seidel(const gauss_seidel_info<type>& info, const ::points<type>& points, const ::function<type>& func, ::field<type>& u_field);
}

#endif // CUDA_GAUSS_SEIDEL_H
