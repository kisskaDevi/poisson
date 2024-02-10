#ifndef CUDA_FIELD_H
#define CUDA_FIELD_H

#include "cuda_point.h"

namespace cuda {

enum boundary_type{
    Non,
    Dirichlet,
    Neumann
};

template<typename type>
struct potential{
    boundary_type b_type{Non};
    type value{type(0)};
    point<type> derivative;

    __device__ bool isDirichlet(){
        return b_type == boundary_type::Dirichlet;
    }
    __device__ bool isNeumann(){
        return b_type == boundary_type::Neumann;
    }
};

}

#endif // CUDA_FIELD_H
