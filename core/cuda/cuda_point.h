#ifndef CUDA_POINT_H
#define CUDA_POINT_H

namespace cuda {

template<typename type>
struct point{
    type x{0};
    type y{0};

    __device__ point() {}
    __device__ point(type x, type y): x(x), y(y) {}
    __device__ point operator+(const point& p){
        return point(x + p.x, y + p.y);
    }
    __device__ point operator*(const point<type>& b){
        return point(x * b.x, y * b.y);
    }
};

template<typename type>
__device__ type dot(const point<type>& a, const point<type>& b){
    return a.x * b.x + a.y * b.y;
}

}

#endif // CUDA_POINT_H
