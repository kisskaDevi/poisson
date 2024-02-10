#include "cuda_gauss_seidel.h"

#include "cuda_field.h"

#include <iostream>

#define DEBUG_INFO

#define checkCudaErrors(val) check_cuda( (val), #val, __FILE__, __LINE__)

void check_cuda(cudaError_t result, char const* const func, const char* const file, int const line) {
    if (result) {
        std::cerr << "CUDA error = " << static_cast<unsigned int>(result) << " at " << file << ":" << line << " '" << func << "' \n";
        cudaDeviceReset();
        exit(99);
    }
}

namespace cuda {

template<typename type>
__device__ type get_u(cuda::point<type> delta, const cuda::potential<type>& u_field){
    switch (u_field.b_type) {
    case boundary_type::Non:
    case boundary_type::Dirichlet:
        break;
    case boundary_type::Neumann:
        return dot(delta, u_field.derivative);
    }
    return u_field.value;
}

template<typename type>
type* alloc(const type& value){
    type* gpu_p = nullptr;
    checkCudaErrors(cudaMalloc((void**)&gpu_p, sizeof(type)));
    checkCudaErrors(cudaMemcpy((void*)gpu_p, (void*)&value, sizeof(type), cudaMemcpyHostToDevice));
    return gpu_p;
}

template<typename type_gpu, typename type_cpu>
type_gpu* alloc_data(type_cpu** data, size_t n_x, size_t n_y){
    type_gpu* grid_data = nullptr;
    checkCudaErrors(cudaMalloc((void**)&grid_data, n_x * n_y * sizeof(type_gpu)));

    type_cpu* cpu_data =  new type_cpu[n_x * n_y];
    for(size_t i = 0; i < n_x; i++){
        for(size_t j = 0; j < n_y; j++){
            cpu_data[i * n_y + j] = data[i][j];
        }
    }
    checkCudaErrors(cudaMemcpy((void*)grid_data, (void*)cpu_data, n_x * n_y * sizeof(type_gpu), cudaMemcpyHostToDevice));

    return grid_data;
}

template<typename type>
__global__ void calculate_wave(type* dmax, size_t p, size_t q_0, size_t q_n, size_t n_y, const cuda::point<type>* points, const type* func, cuda::potential<type>* u_field){
    int q = q_0 + threadIdx.x + blockIdx.x * blockDim.x;
    if(q > q_n){
        return;
    }

    size_t i = 1 + p - q, j = 1 + q;

    if(u_field[i * n_y + j].b_type == boundary_type::Non){
        type dx = type(0.5) * (points[(i + 1) * n_y + j].x - points[(i - 1) * n_y + j].x);
        type dy = type(0.5) * (points[i * n_y + (j + 1)].y - points[i * n_y + (j - 1)].y);
        type a = type(1) / dx / dx, b = type(1) / dy / dy, c = type(2) * (a + b);

        type temp = u_field[i * n_y + j].value;
        u_field[i * n_y + j].value = (
                                  a * (get_u({ type(2) * dx, 0}, u_field[(i + 1) * n_y + j]) + (u_field[(i + 1) * n_y + j].isNeumann() ? u_field[(i - 1) * n_y + j].value : type(0))) +
                                  a * (get_u({-type(2) * dx, 0}, u_field[(i - 1) * n_y + j]) + (u_field[(i - 1) * n_y + j].isNeumann() ? u_field[(i + 1) * n_y + j].value : type(0))) +
                                  b * (get_u({0,  type(2) * dy}, u_field[i * n_y + (j + 1)]) + (u_field[i * n_y + (j + 1)].isNeumann() ? u_field[i * n_y + (j - 1)].value : type(0))) +
                                  b * (get_u({0, -type(2) * dy}, u_field[i * n_y + (j - 1)]) + (u_field[i * n_y + (j - 1)].isNeumann() ? u_field[i * n_y + (j + 1)].value : type(0))) -
                                  func[i * n_y + j]) / c;

        if(*dmax < std::abs(temp - u_field[i * n_y + j].value)){
            *dmax = std::abs(temp - u_field[i * n_y + j].value);
        }
    }
}

template<typename type>
void gauss_seidel(const gauss_seidel_info<type>& info, const ::points<type>& points, const ::function<type>& func, ::field<type>& u_field){
    cuda::point<type>* grid_data = alloc_data<cuda::point<type>>(points.data, points.n_x, points.n_y);
    type* func_data = alloc_data<type>(func.data, func.n_x, func.n_y);
    cuda::potential<type>* field_data = alloc_data<cuda::potential<type>>(u_field.data, u_field.n_x, u_field.n_y);

    bool convergence = false;
    size_t it = 0;
    for (; !convergence && it < info.max_it; it++ ){
        type* dmax = alloc(type(0));
        for (size_t p = 0; p <= points.n_x - 3 + points.n_y - 3; p++ ){
            calculate_wave<<<10, std::max(points.n_x, points.n_y) / 10 + 1>>>(
                dmax,
                p,
                (p <= points.n_x - 3) ? 0 : p - (points.n_x - 3),
                std::min(p, std::min(points.n_x - 3, points.n_y - 3)),
                points.n_y,
                grid_data,
                func_data,
                field_data);
        }
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
        type dmax_cpu = info.eps;
        checkCudaErrors(cudaMemcpy((void*)&dmax_cpu, (void*)dmax, sizeof(type), cudaMemcpyDeviceToHost));
        convergence = dmax_cpu < info.eps;
    }
#ifdef DEBUG_INFO
    std::cout << "   convergence = " << (convergence ? "true" : "false") << std::endl;
    std::cout << "   last_it = " << it << std::endl;
#endif

    ::potential<type>* cpu_data = new ::potential<type>[points.n_x * points.n_y];
    checkCudaErrors(cudaMemcpy((void*)cpu_data, (void*)field_data, points.n_x * points.n_y * sizeof(cuda::potential<type>), cudaMemcpyDeviceToHost));
    for(size_t i = 0; i < points.n_x; i++){
        for(size_t j = 0; j < points.n_y; j++){
            u_field[i][j] = cpu_data[i * points.n_y + j];
        }
    }
}

template void gauss_seidel<float>(const gauss_seidel_info<float>& info, const ::points<float>& points, const ::function<float>& func, ::field<float>& u_field);
template void gauss_seidel<double>(const gauss_seidel_info<double>& info, const ::points<double>& points, const ::function<double>& func, ::field<double>& u_field);
}
