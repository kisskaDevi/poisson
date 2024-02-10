#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <fstream>
#include <filesystem>

template<typename type>
void destroy(type** pointer, size_t n){
    if(pointer){
        for(size_t i = 0; i < n; i++){
            delete[] pointer[i];
        }
        delete[] pointer;
    }
}

template<typename type>
void out_if_file(const std::filesystem::path& file_name, type** data, size_t n_x, size_t n_y){
    if(data){
        std::ofstream file(file_name, std::ofstream::out);
        for(size_t i = 0; i < n_x; i++){
            for(size_t j = 0; j < n_y; j++){
                file << data[i][j] << '\t';
            }
            file << '\n';
        }
    }
#ifdef DEBUG_INFO
    else {
        std::cout << "out_if_file " << file_name << ": data is nullptr" << std::endl;
    }
#endif
}

template<typename type>
type** generate(size_t n_x, size_t n_y){
    type** mem = new type*[n_x];
    for(size_t i = 0; i < n_x; i++){
        mem[i] = new type[n_y];
        for(size_t j = 0; j < n_y; j++){
            mem[i][j] = type();
        }
    }
    return mem;
}

#endif // OPERATIONS_H
