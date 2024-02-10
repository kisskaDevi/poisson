TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

CUDA_DIR = $$(CUDA_PATH)
CUDA_INCLUDE_DIR = $$CUDA_DIR/include
CUDA_LIBS_DIR = $$CUDA_DIR/lib/x64

INCLUDEPATH += \
        -I$$PWD/core \
        -I$$PWD/core/cuda \
        -I$$CUDA_INCLUDE_DIR

SOURCES += \
        main.cpp \
        core/field.cpp \
        core/function.cpp \
        core/gauss_seidel.cpp \
        core/point.cpp

HEADERS += \
        core/cuda/cuda_field.h \
        core/cuda/cuda_gauss_seidel.h \
        core/cuda/cuda_point.h \
        core/field.h \
        core/function.h \
        core/gauss_seidel.h \
        core/operations.h \
        core/point.h

DISTFILES += \
    CMakelists.txt \
    plot_script.py

CUDA_SOURCES += \
        core/cuda/cuda_gauss_seidel.cu

CUDA_ARCH = compute_86
CUDA_SM = sm_86
SYSTEM_TYPE = 64
Release:MSVCRT_LINK_FLAG = "/MD"
Debug:MSVCRT_LINK_FLAG = "/MDd"
XFLAGS = -Xcompiler "/EHsc,/W3,/nologo,/O2,/FS,$$MSVCRT_LINK_FLAG"

QMAKE_LIBDIR += -L$$CUDA_LIBS_DIR
LIBS += $$QMAKE_LIBDIR -lcuda -lcudart -lcudadevrt -lcudart_static

Release: config = release
Debug: config = debug

cuda.input = CUDA_SOURCES
cuda.output =  $$OUT_PWD/$$config/cuda_obj/${QMAKE_FILE_BASE}.obj
cuda.commands =         $${CUDA_DIR}/bin/nvcc.exe \
                        -gencode arch=$$CUDA_ARCH,code=$$CUDA_SM \
                        -x cu \
                        -rdc=true \
                        -std=c++17 \
                        $$CUDA_INCLUDEPATH \
                        --machine $$SYSTEM_TYPE \
                        --compile -cudart static \
                        -DNDEBUG -D_CONSOLE -D_UNICODE -DUNICODE \
                        $$XFLAGS \
                        -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
QMAKE_EXTRA_COMPILERS += cuda

CUDA_OBJS += \
        $$OUT_PWD/release/cuda_obj/cuda_gauss_seidel.obj

cuda_l.input =      CUDA_OBJS
cuda_l.output =     $$OUT_PWD/$$config/cuda_obj/device_link.o
cuda_l.commands =   $$CUDA_DIR/bin/nvcc \
                    -gencode arch=$$CUDA_ARCH,code=$$CUDA_SM \
                    $$NVCCFLAGS \
                    -std=c++17 \
                    -dlink \
                    -o ${QMAKE_FILE_OUT} \
                    $$LIBS $$CUDA_OBJS
cuda_l.dependency_type = TYPE_C
cuda_l.clean = $$OUT_PWD/$$config/cuda_obj/device_link.o
QMAKE_EXTRA_UNIX_COMPILERS += cuda_l
