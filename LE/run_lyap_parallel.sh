#!/bin/bash
#===============================================================================
# run_lyap_parallel.sh
# 
# 编译并运行 OpenMP 并行版本的 Lyapunov 指数计算程序
#
# 使用方法: 
#   chmod +x run_lyap_parallel.sh
#   ./run_lyap_parallel.sh
#===============================================================================

echo "=========================================="
echo "Compiling OpenMP Parallel Lyapunov Code"
echo "=========================================="

# 编译选项说明：
# -O3: 高级优化
# -fopenmp: 启用 OpenMP 支持
# -o lyap: 输出可执行文件名为 lyap

gfortran -O3 -fopenmp -o lyap fieldline_lyap_omp.f90

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Compilation successful!"
echo ""

# 设置 OpenMP 线程数为 8
export OMP_NUM_THREADS=8

# 可选：设置线程亲和性，提高缓存性能
# export OMP_PROC_BIND=close
# export OMP_PLACES=cores

echo "=========================================="
echo "Running Lyapunov Calculation with 8 threads"
echo "=========================================="
echo ""

# 运行程序
time ./lyap

echo ""
echo "=========================================="
echo "Execution completed"
echo "=========================================="
