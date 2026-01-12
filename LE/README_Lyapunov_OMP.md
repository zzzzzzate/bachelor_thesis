# Lyapunov 指数计算程序 - OpenMP 并行版本

## 文件说明

- `fieldline_lyap_omp.f90` - OpenMP 并行化的 Fortran 源代码
- `run_lyap_parallel.sh` - 编译和运行脚本（Linux）

## 并行化特点

### 主要改进

1. **OpenMP 并行化**：使用 `!$OMP PARALLEL DO` 指令并行化计算不同初始点的 Lyapunov 指数
2. **动态调度**：使用 `SCHEDULE(dynamic, 1)` 实现负载均衡
3. **线程安全输出**：使用 `!$OMP CRITICAL` 保护打印语句
4. **性能统计**：输出总计算时间和加速比

### 并行化策略

```fortran
!$OMP PARALLEL DO PRIVATE(i_init, ix, iy, R_init, Z_init, phi_init, lambda, thread_id) &
!$OMP SHARED(n_init_points, nx, R_grid, Z_grid, rho_grid, psi_grid, lyap_results) &
!$OMP SHARED(n_turns, dphi_step, n_renorm, delta0) &
!$OMP SCHEDULE(dynamic, 1)
```

- **PRIVATE 变量**：每个线程独立的局部变量
- **SHARED 变量**：所有线程共享的数据（只读）
- **动态调度**：避免负载不均衡（因为不同初始点的计算时间可能不同）

## 编译和运行

### 方法 1：使用脚本（推荐）

```bash
# 在 Linux 服务器上
chmod +x run_lyap_parallel.sh
./run_lyap_parallel.sh
```

### 方法 2：手动编译和运行

```bash
# 编译
gfortran -O3 -fopenmp -o lyap fieldline_lyap_omp.f90

# 设置线程数（使用 8 核）
export OMP_NUM_THREADS=8

# 运行
./lyap
```

### 调整线程数

如果服务器有不同的核心数，可以调整：

```bash
export OMP_NUM_THREADS=16  # 使用 16 核
./lyap
```

或者在运行时指定：

```bash
OMP_NUM_THREADS=4 ./lyap  # 使用 4 核
```

## 性能优化选项

### 编译器优化标志

```bash
# 基础优化（已包含在脚本中）
gfortran -O3 -fopenmp -o lyap fieldline_lyap_omp.f90

# 更激进的优化（可能带来数值误差）
gfortran -O3 -fopenmp -march=native -ffast-math -o lyap fieldline_lyap_omp.f90

# 启用向量化报告（调试用）
gfortran -O3 -fopenmp -fopt-info-vec -o lyap fieldline_lyap_omp.f90
```

### OpenMP 环境变量

```bash
# 设置线程数
export OMP_NUM_THREADS=8

# 线程绑定到特定核心（提高缓存性能）
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# 设置调度类型（如果需要覆盖代码中的 SCHEDULE）
export OMP_SCHEDULE="dynamic,1"
```

## 预期性能提升

理论上：
- **单核**：计算 50 个点约需 T 秒
- **8 核**：约需 T/6 ~ T/7 秒（考虑并行开销和负载均衡）
- **加速比**：约 6-7 倍

实际性能取决于：
1. CPU 架构和缓存大小
2. 内存带宽
3. 磁场数据的访问模式
4. 不同初始点的计算时间差异

## 输出说明

程序会输出：

```
==========================================
 Lyapunov Exponent Calculator for
 Tokamak Magnetic Field Lines
 (OpenMP Parallel Version)
==========================================

 Using  8 OpenMP threads

Grid loaded: nx= ...

Computing Lyapunov exponents...

  [Thread  3] Point   5/ 50  rho=  0.0500  lambda= 1.23E-02
  [Thread  1] Point  10/ 50  rho=  0.1000  lambda= 2.45E-02
  ...

Results saved to data_mat/lyapunov_results.dat
Total time:      45.23 seconds
Average time per point:      0.9046 seconds
Speedup factor (approx):       6.40
Done!
```

## 注意事项

1. **内存需求**：所有线程共享磁场数据，确保服务器有足够内存
2. **文件 I/O**：磁场数据读取是串行的，这部分不会被加速
3. **数值稳定性**：OpenMP 并行不应影响计算结果，但浮点运算顺序可能略有不同
4. **线程数选择**：通常设置为 CPU 物理核心数，不建议超过核心数

## 故障排查

### 编译错误

```bash
# 检查 gfortran 版本（需要支持 OpenMP）
gfortran --version

# 如果 gfortran 不支持 OpenMP，尝试
gcc --version  # 需要 >= 4.2
```

### 运行时错误

```bash
# 检查数据文件是否存在
ls data_mat/*.dat

# 查看 OpenMP 设置
echo $OMP_NUM_THREADS
```

### 性能问题

```bash
# 使用 time 命令查看详细时间
time ./lyap

# 查看 CPU 使用情况（另一个终端）
top -H -p $(pgrep lyap)
```

## 与串行版本对比

| 特性 | 串行版本 | OpenMP 并行版本 |
|------|---------|----------------|
| 计算时间 | T | T/6 ~ T/7 |
| 内存使用 | M | M (共享数据) |
| CPU 核心 | 1 | 8 |
| 代码复杂度 | 简单 | 中等 |
| 可扩展性 | 无 | 线性（至核心数）|

## 进一步优化建议

如果需要更高性能：

1. **MPI + OpenMP 混合**：跨节点并行
2. **GPU 加速**：使用 CUDA 或 OpenACC
3. **算法优化**：
   - 减少插值次数
   - 使用更高效的网格搜索算法
   - 优化 RK4 积分器

## 相关文件

- `main.f90` - 主程序（生成磁场数据）
- `data_mat/` - 输入/输出数据目录
- `lyapunov_results.dat` - Lyapunov 指数结果
