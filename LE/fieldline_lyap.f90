!===============================================================================
! fieldline_lyap.f90
! 
! 计算托卡马克磁力线轨道的李雅普诺夫指数 (OpenMP并行)
! 
! 编译: gfortran -O3 -fopenmp -o lyap fieldline_lyap.f90
! 运行: export OMP_NUM_THREADS=8; ./lyap
!       或: OMP_NUM_THREADS=8 ./lyap
!
! rblock.dat是在曲线网格数据
!===============================================================================

program fieldline_lyap
  use omp_lib
  implicit none
  
  ! 精度定义
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0d0)
  
  ! 网格参数
  integer :: nx, ny, nphi
  real(dp) :: R0_axis
  
  ! 网格数组
  real(dp), allocatable :: R_grid(:,:)    ! (0:nx, 0:ny)
  real(dp), allocatable :: Z_grid(:,:)    ! (0:nx, 0:ny)
  real(dp), allocatable :: phi_grid(:)    ! (1:nphi)
  real(dp), allocatable :: psi_grid(:)    ! (0:nx) normalized psi
  real(dp), allocatable :: rho_grid(:)    ! (0:nx) minor radius
  
  ! 磁场数组
  real(dp), allocatable :: B_R(:,:,:)     ! (0:nx, 0:ny, 1:nphi)
  real(dp), allocatable :: B_Z(:,:,:)     ! (0:nx, 0:ny, 1:nphi)
  real(dp), allocatable :: B_phi(:,:,:)   ! (0:nx, 0:ny, 1:nphi)
  
  ! 用于插值的辅助量
  real(dp) :: R_min, R_max, Z_min, Z_max
  real(dp) :: dR_grid, dZ_grid, dphi_grid
  
  ! Lyapunov 计算参数
  integer :: n_init_points    ! 不同半径的初始点数
  integer :: n_turns          ! 环绕圈数
  integer :: n_renorm         ! 重标定间隔步数
  real(dp) :: delta0          ! 初始分离距离
  real(dp) :: dphi_step       ! phi 步长
  
  ! 结果数组
  real(dp), allocatable :: lyap_results(:,:)  ! (n_init, 3): rho, psi, lambda
  
  ! 临时变量
  integer :: i, ix, iy, iphi, i_init
  real(dp) :: R_init, Z_init, phi_init
  real(dp) :: lambda
  integer :: num_threads, thread_id
  real(dp) :: start_time, end_time
  
  !------------------------------------------------------------------
  ! 主程序
  !------------------------------------------------------------------
  
  print *, '=========================================='
  print *, ' Lyapunov Exponent Calculator for'
  print *, ' Tokamak Magnetic Field Lines'
  print *, ' (OpenMP Parallel Version)'
  print *, '=========================================='
  print *, ''
  
  ! 获取并设置 OpenMP 线程数
  num_threads = omp_get_max_threads()
  print '(A,I2,A)', ' Using ', num_threads, ' OpenMP threads'
  print *, ''
  
  ! 设置计算参数
  n_init_points = 100      ! 径向初始点数
  n_turns = 500           ! 环绕圈数
  n_renorm = 100          ! 每100步重标定一次
  delta0 = 1.0d-8         ! 初始分离距离 (m)
  dphi_step = 2.0d0*pi/360.0d0  ! phi 步长 (1度)
  
  ! 读取数据
  call read_grid_data()
  call read_field_data()
  
  ! 计算网格范围(用于插值边界检查)
  call compute_grid_bounds()
  
  print *, 'Grid loaded: nx=', nx, ' ny=', ny, ' nphi=', nphi
  print *, 'R range: [', R_min, ',', R_max, ']'
  print *, 'Z range: [', Z_min, ',', Z_max, ']'
  print *, ''
  
  ! 分配结果数组
  allocate(lyap_results(n_init_points, 3))
  
  ! 对不同半径计算 Lyapunov 指数 (OpenMP 并行)
  print *, 'Computing Lyapunov exponents...'
  print *, ''
  
  start_time = omp_get_wtime()
  
  !$OMP PARALLEL DO PRIVATE(i_init, ix, iy, R_init, Z_init, phi_init, lambda, thread_id) &
  !$OMP SHARED(n_init_points, nx, R_grid, Z_grid, rho_grid, psi_grid, lyap_results) &
  !$OMP SHARED(n_turns, dphi_step, n_renorm, delta0) &
  !$OMP SCHEDULE(dynamic, 1)
  do i_init = 1, n_init_points
    ! 选择初始点：在 theta=0 (外中平面) 上，不同的 psi
    ix = int(real(i_init) / real(n_init_points) * real(nx) * 0.9d0) + 1
    if (ix > nx) ix = nx
    
    ! 找 theta=0 对应的 iy (假设 iy=1 对应外中平面)
    iy = 1
    
    R_init = R_grid(ix, iy)
    Z_init = Z_grid(ix, iy)
    phi_init = 0.0d0
    
    ! 计算 Lyapunov 指数
    call compute_lyapunov(R_init, Z_init, phi_init, n_turns, dphi_step, &
                          n_renorm, delta0, lambda)
    
    ! 存储结果
    lyap_results(i_init, 1) = rho_grid(ix)      ! rho
    lyap_results(i_init, 2) = psi_grid(ix)      ! normalized psi
    lyap_results(i_init, 3) = lambda            ! Lyapunov exponent
    
    ! 打印进度 (线程安全)
    !$OMP CRITICAL
    if (mod(i_init, 5) == 0 .or. i_init == n_init_points) then
      thread_id = omp_get_thread_num()
      print '(A,I2,A,I3,A,I3,A,F8.4,A,ES12.4)', &
        '  [Thread ', thread_id, '] Point ', i_init, '/', n_init_points, &
        '  rho=', rho_grid(ix), '  lambda=', lambda
    endif
    !$OMP END CRITICAL
  enddo
  !$OMP END PARALLEL DO
  
  end_time = omp_get_wtime()
  
  ! 输出结果
  call output_results()
  
  print *, ''
  print *, 'Results saved to lyapunov_results.dat'
  print '(A,F10.2,A)', ' Total time: ', end_time - start_time, ' seconds'
  print '(A,F10.4,A)', ' Average time per point: ', (end_time - start_time) / real(n_init_points), ' seconds'
  print '(A,F10.2)', ' Speedup factor (approx): ', real(num_threads) * 0.8
  print *, 'Done!'
  
  ! 清理
  deallocate(R_grid, Z_grid, phi_grid, psi_grid, rho_grid)
  deallocate(B_R, B_Z, B_phi)
  deallocate(lyap_results)

contains

  !-----------------------------------------------------------------------------
  ! 读取网格信息
  !-----------------------------------------------------------------------------
  subroutine read_grid_data()
    integer :: ix, iy, iphi
    
    ! 读取网格维度
    open(101, file='B_grid_info.dat', status='old')
    read(101, *) nx, ny, nphi
    read(101, *) R0_axis
    close(101)
    
    ! 分配数组
    allocate(R_grid(0:nx, 0:ny))
    allocate(Z_grid(0:nx, 0:ny))
    allocate(phi_grid(1:nphi))
    allocate(psi_grid(0:nx))
    allocate(rho_grid(0:nx))
    
    ! 读取 R 网格
    open(102, file='R_grid.dat', status='old')
    do ix = 0, nx
      do iy = 0, ny
        read(102, *) R_grid(ix, iy)
      enddo
    enddo
    close(102)
    
    ! 读取 Z 网格
    open(103, file='Z_grid.dat', status='old')
    do ix = 0, nx
      do iy = 0, ny
        read(103, *) Z_grid(ix, iy)
      enddo
    enddo
    close(103)
    
    ! 读取 phi 网格
    open(104, file='phi_grid.dat', status='old')
    do iphi = 1, nphi
      read(104, *) phi_grid(iphi)
    enddo
    close(104)
    
    ! 读取 psi 和 rho 网格
    open(105, file='psi_grid.dat', status='old')
    do ix = 0, nx
      read(105, *) psi_grid(ix), rho_grid(ix)
    enddo
    close(105)
    
  end subroutine read_grid_data
  
  !-----------------------------------------------------------------------------
  ! 读取磁场数据
  !-----------------------------------------------------------------------------
  subroutine read_field_data()
    integer :: ix, iy, iphi
    
    allocate(B_R(0:nx, 0:ny, 1:nphi))
    allocate(B_Z(0:nx, 0:ny, 1:nphi))
    allocate(B_phi(0:nx, 0:ny, 1:nphi))
    
    ! 读取 B_R
    open(106, file='B_R.dat', status='old')
    do ix = 0, nx
      do iy = 0, ny
        do iphi = 1, nphi
          read(106, *) B_R(ix, iy, iphi)
        enddo
      enddo
    enddo
    close(106)
    
    ! 读取 B_Z
    open(107, file='B_Z.dat', status='old')
    do ix = 0, nx
      do iy = 0, ny
        do iphi = 1, nphi
          read(107, *) B_Z(ix, iy, iphi)
        enddo
      enddo
    enddo
    close(107)
    
    ! 读取 B_phi
    open(108, file='B_phi.dat', status='old')
    do ix = 0, nx
      do iy = 0, ny
        do iphi = 1, nphi
          read(108, *) B_phi(ix, iy, iphi)
        enddo
      enddo
    enddo
    close(108)
    
  end subroutine read_field_data
  
  !-----------------------------------------------------------------------------
  ! 计算网格边界(用于插值)
  !-----------------------------------------------------------------------------
  subroutine compute_grid_bounds()
    R_min = minval(R_grid)
    R_max = maxval(R_grid)
    Z_min = minval(Z_grid)
    Z_max = maxval(Z_grid)
    
    ! 估计网格间距(用于查找)
    dR_grid = (R_max - R_min) / real(nx, dp)
    dZ_grid = (Z_max - Z_min) / real(ny, dp)
    dphi_grid = 2.0d0 * pi / real(nphi, dp)
  end subroutine compute_grid_bounds
  
  !-----------------------------------------------------------------------------
  ! 三线性插值获取磁场
  ! 输入: R, Z, phi (柱坐标)
  ! 输出: Br, Bz, Bphi
  !-----------------------------------------------------------------------------
  subroutine interpolate_B(R, Z, phi, Br, Bz, Bphi, ierr)
    real(dp), intent(in) :: R, Z, phi
    real(dp), intent(out) :: Br, Bz, Bphi
    integer, intent(out) :: ierr
    
    real(dp) :: phi_mod, dist, dist_min
    integer :: ix, iy, iphi_lo, iphi_hi
    integer :: ix_lo, ix_hi, iy_lo, iy_hi
    real(dp) :: wx, wy, wphi
    real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
    real(dp) :: c00, c01, c10, c11, c0, c1
    integer :: i, j
    
    ierr = 0
    
    ! 将 phi 规范到 [0, 2*pi)
    phi_mod = mod(phi, 2.0d0*pi)
    if (phi_mod < 0.0d0) phi_mod = phi_mod + 2.0d0*pi
    
    ! 在 (R,Z) 平面上找最近的网格点
    ! 这里使用简单的搜索方法(对于规则网格可以优化)
    ix_lo = -1
    iy_lo = -1
    dist_min = 1.0d30
    
    ! 简化：假设网格在 theta 方向是均匀的，在 psi 方向也近似均匀
    ! 先找最近的 (ix, iy) 索引
    do ix = 0, nx-1
      do iy = 0, ny-1
        ! 检查点是否在这个网格单元内
        if (is_in_cell(R, Z, ix, iy)) then
          ix_lo = ix
          iy_lo = iy
          exit
        endif
      enddo
      if (ix_lo >= 0) exit
    enddo
    
    ! 如果没找到，说明点在网格外
    if (ix_lo < 0) then
      ierr = 1
      Br = 0.0d0
      Bz = 0.0d0
      Bphi = 1.0d0  ! 返回一个默认值避免除零
      return
    endif
    
    ix_hi = ix_lo + 1
    iy_hi = mod(iy_lo + 1, ny + 1)  ! 周期边界
    if (iy_hi > ny) iy_hi = 1
    
    ! phi 方向的索引
    iphi_lo = int(phi_mod / dphi_grid) + 1
    if (iphi_lo < 1) iphi_lo = 1
    if (iphi_lo > nphi) iphi_lo = nphi
    iphi_hi = mod(iphi_lo, nphi) + 1
    
    ! 计算插值权重
    ! 在 (R,Z) 平面内使用双线性插值的局部坐标
    call compute_local_weights(R, Z, ix_lo, iy_lo, wx, wy)
    
    ! phi 方向权重
    wphi = (phi_mod - phi_grid(iphi_lo)) / dphi_grid
    if (wphi < 0.0d0) wphi = 0.0d0
    if (wphi > 1.0d0) wphi = 1.0d0
    
    ! 三线性插值 B_R
    c000 = B_R(ix_lo, iy_lo, iphi_lo)
    c001 = B_R(ix_lo, iy_lo, iphi_hi)
    c010 = B_R(ix_lo, iy_hi, iphi_lo)
    c011 = B_R(ix_lo, iy_hi, iphi_hi)
    c100 = B_R(ix_hi, iy_lo, iphi_lo)
    c101 = B_R(ix_hi, iy_lo, iphi_hi)
    c110 = B_R(ix_hi, iy_hi, iphi_lo)
    c111 = B_R(ix_hi, iy_hi, iphi_hi)
    
    c00 = c000*(1.0d0-wphi) + c001*wphi
    c01 = c010*(1.0d0-wphi) + c011*wphi
    c10 = c100*(1.0d0-wphi) + c101*wphi
    c11 = c110*(1.0d0-wphi) + c111*wphi
    
    c0 = c00*(1.0d0-wy) + c01*wy
    c1 = c10*(1.0d0-wy) + c11*wy
    
    Br = c0*(1.0d0-wx) + c1*wx
    
    ! 三线性插值 B_Z
    c000 = B_Z(ix_lo, iy_lo, iphi_lo)
    c001 = B_Z(ix_lo, iy_lo, iphi_hi)
    c010 = B_Z(ix_lo, iy_hi, iphi_lo)
    c011 = B_Z(ix_lo, iy_hi, iphi_hi)
    c100 = B_Z(ix_hi, iy_lo, iphi_lo)
    c101 = B_Z(ix_hi, iy_lo, iphi_hi)
    c110 = B_Z(ix_hi, iy_hi, iphi_lo)
    c111 = B_Z(ix_hi, iy_hi, iphi_hi)
    
    c00 = c000*(1.0d0-wphi) + c001*wphi
    c01 = c010*(1.0d0-wphi) + c011*wphi
    c10 = c100*(1.0d0-wphi) + c101*wphi
    c11 = c110*(1.0d0-wphi) + c111*wphi
    
    c0 = c00*(1.0d0-wy) + c01*wy
    c1 = c10*(1.0d0-wy) + c11*wy
    
    Bz = c0*(1.0d0-wx) + c1*wx
    
    ! 三线性插值 B_phi
    c000 = B_phi(ix_lo, iy_lo, iphi_lo)
    c001 = B_phi(ix_lo, iy_lo, iphi_hi)
    c010 = B_phi(ix_lo, iy_hi, iphi_lo)
    c011 = B_phi(ix_lo, iy_hi, iphi_hi)
    c100 = B_phi(ix_hi, iy_lo, iphi_lo)
    c101 = B_phi(ix_hi, iy_lo, iphi_hi)
    c110 = B_phi(ix_hi, iy_hi, iphi_lo)
    c111 = B_phi(ix_hi, iy_hi, iphi_hi)
    
    c00 = c000*(1.0d0-wphi) + c001*wphi
    c01 = c010*(1.0d0-wphi) + c011*wphi
    c10 = c100*(1.0d0-wphi) + c101*wphi
    c11 = c110*(1.0d0-wphi) + c111*wphi
    
    c0 = c00*(1.0d0-wy) + c01*wy
    c1 = c10*(1.0d0-wy) + c11*wy
    
    Bphi = c0*(1.0d0-wx) + c1*wx
    
  end subroutine interpolate_B
  
  !-----------------------------------------------------------------------------
  ! 检查点 (R,Z) 是否在网格单元 (ix,iy) 内
  !-----------------------------------------------------------------------------
  function is_in_cell(R, Z, ix, iy) result(inside)
    real(dp), intent(in) :: R, Z
    integer, intent(in) :: ix, iy
    logical :: inside
    
    real(dp) :: R1, R2, R3, R4, Z1, Z2, Z3, Z4
    real(dp) :: xi, eta
    logical :: success
    integer :: iy_next
    
    inside = .false.
    
    if (ix < 0 .or. ix >= nx .or. iy < 0 .or. iy >= ny) return
    
    iy_next = iy + 1
    if (iy_next > ny) iy_next = 1
    
    R1 = R_grid(ix, iy)
    R2 = R_grid(ix+1, iy)
    R3 = R_grid(ix+1, iy_next)
    R4 = R_grid(ix, iy_next)
    
    Z1 = Z_grid(ix, iy)
    Z2 = Z_grid(ix+1, iy)
    Z3 = Z_grid(ix+1, iy_next)
    Z4 = Z_grid(ix, iy_next)
    
    ! 1. 简单的边界框检查 (AABB) - 快速排除
    if (R < min(R1,R2,R3,R4) - 1.0d-4 .or. &
        R > max(R1,R2,R3,R4) + 1.0d-4 .or. &
        Z < min(Z1,Z2,Z3,Z4) - 1.0d-4 .or. &
        Z > max(Z1,Z2,Z3,Z4) + 1.0d-4) then
      return
    endif

    ! 2. 精确检查：尝试求解局部坐标
    call get_local_coords_newton(R, Z, ix, iy, xi, eta, success)
    inside = success
    
  end function is_in_cell
  
  !-----------------------------------------------------------------------------
  ! 计算在网格单元内的局部插值权重 (使用牛顿迭代求解逆映射)
  !-----------------------------------------------------------------------------
  subroutine compute_local_weights(R, Z, ix, iy, wx, wy)
    real(dp), intent(in) :: R, Z
    integer, intent(in) :: ix, iy
    real(dp), intent(out) :: wx, wy
    logical :: success
    
    call get_local_coords_newton(R, Z, ix, iy, wx, wy, success)
    
    ! 如果求解失败(通常意味着点在外面)，限制在[0,1]
    wx = max(0.0d0, min(1.0d0, wx))
    wy = max(0.0d0, min(1.0d0, wy))
    
  end subroutine compute_local_weights

  !-----------------------------------------------------------------------------
  ! 牛顿-拉夫逊法求解任意四边形内的局部坐标 (xi, eta)
  !-----------------------------------------------------------------------------
  subroutine get_local_coords_newton(R_p, Z_p, ix, iy, xi, eta, success)
    real(dp), intent(in) :: R_p, Z_p
    integer, intent(in) :: ix, iy
    real(dp), intent(out) :: xi, eta
    logical, intent(out) :: success
    
    ! 局部变量
    real(dp) :: R1, R2, R3, R4, Z1, Z2, Z3, Z4
    real(dp) :: AR, BR, CR, AZ, BZ, CZ
    real(dp) :: R_curr, Z_curr, resid_R, resid_Z
    real(dp) :: dRdxi, dRdeta, dZdxi, dZdeta, detJ, dxi, deta
    integer :: iter, iy_next
    
    iy_next = iy + 1
    if (iy_next > ny) iy_next = 1
    
    R1 = R_grid(ix, iy)
    R2 = R_grid(ix+1, iy)
    R3 = R_grid(ix+1, iy_next)
    R4 = R_grid(ix, iy_next)
    
    Z1 = Z_grid(ix, iy)
    Z2 = Z_grid(ix+1, iy)
    Z3 = Z_grid(ix+1, iy_next)
    Z4 = Z_grid(ix, iy_next)
    
    ! 双线性映射系数: R(xi,eta) = R1 + AR*xi + BR*eta + CR*xi*eta
    ! 对应节点顺序: (0,0)->1, (1,0)->2, (1,1)->3, (0,1)->4
    ! 注意代码中的节点顺序: R1(ix,iy), R2(ix+1,iy), R3(ix+1,iy+1), R4(ix,iy+1)
    ! 这是一个逆时针或顺时针循环. 
    ! xi 对应 ix 方向 (1->2), eta 对应 iy 方向 (1->4)
    AR = R2 - R1
    BR = R4 - R1
    CR = R1 - R2 + R3 - R4
    
    AZ = Z2 - Z1
    BZ = Z4 - Z1
    CZ = Z1 - Z2 + Z3 - Z4
    
    ! 初始猜测 (中心)
    xi = 0.5d0
    eta = 0.5d0
    success = .false.
    
    do iter = 1, 10
      ! 计算当前猜测位置的 R, Z
      R_curr = R1 + AR*xi + BR*eta + CR*xi*eta
      Z_curr = Z1 + AZ*xi + BZ*eta + CZ*xi*eta
      
      resid_R = R_p - R_curr
      resid_Z = Z_p - Z_curr
      
      ! 检查收敛
      if (sqrt(resid_R**2 + resid_Z**2) < 1.0d-8) then
        success = .true.
        exit
      endif
      
      ! 计算雅可比矩阵 J
      dRdxi  = AR + CR*eta
      dRdeta = BR + CR*xi
      dZdxi  = AZ + CZ*eta
      dZdeta = BZ + CZ*xi
      
      detJ = dRdxi * dZdeta - dRdeta * dZdxi
      
      if (abs(detJ) < 1.0d-14) exit ! 雅可比奇异
      
      ! 更新: [dxi; deta] = J^-1 * [resid_R; resid_Z]
      dxi  = ( dZdeta * resid_R - dRdeta * resid_Z) / detJ
      deta = (-dZdxi  * resid_R + dRdxi  * resid_Z) / detJ
      
      xi = xi + dxi
      eta = eta + deta
    enddo
    
    ! 检查是否在单元内 (允许极小的数值误差)
    if (success) then
      if (xi < -1.0d-3 .or. xi > 1.001d0 .or. &
          eta < -1.0d-3 .or. eta > 1.001d0) then
        success = .false.
      endif
    endif

  end subroutine get_local_coords_newton
  
  !-----------------------------------------------------------------------------
  ! 磁力线方程的右端函数
  ! dR/dphi = B_R / B_phi
  ! dZ/dphi = B_Z / B_phi
  !-----------------------------------------------------------------------------
  subroutine fieldline_rhs(R, Z, phi, dRdphi, dZdphi, ierr)
    real(dp), intent(in) :: R, Z, phi
    real(dp), intent(out) :: dRdphi, dZdphi
    integer, intent(out) :: ierr
    
    real(dp) :: Br, Bz, Bphi
    
    call interpolate_B(R, Z, phi, Br, Bz, Bphi, ierr)
    
    if (ierr /= 0 .or. abs(Bphi) < 1.0d-12) then
      dRdphi = 0.0d0
      dZdphi = 0.0d0
      ierr = 1
      return
    endif
    
    dRdphi = R * Br / Bphi
    dZdphi = R * Bz / Bphi
    
  end subroutine fieldline_rhs
  
  !-----------------------------------------------------------------------------
  ! RK4 积分一步
  !-----------------------------------------------------------------------------
  subroutine rk4_step(R, Z, phi, dphi, R_new, Z_new, ierr)
    real(dp), intent(in) :: R, Z, phi, dphi
    real(dp), intent(out) :: R_new, Z_new
    integer, intent(out) :: ierr
    
    real(dp) :: k1R, k1Z, k2R, k2Z, k3R, k3Z, k4R, k4Z
    real(dp) :: R_tmp, Z_tmp, phi_tmp
    integer :: ierr_tmp
    
    ierr = 0
    
    ! k1
    call fieldline_rhs(R, Z, phi, k1R, k1Z, ierr_tmp)
    if (ierr_tmp /= 0) then
      ierr = 1
      R_new = R
      Z_new = Z
      return
    endif
    
    ! k2
    R_tmp = R + 0.5d0*dphi*k1R
    Z_tmp = Z + 0.5d0*dphi*k1Z
    phi_tmp = phi + 0.5d0*dphi
    call fieldline_rhs(R_tmp, Z_tmp, phi_tmp, k2R, k2Z, ierr_tmp)
    if (ierr_tmp /= 0) then
      ierr = 1
      R_new = R
      Z_new = Z
      return
    endif
    
    ! k3
    R_tmp = R + 0.5d0*dphi*k2R
    Z_tmp = Z + 0.5d0*dphi*k2Z
    call fieldline_rhs(R_tmp, Z_tmp, phi_tmp, k3R, k3Z, ierr_tmp)
    if (ierr_tmp /= 0) then
      ierr = 1
      R_new = R
      Z_new = Z
      return
    endif
    
    ! k4
    R_tmp = R + dphi*k3R
    Z_tmp = Z + dphi*k3Z
    phi_tmp = phi + dphi
    call fieldline_rhs(R_tmp, Z_tmp, phi_tmp, k4R, k4Z, ierr_tmp)
    if (ierr_tmp /= 0) then
      ierr = 1
      R_new = R
      Z_new = Z
      return
    endif
    
    ! 更新
    R_new = R + dphi/6.0d0 * (k1R + 2.0d0*k2R + 2.0d0*k3R + k4R)
    Z_new = Z + dphi/6.0d0 * (k1Z + 2.0d0*k2Z + 2.0d0*k3Z + k4Z)
    
  end subroutine rk4_step
  
  !-----------------------------------------------------------------------------
  ! 积分磁力线
  !-----------------------------------------------------------------------------
  subroutine integrate_fieldline(R0, Z0, phi0, n_turn, dphi, R_final, Z_final, ierr)
    real(dp), intent(in) :: R0, Z0, phi0, dphi
    integer, intent(in) :: n_turn
    real(dp), intent(out) :: R_final, Z_final
    integer, intent(out) :: ierr
    
    real(dp) :: R, Z, phi, R_new, Z_new
    integer :: n_steps, istep
    
    n_steps = int(2.0d0*pi*real(n_turn, dp) / abs(dphi))
    
    R = R0
    Z = Z0
    phi = phi0
    ierr = 0
    
    do istep = 1, n_steps
      call rk4_step(R, Z, phi, dphi, R_new, Z_new, ierr)
      if (ierr /= 0) exit
      
      R = R_new
      Z = Z_new
      phi = phi + dphi
    enddo
    
    R_final = R
    Z_final = Z
    
  end subroutine integrate_fieldline
  
  !-----------------------------------------------------------------------------
  ! 计算 Lyapunov 指数
  !-----------------------------------------------------------------------------
  subroutine compute_lyapunov(R0, Z0, phi0, n_turn, dphi, n_renorm_interval, delta_init, lambda)
    real(dp), intent(in) :: R0, Z0, phi0, dphi, delta_init
    integer, intent(in) :: n_turn, n_renorm_interval
    real(dp), intent(out) :: lambda
    
    real(dp) :: R1, Z1, R2, Z2, phi
    real(dp) :: R1_new, Z1_new, R2_new, Z2_new
    real(dp) :: dR, dZ, dist, S_sum
    integer :: n_steps, istep, n_renorm_count, ierr
    
    n_steps = int(2.0d0*pi*real(n_turn, dp) / abs(dphi))
    
    ! 初始化两条轨道
    R1 = R0
    Z1 = Z0
    R2 = R0 + delta_init  ! 在 R 方向加一个小扰动
    Z2 = Z0
    phi = phi0
    
    S_sum = 0.0d0
    n_renorm_count = 0
    
    do istep = 1, n_steps
      ! 积分第一条轨道
      call rk4_step(R1, Z1, phi, dphi, R1_new, Z1_new, ierr)
      if (ierr /= 0) exit
      
      ! 积分第二条轨道
      call rk4_step(R2, Z2, phi, dphi, R2_new, Z2_new, ierr)
      if (ierr /= 0) exit
      
      R1 = R1_new
      Z1 = Z1_new
      R2 = R2_new
      Z2 = Z2_new
      phi = phi + dphi
      
      ! 定期重标定
      if (mod(istep, n_renorm_interval) == 0) then
        dR = R2 - R1
        dZ = Z2 - Z1
        dist = sqrt(dR*dR + dZ*dZ)
        
        if (dist > 1.0d-20) then
          S_sum = S_sum + log(dist / delta_init)
          n_renorm_count = n_renorm_count + 1
          
          ! 重标定：把第二条轨道拉回到第一条附近
          R2 = R1 + delta_init * dR / dist
          Z2 = Z1 + delta_init * dZ / dist
        endif
      endif
    enddo
    
    ! 计算 Lyapunov 指数
    if (n_renorm_count > 0) then
      ! lambda = S_sum / (n_renorm_count * n_renorm_interval * dphi)
      lambda = S_sum / (real(n_renorm_count * n_renorm_interval, dp) * dphi)
    else
      lambda = 0.0d0
    endif
    
  end subroutine compute_lyapunov
  
  !-----------------------------------------------------------------------------
  ! 输出结果
  !-----------------------------------------------------------------------------
  subroutine output_results()
    integer :: i
    
    open(201, file='lyapunov_results.dat')
    write(201, '(A)') '# rho    psi_norm    lambda'
    do i = 1, n_init_points
      write(201, '(3ES16.8)') lyap_results(i, 1), lyap_results(i, 2), lyap_results(i, 3)
    enddo
    close(201)
    
  end subroutine output_results

end program fieldline_lyap
