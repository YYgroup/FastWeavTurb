! Generation of knotted fields
! Shiying Xiong, 2019
! Weiyu Shen, 2021, ver-1.1 fix construction method of twisted cases
! Weiyu Shen, 2021 Dec, ver-2.0 add sigmaFunc/etaFunc for bursting cases and fix bugs
! Weiyu Shen, 2022 Jan, ver-2.1 add twist of different vortex surfaces
! Weiyu Shen, 2022 Apr, ver-2.2 unify the output format
! Weiyu Shen, 2022 Jul, ver-3.0 update construction algorithm to closed flux tubes with arbitrary centerlines, vortex surface field, tube thickness distribution and local twist rate distribution
! Weiyu Shen, 2022 Nov, ver-4.0 update construction algorithm to accommodate centerline across boundaries
! Weiyu Shen, 2022 Nov, ver-4.1 use cubic splines to fit the scatter centerline
! Weiyu Shen, 2023 Feb, ver-4.2 use cubic splines to fit multiple scatter centerlines
! Weiyu Shen, 2023 Mar, ver-4.3 fix bugs
! Weiyu Shen, 2023 Mar, ver-4.4 update construction algorithm for stop-points and straight lines
! Weiyu Shen, 2023 Mar, ver-4.5 fix bugs
! Weiyu Shen, 2023 Mar, ver-4.6 add the output for structure functions
! Weiyu Shen, 2023 Mar, ver-4.7 fix bugs; add the calculation of the total length, total writhe, average vortex length density and average vortex volume ratio; update sigmaFunc/etaFunc for arc-length parameter
! Weiyu Shen, 2023 Apr, ver-4.8 fix bugs
! Weiyu Shen, 2023 Apr, ver-4.9 fix bugs; adaptive nt, ntc
! Weiyu Shen, 2023 Apr, ver-4.10 change output field from dat file to plt file
! Weiyu Shen, 2023 June, ver-4.11 add the 7box mode
! Zishuo Han, 2024 April, ver-5.0 comprehensive improvement; efficiency niubi
! Weiyu Shen, 2024 April, ver-5.1 fix bugs
! Zishuo Han, 2024 April, ver-5.2 change the handling of straight line segments (see mark v5_3_1);deal with the bug when ci is on the mesh (see mark v5_3_2)
! Zishuo Han, 2024 April, ver-5.3 use the theoretical derivative instead of the central difference
! Zishuo Han, 2024 April, ver-5.4 employ Len_Disc(:) to avoid numerical issues caused by small cze(i+1)-cze(i)
! Zishuo Han, 2024 April, ver-5.5 employ 5th spline to ensure the continuity of Frenet vector
! Zishuo Han, 2024 April, ver-5.6 employ dynamic dtf
! Zishuo Han, 2024 Oct, ver-5.7 adopt adjustable multi-scale vortices
! Zishuo Han, 2025 Apr, ver-5.8 modify input and output

program Weav_ver5_8_2
   implicit none
   include 'mpif.h'
   include "fftw3.f"

   !==parameters==
   integer nx, ny, nz
   real*8 sigma_min
   real*8 sigma_max, sigma_ratio
   real*8 Gamma
   real*8 Gamma_line_ratio, sigma_line_ratio
   real*8, allocatable :: Gamma_line(:), sigma_line(:)

   integer ntf_times
   integer ntf_pre
   real*8 Rtube, Rtube_ratio
   real*8 dzeta_sup

   integer, parameter :: igst = 3
   !==parameters==

   !==switch==
   !=0:close,=1:open
   !integer, parameter :: if_filter = 1
   integer, parameter :: if_spline_out = 0
   integer, parameter :: vorticity_output = 1, u_out = 0, d_out = 0, multi_out = 0
   integer, parameter :: field2d_output = 0
   integer, parameter :: if_pdf = 1
   integer, parameter :: if_Sp = 1, Sp_order = 10
   !==switch==

   !==centerline==
   integer line_s, line_e
   integer nline, npoint, npointall
   integer, allocatable :: npointlist(:), npointpast(:)
   integer n_sec ! total number of sections
   integer i_line, nt, t_i, t_j
   real*8 Rkappa
   real*8 writhe, twist, slk, dist, cdx, cdy, cdz, dist2, writhetot, lengthtot
   real*8, allocatable :: length(:), Len_Disc(:) !temp
   real*8, allocatable :: c(:, :), dc(:, :), ddc(:, :), dddc(:, :), dcddc(:, :)
   real*8, allocatable :: kappa(:), tau(:), ndc(:), ndcddc(:), vc(:)
   real*8, allocatable :: cx(:), cy(:), cz(:), cze(:), ckx(:), cky(:), ckz(:), ckkx(:), ckky(:), ckkz(:)
 real*8, allocatable :: cxall(:), cyall(:), czall(:), czeall(:), ckxall(:), ckyall(:), ckzall(:), ckkxall(:), ckkyall(:), ckkzall(:)
 real*8, allocatable :: cx_pl(:), cy_pl(:), cz_pl(:), cze_pl(:), ckx_pl(:), cky_pl(:), ckz_pl(:), ckkx_pl(:), ckky_pl(:), ckkz_pl(:)
   real*8, allocatable, dimension(:, :, :) :: coe_list
   real*8, dimension(0:5, 3) :: coe_3
   real*8 dtf_min, dtf_length_max
   integer ntf_max
   real*8, allocatable, dimension(:)  :: dtf_list
   integer, allocatable, dimension(:)  :: ntf_list
   real*8, allocatable, dimension(:)  :: kappa_list
   real*8, allocatable, dimension(:, :, :)  :: ci_list2, dci_list2, tci_list2
   real*8, allocatable, dimension(:, :)  :: ndci_list2, s_list2
   real*8, allocatable, dimension(:, :, :)  :: ci_list_pre, dci_list_pre
   real*8, allocatable, dimension(:)  :: dtf_list_pre
   real*8, allocatable, dimension(:, :)  :: s_list_pre, ndci_list_pre !think later: theory cum?
   real*8, dimension(3) :: ci0, ci1, ci2, dci0, dci1, tci0, tci1
   real*8 s0, s1, s2, ndci0, ndci1
   real*8, dimension(3) :: xczeta, Tzeta, Nzeta, Bzeta, e_th, e_rho
   real*8, dimension(3) :: czeta, dczeta, ddczeta, dddczeta, dcddczeta
   real*8 zeta, ndczeta, ndcddczeta, kappazeta, tauzeta, nxczeta
   real*8 kappa_sup, theta
   real*8 dzeta, zeta_sec
   !==centerline==

   !==vortex tube==
   real*8 eta, eta0
   real*8 pdx ! ??
   real*8 sigma, dsigmads
   real*8 rho, costh, sinth
   real*8 ft_delta
   integer n_T, n_T_line_ratio
   integer, allocatable :: n_T_line(:), if_filter_line(:)
   real*4 sigma_out, Gamma_out
   !==vortex tube==

   !==field==
   real*8 xstart, ystart, zstart, lx, ly, lz, lt
   real*8 dx, dy, dz, dt, dv
   real*8, dimension(3) :: meshxyz
   real*8, allocatable :: meshx_plus(:), meshy_plus(:), meshz_plus(:)
   integer nx_plus, ny_plus, nz_plus, i_plus, j_plus, k_plus
   real*8, allocatable :: kx(:), ky(:), kz(:), k2(:, :, :)
   !real*8, allocatable :: kx_d1(:), ky_d1(:), kz_d1(:), k2_d1(:, :, :)
   real*8, allocatable :: vorx(:, :, :), vory(:, :, :), vorz(:, :, :), phiv(:, :, :)
   !real*8, allocatable :: vorx_d1(:, :, :), vory_d1(:, :, :), vorz_d1(:, :, :)
   real*4, allocatable :: tmp_real4(:, :, :)
   real*8, allocatable :: velx(:, :, :), vely(:, :, :), velz(:, :, :)
   real*8, allocatable :: disp_field(:, :, :)
   complex(8), allocatable :: spec_ux(:, :, :), spec_uy(:, :, :), spec_uz(:, :, :)
   complex(8), allocatable :: spec_wx(:, :, :), spec_wy(:, :, :), spec_wz(:, :, :)
   complex(8), allocatable :: spec(:, :, :)
   !==field==

   !==statistics==
   real*8 deviation, helicity, Etot, Diss, Enstrophy, up
   integer norder
   real*8, allocatable :: Sn(:), r(:)
   integer, parameter :: bin_num = 1000
   !==statistics==

   !==basic==
   real*8, parameter :: pi = 4.0d0*datan(1.0d0), pi2 = 8.0d0*datan(1.0d0)

   integer i, j, k, ii
   integer i_min, i_max, j_min, j_max, k_min, k_max
   integer ijk_co
   integer k_c
   integer id1, id2, id3, id4, id5

   integer r_d
   integer n_d

   real*8 tp, t0, t1, t2, t3, t4, t5, t6, t7, t8

   character*200 name
   real, allocatable :: data_box(:, :, :, :)
   character*40, allocatable :: varname(:)

   integer status(MPI_STATUS_SIZE)
   integer ierr
   integer id, id_l, id_r, nproc
   integer nzp

   integer*8 planxf, planxb, planyf, planyb, planzf, planzb !!! fft plans
   !integer*8 planxf_d1, planxb_d1, planyf_d1, planyb_d1, planzf_d1, planzb_d1 !!! fft plans

   double precision starttime, endtime, start_time_total, end_time_total, starttime_part, endtime_part
   !==basic==

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
   !desk
   ! nproc = 1
   ! id = 0

   !read nline
   if (id == 0) then
      open (90, file='./centerline_input.dat', status='unknown')
      open (31, file='./output/centerline_output.dat', status='unknown')
      read (90, *) nline
      write (31, *) nline
      print *, '=============================='
      print *, 'number of lines:', nline
   end if
   call mpi_bcast(nline, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)

   allocate (Gamma_line(nline))
   allocate (sigma_line(nline))
   allocate (n_T_line(nline))
   allocate (if_filter_line(nline))

   ! do i = 1, nline
   !    Gamma_line(i) = Gamma*Gamma_line_ratio**(i - 1)
   !    sigma_line(i) = sigma_min*sigma_line_ratio**(i - 1)
   !    n_T_line(i) = n_T*n_T_line_ratio**(i - 1)
   ! end do

   !parameters input
   open (512, file='./parameters_input.dat', status='unknown')
   read (512, *)
   read (512, *) nx
   read (512, *)
   read (512, *) n_d
   read (512, *)
   read (512, *) line_s
   read (512, *)
   read (512, *) line_e
   ny = nx
   nz = nx
   read (512, *)
   read (512, *) sigma_ratio
   !sigma_max = sigma_ratio*sigma_min
   read (512, *)
   read (512, *) eta0
   !----------
   read (512, *)
   read (512, *)
   read (512, *) ntf_times
   read (512, *)
   read (512, *) ntf_pre
   read (512, *)
   read (512, *) Rtube_ratio
   !Rtube = Rtube_ratio*sigma_max
   read (512, *)
   read (512, *) dzeta_sup
   !----------
   read (512, *)
   read (512, *)
   do i_line = 1, nline
      read (512, *) i, sigma_line(i), Gamma_line(i), n_T_line(i), if_filter_line(i)
   end do
   !parameters input

   !parameters output
   if (id == 0) then
      open (19, file='./output/parameters_output.dat', status='unknown')
      write (19, *) 'N'
      write (19, *) nx
      write (19, *) 'N_downsample'
      write (19, *) n_d
      write (19, *) 'line_s'
      write (19, *) line_s
      write (19, *) 'line_e'
      write (19, *) line_e
      write (19, *) 'sigma_ratio'
      write (19, *) sigma_ratio
      write (19, *) 'eta_t'
      write (19, *) eta0
      write (19, *) '---------------'
      write (19, *) 'ntf_times'
      write (19, *) ntf_times
      write (19, *) 'ntf_pre'
      write (19, *) ntf_pre
      write (19, *) 'Rtube_ratio'
      write (19, *) Rtube_ratio
      write (19, *) 'dzeta_sup'
      write (19, *) dzeta_sup
      write (19, *) '---------------'
      write (19, *) 'variables="i_line","sigma","Gamma","n_T","if_filter"'
      do i_line = 1, nline
         sigma_out = sigma_line(i_line)
         Gamma_out = Gamma_line(i_line)
         write (19, *) i_line, sigma_out, Gamma_out, n_T_line(i_line), if_filter_line(i_line)
      end do
      close (19)
   end if
   !parameters output

   id_l = mod(nproc + id - 1, nproc)
   id_r = mod(id + 1, nproc)
   nzp = nz/nproc
   if (id == 0) print *, 'nproc=', nproc
   if (id == 0) print *, 'nzp=', nzp

   xstart = -pi
   ystart = -pi
   zstart = -pi
   lx = pi2
   ly = pi2
   lz = pi2
   lt = pi2
   dx = lx/nx
   dy = ly/ny
   dz = lz/nz
   dv = dx*dy*dz
   dtf_length_max = 0.5d0*max(lx, ly, lz)!3.0d0*max(dx, dy, dz) !study it later !length
   if (id == 0) then
      print *, 'dtf_length_max=', dtf_length_max
      !print *, 'Rtube = ', Rtube
   end if

   !consider dtf_length_max, see mark1
   nx_plus = ceiling((Rtube + dtf_length_max)/dx)
   ny_plus = ceiling((Rtube + dtf_length_max)/dy)
   nz_plus = ceiling((Rtube + dtf_length_max)/dz)

   allocate (k2(nx, ny, nzp))
   allocate (kx(nx))
   allocate (ky(ny))
   allocate (kz(nzp))
   allocate (Sn(nx/2 + 1))
   allocate (r(nx/2 + 1))
   allocate (meshx_plus(1 - nx_plus:nx + nx_plus))
   allocate (meshy_plus(1 - ny_plus:ny + ny_plus))
   allocate (meshz_plus(1 - nz_plus:nz + nz_plus)) ! nz, not nzp
   allocate (vorx(nx, ny, nzp))
   allocate (vory(nx, ny, nzp))
   allocate (vorz(nx, ny, nzp))
   !allocate (specx(nx, ny, nzp))
   !allocate (specy(nx, ny, nzp))
   !allocate (specz(nx, ny, nzp))
   !allocate (phiv(1 - igst:nx + igst, 1 - igst:ny + igst, 1 - igst:nzp + igst))

   call set_fft_plans(nx, ny, nz, planxf, planxb, planyf, planyb, planzf, planzb)
   call wavenumber(nx, ny, nz, nzp, lx, ly, lz, kx, ky, kz, k2, id) !mind gyc

   !test
   if (id == 0) then
      print *, 'try to go through'
   end if
   allocate (spec_ux(nx, ny, nzp))
   allocate (spec_uy(nx, ny, nzp))
   allocate (spec_uz(nx, ny, nzp))
   !allocate (spec_wx(nx, ny, nzp))
   !allocate (spec_wy(nx, ny, nzp))
   !allocate (tmp_real4(nx, ny, nzp))
   deallocate (spec_ux)
   deallocate (spec_uy)
   deallocate (spec_uz)
   !deallocate (spec_wx)
   !deallocate (spec_wy)
   !deallocate (tmp_real4)
   if (id == 0) then
      print *, 'go through successfully'
   end if
   !test

   !call initialize_mesh(nx, ny, nzp, xstart, ystart, zstart, dx, dy, dz, meshx, meshy, meshz, id)
   do i_plus = 1 - nx_plus, nx + nx_plus !mind
      meshx_plus(i_plus) = (i_plus - 1.d0)*dx + xstart
   end do
   do j_plus = 1 - ny_plus, ny + ny_plus !mind
      meshy_plus(j_plus) = (j_plus - 1.d0)*dy + ystart
   end do
   do k_plus = 1 - nz_plus, nz + nz_plus !mind
      meshz_plus(k_plus) = (k_plus - 1.d0)*dz + zstart
   end do

   allocate (npointlist(nline))
   allocate (npointpast(nline))
   allocate (length(nline))
   allocate (Len_Disc(nline))
   if (id == 0) then
      do i = 1, nline
         read (90, *) npointlist(i)
         write (31, *) npointlist(i)
         !print *, '===================='
         print *, 'line:', i, 'npoint:', npointlist(i)
      end do
   end if
   call mpi_bcast(npointlist, nline, MPI_INTEGER, 0, mpi_comm_world, ierr)
   npointall = 0
   do i = 1, nline
      npointall = npointall + npointlist(i)
   end do
   allocate (cxall(npointall))
   allocate (cyall(npointall))
   allocate (czall(npointall))
   allocate (czeall(npointall))
   allocate (ckxall(npointall))
   allocate (ckyall(npointall))
   allocate (ckzall(npointall))
   allocate (ckkxall(npointall))
   allocate (ckkyall(npointall))
   allocate (ckkzall(npointall))
   if (id == 0) then
      do i = 1, npointall
         read (90, *) cxall(i), cyall(i), czall(i)
         write (31, *) cxall(i), cyall(i), czall(i)
      end do
      print *, '==================================='
      print *, ' '
      close (90)
   end if
   call mpi_bcast(cxall, npointall, MPI_REAL8, 0, mpi_comm_world, ierr)
   call mpi_bcast(cyall, npointall, MPI_REAL8, 0, mpi_comm_world, ierr)
   call mpi_bcast(czall, npointall, MPI_REAL8, 0, mpi_comm_world, ierr)

   vorx = 0.d0
   vory = 0.d0
   vorz = 0.d0
   !phiv = 0.d0
   lengthtot = 0.d0
   writhetot = 0.d0

   start_time_total = MPI_WTIME()

   starttime = MPI_WTIME()
   !do this in id=0 ??
   npointpast(1) = 0
   do i_line = 1, nline

      npoint = npointlist(i_line)

      allocate (cx(npoint))
      allocate (cy(npoint))
      allocate (cz(npoint))
      allocate (cze(npoint))
      allocate (ckx(npoint))
      allocate (cky(npoint))
      allocate (ckz(npoint))
      allocate (ckkx(npoint))
      allocate (ckky(npoint))
      allocate (ckkz(npoint))

      if (i_line > 1) then
         npointpast(i_line) = npointpast(i_line - 1) + npointlist(i_line - 1)
      end if
      do i = 1, npoint
         cx(i) = cxall(npointpast(i_line) + i)
         cy(i) = cyall(npointpast(i_line) + i)
         cz(i) = czall(npointpast(i_line) + i)
      end do

        !!!!!!get cze: initial zeta
      cze(1) = 0.0d0
      do i = 2, npoint
         cdx = abs(cx(i) - cx(i - 1))
         cdy = abs(cy(i) - cy(i - 1))
         cdz = abs(cz(i) - cz(i - 1))
         if (cdx > pi) cdx = cdx - pi2
         if (cdy > pi) cdy = cdy - pi2
         if (cdz > pi) cdz = cdz - pi2
         dist = dsqrt(cdx**2.0d0 + cdy**2.0d0 + cdz**2.0d0)
         cze(i) = cze(i - 1) + dist
      end do
      cdx = abs(cx(1) - cx(npoint))
      cdy = abs(cy(1) - cy(npoint))
      cdz = abs(cz(1) - cz(npoint))
      if (cdx > pi) cdx = cdx - pi2
      if (cdy > pi) cdy = cdy - pi2
      if (cdz > pi) cdz = cdz - pi2
      dist = dsqrt(cdx**2.0d0 + cdy**2.0d0 + cdz**2.0d0)
      Len_Disc(i_line) = (cze(npoint) + dist)/pi2
      cze = cze*pi2/(cze(npoint) + dist)

      !!!!!!!!get ckx cky ckz
      do i = 1, npoint
         if (i == 1) then
            cdx = cx(2) - cx(npoint)
            cdy = cy(2) - cy(npoint)
            cdz = cz(2) - cz(npoint)
            dist = cze(2) - cze(npoint) + pi2
         else if (i == npoint) then
            cdx = cx(1) - cx(npoint - 1)
            cdy = cy(1) - cy(npoint - 1)
            cdz = cz(1) - cz(npoint - 1)
            dist = cze(1) - cze(npoint - 1) + pi2
         else
            cdx = cx(i + 1) - cx(i - 1)
            cdy = cy(i + 1) - cy(i - 1)
            cdz = cz(i + 1) - cz(i - 1)
            dist = cze(i + 1) - cze(i - 1)
         end if
         if (cdx > pi) cdx = cdx - pi2
         if (cdy > pi) cdy = cdy - pi2
         if (cdz > pi) cdz = cdz - pi2
         if (cdx < -pi) cdx = cdx + pi2
         if (cdy < -pi) cdy = cdy + pi2
         if (cdz < -pi) cdz = cdz + pi2
         ckx(i) = cdx/(dist*Len_Disc(i_line))!temp
         cky(i) = cdy/(dist*Len_Disc(i_line))
         ckz(i) = cdz/(dist*Len_Disc(i_line))
      end do

      !!!!!!!!get ckkx ckky ckkz
      do i = 1, npoint
         if (i == 1) then
            cdx = ckx(2) - ckx(npoint)
            cdy = cky(2) - cky(npoint)
            cdz = ckz(2) - ckz(npoint)
            dist = cze(2) - cze(npoint) + pi2
         else if (i == npoint) then
            cdx = ckx(1) - ckx(npoint - 1)
            cdy = cky(1) - cky(npoint - 1)
            cdz = ckz(1) - ckz(npoint - 1)
            dist = cze(1) - cze(npoint - 1) + pi2
         else
            cdx = ckx(i + 1) - ckx(i - 1)
            cdy = cky(i + 1) - cky(i - 1)
            cdz = ckz(i + 1) - ckz(i - 1)
            dist = cze(i + 1) - cze(i - 1)
         end if
         ckkx(i) = cdx/(dist*Len_Disc(i_line))!temp
         ckky(i) = cdy/(dist*Len_Disc(i_line))
         ckkz(i) = cdz/(dist*Len_Disc(i_line))
      end do

      do i = 1, npoint
         czeall(npointpast(i_line) + i) = cze(i)
         ckxall(npointpast(i_line) + i) = ckx(i)
         ckyall(npointpast(i_line) + i) = cky(i)
         ckzall(npointpast(i_line) + i) = ckz(i)
         ckkxall(npointpast(i_line) + i) = ckkx(i)
         ckkyall(npointpast(i_line) + i) = ckky(i)
         ckkzall(npointpast(i_line) + i) = ckkz(i)
      end do

      length(i_line) = Len_Disc(i_line)*pi2 !temp !bug
      lengthtot = lengthtot + length(i_line)

      deallocate (cx)
      deallocate (cy)
      deallocate (cz)
      deallocate (cze)
      deallocate (ckx)
      deallocate (cky)
      deallocate (ckz)
      deallocate (ckkx)
      deallocate (ckky)
      deallocate (ckkz)

   end do
   endtime = MPI_WTIME()
   if (id == 0) then
      print *, 'Initialization time is', endtime - starttime
      open (192, file='./output/info.dat ', status='unknown')
      write (192, *) 'Initialization time is', endtime - starttime
   end if

   if (id == 0) then
      write (192, *) '=============================='
      write (192, *) '   [i] [Section count] [Spline time] [Construction time]'
   end if
   !mark0
   do i_line = line_s, line_e

      if (id == 0) then
         print *, '-----[ i_line = ', i_line, ' ]-----'
         write (192, '(I6)', advance='no') i_line
      end if
      starttime = MPI_WTIME()

      Gamma = Gamma_line(i_line)
      sigma_min = sigma_line(i_line)
      sigma_max = sigma_ratio*sigma_min
      n_T = n_T_line(i_line)

      !if (id == 0) print *, i_line, Gamma, sigma_min, sigma_max

      n_sec = npointlist(i_line)

      !mind cze_pl when parallel, 2pi?
      allocate (cze_pl(n_sec + 2))
      allocate (cx_pl(n_sec + 2))
      allocate (cy_pl(n_sec + 2))
      allocate (cz_pl(n_sec + 2))
      allocate (ckx_pl(n_sec + 2))
      allocate (cky_pl(n_sec + 2))
      allocate (ckz_pl(n_sec + 2))
      allocate (ckkx_pl(n_sec + 2))
      allocate (ckky_pl(n_sec + 2))
      allocate (ckkz_pl(n_sec + 2))
      allocate (kappa_list(n_sec + 2))

      !mind
      allocate (coe_list(n_sec + 1, 0:5, 3))

      do i = 1, n_sec
         cx_pl(i) = cxall(npointpast(i_line) + i)
         cy_pl(i) = cyall(npointpast(i_line) + i)
         cz_pl(i) = czall(npointpast(i_line) + i)
         cze_pl(i) = czeall(npointpast(i_line) + i)
         ckx_pl(i) = ckxall(npointpast(i_line) + i)
         cky_pl(i) = ckyall(npointpast(i_line) + i)
         ckz_pl(i) = ckzall(npointpast(i_line) + i)
         ckkx_pl(i) = ckkxall(npointpast(i_line) + i)
         ckky_pl(i) = ckkyall(npointpast(i_line) + i)
         ckkz_pl(i) = ckkzall(npointpast(i_line) + i)
      end do

      !Periodic Boundary
      cx_pl(n_sec + 1) = cxall(npointpast(i_line) + 1)
      cy_pl(n_sec + 1) = cyall(npointpast(i_line) + 1)
      cz_pl(n_sec + 1) = czall(npointpast(i_line) + 1)
      cze_pl(n_sec + 1) = pi2 !mind when parallel, consider sigmaFunc
      ckx_pl(n_sec + 1) = ckxall(npointpast(i_line) + 1)
      cky_pl(n_sec + 1) = ckyall(npointpast(i_line) + 1)
      ckz_pl(n_sec + 1) = ckzall(npointpast(i_line) + 1)
      ckkx_pl(n_sec + 1) = ckkxall(npointpast(i_line) + 1)
      ckky_pl(n_sec + 1) = ckkyall(npointpast(i_line) + 1)
      ckkz_pl(n_sec + 1) = ckkzall(npointpast(i_line) + 1)

      !Periodic Boundary
      cx_pl(n_sec + 2) = cxall(npointpast(i_line) + 2)
      cy_pl(n_sec + 2) = cyall(npointpast(i_line) + 2)
      cz_pl(n_sec + 2) = czall(npointpast(i_line) + 2)
      cze_pl(n_sec + 2) = pi2 + czeall(npointpast(i_line) + 2) !mind when parallel, consider sigmaFunc
      ckx_pl(n_sec + 2) = ckxall(npointpast(i_line) + 2)
      cky_pl(n_sec + 2) = ckyall(npointpast(i_line) + 2)
      ckz_pl(n_sec + 2) = ckzall(npointpast(i_line) + 2)
      ckkx_pl(n_sec + 2) = ckkxall(npointpast(i_line) + 2)
      ckky_pl(n_sec + 2) = ckkyall(npointpast(i_line) + 2)
      ckkz_pl(n_sec + 2) = ckkzall(npointpast(i_line) + 2)

      call curve_spline_5th(coe_list, n_sec, cze_pl*Len_Disc(i_line), cx_pl, cy_pl, cz_pl, ckx_pl, cky_pl, ckz_pl, ckkx_pl, ckky_pl, ckkz_pl)

      allocate (ntf_list(n_sec + 1))
      allocate (dtf_list(n_sec + 1))

      ! do i = 1, n_sec + 2
      !    dczeta(1) = ckx_pl(i)
      !    dczeta(2) = cky_pl(i)
      !    dczeta(3) = ckz_pl(i)
      !    ddczeta(1) = ckkx_pl(i)
      !    ddczeta(2) = ckky_pl(i)
      !    ddczeta(3) = ckkz_pl(i)
      !    call cross_zeta(dczeta, ddczeta, dcddczeta)
      !    call norm_zeta(dczeta, ndczeta)
      !    call norm_zeta(dcddczeta, ndcddczeta)
      !    kappa_list(i) = ndcddczeta/(ndczeta**3.d0)
      ! end do

      allocate (dtf_list_pre(n_sec + 1))
      allocate (s_list_pre(n_sec + 1, ntf_pre + 1))
      allocate (ci_list_pre(n_sec + 1, ntf_pre + 1, 3))
      allocate (dci_list_pre(n_sec, ntf_pre + 1, 3))
      allocate (ndci_list_pre(n_sec, ntf_pre + 1))

      if (id == 0) then
         open (519, file='./output/test_data.dat', status='unknown')
      end if
      do t_i = 1, n_sec + 1
         coe_3(:, :) = coe_list(t_i, :, :)
         dtf_list_pre(t_i) = (cze_pl(t_i + 1) - cze_pl(t_i))/(1.0d0*ntf_pre)
         theta = 0
         do t_j = 1, ntf_pre + 1 !s_list(t_i, ntf + 1) = s_list(t_i+1,1)
            s_list_pre(t_i, t_j) = cze_pl(t_i) + dtf_list_pre(t_i)*(t_j - 1) !not arc length
            dzeta = min(dzeta_sup, 1.0d0*dtf_list_pre(t_i))
            call d1_coo_curve_spline_5th((s_list_pre(t_i, t_j) - cze_pl(t_i))*Len_Disc(i_line), dczeta, coe_3)
            !call d1_coo_curve_ip((s_list_pre(t_i, t_j) - cze_pl(t_i))*Len_Disc(i_line), dczeta, dzeta*Len_Disc(i_line), coe_3)
            call norm_zeta(dczeta, ndczeta)
            call d2_coo_curve_spline_5th((s_list_pre(t_i, t_j) - cze_pl(t_i))*Len_Disc(i_line), ddczeta, coe_3)
            !call d2_coo_curve_ip((s_list_pre(t_i, t_j) - cze_pl(t_i))*Len_Disc(i_line), ddczeta, dzeta*Len_Disc(i_line), coe_3)
            call cross_zeta(dczeta, ddczeta, dcddczeta)
            call norm_zeta(dcddczeta, ndcddczeta)
            theta = theta + ndcddczeta/(ndczeta**2.d0)*dtf_list_pre(t_i)*Len_Disc(i_line)
         end do
         ntf_list(t_i) = ceiling(ntf_times*theta/pi2) + 1
         dtf_list(t_i) = (cze_pl(t_i + 1) - cze_pl(t_i))/(1.0d0*ntf_list(t_i))
         ! if (dtf_list(t_i)*Len_Disc(i_line) < dx/10.0d0) then !temp !think later
         !    dtf_list(t_i) = dx/10.0d0/Len_Disc(i_line)
         !    ntf_list(t_i) = ceiling((cze_pl(t_i + 1) - cze_pl(t_i))/dtf_list(t_i))
         !    dtf_list(t_i) = (cze_pl(t_i + 1) - cze_pl(t_i))/(1.0d0*ntf_list(t_i))
         ! end if

         !ntf_list(t_i) = 10 !temp
         !dtf_list(t_i) = (cze_pl(t_i + 1) - cze_pl(t_i))/(1.0d0*ntf_list(t_i))
         ! if (id == 0) then
         !    write (519, *) '----------[t_i=', t_i, ']----------'
         !    write (519, *) 'theta=', theta
         !    write (519, *) 'ntf=', ntf_list(t_i)
         !    write (519, *) 'dtf*Len=', dtf_list(t_i)*Len_Disc(i_line)
         ! end if
      end do

      if (id == 0) then
         if (if_spline_out == 1) then
            id1 = i_line/10000
            id2 = i_line/1000 - (i_line/10000)*10
            id3 = i_line/100 - (i_line/1000)*10
            id4 = i_line/10 - (i_line/100)*10
            id5 = i_line - (i_line/10)*10
            open (312, file='./output/centerline_spline/centerline_spline'//char(id1 + 48)//char(id2 + 48)//char(id3 + 48)//char(id4 + 48)//char(id5 + 48)//'.dat', status='unknown')
            write (312, *) 1
            write (312, *) sum(ntf_list(1:n_sec))
         end if
      end if

      ntf_max = maxval(ntf_list)
      if (id == 0) then
         write (519, *) 'ntf_average=', sum(ntf_list)/(n_sec + 1)
         write (519, *) 'ntf_max=', ntf_max
      end if

      ! do i = 1, n_sec + 1
      !    kappa_sup = max(kappa_list(i), kappa_list(i + 1)) ! ratio ? think later
      !    dist = (cze_pl(i + 1) - cze_pl(i))*Len_Disc(i_line)
      !    theta = 2.0d0*dasin(dist*kappa_sup/2.0d0)

      !    ntf =
      ! end do

      deallocate (dtf_list_pre)
      deallocate (s_list_pre)
      deallocate (ci_list_pre)
      deallocate (dci_list_pre)
      deallocate (ndci_list_pre)

      pdx = 0.0d0 ! ??

      !think again, mind the size
      allocate (ci_list2(n_sec + 1, ntf_max + 1, 3))
      allocate (s_list2(n_sec + 1, ntf_max + 1))
      allocate (dci_list2(n_sec, ntf_max + 1, 3))
      allocate (tci_list2(n_sec, ntf_max + 1, 3))
      allocate (ndci_list2(n_sec, ntf_max + 1))

      ! do t_i = 1, n_sec + 1
      !    !think again
      !    ntf_list(t_i) = ntf_max
      !    dtf_list(t_i) = (cze_pl(t_i + 1) - cze_pl(t_i))/(1.0d0*ntf_list(t_i))
      ! end do
      do t_i = 1, n_sec + 1
         coe_3(:, :) = coe_list(t_i, :, :)
         do t_j = 1, ntf_list(t_i) + 1 !s_list2(t_i, ntf_list(t_i) + 1) = s_list2(t_i+1,1)
            s_list2(t_i, t_j) = cze_pl(t_i) + dtf_list(t_i)*(t_j - 1) !mind !not arc length, 1 to 1 mapping
            !mind
            call coo_curve_spline_5th((s_list2(t_i, t_j) - cze_pl(t_i))*Len_Disc(i_line), ci_list2(t_i, t_j, :), coe_3)
         end do
      end do
      !think later
      ! t_i = n_sec + 1
      ! do t_j = 1,2
      ! enddo

      !mind !mark v5_3_2
      !==deal with the bug when ci is on the mesh
      do t_i = 1, n_sec
         ci_list2(t_i, ntf_list(t_i) + 1, :) = ci_list2(t_i + 1, 1, :)
      end do
      !==deal with the bug when ci is on the mesh

      if (id == 0) then
         if (if_spline_out == 1) then
            do t_i = 1, n_sec
               do t_j = 1, ntf_list(t_i)
                  write (312, *) ci_list2(t_i, t_j, 1), ci_list2(t_i, t_j, 2), ci_list2(t_i, t_j, 3)
               end do
            end do
            close (312)
         end if
      end if

      dtf_min = minval(dtf_list) !think later
      if (id == 0) print *, 'dtf_min =', dtf_min

      do t_i = 1, n_sec
         do t_j = 1, ntf_list(t_i) !mind ntf_list(t_i)+1
            call cminus(ci_list2(t_i, t_j + 1, :), ci_list2(t_i, t_j, :), dci_list2(t_i, t_j, :))
            call norm_zeta(dci_list2(t_i, t_j, :), ndci_list2(t_i, t_j))
            tci_list2(t_i, t_j, :) = dci_list2(t_i, t_j, :)/ndci_list2(t_i, t_j)
         end do
         !mind mark2
         t_j = ntf_list(t_i) + 1
         call cminus(ci_list2(t_i + 1, 2, :), ci_list2(t_i + 1, 1, :), dci_list2(t_i, t_j, :))
         call norm_zeta(dci_list2(t_i, t_j, :), ndci_list2(t_i, t_j))
         tci_list2(t_i, t_j, :) = dci_list2(t_i, t_j, :)/ndci_list2(t_i, t_j) ! = tci_list2(t_i+1, 1, :)
      end do

      endtime = MPI_WTIME()
      if (id == 0) then
         print *, 'Spline interpolation time =', endtime - starttime
         write (192, '(I16)', advance='no') n_sec
         write (192, '(ES18.3)', advance='no') endtime - starttime
      end if

      !think again
      starttime = MPI_WTIME()
      do t_i = 1, n_sec

         coe_3(:, :) = coe_list(t_i, :, :)

         do t_j = 1, ntf_list(t_i)
            s0 = s_list2(t_i, t_j)
            s1 = s_list2(t_i, t_j + 1)
            ci0 = ci_list2(t_i, t_j, :)
            ci1 = ci_list2(t_i, t_j + 1, :)
            dci0 = dci_list2(t_i, t_j, :)
            dci1 = dci_list2(t_i, t_j + 1, :)
            ndci0 = ndci_list2(t_i, t_j)
            ndci1 = ndci_list2(t_i, t_j + 1)
            tci0 = tci_list2(t_i, t_j, :)
            tci1 = tci_list2(t_i, t_j + 1, :)

            !mark1
            do ijk_co = 1, 3
               if (ci1(ijk_co) - ci0(ijk_co) > pi) then
                  ci0(ijk_co) = ci0(ijk_co) + pi2
               elseif (ci1(ijk_co) - ci0(ijk_co) < -pi) then
                  ci1(ijk_co) = ci1(ijk_co) + pi2
               end if
            end do
            !  0 < (xip- xstart) < 2*pi+dtf/(2*pi)*length < 3*pi
            !mark1

            call sigma_func(s0, length(i_line), sigma, sigma_max, sigma_min, n_T)
            Rtube = Rtube_ratio*sigma
            call sigma_func(s1, length(i_line), sigma, sigma_max, sigma_min, n_T)
            Rtube = max(Rtube, Rtube_ratio*sigma)
            i_min = floor((min(ci0(1), ci1(1)) - xstart - Rtube)/dx + 1)  ! mind the meshx_plus(i, :, :) = (i - 1.d0)*dx + xstart
            i_max = ceiling((max(ci0(1), ci1(1)) - xstart + Rtube)/dx + 1)
            j_min = floor((min(ci0(2), ci1(2)) - ystart - Rtube)/dy + 1)
            j_max = ceiling((max(ci0(2), ci1(2)) - ystart + Rtube)/dy + 1)
            k_min = floor((min(ci0(3), ci1(3)) - zstart - Rtube)/dz + 1)
            k_max = ceiling((max(ci0(3), ci1(3)) - zstart + Rtube)/dz + 1)

            do k_plus = k_min, k_max
               do j_plus = j_min, j_max
                  do i_plus = i_min, i_max
                     k = (mod(k_plus + nz - 1, nz) + 1) - nzp*id !mind nzp
                     if (k > 0 .and. k <= nzp) then
                        !see !mark v5_3_2 !mind there may be wrong when ci is on the mesh
                        t0 = (meshx_plus(i_plus) - ci0(1))*tci0(1) + &
                        &(meshy_plus(j_plus) - ci0(2))*tci0(2) + (meshz_plus(k_plus) - ci0(3))*tci0(3)
                        if (t0 >= 0.d0) then

                           t1 = (meshx_plus(i_plus) - ci1(1))*tci1(1) + &
                              &(meshy_plus(j_plus) - ci1(2))*tci1(2) + &
                              &(meshz_plus(k_plus) - ci1(3))*tci1(3)
                           if (t1 < 0.d0) then

                              rho = dsqrt((meshx_plus(i_plus) - ci0(1))**2 + &
                              &(meshy_plus(j_plus) - ci0(2))**2 + &
                              &(meshz_plus(k_plus) - ci0(3))**2 - t0**2)
                              if (rho < Rtube) then
                                 t1 = (meshx_plus(i_plus) - ci1(1))*tci0(1) + &
                                 &(meshy_plus(j_plus) - ci1(2))*tci0(2) + &
                                 &(meshz_plus(k_plus) - ci1(3))*tci0(3)
                                 zeta = (s1*t0 - s0*t1)/ndci0 ! mind the bug!!!! zeta = (s0*t0 - s1*t1)/ndci0
                                 zeta_sec = zeta - s_list2(t_i, 1) !mind

                                 call coo_curve_spline_5th(zeta_sec*Len_Disc(i_line), czeta, coe_3)
                                 dzeta = min(dzeta_sup, 1.0d0*dtf_list(t_i))
                                 !think later, dzeta can not be too small, otherwise the calculation is difficult
                                 !dzeta : 1.0d-4 ~ 1.0d-5 ?
                                 call d1_coo_curve_spline_5th(zeta_sec*Len_Disc(i_line), dczeta, coe_3)
                                 !call d1_coo_curve_ip(zeta_sec*Len_Disc(i_line), dczeta, dzeta*Len_Disc(i_line), coe_3)
                                 call norm_zeta(dczeta, ndczeta)
                                 Tzeta = dczeta/ndczeta
                                 call d2_coo_curve_spline_5th(zeta_sec*Len_Disc(i_line), ddczeta, coe_3)
                                 !call d2_coo_curve_ip(zeta_sec*Len_Disc(i_line), ddczeta, dzeta*Len_Disc(i_line), coe_3)
                                 call cross_zeta(dczeta, ddczeta, dcddczeta)
                                 call norm_zeta(dcddczeta, ndcddczeta)

                                 meshxyz(1) = meshx_plus(i_plus)
                                 meshxyz(2) = meshy_plus(j_plus)
                                 meshxyz(3) = meshz_plus(k_plus)
                                 call cminus(meshxyz, czeta, xczeta)
                                 call norm_zeta(xczeta, nxczeta)
                                 rho = nxczeta
                                 call sigma_func(zeta, length(i_line), sigma, sigma_max, sigma_min, n_T)
                                 tp = dexp(-((rho - pdx)**2.d0)/2.d0/(sigma**2.0d0))/(pi2*(sigma**2.0d0))
                                 ft_delta = -dexp(-((Rtube)**2.d0)/2.d0/(sigma**2.0d0))/(pi2*(sigma**2.0d0))
                                 tp = tp + ft_delta
                                 !think later in parallel
                                 !call dsigmads_func_spline_5th(zeta, zeta_sec, Len_Disc(i_line), length(i_line), dsigmads, sigma_max, sigma_min, dzeta, coe_3)

                    call dsigmadzeta(zeta, zeta_sec, Len_Disc(i_line), length(i_line), dsigmads, sigma_max, sigma_min, dzeta, coe_3)
                    !!!bug zeta is not arc length!!!
                                 dsigmads = dsigmads/ndczeta/Len_Disc(i_line)
                                 if (ndcddczeta < 1.0d-8) then  !mark v5_3_1
                                    e_rho = xczeta/(nxczeta + 1.0d-15) ! = 0 is acceptable, which means rho = 0
                                    call cross_zeta(e_rho, Tzeta, e_th)

                                    i = mod(i_plus + nx - 1, nx) + 1
                                    j = mod(j_plus + ny - 1, ny) + 1
                                    ! k has already been calculated
                                    !phiv(i, j, k) = phiv(i, j, k) + (dexp(-((rho - pdx)**2.d0)/2.d0/(sigma**2.0d0)))
                                    !call etaFunc(zeta, length(i_line), phiv(i, j, k), eta) !temp
                                    eta = eta0 !temp
                               vorx(i, j, k) = vorx(i, j, k) + tp*Gamma*(Tzeta(1) + eta*rho*e_th(1) + (rho*dsigmads/sigma)*e_rho(1))
                               vory(i, j, k) = vory(i, j, k) + tp*Gamma*(Tzeta(2) + eta*rho*e_th(2) + (rho*dsigmads/sigma)*e_rho(2))
                               vorz(i, j, k) = vorz(i, j, k) + tp*Gamma*(Tzeta(3) + eta*rho*e_th(3) + (rho*dsigmads/sigma)*e_rho(3))

                                 else
                                    Bzeta = dcddczeta/ndcddczeta
                                    call cross_zeta(Bzeta, Tzeta, Nzeta)
                                    kappazeta = ndcddczeta/(ndczeta**3.d0)
                                    call dot_zeta(xczeta, Nzeta, t4)
                                    costh = t4/(nxczeta + 1.0d-15)
                                    call dot_zeta(xczeta, Bzeta, t4)
                                    sinth = t4/(nxczeta + 1.0d-15)
                                    t5 = -sinth*rho/(1.d0 - kappazeta*rho*costh)
                                    t6 = costh*rho/(1.d0 - kappazeta*rho*costh)
                                    t7 = costh*rho*dsigmads/sigma/(1.d0 - kappazeta*rho*costh)
                                    t8 = sinth*rho*dsigmads/sigma/(1.d0 - kappazeta*rho*costh)

                                    i = mod(i_plus + nx - 1, nx) + 1
                                    j = mod(j_plus + ny - 1, ny) + 1
                                    ! k has already been calculated
                                    !phiv(i, j, k) = phiv(i, j, k) + (dexp(-((rho - pdx)**2.d0)/2.d0/(sigma**2.0d0)))
                                    !call etaFunc(zeta, length(i_line), phiv(i, j, k), eta)
                                    eta = eta0 !temp
                                    vorx(i, j, k) = vorx(i, j, k) + tp*Gamma*(Tzeta(1) + &
                                    &eta*(t5*Nzeta(1) + t6*Bzeta(1)) + (t7*Nzeta(1) + t8*Bzeta(1)))
                                    vory(i, j, k) = vory(i, j, k) + tp*Gamma*(Tzeta(2) + &
                                    &eta*(t5*Nzeta(2) + t6*Bzeta(2)) + (t7*Nzeta(2) + t8*Bzeta(2)))
                                    vorz(i, j, k) = vorz(i, j, k) + tp*Gamma*(Tzeta(3) + &
                                    &eta*(t5*Nzeta(3) + t6*Bzeta(3)) + (t7*Nzeta(3) + t8*Bzeta(3)))
                                 end if

                              end if
                           end if
                        end if
                     end if
                  end do
               end do
            end do

         end do
      end do

      deallocate (ntf_list)
      deallocate (dtf_list)
      deallocate (kappa_list)
      deallocate (cx_pl)
      deallocate (cy_pl)
      deallocate (cz_pl)
      deallocate (cze_pl)
      deallocate (ckx_pl)
      deallocate (cky_pl)
      deallocate (ckz_pl)
      deallocate (ckkx_pl)
      deallocate (ckky_pl)
      deallocate (ckkz_pl)
      deallocate (coe_list)
      deallocate (ndci_list2)
      deallocate (s_list2)
      deallocate (ci_list2)
      deallocate (dci_list2)
      deallocate (tci_list2)

      if (if_filter_line(i_line)) then
         starttime_part = MPI_WTIME()
         k_c = ceiling(4.0d0/sigma_min) !temp !consider later
     call filter_spec(k_c, vorx, vory, vorz, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, id, nproc)
         !allocate (spec(nx, ny, nzp))
         endtime_part = MPI_WTIME()
         if (id == 0) then
            print *, 'Filtering time =', endtime_part - starttime_part
            write (192, *) 'Filtering time =', endtime_part - starttime_part
         end if
      end if

      endtime = MPI_WTIME()
      if (id == 0) then
         print *, 'Number of section =', n_sec
         print *, 'Construction time =', endtime - starttime
         !write (192, '(I12)', advance='no') n_sec
         write (192, '(ES18.3)', advance='no') endtime - starttime
         write (192, *)
      end if
   end do
   end_time_total = MPI_WTIME()
   if (id == 0) then
      print *, '=============================='
      print *, 'Total construction time =', end_time_total - start_time_total
      write (192, *) '=============================='
      write (192, *) 'Total construction time =', end_time_total - start_time_total
   end if

   if (id == 0) then
      close (31)
   end if

   deallocate (cxall)
   deallocate (cyall)
   deallocate (czall)
   deallocate (czeall)
   deallocate (ckxall)
   deallocate (ckyall)
   deallocate (ckzall)
   deallocate (ckkxall)
   deallocate (ckkyall)
   deallocate (ckkzall)
   deallocate (npointlist)
   deallocate (npointpast)

   if (id == 0) then
      open (499, file='./output/stat/spectrum_vor_temp.dat', status='unknown')
   end if
   call spectrum(nx, ny, nzp, vorx, vory, vorz, 499, k2, id, nproc, planxf, planyf, planzf)
   if (id == 0) then
      close (499)
   end if

   starttime = MPI_WTIME()
   !==vel==
   if (id == 0) then
      write (192, *) 'velocity calculation begin'
   end if
   allocate (spec_wx(nx, ny, nzp))
   allocate (spec_wy(nx, ny, nzp))
   allocate (spec_wz(nx, ny, nzp))
   allocate (spec_ux(nx, ny, nzp))
   allocate (spec_uy(nx, ny, nzp))
   allocate (spec_uz(nx, ny, nzp))
   call fourier_forward(vorx, spec_wx, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(vory, spec_wy, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(vorz, spec_wz, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            spec_ux(i, j, k) = dcmplx(0.0d0, 1.0d0)*(ky(j)*spec_wz(i, j, k) - kz(k)*spec_wy(i, j, k))/k2(i, j, k)
            spec_uy(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kz(k)*spec_wx(i, j, k) - kx(i)*spec_wz(i, j, k))/k2(i, j, k)
            spec_uz(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kx(i)*spec_wy(i, j, k) - ky(j)*spec_wx(i, j, k))/k2(i, j, k)
         end do
      end do
   end do
   if (id == 0) then
      spec_ux(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
      spec_uy(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
      spec_uz(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
   end if
   deallocate (spec_wx)
   deallocate (spec_wy)
   deallocate (spec_wz)
   allocate (velx(nx, ny, nzp))
   allocate (vely(nx, ny, nzp))
   allocate (velz(nx, ny, nzp))
   allocate (disp_field(nx, ny, nzp))
   call fourier_backward(velx, spec_ux, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(vely, spec_uy, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(velz, spec_uz, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   deallocate (spec_ux)
   deallocate (spec_uy)
   deallocate (spec_uz)
   if (id == 0) then
      write (192, *) 'velocity calculation done'
   end if
   !==vel==

   !==up->1==
   !call cal_helicity(vorx, vory, vorz, velx, vely, velz, nx, ny, nz, nzp, helicity, dv, id)
   call cal_helicity(velx, vely, velz, velx, vely, velz, nx, ny, nz, nzp, Etot, dv, id) !cdot
   Etot = Etot/2.d0/lx/ly/lz
   up = sqrt(2*Etot/3)
   call mpi_bcast(up, 1, MPI_REAL8, 0, mpi_comm_world, ierr)
   if (id == 0) then
      write (192, *) 'Etot_temp =', Etot
      write (192, *) 'up_temp =', up
   end if
   vorx = vorx/up
   vory = vory/up
   vorz = vorz/up
   velx = velx/up
   vely = vely/up
   velz = velz/up
   call cal_helicity(vorx, vory, vorz, vorx, vory, vorz, nx, ny, nz, nzp, Diss, dv, id)
   Enstrophy = Diss/2.0d0/lx/ly/lz
   ! Diss = Diss/pi2/pi2/pi2
   !==up->1==

   ! call cross_vector(vorx, vory, vorz, velx, vely, velz, -1, nx, ny, nzp,&
   !     &kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, id, nproc)

   ! call wall_phiv(nx, ny, nzp, phiv, igst, id_l, id_r)
   ! call wall_phiv(nx, ny, nzp, phiv, igst, id_l, id_r)
   ! call wall_phiv(nx, ny, nzp, phiv, igst, id_l, id_r)
   ! call cal_deviation(vorx, vory, vorz, phiv, nx, ny, nz, nzp, deviation, dx, dy, dz, id, igst)

   if (vorticity_output == 1) then
      r_d = nx/n_d !temp
      !---
      allocate (data_box(nx/r_d, ny/r_d, nzp/r_d, 2))
      do k = 1, nzp/r_d
         do j = 1, ny/r_d
            do i = 1, nx/r_d
               data_box(i, j, k, 1) = sqrt(vorx(i*r_d, j*r_d, k*r_d)*vorx(i*r_d, j*r_d, k*r_d) + vory(i*r_d, j*r_d, k*r_d)*vory(i*r_d, j*r_d, k*r_d) + vorz(i*r_d, j*r_d, k*r_d)*vorz(i*r_d, j*r_d, k*r_d))
               data_box(i, j, k, 2) = vorx(i*r_d, j*r_d, k*r_d)*velx(i*r_d, j*r_d, k*r_d) + vory(i*r_d, j*r_d, k*r_d)*vely(i*r_d, j*r_d, k*r_d) + vorz(i*r_d, j*r_d, k*r_d)*velz(i*r_d, j*r_d, k*r_d)            
            end do
         end do
      end do
      name = './output/field_wh.dat'
      allocate (varname(2))
      varname(1) = 'w'
      varname(2) = 'h'
      call output_v(name, varname, nx/r_d, ny/r_d, nzp/r_d, data_box, 2, id, nproc)
      deallocate (data_box)
      deallocate (varname)

      if (u_out == 1) then
         allocate (data_box(nx/r_d, ny/r_d, nzp/r_d, 6))
         do k = 1, nzp/r_d
            do j = 1, ny/r_d
               do i = 1, nx/r_d
                  data_box(i, j, k, 1) = velx(i*r_d, j*r_d, k*r_d)
                  data_box(i, j, k, 2) = vely(i*r_d, j*r_d, k*r_d)
                  data_box(i, j, k, 3) = velz(i*r_d, j*r_d, k*r_d)
                  data_box(i, j, k, 4) = vorx(i*r_d, j*r_d, k*r_d)
                  data_box(i, j, k, 5) = vory(i*r_d, j*r_d, k*r_d)
                  data_box(i, j, k, 6) = vorz(i*r_d, j*r_d, k*r_d)
               end do
            end do
         end do
         name = './output/field_uw.dat'
         allocate (varname(6))
         varname(1) = 'ux'
         varname(2) = 'uy'
         varname(3) = 'uz'
         varname(4) = 'wx'
         varname(5) = 'wy'
         varname(6) = 'wz'
         call output_v(name, varname, nx/r_d, ny/r_d, nzp/r_d, data_box, 6, id, nproc)
         deallocate (data_box)
         deallocate (varname)
      end if

      if (d_out == 1) then
         call dissipation_field(velx,vely,velz,disp_field, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, id, nproc)
         allocate (data_box(nx/r_d, ny/r_d, nzp/r_d, 2))
         do k = 1, nzp/r_d
            do j = 1, ny/r_d
               do i = 1, nx/r_d
                  data_box(i, j, k, 1) = sqrt(vorx(i*r_d, j*r_d, k*r_d)*vorx(i*r_d, j*r_d, k*r_d) + vory(i*r_d, j*r_d, k*r_d)*vory(i*r_d, j*r_d, k*r_d) + vorz(i*r_d, j*r_d, k*r_d)*vorz(i*r_d, j*r_d, k*r_d))
                  data_box(i, j, k, 2) = disp_field(i*r_d, j*r_d, k*r_d)
               end do
            end do
         end do
         name = './output/field_wd.dat'
         allocate (varname(2))
         varname(1) = 'w'
         varname(2) = 'd'
         call output_v(name, varname, nx/r_d, ny/r_d, nzp/r_d, data_box, 2, id, nproc)
         deallocate (data_box)
         deallocate (varname)
      end if

      if (multi_out == 1) then
         r_d = 1
         allocate (data_box(nx/r_d, ny/r_d, nzp/r_d, 2))
         do k = 1, nzp/r_d
            do j = 1, ny/r_d
               do i = 1, nx/r_d
                  data_box(i, j, k, 1) = sqrt(vorx(i*r_d, j*r_d, k*r_d)*vorx(i*r_d, j*r_d, k*r_d) + vory(i*r_d, j*r_d, k*r_d)*vory(i*r_d, j*r_d, k*r_d) + vorz(i*r_d, j*r_d, k*r_d)*vorz(i*r_d, j*r_d, k*r_d))
                  data_box(i, j, k, 2) = vorx(i*r_d, j*r_d, k*r_d)*velx(i*r_d, j*r_d, k*r_d) + vory(i*r_d, j*r_d, k*r_d)*vely(i*r_d, j*r_d, k*r_d) + vorz(i*r_d, j*r_d, k*r_d)*velz(i*r_d, j*r_d, k*r_d)            
               end do
            end do
         end do
         name = './output/w_h_downsample.dat'
         allocate (varname(2))
         varname(1) = 'w'
         varname(2) = 'h'
         call output_v_multi(name, varname, nx/r_d, ny/r_d, nzp/r_d, data_box, 2, id, nproc)
         deallocate (data_box)
         deallocate (varname)
      end if

   end if

   if (field2d_output == 1) then
      if (id == 0) then
         open (1988, file='./output/2d_field.dat ', status='unknown')
         write (1988, *) 'variables="x","y","wx","wy","wz"'
         write (1988, *) 'ZONE I =', nx, ', J =', ny, ', F = POINT'
         do j = 1, ny
            do i = 1, nx
               write (1988, "(7F12.5)") meshx_plus(i), meshy_plus(j), vorx(i, j, 1), vory(i, j, 1), vorz(i, j, 1)
            end do
         end do
         close (1988)
      end if
   end if

   if (id == 0) then
      open (49, file='./output/stat/spectrum.dat', status='unknown')
   end if
   call spectrum(nx, ny, nzp, velx, vely, velz, 49, k2, id, nproc, planxf, planyf, planzf)
   if (id == 0) then
      close (49)
   end if

   if (if_pdf == 1) then
      open (3201, file='./output/stat/pdf_um.dat', status='unknown')
      open (3204, file='./output/stat/pdf_wm.dat', status='unknown')
      allocate (tmp_real4(nx, ny, 3*nzp))
      tmp_real4(:, :, 1:nzp) = velx
      tmp_real4(:, :, nzp + 1:2*nzp) = vely
      tmp_real4(:, :, 2*nzp + 1:3*nzp) = velz
      call export_pdf(0, tmp_real4, nx, ny, 3*nzp, id, nproc, bin_num, 3201)
      tmp_real4(:, :, 1:nzp) = vorx
      tmp_real4(:, :, nzp + 1:2*nzp) = vory
      tmp_real4(:, :, 2*nzp + 1:3*nzp) = vorz
      call export_pdf(0, tmp_real4, nx, ny, 3*nzp, id, nproc, bin_num, 3204)
      deallocate (tmp_real4)
      close (3201)
      close (3202)
      close (3203)
      close (3204)
      close (3205)
      close (3206)
   end if
   if (if_Sp == 1) then
      call Sp_calculate(Sp_order, nx, ny, nzp, velx, vely, velz, id, nproc)
   end if
   endtime = MPI_WTIME()
   if (id == 0) then
      print *, 'Post processing time =', endtime - starttime
      print *, '=============================='
      write (192, *) 'Post processing time =', endtime - starttime
      write (192, *) '=============================='
   end if

   if (id == 0) then
      write (192, *) 'N = ', nx
      write (192, *) 'Total length of curve = ', lengthtot
      !write (192, *) 'Minimum radius of curvature = ', Rkappa
      !write (192, *) 'VSF deviation of vorticity magnitude = ', deviation
      !write (192, *) 'Total helicity = ', helicity
      !write (192, *) 'Total writhe = ', writhetot
      write (192, *) 'Gamma_1 = ', Gamma_line(1)/up
      write (192, *) 'twist eta = ', eta
      write (192, *) 'Etot = ', Etot/up/up
      write (192, *) 'Enstrophy = ', Enstrophy
      ! write (192, *) '(helicity - (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)) / helicity  = ', &
      !     &(helicity - (Gamma**2)*(writhetot + eta*lengthtot/pi2))/helicity
      close (192)
   end if

   call destroy_fft_plans(planxf, planxb, planyf, planyb, planzf, planzb)

   deallocate (Gamma_line)
   deallocate (sigma_line)
   deallocate (n_T_line)
   deallocate (if_filter_line)
   deallocate (length)
   deallocate (Len_Disc)

   deallocate (k2)
   deallocate (kx)
   deallocate (ky)
   deallocate (kz)
   deallocate (Sn)
   deallocate (r)
   deallocate (meshx_plus)
   deallocate (meshy_plus)
   deallocate (meshz_plus)
   deallocate (vorx)
   deallocate (vory)
   deallocate (vorz)
   deallocate (velx)
   deallocate (vely)
   deallocate (velz)
   !deallocate (phiv)
   call MPI_Finalize(ierr)
   !read (*, *)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
