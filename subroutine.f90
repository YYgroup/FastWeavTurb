subroutine cal_pdf(variable, val, pdf_v, nx, ny, nzp, id, nproc, bin_num)
   implicit none
   include 'mpif.h'
   integer(4) :: nx, ny, nzp
   integer(4) :: id, nproc, ierr
   integer(4) :: bin_num
   real(4), dimension(nx, ny, nzp) :: variable
   real(4), dimension(bin_num) :: val, pdf_v, num0
   real(4) :: bin_size, max_varb0, min_varb0, max_varb, min_varb
   integer(4) :: i, j, k, l

   max_varb0 = maxval(variable)
   min_varb0 = minval(variable)
   call MPI_ALLREDUCE(max_varb0, max_varb, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(min_varb0, min_varb, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
   bin_size = (max_varb - min_varb)/real(bin_num)

   do l = 1, bin_num
      val(l) = min_varb + 0.5d0*bin_size + real(l - 1)*bin_size
   end do

   num0 = 0d0
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            do l = 1, bin_num
               if (variable(i, j, k) >= val(l) - 0.5*bin_size .and. variable(i, j, k) < val(l) + 0.5*bin_size) then
                  num0(l) = num0(l) + 1d0
                  exit
               end if
            end do
         end do
      end do
   end do
   call MPI_ALLREDUCE(num0, pdf_v, bin_num, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

   do l = 1, bin_num
      pdf_v(l) = pdf_v(l)/(bin_size*nx*ny*nzp*nproc)
   end do
end subroutine cal_pdf

subroutine export_pdf(time_t, variable, nx, ny, nzp, id, nproc, bin_num, file_num)
   !use temp_array
   implicit none
   include 'mpif.h'
   integer(4) :: nx, ny, nzp
   integer(4) :: id, nproc, ierr
   integer(4) :: bin_num, file_num
   real(4), dimension(bin_num) :: vari, pdf
   real(4), dimension(nx, ny, nzp) :: variable
   real(4) :: time_t
   integer(4) :: i

   call cal_pdf(variable, vari, pdf, nx, ny, nzp, id, nproc, bin_num)
   if (id == 0) then
      ! write (file_num, *) 'variables = "val", "pdf"'
      ! write (file_num, *) 'zone t = "', time_t, '", I = ', bin_num, ', F = point'
      do i = 1, bin_num
         write (file_num, *) vari(i), pdf(i)
      end do
   end if
end subroutine export_pdf

subroutine Sp_calculate(Sp_order, nx, ny, nzp, velx, vely, velz, id, nproc)
   implicit none
   include 'mpif.h'
   real*8 dr, pi, pi2
   integer nx, ny, nzp, nz, id, nproc, ierr, i, j, k
   integer lg2r, lg2L, Sp_order, p, r
   integer, dimension(nproc) :: counts, displs
   real*8, dimension(nx, ny, nzp) :: velx, vely, velz
   real*8, allocatable :: ux(:, :, :), uy(:, :, :), uz(:, :, :)
   real*8, allocatable :: Sp(:, :)
   real*8 Sp_plus, dux, duy, duz
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   dr = pi2/nx
   lg2L = floor(log(real(nx))/log(2.0) + 0.01)
   nz = nzp*nproc
   if (id == 0) then
      allocate (ux(nx, ny, nz))
      allocate (uy(nx, ny, nz))
      allocate (uz(nx, ny, nz))
      allocate (Sp(Sp_order, 0:2*lg2L))
      Sp = 0
   end if
   do i = 1, nproc
      counts(i) = nx*ny*nzp
      displs(i) = (i - 1)*counts(i)
   end do
   call mpi_gatherv(velx, counts(id + 1), MPI_DOUBLE, ux, counts, displs, MPI_DOUBLE, 0, mpi_comm_world, ierr)
   call mpi_gatherv(vely, counts(id + 1), MPI_DOUBLE, uy, counts, displs, MPI_DOUBLE, 0, mpi_comm_world, ierr)
   call mpi_gatherv(velz, counts(id + 1), MPI_DOUBLE, uz, counts, displs, MPI_DOUBLE, 0, mpi_comm_world, ierr)
   if (id == 0) then
      do p = 1, Sp_order

         do lg2r = 0, 2*lg2L
            r = floor(sqrt(2.0)**lg2r + 0.5)
            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx

                     dux = ux(mod(i + r - 1, nx) + 1, j, k) - ux(i, j, k)
                     duy = uy(i, mod(j + r - 1, ny) + 1, k) - uy(i, j, k)
                     duz = uz(i, j, mod(k + r - 1, nz) + 1) - uz(i, j, k)
                     Sp_plus = dux**p + duy**p + duz**p
                     Sp_plus = Sp_plus/nx/ny/nz/3.0d0
                     Sp(p, lg2r) = Sp(p, lg2r) + Sp_plus
                  end do
               end do
            end do
         end do
      end do
      open (3301, file='./output/stat/Sp.dat', status='unknown')
      do lg2r = 0, 2*lg2L
         r = floor(sqrt(2.0)**lg2r + 0.5)
         write (3301, "(ES16.5)", ADVANCE='NO') r*dr
         do p = 1, Sp_order
            write (3301, "(ES16.5)", ADVANCE='NO') Sp(p, lg2r)
         end do
         write (3301, *)
      end do
      close (3301)
   end if
   if (id == 0) then
      deallocate (ux, uy, uz)
      deallocate (Sp)
   end if
end subroutine Sp_calculate

subroutine getStructureFunc_x(Sn, r, norder, velx, nx, ny, nz, nzp, id)
   implicit none
   include 'mpif.h'
   real*8 dr, pi, pi2
   integer i, j, k, m, p, norder
   integer nx, ny, nzp, nz, ierr, id
   real*8, dimension(nx, ny, nzp) :: velx
   real*8, dimension(nx/2 + 1) :: Sn1, Sn, r
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   dr = pi2/nx
   do m = 1, nx/2 + 1
      r(m) = (m - 1)*dr
   end do
   Sn1 = 0.0d0
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            do m = 1, nx/2 + 1
               p = i + (m - 1)
               if (p > nx) p = p - nx
               Sn1(m) = Sn1(m) + (velx(p, j, k) - velx(i, j, k))**norder
            end do
         end do
      end do
   end do

   call MPI_REDUCE(Sn1, Sn, nx/2 + 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

   Sn = Sn/(nx*ny*nz)

end subroutine

subroutine sigma_func(zeta, length, sigma, sigma_max, sigma_min, n_T)
   integer n_T
   real*8 zeta, sigma, sigma_max, sigma_min, pi, pi2, beta, length, omega
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   beta = 1.0d0
   omega = floor(length/pi + 0.5d0)*n_T
   ! omega = n_T
   sigma = sigma_min + (sigma_max - sigma_min)*(dsin(omega*(zeta)) + 1.0d0)/2.0d0
end subroutine

subroutine etaFunc(zeta, length, phivijk, eta)
   real*8 eta, phivijk, zeta, sigma, Tw, pi, pi2, length
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   Tw = 0.0d0
   !sigma = 1.0d0 / 1.6d0 / dsqrt(pi2)
   !if(zeta<0.0d0) then
   !    eta=(dsqrt(pi2)*Tw)/sigma * dexp( - (zeta ** 2.d0) / 2.d0 / (sigma ** 2.0d0))
   !else
   !    eta=(dsqrt(pi2)*Tw)/sigma * dexp( - ((pi2-zeta) ** 2.d0) / 2.d0 / (sigma ** 2.0d0))
   !endif
   eta = 0.0d0
   !eta=10.0d0*(3.0d0*phivijk-1.0d0)*dsin((zeta+pi)/2.0d0)
   !eta = Tw*phivijk
   !eta = Tw*phivijk*dcos(zeta/2.0d0)
   !eta = Tw*dsin(pi*phivijk)*dcos(zeta/2.0d0)
end subroutine

subroutine getsplinecoe(a, b, c, d, zeta1, zeta2, x1, x2, dx1, dx2)
   real*8 a, b, c, d, zeta1, zeta2, x1, x2, x1p, x2p, dx1, dx2, pi, pi2
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   x1p = x1
   x2p = x2
   if (x1 - x2 > pi) then
      x2p = x2 + pi2
   elseif (x1 - x2 < -pi) then
      x1p = x1 + pi2
   end if
   !  0 < (xip- xstart) < 3*pi
   a = -2.0d0/(zeta1 - zeta2)*x1p + 2.0d0/(zeta1 - zeta2)*x2p + dx1 + dx2
   b = 3.0d0*(zeta1+zeta2)/(zeta1-zeta2)*x1p - 3.0d0*(zeta1+zeta2)/(zeta1-zeta2)*x2p + (-zeta1-2.0d0*zeta2)*dx1 + (-2.0d0*zeta1-zeta2)*dx2
   c = -6.0d0*zeta1*zeta2/(zeta1-zeta2)*x1p + 6.0d0*zeta1*zeta2/(zeta1-zeta2)*x2p + (2.0d0*zeta1*zeta2+zeta2**2.0d0)*dx1 + (zeta1**2.0d0+2.0d0*zeta1*zeta2)*dx2
   d = (3.0d0*zeta1*zeta2**2.0d0-zeta2**3.0d0)/(zeta1-zeta2)*x1p + (zeta1**3.0d0-3.0d0*zeta1**2.0d0*zeta2)/(zeta1-zeta2)*x2p - zeta1*zeta2**2.0d0*dx1 - zeta1**2.0d0*zeta2*dx2
   a = a/(zeta1 - zeta2)**2.0d0
   b = b/(zeta1 - zeta2)**2.0d0
   c = c/(zeta1 - zeta2)**2.0d0
   d = d/(zeta1 - zeta2)**2.0d0
end subroutine

subroutine linetocurve(a, b, c, zeta1, zeta2, x1, x2)
   real*8 a, b, c, d, zeta1, zeta2, x1, x2, pi, pi2
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   a = (zeta2 - zeta1)/1.0d10
   b = pi2/(zeta2 - zeta1)
   c = -pi*(3.0d0*zeta1 + zeta2)/2.0d0/(zeta2 - zeta1)
end subroutine

subroutine curve(m, zeta0, c3, npoint, cx, cy, cz, cze, ckx, cky, ckz)
   real*8 zeta, zeta0, pi, pi2, a, b, c, d, flag
   real*8, dimension(3) :: c3
   integer npoint, i, m
   real*8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   zeta = dmod(zeta0, pi2)
   if (zeta < 0.0d0) zeta = zeta + pi2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (cze(npoint) <= zeta) then
      flag = 0.0d0
      call getsplinecoe(a, b, c, d, cze(npoint), pi2, cx(npoint), cx(1), ckx(npoint), ckx(1))
      if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10)) flag = 1.0d0
      c3(1) = a*zeta**3 + b*zeta**2 + c*zeta + d
      call getsplinecoe(a, b, c, d, cze(npoint), pi2, cy(npoint), cy(1), cky(npoint), cky(1))
      if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10) .and. (flag == 1.0d0)) flag = 2.0d0
      c3(2) = a*zeta**3 + b*zeta**2 + c*zeta + d
      call getsplinecoe(a, b, c, d, cze(npoint), pi2, cz(npoint), cz(1), ckz(npoint), ckz(1))
      if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10) .and. (flag == 2.0d0)) flag = 3.0d0
      c3(3) = a*zeta**3 + b*zeta**2 + c*zeta + d
      if (flag == 3.0d0) then
         ! print *,'straight line'
         call linetocurve(a, b, c, cze(npoint), pi2, cx(npoint), cx(1))
         c3(1) = c3(1) + a*dsin(b*zeta + c) + a
         call linetocurve(a, b, c, cze(npoint), pi2, cy(npoint), cy(1))
         c3(2) = c3(2) + a*dsin(b*zeta + c) + a
         call linetocurve(a, b, c, cze(npoint), pi2, cz(npoint), cz(1))
         c3(3) = c3(3) + a*dsin(b*zeta + c) + a
      end if
   else
      do i = 1, npoint - 1
      if ((cze(i) <= zeta) .and. (cze(i + 1) > zeta)) then
         flag = 0.0d0
         call getsplinecoe(a, b, c, d, cze(i), cze(i + 1), cx(i), cx(i + 1), ckx(i), ckx(i + 1))
         if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10)) flag = 1.0d0
         c3(1) = a*zeta**3 + b*zeta**2 + c*zeta + d
         call getsplinecoe(a, b, c, d, cze(i), cze(i + 1), cy(i), cy(i + 1), cky(i), cky(i + 1))
         if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10) .and. (flag == 1.0d0)) flag = 2.0d0
         c3(2) = a*zeta**3 + b*zeta**2 + c*zeta + d
         call getsplinecoe(a, b, c, d, cze(i), cze(i + 1), cz(i), cz(i + 1), ckz(i), ckz(i + 1))
         if ((abs(a) < 1.0d-10) .and. (abs(b) < 1.0d-10) .and. (flag == 2.0d0)) flag = 3.0d0
         c3(3) = a*zeta**3 + b*zeta**2 + c*zeta + d
         if (flag == 3.0d0) then
            ! print *,'straight line'
            call linetocurve(a, b, c, cze(i), cze(i + 1), cx(i), cx(i + 1))
            c3(1) = c3(1) + a*dsin(b*zeta + c) + a
            call linetocurve(a, b, c, cze(i), cze(i + 1), cy(i), cy(i + 1))
            c3(2) = c3(2) + a*dsin(b*zeta + c) + a
            call linetocurve(a, b, c, cze(i), cze(i + 1), cz(i), cz(i + 1))
            c3(3) = c3(3) + a*dsin(b*zeta + c) + a
         end if
      end if
      end do
   end if
   ! if (flag<3.0d0) print *,flag
   if (c3(1) >= pi) c3(1) = c3(1) - pi2
   if (c3(2) >= pi) c3(2) = c3(2) - pi2
   if (c3(3) >= pi) c3(3) = c3(3) - pi2
   if (c3(1) < -pi) c3(1) = c3(1) + pi2
   if (c3(2) < -pi) c3(2) = c3(2) + pi2
   if (c3(3) < -pi) c3(3) = c3(3) + pi2
        !!!!!!!!!!!!!!!!!!!
end subroutine

subroutine curve_spline_pl(a_list, b_list, c_list, d_list, n_sec, cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl)
   real*8 pi, pi2, a, b, c, d
   real*8, dimension(3) :: c3
   integer n_sec, i
   integer flag
   real*8, dimension(n_sec + 2) :: cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl
   real*8, dimension(n_sec + 1, 3) :: a_list, b_list, c_list, d_list
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   do i = 1, n_sec + 1
      flag = 0
      call getsplinecoe(a, b, c, d, cze_pl(i), cze_pl(i + 1), cx_pl(i), cx_pl(i + 1), ckx_pl(i), ckx_pl(i + 1))
      a_list(i, 1) = a
      b_list(i, 1) = b
      c_list(i, 1) = c
      d_list(i, 1) = d
      ! print *, 'i_sec=', i
      ! print *, 'a=', a, 'b=', b
      ! print *, 'c=', c, 'd=', d
      ! print *, 'dcze=', cze_pl(i + 1) - cze_pl(i)
      call getsplinecoe(a, b, c, d, cze_pl(i), cze_pl(i + 1), cy_pl(i), cy_pl(i + 1), cky_pl(i), cky_pl(i + 1))
      a_list(i, 2) = a
      b_list(i, 2) = b
      c_list(i, 2) = c
      d_list(i, 2) = d
      call getsplinecoe(a, b, c, d, cze_pl(i), cze_pl(i + 1), cz_pl(i), cz_pl(i + 1), ckz_pl(i), ckz_pl(i + 1))
      a_list(i, 3) = a
      b_list(i, 3) = b
      c_list(i, 3) = c
      d_list(i, 3) = d
   end do

end subroutine

subroutine curve_spline_pl_v2(a_list, b_list, c_list, d_list, n_sec, cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl)
   real*8 pi, pi2, a, b, c, d
   real*8, dimension(3) :: c3
   integer n_sec, i
   integer flag
   real*8, dimension(n_sec + 2) :: cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl
   real*8, dimension(n_sec + 1, 3) :: a_list, b_list, c_list, d_list
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   do i = 1, n_sec + 1
      flag = 0
      call getsplinecoe(a, b, c, d, 0.0d0, cze_pl(i + 1) - cze_pl(i), cx_pl(i), cx_pl(i + 1), ckx_pl(i), ckx_pl(i + 1))
      a_list(i, 1) = a
      b_list(i, 1) = b
      c_list(i, 1) = c
      d_list(i, 1) = d
      call getsplinecoe(a, b, c, d, 0.0d0, cze_pl(i + 1) - cze_pl(i), cy_pl(i), cy_pl(i + 1), cky_pl(i), cky_pl(i + 1))
      a_list(i, 2) = a
      b_list(i, 2) = b
      c_list(i, 2) = c
      d_list(i, 2) = d
      call getsplinecoe(a, b, c, d, 0.0d0, cze_pl(i + 1) - cze_pl(i), cz_pl(i), cz_pl(i + 1), ckz_pl(i), ckz_pl(i + 1))
      a_list(i, 3) = a
      b_list(i, 3) = b
      c_list(i, 3) = c
      d_list(i, 3) = d
   end do

end subroutine

subroutine coo_curve_ip(zeta0, c3, a, b, c, d)
   real*8 zeta, zeta0, pi, pi2
   real*8, dimension(3) :: c3
   integer ii, i, i_co
   real*8, dimension(3) :: a, b, c, d
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   ! think again
   zeta = zeta0 !bug
   ! zeta = dmod(zeta0, pi2)
   ! if (zeta < 0.0d0) zeta = zeta + pi2

   do i_co = 1, 3
      c3(i_co) = a(i_co)*zeta**3.0d0 + b(i_co)*zeta**2.0d0 + c(i_co)*zeta + d(i_co)
   end do

   !mind
   if (c3(1) >= pi) c3(1) = c3(1) - pi2
   if (c3(2) >= pi) c3(2) = c3(2) - pi2
   if (c3(3) >= pi) c3(3) = c3(3) - pi2
   if (c3(1) < -pi) c3(1) = c3(1) + pi2
   if (c3(2) < -pi) c3(2) = c3(2) + pi2
   if (c3(3) < -pi) c3(3) = c3(3) + pi2

end subroutine

subroutine d1_coo_curve_ip(zeta, dc3, dzeta, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3

   call coo_curve_spline_5th(zeta + dzeta, c3r, coe_3)
   call coo_curve_spline_5th(zeta - dzeta, c3l, coe_3)
   call cminus(c3r, c3l, dc3)
   dc3 = dc3/2.d0/dzeta
end subroutine

subroutine d2_coo_curve_ip(zeta, dc3, dzeta, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3

   call d1_coo_curve_ip(zeta + dzeta, c3r, dzeta, coe_3)
   call d1_coo_curve_ip(zeta - dzeta, c3l, dzeta, coe_3)

   call cminus(c3r, c3l, dc3)
   dc3 = dc3/2.d0/dzeta
end subroutine

subroutine d1_coo_curve_ip_func(zeta, dc3, a, b, c, d)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, i_co
   real*8, dimension(3) :: a, b, c, d
   !think again
   !suppose that zeta +- dzeta belongs to i_p_c, if not, there are minor errors
   !consider the continuity

   do i_co = 1, 3
      dc3(i_co) = 3.0d0*a(i_co)*zeta**2.0d0 + 2.0d0*b(i_co)*zeta + c(i_co)
   end do

end subroutine

subroutine d2_coo_curve_ip_func(zeta, dc3, a, b, c, d)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, i_co
   real*8, dimension(3) :: a, b, c, d
   !think again
   !suppose that zeta +- dzeta belongs to i_p_c, if not, there are minor errors
   !consider the continuity

   do i_co = 1, 3
      dc3(i_co) = 6.0d0*a(i_co)*zeta + 2.0d0*b(i_co)
   end do

end subroutine

subroutine get_spline_5th_coe(coe, zeta_r, x_l, x_r, v_l, v_r, a_l, a_r)
   real*8 pi, pi2
   real*8 zeta_r, x_l, x_r, v_l, v_r, a_l, a_r
   real*8 xp_l, xp_r
   real*8, dimension(0:5) :: coe
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   xp_l = x_l
   xp_r = x_r
   if (x_l - x_r > pi) then
      xp_r = x_r + pi2
   elseif (x_l - x_r < -pi) then
      xp_l = x_l + pi2
   end if
   !  0 < (xp- xstart) < 3*pi
   coe(0) = x_l
   coe(1) = v_l
   coe(2) = a_l/2.0d0
   coe(3) = (20.0d0*(xp_r - xp_l) - (12.0d0*v_l + 8.0d0*v_r)*zeta_r - (3.0d0*a_l - a_r)*zeta_r**2.0d0)/(2.0d0*zeta_r**3.0d0)
   coe(4) = (-30.0d0*(xp_r - xp_l) + (16.0d0*v_l + 14.0d0*v_r)*zeta_r + (3.0d0*a_l - 2.0d0*a_r)*zeta_r**2.0d0)/(2.0d0*zeta_r**4.0d0)
   coe(5) = (12.0d0*(xp_r - xp_l) - 6.0d0*(v_l + v_r)*zeta_r - (a_l - a_r)*zeta_r**2.0d0)/(2.0d0*zeta_r**5.0d0)
end subroutine

subroutine curve_spline_5th(coe_list, n_sec, cze_pl, cx_pl, cy_pl, cz_pl, ckx_pl, cky_pl, ckz_pl, ckkx_pl, ckky_pl, ckkz_pl)
   real*8 pi, pi2
   real*8, dimension(3) :: c3
   integer n_sec, i, j
   real*8, dimension(n_sec + 2) :: cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl, ckkx_pl, ckky_pl, ckkz_pl
   real*8, dimension(n_sec + 1, 0:5, 3) :: coe_list
   real*8, dimension(0:5) :: coe
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   do i = 1, n_sec + 1
      call get_spline_5th_coe(coe, cze_pl(i + 1) - cze_pl(i), cx_pl(i), cx_pl(i + 1), ckx_pl(i), ckx_pl(i + 1), ckkx_pl(i), ckkx_pl(i + 1))
      do j = 0, 5
         coe_list(i, j, 1) = coe(j)
      end do
      call get_spline_5th_coe(coe, cze_pl(i + 1) - cze_pl(i), cy_pl(i), cy_pl(i + 1), cky_pl(i), cky_pl(i + 1), ckky_pl(i), ckky_pl(i + 1))
      do j = 0, 5
         coe_list(i, j, 2) = coe(j)
      end do
      call get_spline_5th_coe(coe, cze_pl(i + 1) - cze_pl(i), cz_pl(i), cz_pl(i + 1), ckz_pl(i), ckz_pl(i + 1), ckkz_pl(i), ckkz_pl(i + 1))
      do j = 0, 5
         coe_list(i, j, 3) = coe(j)
      end do
   end do

end subroutine

subroutine coo_curve_spline_5th(zeta, c3, coe_3)
   real*8 zeta, pi, pi2
   real*8, dimension(3) :: c3
   integer ii, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   ! c3 = 0
   ! do i_co = 1, 3
   !    do j = 0, 5
   !       c3(i_co) = c3(i_co) + coe_3(j, i_co)*zeta**(1.0d0*j)
   !    end do
   ! end do

   do i_co = 1, 3
      c3(i_co) = coe_3(5, i_co)*zeta**5.0d0 + coe_3(4, i_co)*zeta**4.0d0 + &
      &coe_3(3, i_co)*zeta**3.0d0 + coe_3(2, i_co)*zeta**2.0d0 + coe_3(1, i_co)*zeta + coe_3(0, i_co)
   end do

   !mind
   if (c3(1) >= pi) c3(1) = c3(1) - pi2
   if (c3(2) >= pi) c3(2) = c3(2) - pi2
   if (c3(3) >= pi) c3(3) = c3(3) - pi2
   if (c3(1) < -pi) c3(1) = c3(1) + pi2
   if (c3(2) < -pi) c3(2) = c3(2) + pi2
   if (c3(3) < -pi) c3(3) = c3(3) + pi2

end subroutine

subroutine d1_coo_curve_spline_5th(zeta, dc3, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3

   ! dc3 = 0
   ! do i_co = 1, 3
   !    do j = 1, 5
   !       dc3(i_co) = dc3(i_co) + 1.0d0*j*coe_3(j, i_co)*zeta**(1.0d0*j - 1.0d0)
   !    end do
   ! end do

   do i_co = 1, 3
      dc3(i_co) = 5.0d0*coe_3(5, i_co)*zeta**4.0d0 + 4.0d0*coe_3(4, i_co)*zeta**3.0d0 + &
      &3.0d0*coe_3(3, i_co)*zeta**2.0d0 + 2.0d0*coe_3(2, i_co)*zeta + coe_3(1, i_co)
   end do

end subroutine

subroutine d2_coo_curve_spline_5th(zeta, dc3, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3

   ! dc3 = 0
   ! do i_co = 1, 3
   !    do j = 2, 5
   !       dc3(i_co) = dc3(i_co) + 1.0d0*j*(j - 1.0d0)*coe_3(j, i_co)*zeta**(1.0d0*j - 2.0d0)
   !    end do
   ! end do

   do i_co = 1, 3
      dc3(i_co) = 20.0d0*coe_3(5, i_co)*zeta**3.0d0 + 12.0d0*coe_3(4, i_co)*zeta**2.0d0 + &
      &6.0d0*coe_3(3, i_co)*zeta + 2.0d0*coe_3(2, i_co)
   end do

end subroutine

subroutine get_spline_3th_coe(coe, zeta_r, x_l, x_r, v_l, v_r, a_l, a_r)
   real*8 pi, pi2
   real*8 zeta_r, x_l, x_r, v_l, v_r, a_l, a_r
   real*8 xp_l, xp_r
   real*8, dimension(0:5) :: coe
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   xp_l = x_l
   xp_r = x_r
   if (x_l - x_r > pi) then
      xp_r = x_r + pi2
   elseif (x_l - x_r < -pi) then
      xp_l = x_l + pi2
   end if
   !  0 < (xp- xstart) < 3*pi
   coe(0) = x_l
   coe(1) = v_l
   coe(2) = (3.0d0*(xp_r - xp_l) - (2.0d0*v_l + v_r)*zeta_r)/(zeta_r**2.0d0)
   coe(3) = (-2.0d0*(xp_r - xp_l) + (v_l + v_r)*zeta_r)/(zeta_r**3.0d0)
end subroutine

subroutine curve_spline_3th(coe_list, n_sec, cze_pl, cx_pl, cy_pl, cz_pl, ckx_pl, cky_pl, ckz_pl, ckkx_pl, ckky_pl, ckkz_pl)
   real*8 pi, pi2
   real*8, dimension(3) :: c3
   integer n_sec, i, j
   real*8, dimension(n_sec + 2) :: cx_pl, cy_pl, cz_pl, cze_pl, ckx_pl, cky_pl, ckz_pl, ckkx_pl, ckky_pl, ckkz_pl
   real*8, dimension(n_sec + 1, 0:5, 3) :: coe_list
   real*8, dimension(0:5) :: coe
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   do i = 1, n_sec + 1
      call get_spline_3th_coe(coe, cze_pl(i + 1) - cze_pl(i), cx_pl(i), cx_pl(i + 1), ckx_pl(i), ckx_pl(i + 1), ckkx_pl(i), ckkx_pl(i + 1))
      do j = 0, 3
         coe_list(i, j, 1) = coe(j)
      end do
      call get_spline_3th_coe(coe, cze_pl(i + 1) - cze_pl(i), cy_pl(i), cy_pl(i + 1), cky_pl(i), cky_pl(i + 1), ckky_pl(i), ckky_pl(i + 1))
      do j = 0, 3
         coe_list(i, j, 2) = coe(j)
      end do
      call get_spline_3th_coe(coe, cze_pl(i + 1) - cze_pl(i), cz_pl(i), cz_pl(i + 1), ckz_pl(i), ckz_pl(i + 1), ckkz_pl(i), ckkz_pl(i + 1))
      do j = 0, 3
         coe_list(i, j, 3) = coe(j)
      end do
   end do

end subroutine

subroutine coo_curve_spline_3th(zeta, c3, coe_3)
   real*8 zeta, pi, pi2
   real*8, dimension(3) :: c3
   integer ii, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3
   real*8, dimension(0:5) :: zeta_power
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)

   ! c3 = 0
   ! do i_co = 1, 3
   !    do j = 0, 3
   !       c3(i_co) = c3(i_co) + coe_3(j, i_co)*zeta**(1.0d0*j)
   !    end do
   ! end do

   ! do i_co = 1, 3
   !    c3(i_co) = coe_3(3, i_co)*zeta**3.0d0 + coe_3(2, i_co)*zeta**2.0d0 + coe_3(1, i_co)*zeta**1.0d0 + coe_3(0, i_co)
   ! end do

   zeta_power(2) = zeta**2.0d0
   zeta_power(3) = zeta**3.0d0

   do i_co = 1, 3
      c3(i_co) = coe_3(3, i_co)*zeta_power(3) + coe_3(2, i_co)*zeta_power(2) + coe_3(1, i_co)*zeta + coe_3(0, i_co)
   end do

   !mind
   if (c3(1) >= pi) c3(1) = c3(1) - pi2
   if (c3(2) >= pi) c3(2) = c3(2) - pi2
   if (c3(3) >= pi) c3(3) = c3(3) - pi2
   if (c3(1) < -pi) c3(1) = c3(1) + pi2
   if (c3(2) < -pi) c3(2) = c3(2) + pi2
   if (c3(3) < -pi) c3(3) = c3(3) + pi2

end subroutine

subroutine d1_coo_curve_spline_3th(zeta, dc3, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3
   real*8, dimension(0:5) :: zeta_power

   ! dc3 = 0
   ! do i_co = 1, 3
   !    do j = 1, 3
   !       dc3(i_co) = dc3(i_co) + 1.0d0*j*coe_3(j, i_co)*zeta**(1.0d0*j - 1.0d0)
   !    end do
   ! end do

   ! do i_co = 1, 3
   !    dc3(i_co) = 3.0d0*coe_3(3, i_co)*zeta**(2.0d0) + 2.0d0*coe_3(2, i_co)*zeta + coe_3(1, i_co)
   ! end do

   zeta_power(2) = zeta**2.0d0

   do i_co = 1, 3
      dc3(i_co) = 3.0d0*coe_3(3, i_co)*zeta_power(2) + 2.0d0*coe_3(2, i_co)*zeta + coe_3(1, i_co)
   end do

end subroutine

subroutine d2_coo_curve_spline_3th(zeta, dc3, coe_3)
   real*8 zeta, dzeta, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   real*8, dimension(3) :: c3
   integer npoint, i, j, i_co
   real*8, dimension(0:5, 3) :: coe_3

   ! dc3 = 0
   ! do i_co = 1, 3
   !    do j = 2, 3
   !       dc3(i_co) = dc3(i_co) + 1.0d0*j*(j - 1.0d0)*coe_3(j, i_co)*zeta**(1.0d0*j - 2.0d0)
   !    end do
   ! end do

   do i_co = 1, 3
      dc3(i_co) = 6.0d0*coe_3(3, i_co)*zeta + 2.0d0*coe_3(2, i_co)
   end do

end subroutine

!think again
subroutine cminus(c3r, c3l, dc3)
   real*8 cdx, cdy, cdz, pi, pi2
   real*8, dimension(3) :: dc3, c3r, c3l
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   cdx = c3r(1) - c3l(1)
   cdy = c3r(2) - c3l(2)
   cdz = c3r(3) - c3l(3)
   if (cdx > pi) cdx = cdx - pi2
   if (cdy > pi) cdy = cdy - pi2
   if (cdz > pi) cdz = cdz - pi2
   if (cdx < -pi) cdx = cdx + pi2
   if (cdy < -pi) cdy = cdy + pi2
   if (cdz < -pi) cdz = cdz + pi2
   dc3(1) = cdx
   dc3(2) = cdy
   dc3(3) = cdz
end subroutine

subroutine distflag(c3r, c3l, flag)
   real*8 cdx, cdy, cdz, pi, pi2
   real*8, dimension(3) :: c3r, c3l, flag
   pi = 4.0d0*datan(1.0d0)
   pi2 = 8.0d0*datan(1.0d0)
   flag(1) = 0.0d0
   flag(2) = 0.0d0
   flag(3) = 0.0d0
   cdx = c3r(1) - c3l(1)
   cdy = c3r(2) - c3l(2)
   cdz = c3r(3) - c3l(3)
   if (cdx > pi) flag(1) = pi2
   if (cdy > pi) flag(2) = pi2
   if (cdz > pi) flag(3) = pi2
   if (cdx < -pi) flag(1) = -pi2
   if (cdy < -pi) flag(2) = -pi2
   if (cdz < -pi) flag(3) = -pi2
end subroutine

subroutine dsigmads_func_ip_v2(zeta, zeta_sec, Len_Disc, length, dsigmads, sigma_max, sigma_min, dzeta, a, b, c, d)
   real*8 zeta, dsigmads, sigmal, sigmar, dzeta, ds2, length, sigma_max, sigma_min
   real*8 zeta_sec, Len_Disc
   real*8, dimension(3) :: czeta, czetal, czetar, dczetal, dczetar
   real*8, dimension(3) :: dczeta
   real*8, dimension(3) :: a, b, c, d
   !dzeta = 1.0d-4
   call sigma_func(zeta - dzeta, length, sigmal, sigma_max, sigma_min)
   call sigma_func(zeta + dzeta, length, sigmar, sigma_max, sigma_min)
   call d1_coo_curve_ip_func(zeta_sec*Len_Disc, dczeta, a, b, c, d)
   ds2 = dsqrt(dczeta(1)**2.d0 + dczeta(2)**2.d0 + dczeta(3)**2.d0)*2.0d0*dzeta*Len_Disc
   dsigmads = (sigmar - sigmal)/ds2
   !print *, 'dsigmads_func=', dsigmads
end subroutine

subroutine dsigmadzeta(zeta, zeta_sec, Len_Disc, length, dsigmads, sigma_max, sigma_min, dzeta, coe_3)
   real*8 zeta, dsigmads, sigmal, sigmar, dzeta, ds2, length, sigma_max, sigma_min
   real*8 zeta_sec, Len_Disc
   real*8, dimension(3) :: czeta, czetal, czetar, dczetal, dczetar
   real*8, dimension(3) :: dczeta
   real*8, dimension(0:5, 3) :: coe_3
   !dzeta = 1.0d-4
   call sigma_func(zeta - dzeta, length, sigmal, sigma_max, sigma_min)
   call sigma_func(zeta + dzeta, length, sigmar, sigma_max, sigma_min)
   !call d1_coo_curve_spline_5th(zeta_sec*Len_Disc, dczeta, coe_3)
   !ds2 = dsqrt(dczeta(1)**2.d0 + dczeta(2)**2.d0 + dczeta(3)**2.d0)*2.0d0*dzeta*Len_Disc
   !dsigmads = (sigmar - sigmal)/ds2
   dsigmads = (sigmar - sigmal)/(2.0d0*dzeta)
end subroutine

subroutine dsigmads_func_spline_5th(zeta, zeta_sec, Len_Disc, length, dsigmads, sigma_max, sigma_min, dzeta, coe_3)
   real*8 zeta, dsigmads, sigmal, sigmar, dzeta, ds2, length, sigma_max, sigma_min
   real*8 zeta_sec, Len_Disc
   real*8, dimension(3) :: czeta, czetal, czetar, dczetal, dczetar
   real*8, dimension(3) :: dczeta
   real*8, dimension(0:5, 3) :: coe_3
   !dzeta = 1.0d-4
   call sigma_func(zeta - dzeta, length, sigmal, sigma_max, sigma_min)
   call sigma_func(zeta + dzeta, length, sigmar, sigma_max, sigma_min)
   call d1_coo_curve_spline_5th(zeta_sec*Len_Disc, dczeta, coe_3)
   ds2 = dsqrt(dczeta(1)**2.d0 + dczeta(2)**2.d0 + dczeta(3)**2.d0)*2.0d0*dzeta*Len_Disc
   dsigmads = (sigmar - sigmal)/ds2
   !print *, 'dsigmads_func=', dsigmads
end subroutine

subroutine dsigmads_func_spline_3th(zeta, zeta_sec, Len_Disc, length, dsigmads, sigma_max, sigma_min, dzeta, coe_3)
   real*8 zeta, dsigmads, sigmal, sigmar, dzeta, ds2, length, sigma_max, sigma_min
   real*8 zeta_sec, Len_Disc
   real*8, dimension(3) :: czeta, czetal, czetar, dczetal, dczetar
   real*8, dimension(3) :: dczeta
   real*8, dimension(0:5, 3) :: coe_3
   !dzeta = 1.0d-4
   call sigma_func(zeta - dzeta, length, sigmal, sigma_max, sigma_min)
   call sigma_func(zeta + dzeta, length, sigmar, sigma_max, sigma_min)
   call d1_coo_curve_spline_3th(zeta_sec*Len_Disc, dczeta, coe_3)
   ds2 = dsqrt(dczeta(1)**2.d0 + dczeta(2)**2.d0 + dczeta(3)**2.d0)*2.0d0*dzeta*Len_Disc
   dsigmads = (sigmar - sigmal)/ds2
   !print *, 'dsigmads_func=', dsigmads
end subroutine

subroutine norm_zeta(c3, nc3)
   real*8 nc3
   real*8, dimension(3) :: c3
   nc3 = dsqrt(c3(1)**2 + c3(2)**2 + c3(3)**2)
end subroutine

subroutine cross_zeta(vec1, vec2, vec3)
   real*8, dimension(3) :: vec1, vec2, vec3
   vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
   vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
   vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end subroutine

subroutine dot_zeta(vec1, vec2, sc)
   real*8, dimension(3) :: vec1, vec2
   real*8 sc
   sc = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
end subroutine

subroutine norm_t(n, vec, normv)
   integer n
   real*8, dimension(n, 3) :: vec
   real*8, dimension(n) :: normv
   normv = dsqrt(vec(:, 1)*vec(:, 1) + vec(:, 2)*vec(:, 2) + vec(:, 3)*vec(:, 3))
end subroutine

subroutine cross_t(n, vec1, vec2, vec3)
   integer n
   real*8, dimension(n, 3) :: vec1, vec2, vec3
   vec3(:, 1) = vec1(:, 2)*vec2(:, 3) - vec1(:, 3)*vec2(:, 2)
   vec3(:, 2) = vec1(:, 3)*vec2(:, 1) - vec1(:, 1)*vec2(:, 3)
   vec3(:, 3) = vec1(:, 1)*vec2(:, 2) - vec1(:, 2)*vec2(:, 1)
end subroutine

subroutine dot_t(n, vec1, vec2, sc)
   integer n
   real*8, dimension(n, 3) :: vec1, vec2
   real*8, dimension(n) :: sc
   sc = vec1(:, 1)*vec2(:, 1) + vec1(:, 2)*vec2(:, 2) + vec1(:, 3)*vec2(:, 3)
end subroutine

subroutine pt(n, dt, vec, dvec)
   integer n
   real*8 dt
   real*8, dimension(n, 3) :: vec, dvec, vecr, vecl
   vecr(1:n - 1, :) = vec(2:n, :)
   vecr(n, :) = vec(1, :)
   vecl(2:n, :) = vec(1:n - 1, :)
   vecl(1, :) = vec(n, :)
   dvec = (vecr - vecl)/2.d0/dt
end subroutine

subroutine ptboundary(n, dt, vec, dvec)
   integer n, i
   real*8 dt
   real*8, dimension(n, 3) :: vec, dvec, vecr, vecl
   vecr(1:n - 1, :) = vec(2:n, :)
   vecr(n, :) = vec(1, :)
   vecl(2:n, :) = vec(1:n - 1, :)
   vecl(1, :) = vec(n, :)
   !dvec = (vecr - vecl) / 2.d0 / dt
   do i = 1, n
      call cminus(vecr(i, :), vecl(i, :), dvec(i, :))
   end do
   dvec = dvec/2.d0/dt
end subroutine

subroutine set_fft_plans(nx, ny, nz, planxf, planxb, planyf, planyb, planzf, planzb)
   implicit none
   include "fftw3.f"
   integer*8 :: planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz
   double complex, dimension(nx) :: tempx
   double complex, dimension(ny) :: tempy
   double complex, dimension(nz) :: tempz
   call dfftw_plan_dft_1d(planxf, nx, tempx, tempx, FFTW_FORWARD, FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(planxb, nx, tempx, tempx, FFTW_BACKWARD, FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(planyf, ny, tempy, tempy, FFTW_FORWARD, FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(planyb, ny, tempy, tempy, FFTW_BACKWARD, FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(planzf, nz, tempz, tempz, FFTW_FORWARD, FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(planzb, nz, tempz, tempz, FFTW_BACKWARD, FFTW_ESTIMATE)
end subroutine set_fft_plans

subroutine destroy_fft_plans(planxf, planxb, planyf, planyb, planzf, planzb)
   implicit none
   include "fftw3.f"
   integer*8 :: planxf, planyf, planzf, planxb, planyb, planzb
   call dfftw_destroy_plan(planxf)
   call dfftw_destroy_plan(planxb)
   call dfftw_destroy_plan(planyf)
   call dfftw_destroy_plan(planyb)
   call dfftw_destroy_plan(planzf)
   call dfftw_destroy_plan(planzb)
end subroutine destroy_fft_plans

SUBROUTINE wavenumber(nx, ny, nz, nzp, lx, ly, lz, kx, ky, kz, k2, id)
   implicit none
   integer nx, ny, nz, nzp, i, j, k, id
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   real*8, dimension(nx, ny, nzp) :: k2
   real*8 lx, ly, lz, pi2
   real*8 dkx, dky, dkz

   pi2 = 8.0d0*datan(1.0d0)
   dkx = pi2/lx
   dky = pi2/ly
   dkz = pi2/lz

   do i = 1, nx
      kx(i) = (mod(i - 1 + nx/2, nx) - nx/2)*dkx
   end do
   do j = 1, ny
      ky(j) = (mod(j - 1 + ny/2, ny) - ny/2)*dky
   end do
   do k = 1, nzp
      kz(k) = (mod(k - 1 + nz/2 + id*nzp, nz) - nz/2)*dkz
   end do
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            k2(i, j, k) = kx(i)**2 + ky(j)**2 + kz(k)**2
         end do
      end do
   end do
end subroutine wavenumber

subroutine initialize_mesh(nx, ny, nzp, xstart, ystart, zstart, dx, dy, dz, meshx, meshy, meshz, id)
   implicit none
   integer nx, ny, nzp, i, j, k, id
   real*8 xstart, ystart, zstart, dx, dy, dz
   real*8, dimension(nx, ny, nzp) :: meshx, meshy, meshz
   do i = 1, nx
      meshx(i, :, :) = (i - 1.d0)*dx + xstart
   end do
   do j = 1, ny
      meshy(:, j, :) = (j - 1.d0)*dy + ystart
   end do
   do k = 1, nzp
      meshz(:, :, k) = (k - 1.d0 + nzp*id)*dz + zstart
   end do
end subroutine initialize_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_yz(mx, my, mz, nproc, spec, spectemp, id)
   implicit none
   include 'mpif.h'
   real ierr
   integer mx, my, mz, nproc, i, j, k, imp, mzp, myp, id
   integer, dimension(nproc) :: counts, displs
   double complex, dimension(mx, my, mz/nproc) :: spec
   double complex, dimension(mx, mz, my/nproc) :: spectemp
   double complex, dimension(mx, my/nproc, mz) :: spectemp1
   double complex, dimension(mx, my/nproc, mz/nproc) :: specb
   mzp = mz/nproc
   myp = my/nproc
   do i = 1, nproc
      counts(i) = mx*myp*mzp
      displs(i) = (i - 1)*counts(i)
   end do
   do imp = 1, nproc
      do k = 1, mzp
         do j = 1, myp
            specb(:, j, k) = spec(:, j + (imp - 1)*myp, k)
         end do
      end do
      call mpi_gatherv(specb, counts(id + 1), MPI_DOUBLE_COMPLEX, spectemp1,&
          &counts, displs, MPI_DOUBLE_COMPLEX, imp - 1, mpi_comm_world, ierr)
   end do
   do j = 1, myp
      do k = 1, mz
         spectemp(:, k, j) = spectemp1(:, j, k)
      end do
   end do
end subroutine transpose_yz

!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D forward fast fourier transformation
!! mpi for z direction,
!!    number of blocks: nproc=nz/nzp
!!    id = 0, 1, ... , nproc-1
!! input:
!!    phy: real, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!          fft plans: planxf, planyf, planzf
!! output:
!!          spec: complex,dimension (nx,ny,nzp)
subroutine fourier_forward(phy, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   implicit none
   include 'mpif.h'
   include "fftw3.f"
   integer*8 planxf, planyf, planzf
    !!!MPI
   integer :: id, nproc, ierr
   integer nx, ny, nz, nzp, nyp, i, j, k
   real*8, dimension(nx, ny, nz/nproc) :: phy
   double complex, dimension(nx, ny, nz/nproc) :: spec
   double complex, dimension(nx) :: tempx
   double complex, dimension(ny) :: tempy
   double complex, dimension(nz) :: tempz
   double complex, dimension(nx, nz, ny/nproc) :: spectemp
   nyp = ny/nproc
   nzp = nz/nproc
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            tempx(i) = dcmplx(phy(i, j, k)/nx + 0.0d0, 0.0d0)
         end do
         call dfftw_execute_dft(planxf, tempx, tempx)
         spec(:, j, k) = tempx
      end do
   end do
   do k = 1, nzp
      do i = 1, nx
         tempy = spec(i, :, k)/ny
         call dfftw_execute_dft(planyf, tempy, tempy)
         spec(i, :, k) = tempy
      end do
   end do
   call transpose_yz(nx, ny, nz, nproc, spec, spectemp, id)
   do j = 1, nyp
      do i = 1, nx
         tempz = spectemp(i, :, j)/nz
         call dfftw_execute_dft(planzf, tempz, tempz)
         spectemp(i, :, j) = tempz
      end do
   end do
   call transpose_yz(nx, nz, ny, nproc, spectemp, spec, id)
end subroutine fourier_forward

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D backward fast fourier transformation
!! mpi for z direction,
!!    number of blocks: nproc=nz/nzp
!!    id = 0, 1, ... , nproc-1
!! input:
!!    phy: real, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!          fft plans: planxf, planyf, planzf
!! output:
!!          spec: complex,dimension (nx,ny,nzp)
subroutine fourier_backward(phy, spec, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   implicit none
   include 'mpif.h'
   include "fftw3.f"
   integer*8 planxb, planyb, planzb
   integer :: id, nproc, ierr
   integer nx, ny, nz, nzp, nyp, i, j, k
   real*8, dimension(nx, ny, nz/nproc) :: phy
   double complex, dimension(nx, ny, nz/nproc) :: spec
   double complex, dimension(nx) :: tempx
   double complex, dimension(ny) :: tempy
   double complex, dimension(nz) :: tempz
   double complex, dimension(nx, nz, ny/nproc) :: spectemp
   double complex, allocatable :: specall(:, :, :)
   double complex, allocatable :: specalltemp(:, :, :)
   ! ?????
   ! if (id == 0) then
   !    allocate (specall(nx, ny, nz))
   !    allocate (specalltemp(nx, nz, ny))
   ! end if
   ! ?????
   nyp = ny/nproc
   nzp = nz/nproc
   do k = 1, nzp
      do j = 1, ny
         tempx = spec(:, j, k)
         call dfftw_execute_dft(planxb, tempx, tempx)
         spec(:, j, k) = tempx
      end do
   end do
   do k = 1, nzp
      do i = 1, nx
         tempy = spec(i, :, k)
         call dfftw_execute_dft(planyb, tempy, tempy)
         spec(i, :, k) = tempy
      end do
   end do
   call transpose_yz(nx, ny, nz, nproc, spec, spectemp, id)
   do j = 1, nyp
      do i = 1, nx
         tempz = spectemp(i, :, j)
         call dfftw_execute_dft(planzb, tempz, tempz)
         spectemp(i, :, j) = tempz
      end do
   end do
   call transpose_yz(nx, nz, ny, nproc, spectemp, spec, id)
   phy = dreal(spec)
end subroutine fourier_backward

subroutine dissipation_field(velx, vely, velz, disp_field, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz, nzp, i, j, k, ii, jj, kk
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   real*8, dimension(nx, ny, nzp) :: k2
   integer :: id, nproc, switch_d
   real*8, dimension(nx, ny, nzp) :: velx, vely, velz, disp_field
   real*8, dimension(3, 3, nx, ny, nzp) ::dudx, sij
   double complex, dimension(nx, ny, nzp) :: spec
   nz = nzp*nproc

   do j = 1, 3
call dx_dy_dz_dp_dm(velx, dudx(1,j,:,:,:), j, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,planxb, planyb, planzb, id, nproc)
call dx_dy_dz_dp_dm(vely, dudx(2,j,:,:,:), j, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,planxb, planyb, planzb, id, nproc)
call dx_dy_dz_dp_dm(velz, dudx(3,j,:,:,:), j, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,planxb, planyb, planzb, id, nproc)
   end do

   do i = 1, 3
      do j = 1, 3
         sij(i, j, :, :, :) = (dudx(i, j, :, :, :) + dudx(j, i, :, :, :))/2.0d0
      end do
   end do

   disp_field = 0
   do kk = 1, nzp
      do jj = 1, ny
         do ii = 1, nx
            do i = 1, 3
               do j = 1, 3
                  disp_field(ii, jj, kk) = disp_field(ii, jj, kk) + sij(i, j, ii, jj, kk)*sij(i, j, ii, jj, kk)
               end do
            end do
         end do
      end do
   end do

end subroutine dissipation_field

subroutine dx_dy_dz_dp_dm(phy, dphy, switch_d, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz, nzp, i, j, k
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   real*8, dimension(nx, ny, nzp) :: k2
   integer :: id, nproc, switch_d
   real*8, dimension(nx, ny, nzp) :: phy, dphy
   double complex, dimension(nx, ny, nzp) :: spec
   nz = nzp*nproc
   call fourier_forward(phy, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   if (switch_d == 1) then
      do i = 1, nx
         spec(i, :, :) = spec(i, :, :)*dcmplx(0.0d0, kx(i) + 0.0d0)
      end do
   elseif (switch_d == 2) then
      do j = 1, ny
         spec(:, j, :) = spec(:, j, :)*dcmplx(0.0d0, ky(j) + 0.0d0)
      end do
   elseif (switch_d == 3) then
      do k = 1, nzp
         spec(:, :, k) = spec(:, :, k)*dcmplx(0.0d0, kz(k) + 0.0d0)
      end do
   elseif (switch_d == 6) then
      spec = -k2*spec
   elseif (switch_d == -6) then
      spec = -spec/k2
      if (id == 0) then
         spec(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
      end if
   end if
   call fourier_backward(dphy, spec, nx, ny, nz, planxb, planyb, planzb, id, nproc)
end subroutine dx_dy_dz_dp_dm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!if (switch_d==-1) (dphy1,dphy2,dphy3) is divergence free
subroutine cross_vector(phy1, phy2, phy3, dphy1, dphy2, dphy3, switch_d, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz, nzp, i, j, k
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   real*8, dimension(nx, ny, nzp) :: k2
   integer :: id, nproc, switch_d
   real*8, dimension(nx, ny, nzp) :: phy1, phy2, phy3, dphy1, dphy2, dphy3
   double complex, dimension(nx, ny, nzp) :: spec1, spec2, spec3, spec4, spec5, spec6
   nz = nzp*nproc
   call fourier_forward(phy1, spec1, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(phy2, spec2, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(phy3, spec3, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   if (switch_d == 1) then
      do k = 1, nzp
         do j = 1, ny
            do i = 1, nx
               spec4(i, j, k) = dcmplx(0.0d0, 1.0d0)*(ky(j)*spec3(i, j, k) - kz(k)*spec2(i, j, k))
               spec5(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kz(k)*spec1(i, j, k) - kx(i)*spec3(i, j, k))
               spec6(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kx(i)*spec2(i, j, k) - ky(j)*spec1(i, j, k))
            end do
         end do
      end do
   elseif (switch_d == -1) then
      do k = 1, nzp
         do j = 1, ny
            do i = 1, nx
               spec4(i, j, k) = dcmplx(0.0d0, 1.0d0)*(ky(j)*spec3(i, j, k) - kz(k)*spec2(i, j, k))/k2(i, j, k)
               spec5(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kz(k)*spec1(i, j, k) - kx(i)*spec3(i, j, k))/k2(i, j, k)
               spec6(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kx(i)*spec2(i, j, k) - ky(j)*spec1(i, j, k))/k2(i, j, k)
            end do
         end do
      end do
      if (id == 0) then
         spec4(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
         spec5(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
         spec6(1, 1, 1) = dcmplx(0.0d0, 0.0d0)
      end if
   end if
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            if (k2(i, j, k) > (nx**2 + ny**2 + nz**2)/12) then !kmax = N/2
               spec4(i, j, k) = dcmplx(0.0d0, 0.0d0)
               spec5(i, j, k) = dcmplx(0.0d0, 0.0d0)
               spec6(i, j, k) = dcmplx(0.0d0, 0.0d0)
            end if
         end do
      end do
   end do
   call fourier_backward(dphy1, spec4, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(dphy2, spec5, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(dphy3, spec6, nx, ny, nz, planxb, planyb, planzb, id, nproc)
end subroutine cross_vector

subroutine filter_spec(k_c, phy1, phy2, phy3, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz, nzp, i, j, k, k_c
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   real*8, dimension(nx, ny, nzp) :: k2
   integer :: id, nproc, switch_d
   real*8, dimension(nx, ny, nzp) :: phy1, phy2, phy3
   double complex, dimension(nx, ny, nzp) :: spec
   nz = nzp*nproc

   call fourier_forward(phy1, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)

   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            if (k2(i, j, k) > k_c**2) then
               spec(i, j, k) = dcmplx(0.0d0, 0.0d0)
            end if
         end do
      end do
   end do

   call fourier_backward(phy1, spec, nx, ny, nz, planxb, planyb, planzb, id, nproc)

   call fourier_forward(phy2, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            if (k2(i, j, k) > k_c**2) then
               spec(i, j, k) = dcmplx(0.0d0, 0.0d0)
            end if
         end do
      end do
   end do
   call fourier_backward(phy2, spec, nx, ny, nz, planxb, planyb, planzb, id, nproc)

   call fourier_forward(phy3, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            if (k2(i, j, k) > k_c**2) then
               spec(i, j, k) = dcmplx(0.0d0, 0.0d0)
            end if
         end do
      end do
   end do
   call fourier_backward(phy3, spec, nx, ny, nz, planxb, planyb, planzb, id, nproc)

end subroutine filter_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine divergence(phy1, phy2, phy3, dphy, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   integer nx, ny, nz, nzp, i, j, k
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   integer :: id, nproc
   real*8, dimension(nx, ny, nzp) :: phy1, phy2, phy3, dphy
   double complex, dimension(nx, ny, nzp) :: spec1, spec2, spec3
   nz = nzp*nproc
   call fourier_forward(phy1, spec1, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(phy2, spec2, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   call fourier_forward(phy3, spec3, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            spec1(i, j, k) = dcmplx(0.0d0, 1.0d0)*(kx(i)*spec1(i, j, k) + ky(j)*spec2(i, j, k) + kz(k)*spec3(i, j, k))
         end do
      end do
   end do
   call fourier_backward(dphy, spec1, nx, ny, nz, planxb, planyb, planzb, id, nproc)
end subroutine divergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gradient(phy, dphy1, dphy2, dphy3, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf,&
    &planxb, planyb, planzb, id, nproc)
   implicit none
   integer nx, ny, nz, nzp, i, j, k
   integer*8 planxf, planyf, planzf, planxb, planyb, planzb
   real*8, DIMENSION(nx) :: kx
   real*8, DIMENSION(ny) :: ky
   real*8, DIMENSION(nzp) :: kz
   integer :: id, nproc
   real*8, dimension(nx, ny, nzp) :: dphy1, dphy2, dphy3, phy
   double complex, dimension(nx, ny, nzp) :: spec, spec1, spec2, spec3
   nz = nzp*nproc
   call fourier_forward(phy, spec, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            spec1(i, j, k) = dcmplx(0.0d0, 1.0d0)*kx(i)*spec(i, j, k)
            spec2(i, j, k) = dcmplx(0.0d0, 1.0d0)*ky(j)*spec(i, j, k)
            spec3(i, j, k) = dcmplx(0.0d0, 1.0d0)*kz(k)*spec(i, j, k)
         end do
      end do
   end do
   call fourier_backward(dphy1, spec1, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(dphy2, spec2, nx, ny, nz, planxb, planyb, planzb, id, nproc)
   call fourier_backward(dphy3, spec3, nx, ny, nz, planxb, planyb, planzb, id, nproc)
end subroutine gradient

subroutine cross_product(nx, ny, nzp, phy1, phy2, phy3, phy11, phy12, phy13, phy21, phy22, phy23)
   implicit none
   integer nx, ny, nzp
   real*8, dimension(nx, ny, nzp) :: phy1, phy2, phy3, phy11, phy12, phy13, phy21, phy22, phy23
   phy21 = phy2*phy13 - phy3*phy12
   phy22 = phy3*phy11 - phy1*phy13
   phy23 = phy1*phy12 - phy2*phy11
end subroutine cross_product

subroutine output_v(name, varname, nx, ny, nzp, data_box, nbox, id, nproc)
   implicit none
   include 'mpif.h'
   character*200 name
   integer nx, ny, nzp, nz, id, nproc, ierr, i, nbox
   integer, dimension(nproc) :: counts, displs
   real, dimension(nx, ny, nzp, nbox) :: data_box
   real, dimension(nx, ny, nzp) :: temp
   real, allocatable :: v_box(:, :, :, :)
   character*40, dimension(nbox) :: varname
   if (id == 0) then
      nz = nzp*nproc
      allocate (v_box(nx, ny, nz, nbox))
   end if
   do i = 1, nproc
      counts(i) = nx*ny*nzp
      displs(i) = (i - 1)*counts(i)
   end do
   do i = 1, nbox
      temp = data_box(:, :, :, i)
      call mpi_gatherv(temp, counts(id + 1), MPI_real, v_box(:, :, :, i), counts, displs,&
          &MPI_real, 0, mpi_comm_world, ierr)
   end do
   if (id == 0) call output_puredata(v_box, nx, ny, nz, nbox, name, varname)
   if (id == 0) deallocate (v_box)
end subroutine output_v

subroutine output_v_multi(name, varname, nx, ny, nzp, data_box, nbox, id, nproc)
   implicit none
   include 'mpif.h'
   character*200 name
   integer nx, ny, nzp, nz, id, nproc, ierr, i, nbox
   integer, dimension(nproc) :: counts, displs
   real, dimension(nx, ny, nzp, nbox) :: data_box
   real, dimension(nx, ny, nzp) :: temp
   real, allocatable :: v_box(:, :, :, :)
   character*40, dimension(nbox) :: varname
   INTEGER    ::  i1dd, i2dd, i3dd, i4dd, i5dd, i6dd

   i3dd = mod(id + 1, 10000)/1000
   i4dd = mod(id + 1, 1000)/100
   i5dd = mod(id + 1, 100)/10
   i6dd = mod(id + 1, 10)/1
   name = './output/field_multi'//char(i3dd + 48)//char(i4dd + 48)//char(i5dd + 48)//char(i6dd + 48)//'.dat'
   call output_puredata(data_box, nx, ny, nzp, nbox, name, varname)
end subroutine output_v_multi

subroutine output_puredata(v_box, nx, ny, nz, nbox, name, varname)
   implicit none
   integer nx, ny, nz, nbox
   real, dimension(nx, ny, nz, nbox) :: v_box
   character*200 name
   character*40, dimension(nbox) :: varname

   open (unit=99, file=name, form="BINARY")
   write (99) v_box
   close (99)
end subroutine output_puredata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output binary data for tecplot
!! input:
!!    v_box: real, dimension (nx,ny,nz,nbox)
!!    nbox: number of output varible
!!    mesh size: nx, ny, nz
!!          name: character*200 name of the output data
!!    varname: character*40, dimension (nbox) names of varible
subroutine output_3d_tecplot_bin(v_box, nx, ny, nz, nbox, name, varname)
   implicit none
   integer nx, ny, nz, nbox
   real, dimension(nx, ny, nz, nbox) :: v_box
   real*8, dimension(nbox) :: min_value, max_value
   real*4 ZONEMARKER, EOHMARKER
   integer len, i
   character*40 Title, var
   character*200 name
   character*40, dimension(nbox) :: varname
   character*40 Zonename
   character(40) instring
   ZONEMARKER = 299.0
   EOHMARKER = 357.0
   do i = 1, nbox
      min_value(i) = minval(v_box(:, :, :, i))
      max_value(i) = maxval(v_box(:, :, :, i))
   end do
   open (unit=99, file=name, form="BINARY")
   !I. The header section.
   !1.1 Magic number, Version number
   write (99) "#!TDV112"
   !1.2. Integer value of 1.
   write (99) 1
   !1.3. Title and variable names.
   !Filetype
   write (99) 0
   !1.3.1. The TITLE.
   Title = ""
   call dumpstring(Title)
   !1.3.2 Number of variables (NumVar) in the datafile.
   write (99) nbox
   !1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]
   do i = 1, nbox
      call dumpstring(varname(i))
   end do
   !1.4. Zones
   !Zone marker. Value = 299.0
   write (99) ZONEMARKER
   !Zone name.
   Zonename = "ZONE 001"
   call dumpstring(Zonename)
   !ParentZone
   write (99) - 1
   !StrandID
   write (99) - 1
   !solution time
   write (99) 0
   write (99) 0
   !not used
   write (99) - 1
   !ZoneType
   write (99) 0
   !DataPacking 0=Block, 1=Point
   write (99) 0
   !Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
   write (99) 0
   !Number of user defined face neighbor connections (value >= 0)
   write (99) 0
   !IMax,JMax,KMax
   write (99) nx
   write (99) ny
   write (99) nz
   !1=Auxiliary name/value pair to follow   0=No more Auxiliar name/value pairs.
   write (99) 0
   !I HEADER OVER
   !EOHMARKER, value=357.0
   write (99) EOHMARKER
   !II. Data section
   !2.1 zone
   write (99) Zonemarker
   !variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
   do i = 1, nbox
      write (99) 1
   end do
   !Has variable sharing 0 = no, 1 = yes.
   write (99) 0
   !Has passive variables 0 = no, 1 = yes.
   write (99) 0
   !Zone number to share connectivity list with (-1 = no sharing).
   write (99) - 1
   !min value
   !max value
   do i = 1, nbox
      write (99) min_value(i)
      write (99) max_value(i)
   end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Zone Data. Each variable is in data format asspecified above.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write (99) v_box
   close (99)
end subroutine output_3d_tecplot_bin

!

!
subroutine dumpstring(instring)
    !!!for binary output
   character(40) instring
   integer len
   len = LEN_TRIM(instring)
   do i = 1, len
      ii = ICHAR(instring(i:i))
      write (99) ii
   end do
   write (99) 0
end

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cal_deviation(vorx, vory, vorz, phiv, nx, ny, nz, nzp, deviation, dx, dy, dz, id, igst)
   implicit none
   include 'mpif.h'
   integer i, j, k
   real*8 dx, dy, dz
   integer nx, ny, nzp, nz, igst, ierr, id
   real*8 deviation, deviation1
   real*8 eps
   real*8 dxphi, dyphi, dzphi, ome_dot_graph, graph, vorti
   real*8, dimension(nx, ny, nzp) :: vorx, vory, vorz
   real*8 phiv(1 - igst:nx + igst, 1 - igst:ny + igst, 1 - igst:nzp + igst)
   eps = 0.1d-10
   deviation1 = 0.d0
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
     dxphi=(9.0d0*phiv(i-2,j,k)-phiv(i-3,j,k)-45.0d0*phiv(i-1,j,k)+45.0d0*phiv(i+1,j,k)-9.0d0*phiv(i+2,j,k)+phiv(i+3,j,k))/60.0d0/dx
     dyphi=(9.0d0*phiv(i,j-2,k)-phiv(i,j-3,k)-45.0d0*phiv(i,j-1,k)+45.0d0*phiv(i,j+1,k)-9.0d0*phiv(i,j+2,k)+phiv(i,j+3,k))/60.0d0/dy
     dzphi=(9.0d0*phiv(i,j,k-2)-phiv(i,j,k-3)-45.0d0*phiv(i,j,k-1)+45.0d0*phiv(i,j,k+1)-9.0d0*phiv(i,j,k+2)+phiv(i,j,k+3))/60.0d0/dz
            graph = dsqrt(dxphi**2.d0 + dyphi**2.d0 + dzphi**2.d0)
            vorti = dsqrt((vorx(i, j, k))**2.d0 + (vory(i, j, k))**2.d0 + (vorz(i, j, k))**2.d0)
            ome_dot_graph = vorx(i, j, k)*dxphi + vory(i, j, k)*dyphi + vorz(i, j, k)*dzphi
            deviation1 = deviation1 + dabs(ome_dot_graph)/(vorti*graph + eps)/nx/ny/nz
         end do
      end do
   end do
   call MPI_REDUCE(deviation1, deviation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
end subroutine cal_deviation
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cal_helicity(vorx, vory, vorz, velx, vely, velz, nx, ny, nz, nzp, helicity, dv, id)
   implicit none
   include 'mpif.h'
   integer i, j, k
   real*8 dv
   integer nx, ny, nzp, nz, ierr, id
   real*8 helicity, helicity1
   real*8, dimension(nx, ny, nzp) :: vorx, vory, vorz, velx, vely, velz
   helicity1 = 0.d0
   do k = 1, nzp
      do j = 1, ny
         do i = 1, nx
            helicity1 = helicity1 + vorx(i, j, k)*velx(i, j, k) + &
            &vory(i, j, k)*vely(i, j, k) + vorz(i, j, k)*velz(i, j, k)
         end do
      end do
   end do
   call MPI_REDUCE(helicity1, helicity, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   helicity = helicity*dv
end subroutine cal_helicity

!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!
!Set the value on the ghost points for phiv
subroutine wall_phiv(nx, ny, nzp, phiv, igst, id_l, id_r)
   implicit none
   include 'mpif.h'
   integer ierr, i, j
   integer nx, ny, nzp, igst, nxy, id_l, id_r
   integer status(MPI_STATUS_SIZE)
   real*8 phiv(1 - igst:nx + igst, 1 - igst:ny + igst, 1 - igst:nzp + igst)
   if (igst /= 0) then
      nxy = (nx + 2*igst)*(ny + 2*igst)*igst
      do i = 1, igst
         phiv(nx + i, :, :) = phiv(i, :, :)
         phiv(-igst + i, :, :) = phiv(nx - igst + i, :, :)
      end do
      do j = 1, igst
         phiv(:, ny + j, :) = phiv(:, j, :)
         phiv(:, -igst + j, :) = phiv(:, ny - igst + j, :)
      end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!z-direction, need communication between nodes
        !!!Communication with each other as on a ring.
      call mpi_sendrecv(phiv(:, :, 1:igst), nxy, MPI_DOUBLE_PRECISION, id_l, 99, &
                        phiv(:, :, nzp + 1:nzp + igst), nxy, MPI_DOUBLE_PRECISION, id_r, 99, mpi_comm_world, status, ierr)
      call mpi_sendrecv(phiv(:, :, nzp - igst + 1:nzp), nxy, MPI_DOUBLE_PRECISION, id_r, 99, &
                        phiv(:, :, 1 - igst:0), nxy, MPI_DOUBLE_PRECISION, id_l, 99, mpi_comm_world, status, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
   end if
end subroutine
!!!!!!!!!

subroutine spectrum(nx, ny, nzp, velx, vely, velz, num_data, k2, id, nproc, planxf, planyf, planzf)
   implicit none
   integer*8 planxf, planyf, planzf
   integer nx, ny, nzp, nz, num_data
   real*8, dimension(nx, ny, nzp) :: velx, vely, velz
   integer id, nproc
   real*8, dimension(nx, ny, nzp) :: k2
   double complex, dimension(nx, ny, nzp) :: spec_v
   real, dimension(nx, ny, nzp) :: spec
   nz = nzp*nproc
   call fourier_forward(velx, spec_v, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   spec = real(spec_v*conjg(spec_v))/2.0
   call fourier_forward(vely, spec_v, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   spec = spec + real(spec_v*conjg(spec_v))/2.0
   call fourier_forward(velz, spec_v, nx, ny, nz, planxf, planyf, planzf, id, nproc)
   spec = spec + real(spec_v*conjg(spec_v))/2.0
   call get_spectrum(nx, ny, nzp, nz, spec, num_data, k2, id)
end subroutine spectrum

subroutine get_spectrum(nx, ny, nzp, nz, spec, num_data, k2, id)
   implicit none
   include 'mpif.h'
   real ierr
   integer nx, ny, nz, nzp, nproc, num_data, nek, i, id
   real, dimension(nx, ny, nzp) :: spec
   real*8, dimension(nx, ny, nzp) :: k2, ik2
   real spectrum_k1, spectrum_k
   nek = int(sqrt((nx**2 + ny**2 + nz**2)/12.0)) !kmax = N/2
   ik2 = int(sqrt(k2 + 0.00000001) + 0.5)
   do i = 1, nek
      spectrum_k1 = sum(spec, mask=(ik2 .eq. i))
      call MPI_REDUCE(spectrum_k1, spectrum_k, 1, MPI_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (id == 0) write (num_data, *) i, spectrum_k
   end do
end subroutine get_spectrum
