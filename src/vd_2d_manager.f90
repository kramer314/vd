module vd_2d_manager

  use numerics, only: numerics_linspace, numerics_trapz
  use precision, only: ip, fp

  use vd_2d, only: vd_2d_obj
  use vd, only: vd_get_indices

  implicit none

  private

  public vd_2d_manager_obj
  type vd_2d_manager_obj
     real(fp), pointer :: vd_p_arr(:,:)
     real(fp) :: j_eps
     logical :: semi_classical
     logical :: vd_disjoint
     integer(ip) :: vd_np_stdev
     integer(ip) :: nxl, nxr, nyl, nyr
     real(fp) :: dx, dy
     integer(ip) :: nx, ny
     integer(ip) :: xl_min, xl_max, xr_min, xr_max
     integer(ip) :: yl_min, yl_max, yr_min, yr_max

     integer(ip) :: npx, npy
     real(fp) :: px_min, px_max, py_min, py_max
     real(fp) :: dpx, dpy

     real(fp), pointer :: px_range(:), py_range(:)

     real(fp) :: dt

     real(fp) :: hbar, m
     real(fp) :: net_flux
     type(vd_2d_obj), allocatable :: vdl_arr(:,:), vdr_arr(:,:), &
          vdt_arr(:,:), vdb_arr(:,:)

   contains
     procedure :: init => vd_2d_manager_init
     procedure :: cleanup => vd_2d_manager_cleanup
     procedure :: finalize => vd_2d_manager_finalize
     procedure :: update => vd_2d_manager_update
  end type vd_2d_manager_obj

contains

  subroutine vd_2d_manager_init(this, nx, ny, nxl_ext, nxr_ext, nxl_vd, &
       nxr_vd, nyl_ext, nyr_ext, nyl_vd, nyr_vd, dx, dy, npx, npy, px_min, &
       px_max, py_min, py_max, dt, j_eps, sc, vd_disjoint, hbar, m,&
       vd_np_stdev)
    class(vd_2d_manager_obj), intent(inout) :: this
    integer(ip), intent(in) :: nx, ny
    integer(ip), intent(in) :: nxl_ext, nxr_ext, nxl_vd, nxr_vd
    integer(ip), intent(in) :: nyl_ext, nyr_ext, nyl_vd, nyr_vd
    real(fp), intent(in) :: dx, dy
    integer(ip), intent(in) :: npx, npy
    real(fp), intent(in) :: px_min, px_max, py_min, py_max
    real(fp), intent(in) :: dt
    real(fp), intent(in) :: j_eps
    logical :: sc
    logical :: vd_disjoint
    real(fp) :: hbar, m
    integer(ip), optional :: vd_np_stdev

    integer(ip) :: i_x, i_y

    this%nx = nx
    this%ny = ny

    this%nxl = nxl_vd
    this%nxr = nxr_vd
    this%nyl = nyl_vd
    this%nyr = nyr_vd

    this%dx = dx
    this%dy = dy

    this%npx = npx
    this%npy = npy
    this%px_min = px_min
    this%px_max = px_max
    this%py_min = py_min
    this%py_max = py_max

    this%dt = dt

    this%j_eps = j_eps
    this%semi_classical = sc
    this%vd_disjoint = vd_disjoint

    if (present(vd_np_stdev)) then
       this%vd_np_stdev = vd_np_stdev
    end if

    this%hbar = hbar
    this%m = m

    this%net_flux = 0.0_fp

    call vd_get_indices(this%nx, nxl_ext, nxr_ext, this%nxl, this%nxr, &
         this%xl_min, this%xl_max, this%xr_min, this%xr_max)
    call vd_get_indices(this%ny, nyl_ext, nyr_ext, this%nyl, this%nyr, &
         this%yl_min, this%yl_max, this%yr_min, this%yr_max)

    allocate(this%vd_p_arr(this%npx, this%npy))
    !$omp parallel do private(i_x, i_y)
    do i_y = 1, this%npx
       do i_x = 1, this%npy
          this%vd_p_arr(i_x, i_y) = 0.0_fp
       end do
    end do
    !$omp end parallel do

    allocate(this%px_range(this%npx))
    allocate(this%py_range(this%npy))
    call numerics_linspace(this%px_min, this%px_max, this%px_range, this%dpx)
    call numerics_linspace(this%py_min, this%py_max, this%py_range, this%dpy)

    ! Top and bottom arrays go all the way across
    ! Left and right array heights are truncated so corners don't overlap
    allocate(this%vdb_arr(this%xl_min:this%xr_max, this%yl_min:this%yl_max))
    allocate(this%vdt_arr(this%xl_min:this%xr_max, this%yr_min:this%yr_max))
    allocate(this%vdl_arr(this%xl_min:this%xl_max, this%yl_max + 1:this%yr_min - 1))
    allocate(this%vdr_arr(this%xr_min:this%xr_max, this%yl_max + 1:this%yr_min - 1))

    do i_y = this%yl_min, this%yl_max
       do i_x = this%xl_min, this%xr_max
          call this%vdb_arr(i_x, i_y)%init(this%dx, this%dy, this%npx, &
               this%npy, this%px_min, this%px_max, this%py_min, this%py_max, &
               this%dpx, this%dpy, this%px_range, this%py_range, this%dt, &
               this%j_eps, this%semi_classical, this%vd_disjoint, this%hbar, &
               this%m, vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
       end do
    end do

    do i_y = this%yr_min, this%yr_max
       do i_x = this%xl_min, this%xr_max
          call this%vdt_arr(i_x, i_y)%init(this%dx, this%dy, this%npx, &
               this%npy, this%px_min, this%px_max, this%py_min, this%py_max, &
               this%dpx, this%dpy, this%px_range, this%py_range, this%dt, &
               this%j_eps, this%semi_classical, this%vd_disjoint, this%hbar, &
               this%m, vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
       end do
    end do

    do i_y = this%yl_max + 1, this%yr_min - 1
       do i_x = this%xl_min, this%xl_max
          call this%vdl_arr(i_x, i_y)%init(this%dx, this%dy, this%npx, &
               this%npy, this%px_min, this%px_max, this%py_min, this%py_max, &
               this%dpx, this%dpy, this%px_range, this%py_range, this%dt, &
               this%j_eps, this%semi_classical, this%vd_disjoint, this%hbar, &
               this%m, vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
       end do

       do i_x = this%xr_min, this%xr_max
          call this%vdr_arr(i_x, i_y)%init(this%dx, this%dy, this%npx, &
               this%npy, this%px_min, this%px_max, this%py_min, this%py_max, &
               this%dpx, this%dpy, this%px_range, this%py_range, this%dt, &
               this%j_eps, this%semi_classical, this%vd_disjoint, this%hbar, &
               this%m, vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
       end do
    end do
  end subroutine vd_2d_manager_init

  subroutine vd_2d_manager_cleanup(this)
    class(vd_2d_manager_obj), intent(inout) :: this

    integer(ip) :: i_x, i_y

    do i_y = this%yl_min, this%yl_max
       do i_x = this%xl_min, this%xr_max
          call this%vdb_arr(i_x, i_y)%cleanup()
       end do
    end do
    deallocate(this%vdb_arr)

    do i_y = this%yl_max + 1, this%yr_min - 1
       do i_x = this%xl_min, this%xl_max
          call this%vdl_arr(i_x, i_y)%cleanup()
       end do

       do i_x = this%xr_min, this%xr_max
          call this%vdr_arr(i_x, i_y)%cleanup()
       end do
    end do
    deallocate(this%vdl_arr)
    deallocate(this%vdr_arr)

    do i_y = this%yr_min, this%yr_max
       do i_x = this%xl_min, this%xr_max
          call this%vdt_arr(i_x, i_y)%cleanup()
       end do
    end do
    deallocate(this%vdt_arr)

    deallocate(this%px_range)
    deallocate(this%py_range)

    deallocate(this%vd_p_arr)

  end subroutine vd_2d_manager_cleanup

  subroutine vd_2d_manager_update(this, psi_arr)
    class(vd_2d_manager_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:,:)

    integer(ip) :: i_x, i_y

    ! Update bottom VD objects
    do i_y = this%yl_min, this%yl_max
       !$omp parallel do private(i_x)
       do i_x = this%xl_min, this%xr_max
          call this%vdb_arr(i_x, i_y)%update( &
               psi_arr(i_x - 1 : i_x + 1, i_y - 1: i_y + 1) )
       end do
       !#omp end parallel do
    end do

    ! Update left / right VD objects
    do i_y = this%yl_max + 1, this%yr_min - 1
       !$omp parallel do private(i_x)
       do i_x = this%xl_min, this%xl_max
          call this%vdl_arr(i_x, i_y)%update( &
               psi_arr(i_x - 1 : i_x + 1, i_y - 1: i_y + 1) )
       end do
       !$omp end parallel do

       !$omp parallel do private(i_x)
       do i_x = this%xr_min, this%xr_max
          call this%vdr_arr(i_x, i_y)%update( &
               psi_arr(i_x - 1 : i_x + 1, i_y - 1: i_y + 1) )
       end do
       !$omp end parallel do
    end do

    ! Update top VD objects
    do i_y = this%yr_min, this%yr_max
       !$omp parallel do private(i_x)
       do i_x = this%xl_min, this%xr_max
          call this%vdt_arr(i_x, i_y)%update( &
               psi_arr(i_x - 1 : i_x + 1, i_y - 1: i_y + 1) )
       end do
       !$omp end parallel do
    end do

  end subroutine vd_2d_manager_update

  subroutine vd_2d_manager_finalize(this)
    class(vd_2d_manager_obj), intent(inout) :: this

    this%vd_p_arr = this%vd_p_arr * 1

    ! integer(ip) :: i_x, i_y

    ! real(fp) :: norm
    ! real(fp) :: norm_integrand(this%npx)

    ! Right now, we assume the VDs are disjoint. The commented code is for
    ! non-disjoint VDs.

    ! UNCOMMENT THIS TO DO NORMALIZATION
    ! !$omp parallel do private(i_x)
    ! do i_x = 1, this%npx
    !    norm_integrand(i_x) = numerics_trapz(this%vd_p_arr(i_x, :), this%dpy)
    ! end do
    ! !$omp end parallel do
    ! norm = numerics_trapz(norm_integrand(:), this%dpx)

    ! !$omp parallel do private(i_x, i_y)
    ! do i_y = 1, this%npy
    !    do i_x = 1, this%npx
    !       this%vd_p_arr(i_x, i_y) = this%vd_p_arr(i_x, i_y) / norm
    !    end do
    ! end do


    ! ! Get total flux
    ! do i_x = 1, this%nx
    !    do i_y = 1, this%nyl
    !       this%net_flux = this%net_flux + this%vdb_arr(i_x, i_y)%net_flux
    !    end do
    !    do i_y =1 , this%nyr
    !       this%net_flux = this%net_flux + this%vdt_arr(i_x, i_y)%net_flux
    !    end do
    ! end do
    ! do i_y = 1, this%ny - (this%nyl + this%nyr)
    !    do i_x = 1, this%nxl
    !       this%net_flux = this%net_flux + this%vdl_arr(i_x, i_y)%net_flux
    !    end do
    !    do i_x = 1, this%nxr
    !       this%net_flux = this%net_flux + this%vdr_arr(i_x, i_y)%net_flux
    !    end do
    ! end do

    ! ! Combine results
    ! do i_x = 1, this%nx
    !    do i_y = 1, this%nyl
    !       weight = this%vdb_arr(i_x, i_y)%net_flux / this%net_flux
    !       this%vd_p_arr(:,:) = this%vd_p_arr(:,:) + &
    !            weight * this%vdb_arr(i_x, i_y)%vd_p_arr(:,:)
    !    end do
    !    do i_y = 1, this%nyr
    !       weight = this%vdt_arr(i_x, i_y)%net_flux / this%net_flux
    !       this%vd_p_arr(:,:) = this%vd_p_arr(:,:) + &
    !            weight * this%vdt_arr(i_x, i_y)%vd_p_arr(:,:)
    !    end do
    ! end do

    ! do i_y = 1, this%ny - (this%nyl + this%nyr)
    !    do i_x = 1, this%nxl
    !       weight = this%vdl_arr(i_x, i_y)%net_flux / this%net_flux
    !       this%vd_p_arr(:,:) = this%vd_p_arr(:,:) + &
    !            weight * this%vdl_arr(i_x, i_y)%vd_p_arr(:,:)
    !    end do
    !    do i_x = 1, this%nxr
    !       weight = this%vdr_arr(i_x, i_y)%net_flux / this%net_flux
    !       this%vd_p_arr(:,:) = this%vd_p_arr(:,:) + &
    !            weight * this%vdr_arr(i_x, i_y)%vd_p_arr(:,:)
    !    end do
    ! end do

    ! Normalize distribution


  end subroutine vd_2d_manager_finalize

end module vd_2d_manager
