! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module vd_1d_manager

  use numerics, only: numerics_linspace, numerics_trapz
  use precision, only: ip, fp

  use vd_1d, only: vd_1d_obj
  use vd, only: vd_get_indices

  implicit none

  private

  public vd_1d_manager_obj
  type vd_1d_manager_obj
     ! 1-dimensional VD manager

     ! Virtual detector result; array of momentum counts
     real(fp), pointer :: vd_p_arr(:)

     ! Semi-classical binning flag
     logical :: semi_classical
     ! Number of standard deviations to include in quantum VD binning
     integer(ip) :: vd_np_stdev

     logical :: vd_disjoint

     ! Number of left / right virtual detector points
     integer(ip) :: nxl, nxr

     ! Spatial grid parameters
     real(fp) :: dx
     integer(ip) :: nx

     ! Spatial grid index bounds for left / right virtual detector points
     integer(ip) :: xl_min, xl_max
     integer(ip) :: xr_min, xr_max

     ! Momentum grid parameters
     integer(ip) :: np
     real(fp) :: p_min, p_max
     real(fp) :: dp

     ! Momentum grid; we use an allocatable pointer here since the VD objects
     ! this class managers need this array as well, but we don't want to
     ! allocate a new array for each VD object.
     real(fp), pointer :: p_range(:)

     ! Temporal grid step
     real(fp) :: dt

     ! Flux threshold
     real(fp) :: j_eps

     ! Units
     real(fp) :: hbar
     real(fp) :: m

     ! Total probability flux through all virtual detectors
     real(fp) :: net_flux

     ! Array of left/right virtual detector point objects
     type(vd_1d_obj), allocatable :: vdl_arr(:), vdr_arr(:)

   contains
     ! Initialize object
     procedure :: init => vd_1d_manager_init
     ! Cleanup object, deallocating memory
     procedure :: cleanup => vd_1d_manager_cleanup
     ! Update virtual detector point counts
     procedure :: update => vd_1d_manager_update
     ! Combine virtual detector point counts and obtain final result
     procedure :: finalize => vd_1d_manager_finalize

  end type vd_1d_manager_obj

contains

  subroutine vd_1d_manager_init(this, nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd, dx, np, &
       p_min, p_max, dt, j_eps, sc, vd_disjoint, hbar, m, vd_np_stdev)
    ! Initialize 1D VD manager object
    !
    ! This method is exposed as vd_1d_manager%init
    !
    ! this :: vd_1d_manager instance
    ! nx :: length of spatial grid
    ! nxl_ext :: number of left external spatial points
    ! nxr_ext :: number of right external spatial points
    ! nxl_vd :: number of left VD points
    ! nxr_vd :: number of right VD points
    ! dx :: spatial grid step
    ! np :: length of momentum grid
    ! p_min :: lower bound of momentum grid
    ! p_max :: upper bound of momentum grid
    ! dt :: temporal grid step
    ! sc :: semi-classical binning flag
    ! hbar :: hbar units
    ! m :: particle mass
    ! vd_disjoint :: whether to treat detectors as disjoint or not
    !   treating detectors as non-disjoint means each VD allocates their own
    !   momentum distribution counts which are recombined in this%finalize()
    !   so it uses (potentially drastically) more memory.
    ! vd_np_stdev :: REQUIRED if sc = False, number of standard deviations to
    !   include in quantum VD binning
    class(vd_1d_manager_obj), intent(inout) :: this
    integer(ip), intent(in) :: nx
    integer(ip), intent(in) :: nxl_ext, nxr_ext
    integer(ip), intent(in) :: nxl_vd, nxr_vd
    real(fp), intent(in) :: dx
    integer(ip), intent(in) :: np
    real(fp), intent(in) :: p_min, p_max
    real(fp), intent(in) :: dt
    real(fp), intent(in) :: j_eps
    logical, intent(in) :: sc
    logical, intent(in) :: vd_disjoint
    real(fp), intent(in) :: hbar, m
    integer(ip), intent(in), optional :: vd_np_stdev

    integer(ip) :: i_x

    this%nx = nx
    this%dx = dx
    this%nxl = nxl_vd
    this%nxr = nxr_vd

    this%np = np
    this%p_min = p_min
    this%p_max = p_max

    this%dt=  dt
    this%j_eps = j_eps

    this%semi_classical = sc

    this%vd_disjoint = vd_disjoint

    if (present(vd_np_stdev)) then
       this%vd_np_stdev = vd_np_stdev
    end if

    this %hbar = hbar
    this%m = m

    this%net_flux = 0.0_fp

    ! Verify and calculate VD grid edge indices on the total spatial grid
    call vd_get_indices(this%nx, nxl_ext, nxr_ext, this%nxl, this%nxr, &
         this%xl_min, this%xl_max, this%xr_min, this%xr_max)

    ! Construct momentum grid and count arrays
    allocate(this%vd_p_arr(this%np))
    this%vd_p_arr(:) = 0.0_fp

    allocate(this%p_range(this%np))
    call numerics_linspace(this%p_min, this%p_max, this%p_range, this%dp)

    ! Construct VD arrays
    allocate(this%vdl_arr(this%nxl))
    allocate(this%vdr_arr(this%nxr))

    ! Initialize individual VD objects
    do i_x = 1, this%nxl
       call this%vdl_arr(i_x)%init(this%dx, this%np, this%p_min, this%p_max, &
            this%dp, this%p_range, this%dt, this%j_eps, this%semi_classical, &
            this%vd_disjoint, this%hbar, this%m, &
            vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
    end do

    do i_x = 1, this%nxr
       call this%vdr_arr(i_x)%init(this%dx, this%np, this%p_min, this%p_max, &
            this%dp, this%p_range, this%dt, this%j_eps, this%semi_classical, &
            this%vd_disjoint, this%hbar, this%m, &
            vd_np_stdev=this%vd_np_stdev, vd_p_arr=this%vd_p_arr)
    end do
  end subroutine vd_1d_manager_init

  subroutine vd_1d_manager_update(this, psi_arr)
    ! Update VD counts
    !
    ! This method is exposed as vd_1d_manager%update
    !
    ! this :: vd_1d_manager instancex
    ! psi_arr :: wavefunction
    class(vd_1d_manager_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:)

    integer(ip) :: i_x
    integer(ip) :: psi_i_x

    ! Update left VD objects
    !$omp parallel do private(i_x, psi_i_x)
    do i_x = 1, this%nxl
       psi_i_x = this%xl_min + (i_x - 1)
       call this%vdl_arr(i_x)%update(psi_arr(psi_i_x - 1 : psi_i_x + 1))
    end do
    !$omp end parallel do

    ! Update right VD objects
    !$omp parallel do private(i_x, psi_i_x)
    do i_x = 1, this%nxr
       psi_i_x = this%xr_min + (i_x - 1)
       call this%vdr_arr(i_x)%update(psi_arr(psi_i_x - 1 : psi_i_x + 1))
    end do
    !$omp end parallel do

  end subroutine vd_1d_manager_update

  subroutine vd_1d_manager_cleanup(this)
    ! Cleanup / deallocate internal arrays
    !
    ! This method is exposed as vd_1d_manager%cleanup
    !
    ! this :: vd_1d_manager instance
    class(vd_1d_manager_obj), intent(inout) :: this

    integer(ip) :: i_x

    ! Cleanup left VDs
    do i_x = 1, this%nxl
       call this%vdl_arr(i_x)%cleanup()
    end do
    deallocate(this%vdl_arr)

    ! Cleanup right VDs
    do i_x = 1, this%nxr
       call this%vdr_arr(i_x)%cleanup()
    end do
    deallocate(this%vdr_arr)

    deallocate(this%p_range)
    deallocate(this%vd_p_arr)

  end subroutine vd_1d_manager_cleanup

  subroutine vd_1d_manager_finalize(this)
    ! Combine VD point counts and calculate normalized total momentum
    ! distribution counts
    !
    ! this :: vd_1d_manager instance
    class(vd_1d_manager_obj), intent(inout) :: this

    integer(ip) :: i_x

    real(fp) :: weight_i_x
    real(fp) :: norm

    if (.not. this%vd_disjoint) then
       ! Get total flux
       do i_x = 1, this%nxl
          this%net_flux = this%net_flux + this%vdl_arr(i_x)%net_flux
       end do

       do i_x = 1, this%nxr
          this%net_flux = this%net_flux + this%vdr_arr(i_x)%net_flux
       end do

       ! Combine results
       do i_x = 1, this%nxl
          weight_i_x = this%vdl_arr(i_x)%net_flux / this%net_flux
          this%vd_p_arr(:) = this%vd_p_arr(:) + &
               weight_i_x * this%vdl_arr(i_x)%vd_p_arr(:)
       end do

       do i_x = 1, this%nxr
          weight_i_x = this%vdr_arr(i_x)%net_flux / this%net_flux
          this%vd_p_arr(:) = this%vd_p_arr(:) + &
               weight_i_x * this%vdr_arr(i_x)%vd_p_arr(:)
       end do
    end if

    ! Normalize distribution
    norm = numerics_trapz(this%vd_p_arr, this%dp)
    this%vd_p_arr(:) = this%vd_p_arr(:) / norm

  end subroutine vd_1d_manager_finalize

end module vd_1d_manager
