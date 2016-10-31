module vd_2d
  use numerics, only: numerics_linspace, numerics_linspace_index, &
       numerics_trapz
  use dists, only: dists_gaussian
  use precision, only: ip, fp

  use vd, only: vd_get_local_quantities, vd_get_indices, &
       vd_validate_quantum_update, vd_obj

  implicit none

  private

  public vd_2d_obj
  type, extends(vd_obj) :: vd_2d_obj
     real(fp), pointer :: vd_p_arr(:,:)
     real(fp) :: dx, dy
     integer(ip) :: npx, npy
     real(fp) :: px_min, px_max, py_min, py_max
     real(fp) :: dpx, dpy
     real(fp), pointer :: px_range(:), py_range(:)
   contains
     procedure :: init => vd_2d_obj_init
     procedure :: cleanup => vd_2d_obj_cleanup
     procedure :: update => vd_2d_obj_update
  end type vd_2d_obj

contains

  subroutine vd_2d_obj_init(this, dx, dy, npx, npy, px_min, px_max, py_min, &
       py_max, dpx, dpy, px_range, py_range, dt, j_eps, sc, vd_disjoint, &
       hbar, m, vd_np_stdev, vd_p_arr)
    class(vd_2d_obj), intent(inout) :: this
    real(fp), intent(in) :: dx, dy
    integer(ip), intent(in) :: npx, npy
    real(fp), intent(in) :: px_min, px_max, py_min, py_max
    real(fp), intent(in) :: dpx, dpy
    real(fp), pointer, intent(in) :: px_range(:), py_range(:)
    real(fp), intent(in) :: dt
    real(fp), intent(in) :: j_eps
    logical, intent(in) :: sc
    logical, intent(in) :: vd_disjoint
    real(fp), intent(in) :: hbar, m
    integer(ip), intent(in), optional :: vd_np_stdev
    real(fp), pointer, intent(in), optional :: vd_p_arr(:, :)

    call this%vd_obj%init_vd_base(dt, j_eps, sc, vd_disjoint, hbar, m, &
         vd_np_stdev=vd_np_stdev)

    this%dx = dx
    this%dy = dy

    this%npx = npx
    this%npy = npy
    this%px_min = px_min
    this%px_max = px_max
    this%py_min = py_min
    this%py_max = py_max
    this%dpx = dpx
    this%dpy = dpy

    this%px_range => px_range
    this%py_range => py_range

    if (this%vd_disjoint) then
       this%vd_p_arr => vd_p_arr
    else
       allocate(this%vd_p_arr(this%npx, this%npy))
       this%vd_p_arr = 0.0_fp
    end if

  end subroutine vd_2d_obj_init

  subroutine vd_2d_obj_cleanup(this)
    class(vd_2d_obj), intent(inout) :: this

    nullify(this%px_range)
    nullify(this%py_range)

    if (this%vd_disjoint) then
       nullify(this%vd_p_arr)
    else
       deallocate(this%vd_p_arr)
    end if

  end subroutine vd_2d_obj_cleanup

  subroutine vd_2d_obj_update(this, psi_arr)
    class(vd_2d_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:,:)

    real(fp) :: px_mu, py_mu
    real(fp) :: px_var, py_var
    real(fp) :: jx, jy, j

    real(fp) :: px_min, px_max, py_min, py_max
    real(fp) :: px_stdev, py_stdev
    integer(ip) :: i_px_min, i_px_max, i_py_min, i_py_max
    logical :: valid_x, valid_y

    integer(ip) :: i_px, i_py
    real(fp) :: px, py
    real(fp) :: scale
    real(fp) :: gx, gy

    call vd_get_local_quantities(psi_arr(:, 2), this%dx, this%m, this%hbar, &
         px_mu, px_var, jx)
    call vd_get_local_quantities(psi_arr(2, :), this%dy, this%m, this%hbar, &
         py_mu, py_var, jy)

    j = sqrt(jx**2 + jy**2)
    scale = this%dt * j

    if (j .gt. this%j_eps .and. this%semi_classical) then
       i_px = numerics_linspace_index(px_mu, this%px_range)
       i_py = numerics_linspace_index(py_mu, this%py_range)

       if (i_px .gt. 0_ip .and. i_py .gt. 0_ip) then
          this%vd_p_arr(i_px, i_py) = this%vd_p_arr(i_px, i_py) + &
               scale / (this%dpx * this%dpy)
       end if

    else if (j .gt. this%j_eps) then

       px_stdev= sqrt(px_var)
       py_stdev = sqrt(py_var)

       px_min = px_mu - this%vd_np_stdev * px_stdev
       py_min = py_mu - this%vd_np_stdev * py_stdev

       px_max = px_mu + this%vd_np_stdev * px_stdev
       py_max = py_mu + this%vd_np_stdev * py_stdev

       i_px_min = numerics_linspace_index(px_min, this%px_range)
       i_px_max = numerics_linspace_index(px_max, this%px_range)

       i_py_min = numerics_linspace_index(py_min, this%py_range)
       i_py_max = numerics_linspace_index(py_max, this%py_range)

       call vd_validate_quantum_update(i_px_min, i_px_max, this%npx, valid_x)
       call vd_validate_quantum_update(i_py_min, i_py_max, this%npy, valid_y)

       if (valid_x .and. valid_y) then
          do i_py = i_py_min, i_py_max
             py = this%py_range(i_py)
             gy = dists_gaussian(py, py_mu, py_var)

             do i_px = i_px_min, i_px_max
                px = this%px_range(i_px)
                gx = dists_gaussian(px, px_mu, px_var)

                this%vd_p_arr(i_px, i_py) = this%vd_p_arr(i_px, i_py) + &
                     scale * gx * gy / (this%dpx * this%dpy)
             end do
          end do
       end if
    end if

    this%net_flux = this%net_flux + scale
  end subroutine vd_2d_obj_update

end module vd_2d
