module vd_1d
  use numerics, only: numerics_linspace, numerics_linspace_index, &
       numerics_trapz
  use dists, only: dists_gaussian
  use precision, only: ip, fp

  use vd, only: vd_get_local_quantities, vd_get_indices, &
       vd_validate_quantum_update, vd_obj

  implicit none

  private

  public vd_1d_obj
  type, extends(vd_obj) :: vd_1d_obj
     real(fp), pointer :: vd_p_arr(:)
     real(fp) :: dx
     integer(ip) :: npx
     real(fp) :: px_min, px_max
     real(fp) :: dpx
     real(fp), pointer :: px_range(:)
   contains
     procedure :: init => vd_1d_obj_init
     procedure :: cleanup => vd_1d_obj_cleanup
     procedure :: update => vd_1d_obj_update
  end type vd_1d_obj

contains

  subroutine vd_1d_obj_init(this, dx, npx, px_min, px_max,dpx, px_range, dt, &
       j_eps, sc, vd_disjoint, hbar, m, vd_np_stdev, vd_p_arr)

    class(vd_1d_obj), intent(inout) :: this
    real(fp), intent(in) :: dx
    integer(ip), intent(in) :: npx
    real(fp), intent(in) :: px_min, px_max
    real(fp), intent(in) :: dpx
    real(fp), pointer, intent(in) :: px_range(:)
    real(fp), intent(in) :: dt
    real(fp), intent(in) :: j_Eps
    logical, intent(in) :: sc
    logical, intent(in) :: vd_disjoint
    real(fp), intent(in) :: hbar, m
    integer(ip), intent(in), optional :: vd_np_stdev
    real(fp), pointer, intent(in), optional :: vd_p_arr(:)

    call this%vd_obj%init_vd_base(dt, j_eps, sc, vd_disjoint, hbar, m, &
         vd_np_stdev=vd_np_stdev)

    this%dx = dx

    this%npx = npx
    this%px_min = px_min
    this%px_max = px_max
    this%dpx = dpx

    this%px_range => px_range

    if (this%vd_disjoint) then
       this%vd_p_arr => vd_p_arr
    else
       allocate(this%vd_p_arr(this%npx))
       this%vd_p_arr = 0.0_fp
    end if

  end subroutine vd_1d_obj_init

  subroutine vd_1d_obj_cleanup(this)
    class(vd_1d_obj), intent(inout) :: this

    nullify(this%px_range)

    if (this%vd_disjoint) then
       nullify(this%vd_p_arr)
    else
       deallocate(this%vd_p_arr)
    end if

  end subroutine vd_1d_obj_cleanup

  subroutine vd_1d_obj_update(this, psi_arr)
    class(vd_1d_obj), intent(inout) :: this

    complex(fp), intent(in) :: psi_arr(:)

    real(fp) :: px_mu, px_var
    real(fp) :: jx

    real(fp) :: px_min, px_max, px_stdev

    integer(ip) :: i_px_min, i_px_max
    logical :: valid

    integer(ip) :: i_px
    real(fp) :: px
    real(fp) :: scale
    real(fp) :: gx

    call vd_get_local_quantities(psi_arr(:), this%dx, this%m, this%hbar, &
         px_mu, px_var, jx)
    scale = this%dt * abs(jx)

    if (jx .gt. this%j_eps .and. this%semi_classical) then

       i_px = numerics_linspace_index(px_mu, this%px_range)
       if (i_px .gt. 0_ip) then
          this%vd_p_arr(i_px) = this%vd_p_arr(i_px) + scale / this%dpx
       end if

    else if (jx .gt. this%j_eps) then

       px_stdev = sqrt(px_var)
       px_min = px_mu - this%vd_np_stdev * px_stdev
       px_max = px_mu + this%vd_np_stdev * px_stdev

       i_px_min = numerics_linspace_index(px_min, this%px_range)
       i_px_max = numerics_linspace_index(px_max, this%px_range)

       call vd_validate_quantum_update(i_px_min, i_px_max, this%npx, valid)

       if (valid) then

          if (i_px_min .eq. i_px_max) then

             this%vd_p_arr(i_px_min) = this%vd_p_arr(i_px_min) + scale / this%dpx

          else

             do i_px = i_px_min, i_px_max
                px = this%px_range(i_px)
                gx = dists_gaussian(px, px_mu, px_var)
                this%vd_p_arr(i_px) = this%vd_p_arr(i_px) + scale * gx
             end do
          end if

          
       end if

    end if
  end subroutine vd_1d_obj_update
end module vd_1d
