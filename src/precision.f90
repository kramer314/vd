! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module precision
  ! Specify floating point precision

  use globvars, only: dp, dp_format, dp_format_raw, sp, sp_format, &
       sp_format_raw, ip, ip_format, ip_format_raw

  implicit none

  ! Real precision kind parameter / formatting
  integer(ip), parameter :: fp = dp
  character(*), parameter :: fp_format = dp_format
  character(*), parameter :: fp_format_raw = dp_format_raw

  real(fp), parameter :: fp_eps = epsilon(1.0_fp)
  
end module precision
