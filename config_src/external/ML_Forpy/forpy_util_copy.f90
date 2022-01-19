module forpy_util

use forpy_mod
use iso_fortran_env, only: real64

implicit none ; private

public :: forpy_run_python

integer :: ierror
type(module_py) :: f2p2f
type(list) :: paths

contains

!> Send a variable to a python script and output the results
subroutine forpy_run_python(diffu,diffv)
  real,dimension(:,:,:), intent(inout):: diffu,diffv
  ! Local Variables
  type(ndarray) :: diffu_py,diffv_py,diffu_arr,diffv_arr
  type(object) :: obj_diffu,obj_diffv
  type(tuple) :: args_u,args_v
  real, dimension(:,:,:), pointer :: diffu_for,diffv_for
  real :: io_diffu_min,io_diffv_min

  ! initialize forpy
  ierror = forpy_initialize()
  ierror = get_sys_path(paths)
  ierror = paths%append(".")
  ierror = import_py(f2p2f, "f2p2f")
  ierror = print_py(f2p2f)
  
  ! Covert input into Forpy Numpy Arrays 
  ierror = ndarray_create(diffu_py, diffu)
  ierror = ndarray_create(diffv_py, diffv)

  ! Create Python Argument 
  ierror = tuple_create(args_u,1)
  ierror = tuple_create(args_v,1)
  ierror = args_u%setitem(0,diffu_py)
  ierror = args_v%setitem(0,diffv_py)

  ! Invoke Python 
  ierror = call_py(obj_diffu, f2p2f, "iou_py", args_u)
  ierror = call_py(obj_diffv, f2p2f, "iov_py", args_v)
  ierror = cast(diffu_arr, obj_diffu)
  ierror = cast(diffv_arr, obj_diffv)
  ierror = diffu_arr%get_data(diffu_for)
  ierror = diffv_arr%get_data(diffv_for)

  ! Destroy Objects
  call diffu_py%destroy
  call diffv_py%destroy
  call diffu_arr%destroy
  call diffv_arr%destroy
  call obj_diffu%destroy
  call obj_diffv%destroy
  call args_u%destroy
  call args_v%destroy

  ! Difference before and after passing into Python
  io_diffu_min  = minval(abs(diffu-diffu_for))
  io_diffv_min  = minval(abs(diffv-diffv_for))
  print*,' diffu before and after passing into Python >',io_diffu_min
  print*,' diffv before and after passing into Python >',io_diffv_min

  ! Output
  diffu = diffu_for
  diffv = diffv_for

  call f2p2f%destroy
  call paths%destroy

  call forpy_finalize

end subroutine forpy_run_python
  
end module forpy_util
 
