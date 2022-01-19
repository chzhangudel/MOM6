module forpy_util

use forpy_mod
use iso_fortran_env, only: real64

implicit none ; private

public :: forpy_run_python_init,forpy_run_python,forpy_run_python_finalize

integer :: ierror
type(module_py) :: f2p2f
type(list) :: paths

contains

!> Send a variable to a python script and output the results
subroutine forpy_run_python_init()

  ierror = forpy_initialize()
  write(*,*) "############ Initialize Forpy ############"
  ierror = get_sys_path(paths)
  ierror = paths%append('/scratch/cimes/cz3321/MOM6/MOM6-examples/src/MOM6/config_src/external/ML_Forpy/')
  ierror = import_py(f2p2f, "f2p2f")
  if (ierror/=0) then; call err_print; endif
  ierror = print_py(f2p2f)
  if (ierror/=0) then; call err_print; endif

end subroutine forpy_run_python_init

!> Send a variable to a python script and output the results
subroutine forpy_run_python

  type(ndarray) :: w0_py,w1_py
  type(object) :: return_obj
  type(tuple) :: args

  real(kind=real64), dimension(:,:), pointer :: w1
  
  real(kind=real64), dimension(6,4) :: w0 
  real(kind=real64) :: wmin 
  w0 = 1.0
  print*,' create a matrix >'
  print*, w0

  ierror = ndarray_create(w0_py, w0)
  if (ierror/=0) then; call err_print; endif
  
  ierror = tuple_create(args,1)
  if (ierror/=0) then; call err_print; endif
  ierror = args%setitem(0,w0_py)
  if (ierror/=0) then; call err_print; endif
  ierror = call_py(return_obj, f2p2f, "wb_out", args)
  if (ierror/=0) then; call err_print; endif
  ierror = cast(w1_py, return_obj)
  if (ierror/=0) then; call err_print; endif
  ierror = w1_py%get_data(w1)
  if (ierror/=0) then; call err_print; endif
  call args%destroy
  call return_obj%destroy
  call w1_py%destroy
  call w0_py%destroy
  
  print*,' transfer back to fortran >'
  print*, w1
  wmin = minval(abs(w1-w0))
  print*,' Difference before and after passing into Python >',wmin

end subroutine forpy_run_python

subroutine forpy_run_python_finalize()
  write(*,*) "############ Finalize Forpy ############"
  call f2p2f%destroy
  call paths%destroy
 
  call forpy_finalize

end subroutine forpy_run_python_finalize 
  
end module forpy_util