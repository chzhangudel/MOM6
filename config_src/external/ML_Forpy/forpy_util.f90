module forpy_util

use forpy_mod,                 only : module_py,list,ndarray,object,tuple
use forpy_mod,                 only : err_print
use forpy_mod,                 only : forpy_initialize,get_sys_path,import_py,print_py
use forpy_mod,                 only : ndarray_create,tuple_create,call_py,cast
use forpy_mod,                 only : forpy_finalize
use iso_fortran_env,           only : real64

implicit none ; private

public :: forpy_run_python_init,forpy_run_python,forpy_run_python_finalize

!> Control structure for Python interface
type, public :: python_interface ; private
  type(module_py) :: pymodule
  type(list) :: paths
end type

contains

!> Send a variable to a python script and output the results
subroutine forpy_run_python_init(CS,python_dir,python_file)
  character(len=*),         intent(in)  :: python_dir   !< The directory in which python scripts are found
  character(len=*),         intent(in)  :: python_file  !< The name of the Python script to read
  type(python_interface), intent(inout) :: CS !< Python interface object
  ! Local Variables
  integer :: ierror ! return code from python interfaces
  ierror = forpy_initialize()
  write(*,*) "############ Initialize Forpy ############"
  ierror = get_sys_path(CS%paths)
  ierror = CS%paths%append(python_dir)
  ierror = import_py(CS%pymodule,python_file)
  if (ierror/=0) then; call err_print; endif
  ierror = print_py(CS%pymodule)
  if (ierror/=0) then; call err_print; endif

end subroutine forpy_run_python_init

!> Send a variable to a python script and output the results
subroutine forpy_run_python(CS,diffu,diffv)
  type(python_interface), intent(in) :: CS !< Python interface object
  real,dimension(:,:,:), intent(inout):: diffu !< u-acceleration due to diffusion [m s-1]
  real,dimension(:,:,:), intent(inout):: diffv !< v-acceleration due to diffusion [m s-1]
  ! Local Variables
  integer :: ierror ! return code from python interfaces
  type(ndarray) :: diffu_py,diffv_py,diffu_arr,diffv_arr !< diffu in the form of numpy array
  type(object)  :: obj_diffu,obj_diffv                   !< return objects
  type(tuple)   :: args_u,args_v                         !< input arguments for the Python module
  real, dimension(:,:,:), pointer :: diffu_for,diffv_for !< outputs from Python module
  real :: io_diffu_max,io_diffv_max

  ! Covert input into Forpy Numpy Arrays 
  ierror = ndarray_create(diffu_py, diffu)
  ierror = ndarray_create(diffv_py, diffv)
  if (ierror/=0) then; call err_print; endif
  
  ! Create Python Argument 
  ierror = tuple_create(args_u,1)
  ierror = tuple_create(args_v,1)
  if (ierror/=0) then; call err_print; endif
  ierror = args_u%setitem(0,diffu_py)
  ierror = args_v%setitem(0,diffv_py)
  if (ierror/=0) then; call err_print; endif

  ! Invoke Python 
  ierror = call_py(obj_diffu, CS%pymodule, "iou_py", args_u)
  ierror = call_py(obj_diffv, CS%pymodule, "iov_py", args_v)
  if (ierror/=0) then; call err_print; endif
  ierror = cast(diffu_arr, obj_diffu)
  ierror = cast(diffv_arr, obj_diffv)
  if (ierror/=0) then; call err_print; endif
  ierror = diffu_arr%get_data(diffu_for)
  ierror = diffv_arr%get_data(diffv_for)
  if (ierror/=0) then; call err_print; endif

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
  io_diffu_max  = maxval(abs(diffu-diffu_for))
  io_diffv_max  = maxval(abs(diffv-diffv_for))
  print*,' diffu before and after passing into Python >',io_diffu_max
  print*,' diffv before and after passing into Python >',io_diffv_max
  
  ! Output
  diffu = diffu_for
  diffv = diffv_for

end subroutine forpy_run_python

subroutine forpy_run_python_finalize(CS)
  type(python_interface), intent(inout) :: CS !< Python interface object
  write(*,*) "############ Finalize Forpy ############"
  call CS%pymodule%destroy
  call CS%paths%destroy
 
  call forpy_finalize

end subroutine forpy_run_python_finalize 
  
end module forpy_util