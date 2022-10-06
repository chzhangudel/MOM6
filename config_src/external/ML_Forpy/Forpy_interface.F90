module Forpy_interface

use MOM_grid,                  only : ocean_grid_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_error_handler,         only : MOM_error, WARNING

implicit none; private

#include <MOM_memory.h>

public :: forpy_run_python_init,forpy_run_python,forpy_run_python_finalize

!> Control structure for Python interface
type, public :: python_interface ; private
end type

contains

!> Initialize Forpy with specify Python script and directory
subroutine forpy_run_python_init(CS,python_dir,python_file)
    character(len=*),         intent(in)  :: python_dir   !< The directory in which python scripts are found
    character(len=*),         intent(in)  :: python_file  !< The name of the Python script to read
    type(python_interface),   intent(inout) :: CS         !< Python interface object
    call MOM_error(WARNING,"Forpy was compiled using the dummy module. If this was &
                            a mistake, please follow the instructions in: &
                            MOM6/config_src/external/ML_Forpy/README.md")

end subroutine forpy_run_python_init

!> !> Send variables to a python script and output the results
subroutine forpy_run_python(WH_u, WH_v, Sx, Sy, G, GV, CS, wh_size_in, BT)
    type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
    type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
    type(python_interface),        intent(in)  :: CS     !< Python interface object
  ! Local Variables
    integer, intent(in) :: wh_size_in(4)  ! Subdomain size with wide halos for input
    logical, intent(in) :: BT             !< If true, momentum forcing from CNN is barotropic.
    real, dimension(wh_size_in(1):wh_size_in(2),&
                    wh_size_in(3):wh_size_in(4),&
                    SZK_(GV)), &
                                    intent(in) :: WH_u     ! The zonal velocity with a wide halo [L T-1 ~> m s-1].
    real, dimension(wh_size_in(1):wh_size_in(2),&
                    wh_size_in(3):wh_size_in(4),&
                    SZK_(GV)), &
                                    intent(in) :: WH_v     ! The meridional velocity with a wide halo [L T-1 ~> m s-1].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(out) :: Sx      ! CNN output Sx
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(out) :: Sy      ! CNN output Sy

end subroutine forpy_run_python 

!> Finalize Forpy
subroutine forpy_run_python_finalize(CS)
    type(python_interface), intent(inout) :: CS !< Python interface object

end subroutine forpy_run_python_finalize

end module Forpy_interface