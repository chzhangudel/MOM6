module forpy_util

use MOM_grid,                  only : ocean_grid_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_domains,               only : clone_MOM_domain,MOM_domain_type
use MOM_domains,               only : pass_vector,pass_var,BGRID_NE
use MOM_domains,               only : create_group_pass,do_group_pass,group_pass_type
use MOM_domains,               only : To_North, To_East
use MOM_diag_mediator,         only : post_data,register_diag_field
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_coms,                  only : PE_here,num_PEs,reproducing_sum
use MOM_unit_scaling,          only : unit_scale_type
use MOM_file_parser,           only : get_param,param_file_type
use forpy_mod,                 only : module_py,list,ndarray,object,tuple
use forpy_mod,                 only : err_print
use forpy_mod,                 only : forpy_initialize,get_sys_path,import_py,print_py
use forpy_mod,                 only : ndarray_create,tuple_create,call_py,cast
use forpy_mod,                 only : forpy_finalize
use iso_fortran_env,           only : real64

implicit none ; private

#include <MOM_memory.h>
#ifdef STATIC_MEMORY_
#  ifndef BTHALO_
#    define BTHALO_ 0
#  endif
#  define WHALOI_ MAX(BTHALO_-NIHALO_,0)
#  define WHALOJ_ MAX(BTHALO_-NJHALO_,0)
#  define NIMEMW_   1-WHALOI_:NIMEM_+WHALOI_
#  define NJMEMW_   1-WHALOJ_:NJMEM_+WHALOJ_
#  define NIMEMBW_  -WHALOI_:NIMEM_+WHALOI_
#  define NJMEMBW_  -WHALOJ_:NJMEM_+WHALOJ_
#  define SZIW_(G)  NIMEMW_
#  define SZJW_(G)  NJMEMW_
#  define SZIBW_(G) NIMEMBW_
#  define SZJBW_(G) NJMEMBW_
#else
#  define NIMEMW_   :
#  define NJMEMW_   :
#  define NIMEMBW_  :
#  define NJMEMBW_  :
#  define SZIW_(G)  G%isdw:G%iedw
#  define SZJW_(G)  G%jsdw:G%jedw
#  define SZIBW_(G) G%isdw-1:G%iedw
#  define SZJBW_(G) G%jsdw-1:G%jedw
#endif

public :: CNN_init,forpy_run_python_init,forpy_run_python,forpy_run_python_finalize

!> Control structure for Python interface
type, public :: python_interface ; private
  type(module_py) :: pymodule
  type(list) :: paths
end type

!> Control structure for CNN
type, public :: CNN_CS ; private
  type(MOM_domain_type), pointer :: CNN_Domain => NULL()  !< Domain for inputs/outputs for CNN
  integer :: isdw !< The lower i-memory limit for the wide halo arrays.
  integer :: iedw !< The upper i-memory limit for the wide halo arrays.
  integer :: jsdw !< The lower j-memory limit for the wide halo arrays.
  integer :: jedw !< The upper j-memory limit for the wide halo arrays.
  integer :: CNN_halo_size  !< Halo size at each side of subdomains

  logical :: CNN_BT    !< If true, momentum forcing from CNN is barotropic.

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_CNNu = -1, id_CNNv = -1, id_KE_CNN = -1
  !>@}
end type CNN_CS

contains

!> Prepare CNN input variables with wide halos
subroutine CNN_init(Time,G,GV,US,param_file,diag,CS)
  type(time_type),               intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),         intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target,       intent(inout) :: diag  !< Diagnostics structure.
  type(CNN_CS),                  intent(inout) :: CS    !< Control structure for CNN
  ! Local Variables
  integer :: wd_halos(2) ! Varies with CNN
  character(len=40)  :: mdl = "MOM_CNN"  ! module name

  ! Register fields for output from this module.
  CS%diag => diag

  CS%id_CNNu = register_diag_field('ocean_model', 'CNNu', diag%axesCuL, Time, &
      'Zonal Acceleration from CNN model', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_CNNv = register_diag_field('ocean_model', 'CNNv', diag%axesCvL, Time, &
      'Meridional Acceleration from CNN model', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_KE_CNN = register_diag_field('ocean_model', 'KE_CNN', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity of CNN model', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
      
  call get_param(param_file, mdl, "CNN_BT", CS%CNN_BT, &
      "If true, momentum forcing from CNN is barotropic, otherwise baroclinic (default).", &
      default=.false.)
  call get_param(param_file, mdl, "CNN_HALO_SIZE", CS%CNN_halo_size, &
      "Halo size at each side of subdomains, depends on CNN architecture.", & 
      units="nondim", default=10)

  wd_halos(1) = CS%CNN_halo_size
  wd_halos(2) = CS%CNN_halo_size
  call clone_MOM_domain(G%Domain, CS%CNN_Domain, min_halo=wd_halos, symmetric=.true.)
  CS%isdw = G%isc-wd_halos(1) ; CS%iedw = G%iec+wd_halos(1)
  CS%jsdw = G%jsc-wd_halos(2) ; CS%jedw = G%jec+wd_halos(2)

end subroutine CNN_init

!> Initialize Forpy with specify Python script and directory
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
subroutine forpy_run_python(u, v, h, diffu, diffv, G, GV, CS, CNN)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(python_interface),        intent(in)  :: CS     !< Python interface object
  type(CNN_CS),                  intent(in)  :: CNN    !< Control structure for CNN
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: diffu  !< Zonal acceleration due to convergence of
                                                         !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(inout) :: diffv  !< Meridional acceleration due to convergence
                                                         !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  ! Local Variables
  real, dimension(SZIW_(CNN),SZJW_(CNN),SZK_(GV)) &
                                    :: WH_u     ! The zonal velocity with a wide halo [L T-1 ~> m s-1].
  real, dimension(SZIW_(CNN),SZJW_(CNN),SZK_(GV)) &
                                    :: WH_v     ! The meridional velocity with a wide halo [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: Sx,  & ! CNN output Sx
                                               Sy     ! CNN output Sy
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)):: fx     ! CNN output Sx at cell faces
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)):: fy     ! CNN output Sy at cell faces
  integer :: i, j, k
  integer :: is, ie, js, je, nz, nztemp
  integer :: isdw, iedw, jsdw, jedw
  ! Local Variables for Forpy
  integer :: ierror ! return code from python interfaces
  type(ndarray) :: u_py,v_py,out_arr      !< u and v in the form of numpy array
  type(object)  :: obj                    !< return objects
  type(tuple)   :: args                   !< input arguments for the Python module
  real, dimension(:,:,:,:), pointer :: out_for  !< outputs from Python module
  integer :: current_pe
!  CHARACTER(LEN=80)::FILE_NAME=' '
!  CHARACTER(LEN=80)::TMP_NAME=' '

  current_pe = PE_here()
!  print*, 'current PE is ',current_pe 

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isdw = CNN%isdw; iedw = CNN%iedw; jsdw = CNN%jsdw; jedw = CNN%jedw
  if (CNN%CNN_BT) then
    nztemp = 1
  else
    nztemp = nz
  endif

  WH_u = 0.0; WH_v = 0.0;
  do k=1,nztemp
    do j=js,je ; do i=is,ie
      WH_u(i,j,k) = 0.5*( u(I-1,j,k) + u(I,j,k) ) ! Copy the computational section from u into cell center
    enddo ; enddo
    do j=js,je ; do i=is,ie
      WH_v(i,j,k) = 0.5*( v(i,J-1,k) + v(i,J,k) ) ! Copy the computational section from v into cell center
    enddo ; enddo
  enddo

!  write(file_name(1:1),'(I1)')current_pe 
!  TMP_NAME = 'WH_u0_'//TRIM(FILE_NAME)
!  open(10,file=TMP_NAME)
!  do j = jsdw,jedw
!    write(10,100) (WH_u(i,j,1),i=isdw,iedw)
!  enddo
!  close(10)

!  call pass_vector(WH_u, CNN_Domain, BGRID_NE) ! Update the wide halos of WH_u WH_v
  call pass_var(WH_u, CNN%CNN_Domain)
  call pass_var(WH_v, CNN%CNN_Domain)

!  TMP_NAME = 'WH_u_'//TRIM(FILE_NAME)
!  open(10,file=TMP_NAME)
!  do j = jsdw,jedw
!    write(10,100) (WH_u(i,j,1),i=isdw,iedw)
!  enddo
!  close(10)
!  100 FORMAT(5000es15.4)

  ! Covert input into Forpy Numpy Arrays 
  if (CNN%CNN_BT) then
    ierror = ndarray_create(u_py, WH_u(:,:,1))
    ierror = ndarray_create(v_py, WH_v(:,:,1))
  else
    ierror = ndarray_create(u_py, WH_u)
    ierror = ndarray_create(v_py, WH_v)
  endif
  if (ierror/=0) then; call err_print; endif

  ! Create Python Argument 
  ierror = tuple_create(args,4)
  if (ierror/=0) then; call err_print; endif
  ierror = args%setitem(0,u_py)
  ierror = args%setitem(1,v_py)
  ierror = args%setitem(2,current_pe)
  ierror = args%setitem(3,num_PEs())
  if (ierror/=0) then; call err_print; endif

  ! Invoke Python 
  ierror = call_py(obj, CS%pymodule, "MOM6_testNN", args)
  if (ierror/=0) then; call err_print; endif
  ierror = cast(out_arr, obj)
  if (ierror/=0) then; call err_print; endif
  ierror = out_arr%get_data(out_for, order='C')
  if (ierror/=0) then; call err_print; endif

  ! Destroy Objects
  call u_py%destroy
  call v_py%destroy
  call out_arr%destroy
  call obj%destroy
  call args%destroy

  ! Output (out_for in C order has index (nk,nj,ni))
                  ! in F order has index (ni,nj,nk)
  Sx = 0.0; Sy = 0.0;
  do k=1,nztemp
    do j=js,je ; do i=is,ie
      Sx(i,j,k) = out_for(k,j-js+1,i-is+1,1) ! if order='C'
      ! Sx(i,j,k) = out_for(1,i-is+1,j-js+1,k) ! if order='F'
      Sy(i,j,k) = out_for(k,j-js+1,i-is+1,2)
    enddo ; enddo
  enddo
  call pass_var(Sx, G%domain)
  call pass_var(Sy, G%domain)

!  open(10,file='Sx')
!  do j = G%jsd,G%jed
!    write(10,100) (Sx(i,j,1),i=G%isd,G%ied)
!  enddo
!  close(10)
 
  fx = 0.0; fy = 0.0;
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      if (CNN%CNN_BT) then
        fx(I,j,k) = 0.5*(Sx(i,j,1) + Sx(i+1,j,1))
      else
        fx(I,j,k) = 0.5*(Sx(i,j,k) + Sx(i+1,j,k))
      endif
      diffu(I,j,k) = diffu(I,j,k) + fx(I,j,k) ! Update diffu with Sx
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      if (CNN%CNN_BT) then
        fy(i,J,k) = 0.5*(Sy(i,j,1) + Sy(i,j+1,1))
      else
        fy(i,J,k) = 0.5*(Sy(i,j,k) + Sy(i,j+1,k))
      endif
      diffv(i,J,k) = diffv(i,J,k) + fy(i,J,k) ! Update diffv with Sy
    enddo ; enddo
  enddo

  if (CNN%id_CNNu>0)   call post_data(CNN%id_CNNu, fx, CNN%diag)
  if (CNN%id_CNNv>0)   call post_data(CNN%id_CNNv, fy, CNN%diag)
  if (CNN%id_KE_CNN>0) call compute_energy_source(u, v, h, fx, fy, G, GV, CNN)

!  stop'debugging!'

end subroutine forpy_run_python

!> Finalize Forpy
subroutine forpy_run_python_finalize(CS)
  type(python_interface), intent(inout) :: CS !< Python interface object
  write(*,*) "############ Finalize Forpy ############"
  call CS%pymodule%destroy
  call CS%paths%destroy
 
  call forpy_finalize

end subroutine forpy_run_python_finalize 

! This is copy-paste from MOM_diagnostics.F90, specifically 1125 line
subroutine compute_energy_source(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(CNN_CS),                  intent(in)  :: CS    !< Control structure for CNN

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: fx     !< Zonal acceleration due to convergence of
                                                      !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in) :: fy     !< Meridional acceleration due to convergence
                                                      !! of along-coordinate stress tensor [L T-2 ~> m s-2]
  
  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real :: tmp(SZI_(G),SZJ_(G),SZK_(GV)) ! temporary array for integration
  real :: KE_u(SZIB_(G),SZJ_(G)) ! The area integral of a KE term in a layer at u-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real :: KE_v(SZI_(G),SZJB_(G)) ! The area integral of a KE term in a layer at v-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k

  real :: uh !< Transport through zonal faces = u*h*dy,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vh !< Transport through meridional faces = v*h*dx,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: global_integral !< Global integral of the energy effect of CNN model [W]

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  
  call create_group_pass(pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  
  KE_term(:,:,:) = 0.
  tmp(:,:,:) = 0.
  ! Calculate the KE source from Zanna-Bolton2020 [H L2 T-3 ~> m3 s-3].
  do k=1,nz
    KE_u(:,:) = 0.
    KE_v(:,:) = 0.
    do j=js,je ; do I=Isq,Ieq
      uh = u(I,j,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) * &
        G%dyCu(I,j)
      KE_u(I,j) = uh * G%dxCu(I,j) * fx(I,j,k)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      vh = v(i,J,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) * &
        G%dxCv(i,J)
      KE_v(i,J) = vh * G%dyCv(i,J) * fy(i,J,k)
    enddo ; enddo
    call do_group_pass(pass_KE_uv, G%domain)
    do j=js,je ; do i=is,ie
      KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
          * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      ! copy-paste from MOM_spatial_means.F90, line 42
      tmp(i,j,k) = KE_term(i,j,k) * G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
  enddo

  global_integral = reproducing_sum(tmp)

  call post_data(CS%id_KE_CNN, KE_term, CS%diag)

end subroutine compute_energy_source
  
end module forpy_util