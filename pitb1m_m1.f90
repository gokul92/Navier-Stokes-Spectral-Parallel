!*********************************************************************
! PARALLEL NAVIER STOKES SOLVER FOR COMPRESSIBEL FLOWS 
! USING SPECTRAL METHODS
!*********************************************************************

!**********************************
! MODULE DEFINITIONS
!**********************************

module gdef

  integer, parameter :: NX=96,NY=96,NZ=96,npy=2,npz=2,ntmax=1700,intrvl=200
  integer ix,iy,iz,ax
  integer fx,lx,fy,ly,fz,lz,dx,dy,dz,nax
  double precision pi,val,val1
  double complex data(NX),trig(NX+15),ic


  double precision, allocatable :: u1(:,:,:)
  double precision, allocatable :: v1(:,:,:),w1(:,:,:)
  !double complex, allocatable :: f(:,:,:),f1(:,:,:),f2(:,:,:)
  double complex, allocatable :: f1(:,:,:),f2(:,:,:)
  double precision, allocatable :: ro1(:,:,:)
  double complex, allocatable :: u(:,:,:),v(:,:,:),ro(:,:,:)
  double complex, allocatable :: w(:,:,:),rou(:,:,:),rov(:,:,:)
  double complex, allocatable :: row(:,:,:),p(:,:,:),T(:,:,:)
  double complex, allocatable :: e(:,:,:),et(:,:,:),roet(:,:,:)
  double complex, allocatable :: rhsc(:,:,:),d(:,:,:),div(:,:,:) 
  double complex, allocatable :: rmx(:,:,:),rmy(:,:,:),rmz(:,:,:)
  double complex, allocatable :: px(:,:,:),py(:,:,:),pz(:,:,:)
  double complex, allocatable :: mvx(:,:,:),mvy(:,:,:),mvz(:,:,:)
  double complex, allocatable :: cx(:,:,:),cy(:,:,:),cz(:,:,:)
  double complex, allocatable :: txx(:,:,:),txy(:,:,:),txz(:,:,:)
  double complex, allocatable :: tyx(:,:,:),tyy(:,:,:),tyz(:,:,:)
  double complex, allocatable :: tzx(:,:,:),tzy(:,:,:),tzz(:,:,:)
  double complex, allocatable :: reng(:,:,:),ep(:,:,:),ec(:,:,:)
  double complex, allocatable :: ev(:,:,:),eq(:,:,:),qx(:,:,:)
  double complex, allocatable :: qy(:,:,:),qz(:,:,:),rod(:,:,:)
  double complex, allocatable :: kc1(:,:,:),kc2(:,:,:),roud(:,:,:)
  double complex, allocatable :: kc3(:,:,:),kc4(:,:,:),rovd(:,:,:)
  double complex, allocatable :: kmx1(:,:,:),kmx2(:,:,:),rowd(:,:,:)
  double complex, allocatable :: kmx3(:,:,:),kmx4(:,:,:),roetd(:,:,:)
  double complex, allocatable :: kmy1(:,:,:),kmy2(:,:,:),roe(:,:,:)
  double complex, allocatable :: kmy3(:,:,:),kmy4(:,:,:)
  double complex, allocatable :: kmz1(:,:,:),kmz2(:,:,:)
  double complex, allocatable :: kmz3(:,:,:),kmz4(:,:,:)
  double complex, allocatable :: ke1(:,:,:),ke2(:,:,:)
  double complex, allocatable :: ke3(:,:,:),ke4(:,:,:)
  double complex, allocatable :: dudx(:,:,:),dudy(:,:,:),dudz(:,:,:)
  double complex, allocatable :: dvdx(:,:,:),dvdy(:,:,:),dvdz(:,:,:)
  double complex, allocatable :: dwdx(:,:,:),dwdy(:,:,:),dwdz(:,:,:)
  double complex, allocatable :: droudx(:,:,:),drovdy(:,:,:),drowdz(:,:,:)
  double precision, allocatable :: Kt(:,:,:)  
  double complex, allocatable :: engi(:,:,:),eng(:,:,:),aeng(:,:,:)
  double precision, allocatable :: mu(:,:,:)
  double precision, allocatable :: roin(:,:,:)
  double precision, allocatable :: q(:,:,:)
  double complex, allocatable :: pin(:,:,:),q1(:,:,:)
  double complex, allocatable :: du(:,:,:),dv(:,:,:),dw(:,:,:)

  ! Modification
  double precision, allocatable :: vort(:,:,:), temp1(:,:,:), avgcalc(:,:,:)
  double complex, allocatable :: vortx(:,:,:), vorty(:,:,:), vortz(:,:,:)
  double complex, allocatable :: tfluc(:,:,:), rhofluc(:,:,:), pfluc(:,:,:)
  double complex, allocatable :: vort_matrix(:,:,:)

  ! Modification for tke equation
  double precision, allocatable :: tke_temp_real(:,:,:), udp(:,:,:), vdp(:,:,:), wdp(:,:,:)
  double complex, allocatable :: tke_temp_cmplx(:,:,:), tke_matrix(:,:,:), rst(:,:,:)
  double complex, allocatable :: convec_1(:,:,:), convec_2(:,:,:), convec_3(:,:,:)
  double complex, allocatable :: diffusion_1(:,:,:), diffusion_2(:,:,:), diffusion_3(:,:,:)
  double precision, allocatable :: convection(:,:,:), production(:,:,:), dissipation(:,:,:), diffusion(:,:,:), pressure_work(:,:,:), pressure_dilatation(:,:,:)

  ! Modification for enstrophy equation
  double complex, allocatable :: vortdp1(:,:,:), vortdp2(:,:,:), vortdp3(:,:,:), A11(:,:,:), A12(:,:,:), A13(:,:,:), A21(:,:,:), A22(:,:,:), A23(:,:,:), A31(:,:,:), A32(:,:,:), A33(:,:,:)
  double complex, allocatable :: crossprodresult(:,:,:), dotprodresult(:,:,:)
  double complex, allocatable :: dudpdx(:,:,:), dudpdy(:,:,:), dudpdz(:,:,:), dvdpdx(:,:,:), dvdpdy(:,:,:), dvdpdz(:,:,:), dwdpdx(:,:,:), dwdpdy(:,:,:), dwdpdz(:,:,:)
  double precision, allocatable :: realtemp(:,:,:)
  double complex, allocatable :: dpdx(:,:,:), dpdy(:,:,:), dpdz(:,:,:), drodx(:,:,:), drody(:,:,:), drodz(:,:,:)
  double complex, allocatable :: dvortdp1dx(:,:,:), dvortdp1dy(:,:,:), dvortdp1dz(:,:,:), dvortdp2dx(:,:,:), dvortdp2dy(:,:,:), dvortdp2dz(:,:,:), dvortdp3dx(:,:,:), dvortdp3dy(:,:,:), dvortdp3dz(:,:,:)
  double complex, allocatable :: uprime(:,:,:), vprime(:,:,:), wprime(:,:,:), pprime(:,:,:), roprime(:,:,:), tprime(:,:,:)
  double complex, allocatable :: complexmu(:,:,:)

  double complex, allocatable :: txxprime(:,:,:), tyyprime(:,:,:), tzzprime(:,:,:), tyxprime(:,:,:), txyprime(:,:,:), tzxprime(:,:,:), txzprime(:,:,:), tzyprime(:,:,:), tyzprime(:,:,:)

  double precision, dimension(:,:,:), allocatable :: ro1_read, p1_read, t1_read, e1_read, et_read

end module gdef

module mpidef

  include 'mpif.h'
!!$  include '/home/cluster/ghosh/FFT/include/fftw3.f'
  include 'fftw3.f'
  integer myrank
  integer p1
  integer tag1
  character*25 messg
  character*3  digit_string
  integer stat(mpi_status_size)
  integer ierr,comm2d
  integer dims(2),ndim,coords(2)
  logical isperiodic(2),reorder

end module mpidef

!*********************************
! MAIN PROGRAM
!*********************************

program main

  use gdef
  use mpidef

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myrank,ierr)
  call mpi_comm_size(mpi_comm_world,p1,ierr)

  dims(1)=npy
  dims(2)=npz
  isperiodic(1)=.true.
  isperiodic(2)=.true.
  reorder=.true.
  ndim=2

  call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,reorder,comm2d,ierr)
  call mpi_cart_get(comm2d,2,dims,isperiodic,coords,ierr)

  iy=coords(1);iz=coords(2)

  dy=NY/npy
  dz=NZ/npz

  fx=1
  lx=NX
  fy=(NY/npy)*iy+1
  ly=(NY/npy)*(iy+1)
  fz=(NZ/npz)*iz+1
  lz=(NZ/npz)*(iz+1)

  ax=(lx-fx+1)*8
  !write(*,*) 'hello'

  allocate(u1(fx:lx,fy:ly,fz:lz))
  allocate(v1(fx:lx,fy:ly,fz:lz),w1(fx:lx,fy:ly,fz:lz))
  !allocate(f(fx:lx,fy:ly,fz:lz),f1(fx:lx,fy:ly,fz:lz),f2(fx:lx,fy:ly,fz:lz))
  allocate(f1(fx:lx,fy:ly,fz:lz),f2(fx:lx,fy:ly,fz:lz))
  allocate(ro1(fx:lx,fy:ly,fz:lz))
  allocate(u(fx:lx,fy:ly,fz:lz),v(fx:lx,fy:ly,fz:lz),ro(fx:lx,fy:ly,fz:lz))
  allocate(w(fx:lx,fy:ly,fz:lz),rou(fx:lx,fy:ly,fz:lz),rov(fx:lx,fy:ly,fz:lz))
  allocate(row(fx:lx,fy:ly,fz:lz),p(fx:lx,fy:ly,fz:lz),T(fx:lx,fy:ly,fz:lz))
  allocate(e(fx:lx,fy:ly,fz:lz),et(fx:lx,fy:ly,fz:lz),roet(fx:lx,fy:ly,fz:lz))
  allocate(rhsc(fx:lx,fy:ly,fz:lz),d(fx:lx,fy:ly,fz:lz),div(fx:lx,fy:ly,fz:lz)) 
  allocate(rmx(fx:lx,fy:ly,fz:lz),rmy(fx:lx,fy:ly,fz:lz),rmz(fx:lx,fy:ly,fz:lz))
  allocate(px(fx:lx,fy:ly,fz:lz),py(fx:lx,fy:ly,fz:lz),pz(fx:lx,fy:ly,fz:lz))
  allocate(mvx(fx:lx,fy:ly,fz:lz),mvy(fx:lx,fy:ly,fz:lz),mvz(fx:lx,fy:ly,fz:lz))
  allocate(cx(fx:lx,fy:ly,fz:lz),cy(fx:lx,fy:ly,fz:lz),cz(fx:lx,fy:ly,fz:lz))
  allocate(txx(fx:lx,fy:ly,fz:lz),txy(fx:lx,fy:ly,fz:lz),txz(fx:lx,fy:ly,fz:lz))
  allocate(tyx(fx:lx,fy:ly,fz:lz),tyy(fx:lx,fy:ly,fz:lz),tyz(fx:lx,fy:ly,fz:lz))
  allocate(tzx(fx:lx,fy:ly,fz:lz),tzy(fx:lx,fy:ly,fz:lz),tzz(fx:lx,fy:ly,fz:lz))
  allocate(reng(fx:lx,fy:ly,fz:lz),ep(fx:lx,fy:ly,fz:lz),ec(fx:lx,fy:ly,fz:lz))
  allocate(ev(fx:lx,fy:ly,fz:lz),eq(fx:lx,fy:ly,fz:lz),qx(fx:lx,fy:ly,fz:lz))
  allocate(qy(fx:lx,fy:ly,fz:lz),qz(fx:lx,fy:ly,fz:lz),rod(fx:lx,fy:ly,fz:lz))
  allocate(kc1(fx:lx,fy:ly,fz:lz),kc2(fx:lx,fy:ly,fz:lz),roud(fx:lx,fy:ly,fz:lz))
  allocate(kc3(fx:lx,fy:ly,fz:lz),kc4(fx:lx,fy:ly,fz:lz),rovd(fx:lx,fy:ly,fz:lz))
  allocate(kmx1(fx:lx,fy:ly,fz:lz),kmx2(fx:lx,fy:ly,fz:lz),rowd(fx:lx,fy:ly,fz:lz))
  allocate(kmx3(fx:lx,fy:ly,fz:lz),kmx4(fx:lx,fy:ly,fz:lz),roetd(fx:lx,fy:ly,fz:lz))
  allocate(kmy1(fx:lx,fy:ly,fz:lz),kmy2(fx:lx,fy:ly,fz:lz),roe(fx:lx,fy:ly,fz:lz))
  allocate(kmy3(fx:lx,fy:ly,fz:lz),kmy4(fx:lx,fy:ly,fz:lz))
  allocate(kmz1(fx:lx,fy:ly,fz:lz),kmz2(fx:lx,fy:ly,fz:lz))
  allocate(kmz3(fx:lx,fy:ly,fz:lz),kmz4(fx:lx,fy:ly,fz:lz))
  allocate(ke1(fx:lx,fy:ly,fz:lz),ke2(fx:lx,fy:ly,fz:lz))
  allocate(ke3(fx:lx,fy:ly,fz:lz),ke4(fx:lx,fy:ly,fz:lz))
  allocate(dudx(fx:lx,fy:ly,fz:lz),dudy(fx:lx,fy:ly,fz:lz),dudz(fx:lx,fy:ly,fz:lz))
  allocate(dvdx(fx:lx,fy:ly,fz:lz),dvdy(fx:lx,fy:ly,fz:lz),dvdz(fx:lx,fy:ly,fz:lz))
  allocate(dwdx(fx:lx,fy:ly,fz:lz),dwdy(fx:lx,fy:ly,fz:lz),dwdz(fx:lx,fy:ly,fz:lz))
  allocate(droudx(fx:lx,fy:ly,fz:lz),drovdy(fx:lx,fy:ly,fz:lz),drowdz(fx:lx,fy:ly,fz:lz))
  allocate(Kt(fx:lx,fy:ly,fz:lz))
  allocate(engi(fx:lx,fy:ly,fz:lz),eng(fx:lx,fy:ly,fz:lz),aeng(fx:lx,fy:ly,fz:lz))
  allocate(mu(fx:lx,fy:ly,fz:lz))
  allocate(roin(fx:lx,fy:ly,fz:lz))
  allocate(q(fx:lx,fy:ly,fz:lz))
  allocate(pin(fx:lx,fy:ly,fz:lz),q1(fx:lx,fy:ly,fz:lz))
  allocate(du(fx:lx,fy:ly,fz:lz),dv(fx:lx,fy:ly,fz:lz),dw(fx:lx,fy:ly,fz:lz))

  ! Modification
  allocate(vort(fx:lx,fy:ly,fz:lz), temp1(fx:lx,fy:ly,fz:lz))
  allocate(vortx(fx:lx,fy:ly,fz:lz), vorty(fx:lx,fy:ly,fz:lz), vortz(fx:lx,fy:ly,fz:lz))
  allocate(tfluc(fx:lx,fy:ly,fz:lz), rhofluc(fx:lx,fy:ly,fz:lz), pfluc(fx:lx,fy:ly,fz:lz))
  allocate(vort_matrix(fx:lx,fy:ly,fz:lz))
  allocate(avgcalc(fx:lx,fy:ly,fz:lz))

  ! Modification for tke equation
  allocate(tke_temp_real(fx:lx,fy:ly,fz:lz), udp(fx:lx,fy:ly,fz:lz), vdp(fx:lx,fy:ly,fz:lz), wdp(fx:lx,fy:ly,fz:lz))
  allocate(tke_temp_cmplx(fx:lx,fy:ly,fz:lz), tke_matrix(fx:lx,fy:ly,fz:lz), rst(fx:lx,fy:ly,fz:lz))
  allocate(convec_1(fx:lx,fy:ly,fz:lz), convec_2(fx:lx,fy:ly,fz:lz), convec_3(fx:lx,fy:ly,fz:lz))
  allocate(diffusion_1(fx:lx,fy:ly,fz:lz), diffusion_2(fx:lx,fy:ly,fz:lz), diffusion_3(fx:lx,fy:ly,fz:lz))
  allocate(convection(fx:lx,fy:ly,fz:lz), production(fx:lx,fy:ly,fz:lz), dissipation(fx:lx,fy:ly,fz:lz), diffusion(fx:lx,fy:ly,fz:lz), pressure_work(fx:lx,fy:ly,fz:lz), pressure_dilatation(fx:lx,fy:ly,fz:lz))

  ! Modification for enstrophy equation
  allocate(vortdp1(fx:lx,fy:ly,fz:lz), vortdp2(fx:lx,fy:ly,fz:lz), vortdp3(fx:lx,fy:ly,fz:lz), A11(fx:lx,fy:ly,fz:lz), A12(fx:lx,fy:ly,fz:lz), A13(fx:lx,fy:ly,fz:lz), A21(fx:lx,fy:ly,fz:lz), A22(fx:lx,fy:ly,fz:lz), A23(fx:lx,fy:ly,fz:lz), A31(fx:lx,fy:ly,fz:lz), A32(fx:lx,fy:ly,fz:lz), A33(fx:lx,fy:ly,fz:lz))
  allocate(crossprodresult(fx:lx,fy:ly,fz:lz), dotprodresult(fx:lx,fy:ly,fz:lz))
  allocate(dudpdx(fx:lx,fy:ly,fz:lz), dudpdy(fx:lx,fy:ly,fz:lz), dudpdz(fx:lx,fy:ly,fz:lz), dvdpdx(fx:lx,fy:ly,fz:lz), dvdpdy(fx:lx,fy:ly,fz:lz), dvdpdz(fx:lx,fy:ly,fz:lz), dwdpdx(fx:lx,fy:ly,fz:lz), dwdpdy(fx:lx,fy:ly,fz:lz), dwdpdz(fx:lx,fy:ly,fz:lz))
  allocate(realtemp(fx:lx,fy:ly,fz:lz))
  allocate(dpdx(fx:lx,fy:ly,fz:lz), dpdy(fx:lx,fy:ly,fz:lz), dpdz(fx:lx,fy:ly,fz:lz), drodx(fx:lx,fy:ly,fz:lz), drody(fx:lx,fy:ly,fz:lz), drodz(fx:lx,fy:ly,fz:lz))
  allocate(dvortdp1dx(fx:lx,fy:ly,fz:lz), dvortdp1dy(fx:lx,fy:ly,fz:lz), dvortdp1dz(fx:lx,fy:ly,fz:lz), dvortdp2dx(fx:lx,fy:ly,fz:lz), dvortdp2dy(fx:lx,fy:ly,fz:lz), dvortdp2dz(fx:lx,fy:ly,fz:lz), dvortdp3dx(fx:lx,fy:ly,fz:lz), dvortdp3dy(fx:lx,fy:ly,fz:lz), dvortdp3dz(fx:lx,fy:ly,fz:lz))
  allocate(uprime(fx:lx,fy:ly,fz:lz),vprime(fx:lx,fy:ly,fz:lz),wprime(fx:lx,fy:ly,fz:lz),roprime(fx:lx,fy:ly,fz:lz),pprime(fx:lx,fy:ly,fz:lz),tprime(fx:lx,fy:ly,fz:lz))
  allocate(complexmu(fx:lx,fy:ly,fz:lz))

  allocate(txxprime(fx:lx,fy:ly,fz:lz), tyyprime(fx:lx,fy:ly,fz:lz), tzzprime(fx:lx,fy:ly,fz:lz), tyxprime(fx:lx,fy:ly,fz:lz), txyprime(fx:lx,fy:ly,fz:lz), tzxprime(fx:lx,fy:ly,fz:lz), txzprime(fx:lx,fy:ly,fz:lz), tzyprime(fx:lx,fy:ly,fz:lz), tyzprime(fx:lx,fy:ly,fz:lz))

  allocate(ro1_read(fx:lx,fy:ly,fz:lz), p1_read(fx:lx,fy:ly,fz:lz), t1_read(fx:lx,fy:ly,fz:lz), e1_read(fx:lx,fy:ly,fz:lz), et_read(fx:lx,fy:ly,fz:lz))

  call isoturb()

  call mpi_finalize(ierr)

end program main

!************************************
! SUBROUTINE TO SOLVE NS EQUATIONS
!************************************

subroutine isoturb()

  use gdef
  use mpidef
  implicit none

  integer numrc3,numrc4                           !modification
  character*40 :: filename                        !modification
  double precision x_un,y_un,z_un                 !modification   
  double precision, parameter :: no=0.670d0
  double precision, parameter :: gamma=1.40d0,PR=0.70d0,Mno=0.30d0
  integer my,mz
  integer i,j,k,tl,vl,isign,flag,co,nt
  integer count,datasize,nx1,ny1,nz1,fh
  integer (kind=mpi_offset_kind) offset
  integer status(mpi_status_size)
  double precision buf(ax)
  double precision x(fx:lx),kx(fy:ly),y(fy:ly),ky(fz:lz),z(fz:lz),kz(fx:lx)
  double precision e1, t1(2),x1,x2,b,pnot,sc,dt
  double precision Ei,Eval,aEval,tval
  double precision mubar,robar,mecher(ntmax),mecher1(ntmax)
  double precision K0,U0,C0,Urms,Relam(ntmax),lam1(ntmax),lam2(ntmax),lam3(ntmax),lam(ntmax),qr,qsnot
  double precision s(ntmax),s1(ntmax),s2(ntmax),s3(ntmax),duc,dus,dvc,dvs,dwc,dws,duc1,dus1,dvc1,dvs1,dwc1,dws1
  double precision us,us1,vs,vs1,ws,ws1
  double complex uin,vin,win,Tin,etin
  double precision ur(ntmax),vr(ntmax),wr(ntmax),prm(ntmax),Tr(ntmax),ror(ntmax),er(ntmax),etr(ntmax),qval,qval1
  double precision initial_divergence, initial_compressible_divergence
  double precision dilatational_dissipation, solenoidal_dissipation, velocity_scale, dummy_var, velocity_scaling
  double precision compressibility(ntmax), returb(ntmax), divergence(ntmax), vorticity(ntmax), initial_vorticity
  double precision :: lam1_init, lam2_init, lam3_init, lam_init, Relam_init, dummytemp1, dummytemp2, dummytemp3, dummytemp4
  double precision :: etotime, total_dissipation

  double complex, dimension(:), allocatable :: in,out
  integer*8 plan,plan1
  integer digit1,digit2, statusint

  !* quantities for dimensionalizing
  double precision, parameter :: t0star = 300.0d0, Rair = 287.150d0, c0star = dsqrt(gamma*Rair*t0star), Cp = Rair*gamma/(gamma-1.0d0), Cv = Rair/(gamma-1.0d0)
  double precision, parameter :: bconstant = 4.34161664d0/(10**(7)), mu0star = bconstant*(t0star)**0.67d0, ro0star = 2.320000092188518E-002, p0star = ro0star*c0star**2
  !double precision, parameter :: l0star = Re*mu0star/(ro0star*c0star)

  double precision, parameter :: t0star_base = 300.0d0, c0star_base = dsqrt(gamma*Rair*t0star_base)
  double precision, parameter :: bconstant_base = 4.34161664d0/(10**(7)), mu0star_base = bconstant_base*(t0star_base)**0.67d0, ro0star_base = 1.20d0, p0star_base = ro0star_base*c0star_base**2
  !double precision, parameter :: l0star_base = Re*mu0star_base/(ro0star_base*c0star_base)

  integer restart, nt_digit1, nt_digit2, nt_digit3, nt_digit4, nt_digit5, nt_digit6, nt_digit7, nt_digit8
  ! calculating various means
  double precision :: tbar(ntmax), pbar(ntmax)
  double precision :: robar0, mubar0
  !modification for tke equation
  double precision :: robar_prev, tke_favreavg_prev, tke_favreavg0, tke_lhs, tke_rhs
  double precision :: ubar, vbar, wbar, favreavg
  double precision :: udpbar, vdpbar, wdpbar, rst_temp, tke_favreavg
  !*Computational Reynolds number is set to obtain appropriate taylor microscale reynolds number.
  double precision :: Re
  double precision, parameter :: relam_m1 = 70.0d0

  !  Modification for enstrophy equation
  double precision :: enstrophy_rhs, enstrophy_rhs1, enstrophy_rhs2, enstrophy_rhs3, enstrophy_rhs4, enstrophy_rhs5, enstrophy_rhs6, enstrophy_rhs7, enstrophy_rhs8, enstrophy_rhs9
  double precision :: mufavg

  ! Modification for internal energy equation
  double precision :: e_favreavg, e_favreavg_prev, e_favreavg0, elhs, erhs
  double precision :: ufavg, vfavg, wfavg

  ! To calculate two point correlations
  double precision, dimension(NZ) :: uu_corr
  double precision, dimension(NZ/2) :: q_array
  double precision :: tempsum
  integer :: counter, num

  ! initialise variables and constants

  statusint = 0

  dt=0.0050d0

  ic=(0.0d0,1.0d0)
  pi=2*asin(1.0d0)

  C0=1.0d0
  K0=12.0d0

  ! initialise fft

  allocate(in(NX))
  allocate(out(NX))

  call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  call dfftw_plan_dft_1d(plan1,NX,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

  ! initialise grid

  dy=NY/npy
  dz=NZ/npz

  do i=fx,lx
     x(i)=(i-1)*2.0d0*pi/NX
     if(i<=NX/2)then
        kz(i)=i-1
     else
        kz(i)=i-NX-1
     end if
  end do

  do j=fy,ly
     y(j)=(j-1)*2.0d0*pi/NY
     if(j<=NY/2)then
        kx(j)=j-1
     else
        kx(j)=j-NY-1
     end if
  end do

  do k=fz,lz
     z(k)=(k-1)*2.0d0*pi/NZ
     if(k<=NZ/2)then
        ky(k)=k-1
     else
        ky(k)=k-NZ-1
     end if
  end do

  !/////////////////////////////////////////////////////////////////////////////
  numrc3 = myrank/10
  numrc4 = mod(myrank,10)
  filename='p'//char(numrc3+48)//char(numrc4+48)//'.dat'
  !write(6,*) myrank,filename
  open(unit=50,file=filename,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
  read(50,*)
  read(50,*)
  do k=fz,lz
     do j=fy,ly
        do i=fx,lx
           read(50,'(6E30.18)') u1(i,j,k),v1(i,j,k),w1(i,j,k), dummytemp1, dummytemp2, dummytemp3
           pfluc(i,j,k) = dummytemp1
           rhofluc(i,j,k) = dummytemp2
           tfluc(i,j,k) = dummytemp3
        end do
     end do
  end do
  close(50)

  do k=fz,lz
     do j=fy,ly
        u(fx:lx,j,k)=u1(fx:lx,j,k)
     end do
  end do

  do k=fz,lz
     do j=fy,ly
        v(fx:lx,j,k)=v1(fx:lx,j,k)
     end do
  end do

  do k=fz,lz
     do j=fy,ly
        w(fx:lx,j,k)=w1(fx:lx,j,k)
     end do
  end do

  !///////////////////////////////////////////////////////////////////////////

  ! Rescaling fluctuations
  u = u*c0star/c0star_base
  v = v*c0star/c0star_base
  w = w*c0star/c0star_base

  pfluc = pfluc*(ro0star*c0star**2)/(ro0star_base*c0star_base**2)
  rhofluc = rhofluc*ro0star/ro0star_base
  tfluc = tfluc*t0star/t0star_base

  du=u
  call deriv(du,1)

  dv=v
  call deriv(dv,2)

  dw=w
  call deriv(dw,3)

  temp1 = real(du) + real(dv) + real(dw)
  temp1 = temp1**2
  call getavg(temp1)

  initial_divergence = val1

  if ( myrank .eq. 0 ) then
     write(*,*) "initial dilatation", initial_divergence
  end if

  !call exit(statusint)

  qsnot=0.0d0
  call getrms(u)
  qsnot=val**2
  call getrms(v)
  qsnot=qsnot+val**2
  call getrms(w)
  qsnot=qsnot+val**2

  if ( myrank .eq. 0 ) then
     write(*,*) "turbulent mach number before velocity scaling", dsqrt(qsnot)
  end if

  Tin=t0star/t0star_base + ic*0.0d0
  roin=ro0star/ro0star_base
  pin=roin*Tin/gamma

  T=Tin+tfluc
  ro=roin+rhofluc
  p=pin+pfluc

  e=T/(gamma*(gamma-1))
  et=e+0.50d0*(u**2+v**2+w**2)

  !modification for internal energy equation
  call favreavgfunc(e, ro, favreavg)
  e_favreavg0 =  favreavg

  mu=T**no
  call getavg(mu)
  mubar=val1
  mubar0=mubar

  us=0.0d0
  us1 = 0.0d0
  vs=0.0d0
  vs1 = 0.0d0
  ws=0.0d0
  ws1 = 0.0d0
  dus = 0.0d0
  dvs = 0.0d0
  dws = 0.0d0
  qval = 0.0d0
  qval1 = 0.0d0

  do k = fz, lz
     do j = fy, ly
        do i = fx, lx

           us=us+(real(u(i,j,k)))**2
           dus=dus+(real(du(i,j,k)))**2

           vs=vs+(real(v(i,j,k)))**2
           dvs=dvs+(real(dv(i,j,k)))**2

           ws=ws+(real(w(i,j,k)))**2
           dws=dws+(real(dw(i,j,k)))**2

           qval=qval+(real(u(i,j,k)))**2+(real(v(i,j,k)))**2+(real(w(i,j,k)))**2

        end do
     end do
  end do

  us=us/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
  dus=dus/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

  vs=vs/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
  dvs=dvs/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

  ws=ws/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
  dws=dws/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

  qval=qval/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

  call mpi_allreduce(us,us1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
  call mpi_allreduce(dus,dus1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

  call mpi_allreduce(vs,vs1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
  call mpi_allreduce(dvs,dvs1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

  call mpi_allreduce(ws,ws1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
  call mpi_allreduce(dws,dws1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

  call mpi_allreduce(qval,qval1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

  us1=us1/npy/npz
  dus1=dus1/npy/npz

  vs1=vs1/npy/npz
  dvs1=dvs1/npy/npz

  ws1=ws1/npy/npz
  dws1=dws1/npy/npz

  qval1=qval1/npy/npz

  qval1=dsqrt(qval1)

  lam1_init=dsqrt(us1/dus1)
  lam2_init=dsqrt(vs1/dvs1)
  lam3_init=dsqrt(ws1/dws1)
  lam_init=(lam1_init+lam2_init+lam3_init)/3.0d0

  call getavg(roin)
  robar=val1
  robar0=robar

  ! Initial Taylor Reynolds Number
  Relam_init=relam_m1
  Re = mubar*relam_m1/(lam_init*dsqrt(qsnot/3.0d0)*robar)
  if ( myrank .eq. 0 ) then
     write(*,*) "Relam", Relam_init
  end if

  !********************************************************************************!

  ! Calculating Dissipation

  ! Calculating Dilatational Dissipation

  du = u
  call deriv(du,1)
  dv = v
  call deriv(dv,2)
  dw = w
  call deriv(dw,3)

  temp1 = real(du) + real(dv) + real(dw)
  temp1 = temp1**2
  call getavg(temp1)
  dilatational_dissipation = val1*mubar*(4.0d0/(Re*3.0d0))

  ! Calculating Solenoidal Dissipation

  ! Calculating Vorticity
  vortx = 0.0d0 + ic*0.0d0
  vorty = 0.0d0 + ic*0.0d0
  vortz = 0.0d0 + ic*0.0d0

  d = w
  call deriv(d, 2)
  vortx = vortx+d
  d = v
  call deriv(d, 3)
  vortx = vortx-d

  d = u
  call deriv(d, 3)
  vorty = vorty+d
  d = w
  call deriv(d, 1)
  vorty = vorty-d

  d = v
  call deriv(d, 1)
  vortz = vortz+d
  d = u
  call deriv(d, 2)
  vortz = vortz-d

  ! vort is real matrix
  vort = real( (vortx)**2 + (vorty)**2 + (vortz)**2 )
  call getavg(vort)
  initial_vorticity = val1
  solenoidal_dissipation = mubar*val1/Re

  ! compressibility is a vector of size ntmax
  total_dissipation = (solenoidal_dissipation + dilatational_dissipation)

  etotime = robar*qsnot/(3.0d0*total_dissipation)

  !********************************************************************************!

  Kt=(mu)/(Re*Pr*(gamma-1.0d0))

  du=u
  call deriv(du,1)

  dv=v
  call deriv(dv,2)

  dw=w
  call deriv(dw,3)

  temp1 = real(du) + real(dv) + real(dw)
  temp1 = temp1**2
  call getavg(temp1)

  initial_divergence = val1

  if (myrank .eq. 0) then
     write(*,*) "initial dilatation", initial_divergence
     write(*,*) "total dissipation", total_dissipation
     write(*,*) "Computational Reynolds number", Re
     write(*,*) "Relam initial", Relam_init
     write(*,*) "lam", lam_init
     write(*,*) "qval1", qval1
     write(*,*) "robar", robar
     write(*,*) "mubar", mubar
     write(*,*) "etotime", etotime
  end if

  !call exit(statusint)

  !*****calculating initial favre averaged turbulent kinetic energy
  call favreavgfunc(u, ro, favreavg)
  udp = u - favreavg

  call favreavgfunc(v, ro, favreavg)
  vdp = v - favreavg

  call favreavgfunc(w, ro, favreavg)
  wdp = w - favreavg

  tke_matrix = 0.50d0*(udp**2 + vdp**2 + wdp**2)
  call favreavgfunc(tke_matrix, ro, favreavg)
  tke_favreavg0 = favreavg

  !call exit(statusint)
  ! assign values for restart
  restart = 0

  if ( restart .eq. 1 ) then
     write(*,*) "restart activated"

     nt = int(1.20d0/dt + 1.0d0)

     numrc3 = myrank/10
     numrc4 = mod(myrank,10)
     !filename='valp_'//char(numrc3+48)//char(numrc4+48)//'.dat'

     nt_digit1 = (nt-1)/10
     nt_digit2 = mod((nt-1),10)

     nt_digit3 = nt_digit1/10
     nt_digit4 = mod(nt_digit1,10)

     nt_digit5 = nt_digit3/10
     nt_digit6 = mod(nt_digit3,10)

     nt_digit7 = nt_digit5/10
     nt_digit8 = mod(nt_digit5,10)

     digit1=myrank/10
     digit2=mod(myrank,10)
     filename='vars_nt_'//char(nt_digit8+48)//char(nt_digit6+48)//char(nt_digit4+48) &
          //char(nt_digit2+48)//'_p_'//char(digit1+48)//char(digit2+48)//'.dat'

     open(unit=50,file=filename,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
     !read(50,*)
     !read(50,*)
     !do i=fx,lx
     !   do j=fy,ly
     !      do k=fz,lz
     do k=fz,lz
        do j=fy,ly
           do i=fx,lx
              !read(50,'(1x,12E15.6)') dummytemp1, dummytemp2, dummytemp3, u1(i,j,k), v1(i,j,k), w1(i,j,k), ro1_read(i,j,k), p1_read(i,j,k), t1_read(i,j,k), e1_read(i,j,k), et_read(i,j,k), dummytemp4
              read(50,'(6E30.18)') u1(i,j,k), v1(i,j,k),w1(i,j,k),p1_read(i,j,k),ro1_read(i,j,k),t1_read(i,j,k)
              u(i,j,k) = u1(i,j,k)
              v(i,j,k) = v1(i,j,k)
              w(i,j,k) = w1(i,j,k)
           end do
        end do
     end do
     close(50)
     ro = ro1_read
     p = p1_read
     T = t1_read
     e=T/(gamma*(gamma-1.0d0))
     et=e+0.50d0*(u**2 + v**2 + w**2)

     mu = T**no
     Kt=(mu)/(Re*Pr*(gamma-1.0d0))

     write(*,*) "nt value for restart is", nt
  end if

  rou=ro*u
  rov=ro*v
  row=ro*w
  roet=ro*et

  nt = 1

  if ( restart .eq. 1 ) then
     open(unit=10, file='stat3.dat',status='old',action='write',position='append',iostat=ierr)
  else if ( restart .eq. 0 ) then
     open(file='stat3.dat',unit=10)
     open(file='two_pt_corr.dat',unit=68)
     write(68,*)'variables="r", "Q"'
     write(75,*)'variables="time","convection","production","dissipation","diffusion","pressure_work","pressure_dilatation","elhs","erhs","tke_lhs","tke_rhs"'
     write(85,*)'variables="time","ens1","ens2","ens3","ens4","ens5","ens6","ens7","ens8","ens9","total_ens"'
     !write(85,*)'variables="time","convection","production","dissipation","diffusion","pressure_work","pressure_dilatation"'
  end if
  do while ( nt .le. ntmax )

     !**********Modification for tke verification************!
     if ( nt .eq. 1) then
        robar_prev = robar0
        tke_favreavg_prev = tke_favreavg0

        e_favreavg_prev = e_favreavg0
     else
        robar_prev = robar
        tke_favreavg_prev = tke_favreavg

        e_favreavg_prev = e_favreavg
     end if
     !**********Modification for tke verification************!

     call rungakutta()

     s(nt)=0.0d0
     s1(nt)=0.0d0
     s2(nt)=0.0d0
     s3(nt)=0.0d0

     du=0.0d0+ic*0.0d0
     dv=0.0d0+ic*0.0d0
     dw=0.0d0+ic*0.0d0

     us=0.0d0
     vs=0.0d0
     ws=0.0d0
     duc=0.0d0
     dvc=0.0d0
     dwc=0.0d0
     dus=0.0d0
     dvs=0.0d0
     dws=0.0d0
     qval=0.0d0
     qval1=0.0d0

     du=u
     call deriv(du,1)

     dv=v
     call deriv(dv,2)

     dw=w
     call deriv(dw,3)

     do i=fx,lx
        do j=fy,ly
           do k=fz,lz

              us=us+(real(u(i,j,k)))**2
              duc=duc+(real(du(i,j,k)))**3
              dus=dus+(real(du(i,j,k)))**2

              vs=vs+(real(v(i,j,k)))**2
              dvc=dvc+(real(dv(i,j,k)))**3
              dvs=dvs+(real(dv(i,j,k)))**2

              ws=ws+(real(w(i,j,k)))**2
              dwc=dwc+(real(dw(i,j,k)))**3
              dws=dws+(real(dw(i,j,k)))**2

              qval=qval+(real(u(i,j,k)))**2+(real(v(i,j,k)))**2+(real(w(i,j,k)))**2

           end do
        end do
     end do

     us=us/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     duc=duc/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     dus=dus/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

     vs=vs/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     dvc=dvc/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     dvs=dvs/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

     ws=ws/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     dwc=dwc/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
     dws=dws/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

     qval=qval/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)

     call mpi_allreduce(us,us1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(duc,duc1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(dus,dus1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

     call mpi_allreduce(vs,vs1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(dvc,dvc1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(dvs,dvs1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

     call mpi_allreduce(ws,ws1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(dwc,dwc1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
     call mpi_allreduce(dws,dws1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

     call mpi_allreduce(qval,qval1,1,mpi_double_precision,mpi_sum,comm2d,ierr)

     us1=us1/npy/npz
     duc1=duc1/npy/npz
     dus1=dus1/npy/npz

     vs1=vs1/npy/npz
     dvc1=dvc1/npy/npz
     dvs1=dvs1/npy/npz

     ws1=ws1/npy/npz
     dwc1=dwc1/npy/npz
     dws1=dws1/npy/npz

     qval1=qval1/npy/npz

     s1(nt)=duc1/((dus1**1.50d0))
     s2(nt)=dvc1/((dvs1**1.50d0))
     s3(nt)=dwc1/((dws1**1.50d0))
     s(nt)=(s1(nt)+s2(nt)+s3(nt))/3.0d0

     lam1(nt)=dsqrt(us1/dus1)
     lam2(nt)=dsqrt(vs1/dvs1)
     lam3(nt)=dsqrt(ws1/dws1)
     lam(nt)=(lam1(nt)+lam2(nt)+lam3(nt))/3.0d0

     qval1=dsqrt(qval1)

     call getavg(mu)
     mubar=val1

     ro1=real(ro)
     call getavg(ro1)
     robar=val1

     Relam(nt)=Re*lam(nt)*(qval1/dsqrt(3.0d0))*(robar/mubar)

     avgcalc = real(T)
     call getavg(avgcalc)
     tbar(nt)=val1

     avgcalc = real(p)
     call getavg(avgcalc)
     pbar(nt)=val1

     tke_temp_real = real(u)
     call getavg(tke_temp_real)
     ubar = val1

     tke_temp_real = real(v)
     call getavg(tke_temp_real)
     vbar = val1

     tke_temp_real = real(w)
     call getavg(tke_temp_real)
     wbar = val1

     call getrms(u)
     ur(nt)=val
     call getrms(v)
     vr(nt)=val
     call getrms(w)
     wr(nt)=val
     call getrms(p)
     prm(nt)=val
     call getrms(T)
     Tr(nt)=val
     call getrms(ro)
     ror(nt)=val
     call getrms(e)
     er(nt)=val
     call getrms(et)
     etr(nt)=val
     q1=real(u)**2+real(v)**2+real(w)**2
     q=real(q1)
     call getavg(q)
     qr=val1
     mecher(nt)=(ur(nt)**2+vr(nt)**2+wr(nt)**2)/qsnot
     mecher1(nt)=qr/qsnot

     call favreavgfunc(u, ro, favreavg)
     ufavg = favreavg
     udp = u - ufavg
     realtemp = real(udp)
     call getavg(realtemp)
     udpbar = val1

     call favreavgfunc(v, ro, favreavg)
     vfavg = favreavg
     vdp = v - vfavg
     realtemp = real(vdp)
     call getavg(realtemp)
     vdpbar = val1

     call favreavgfunc(w, ro, favreavg)
     wfavg = favreavg
     wdp = w - wfavg
     realtemp = real(wdp)
     call getavg(realtemp)
     wdpbar = val1

     uprime = u - ubar
     vprime = v - vbar
     wprime = w - wbar

     ! Reassigning du, dv and dw to derivatives of uprime, vprime, wprime and not u, v and w.
     d = uprime
     call deriv(d,1)
     du = d

     d = vprime
     call deriv(d,2)
     dv = d

     d = wprime
     call deriv(d,3)
     dw = d

     complexmu = real(mu)
     ! Calculating Dissipation
     call favreavgfunc(complexmu, ro, favreavg)
     mufavg = favreavg

     ! Calculating Dilatational Dissipation
     temp1 = real(du) + real(dv) + real(dw)
     temp1 = temp1**2
     call getavg(temp1)
     dilatational_dissipation = val1*mufavg*(4.0d0/(Re*3.0d0))

     ! Calculating divergence squared averaged
     divergence(nt) = val1

     ! Calculating Solenoidal Dissipation

     ! Calculating Vorticity
     vortx = 0.0d0 + ic*0.0d0
     vorty = 0.0d0 + ic*0.0d0
     vortz = 0.0d0 + ic*0.0d0

     d = wprime
     call deriv(d, 2)
     vortx = vortx+d
     d = vprime
     call deriv(d, 3)
     vortx = vortx-d

     d = uprime
     call deriv(d, 3)
     vorty = vorty+d
     d = wprime
     call deriv(d, 1)
     vorty = vorty-d

     d = vprime
     call deriv(d, 1)
     vortz = vortz+d
     d = uprime
     call deriv(d, 2)
     vortz = vortz-d

     ! vort is real matrix
     vort =  real(vortx)**2 + real(vorty)**2 + real(vortz)**2
     call getavg(vort)
     solenoidal_dissipation = mufavg*val1/Re

     ! Calculating vorticity squared averaged
     vorticity(nt) = val1

     ! Calculating total dissipation
     total_dissipation = solenoidal_dissipation + dilatational_dissipation

     ! compressibility is a vector of size ntmax
     compressibility(nt) = dilatational_dissipation/(solenoidal_dissipation + dilatational_dissipation)

     ! writing velocity field to file at regular intervals to calculate energy spectrum
     if ( mod(nt,intrvl) .eq. 0 ) then

        nt_digit1 = nt/10
        nt_digit2 = mod(nt,10)

        nt_digit3 = nt_digit1/10
        nt_digit4 = mod(nt_digit1,10)

        nt_digit5 = nt_digit3/10
        nt_digit6 = mod(nt_digit3,10)

        nt_digit7 = nt_digit5/10
        nt_digit8 = mod(nt_digit5,10)

        digit1=myrank/10
        digit2=mod(myrank,10)
        filename='vars_nt_'//char(nt_digit8+48)//char(nt_digit6+48)//char(nt_digit4+48) &
             //char(nt_digit2+48)//'_p_'//char(digit1+48)//char(digit2+48)//'.dat'

        open(unit = 67+myrank, file=filename, action='write')

        do k = fz, lz
           do j = fy, ly
              do i= fx, lx
                 write(67+myrank, '(6E30.18)') real(u(i,j,k)), real(v(i,j,k)), real(w(i,j,k)), real(p(i,j,k)), real(ro(i,j,k)), real(T(i,j,k))
              end do
           end do
        end do

        !write(67+myrank)real(u),real(v),real(w),real(p),real(ro),real(T),real(div),dsqrt(vort)

        close(67+myrank)

     end if

     !Calculating turbulent kinetic energy equation rhs
     ! subroutines to calculate the individual RHS terms of the turbulent kinetic energy equations
     call favrevelocitygradients()
     call favrevelocityvorticity()
     !call Amatrix()
     call favrevorticitygradients()
     call densitypressuregradients()
     !call touprime()

     call convectioncalc()
     call productioncalc()
     call dissipationcalc()
     call diffusioncalc()
     call pressureworkcalc()
     call pressuredilatationcalc()
     call totaltkeequationcalc()

     ! Computing enstrophy equation rhs
     ! subroutines to calculate the individual RHS terms of the enstrophy equations
     enstrophy_rhs1 = 0.0d0
     enstrophy_rhs2 = 0.0d0
     enstrophy_rhs3 = 0.0d0
     enstrophy_rhs4 = 0.0d0
     enstrophy_rhs5 = 0.0d0
     enstrophy_rhs6 = 0.0d0
     enstrophy_rhs7 = 0.0d0
     enstrophy_rhs8 = 0.0d0
     enstrophy_rhs9 = 0.0d0

     call enstrophyrhscalc1()
     call enstrophyrhscalc2()
     call enstrophyrhscalc3()
     call enstrophyrhscalc4()
     call enstrophyrhscalc5()
     call enstrophyrhscalc6()
     call enstrophyrhscalc7()
     call enstrophyrhscalc8()
     call enstrophyrhscalc9()
     call totalenstrophyrhscalc()

     ! Computing lhs of internal energy equation for verification

     !call favreavgfunc(e, ro, favreavg)
     !e_favreavg = favreavg
     !elhs = ( (robar*e_favreavg) - (robar_prev*e_favreavg_prev) )/dt
     ! To write enstrophy equations for rhs8 and rhs9

     tval=dt*nt

     if ( myrank .eq. 0) then
        write(*,*) nt, mecher(nt), Relam(nt)
     end if
     write(10,'(1x,31E14.6)')tval,tval/etotime,s1(nt),s2(nt),s3(nt),s(nt),ur(nt)**2,vr(nt)**2,wr(nt)**2,prm(nt)**2,Tr(nt)**2,ror(nt)**2,mecher(nt),lam1(nt),lam2(nt),lam3(nt),lam(nt),Relam(nt),compressibility(nt),divergence(nt),vorticity(nt)/initial_vorticity,pbar(nt),robar,tbar(nt),mubar,vorticity(nt),total_dissipation,enstrophy_rhs,ubar,vbar,wbar

     if (myrank .eq. 0 ) then
        write(75,'(11E30.18)') tval, maxval(convection), maxval(production), maxval(dissipation), maxval(diffusion), maxval(pressure_work), maxval(pressure_dilatation), elhs, erhs, tke_lhs, tke_rhs
        write(85,'(11E30.18)') tval, enstrophy_rhs1, enstrophy_rhs2, enstrophy_rhs3, enstrophy_rhs4, enstrophy_rhs5, enstrophy_rhs6, enstrophy_rhs7, enstrophy_rhs8, enstrophy_rhs9, enstrophy_rhs
     end if

     nt = nt + 1
  end do
  close(10)
  !close(45)
  !close(55)
  close(75)
  close(85)
  close(85)
  close(68)

  digit1=myrank/10
  digit2=mod(myrank,10)
  filename='valp_'//char(digit1+48)//char(digit2+48)//'.dat'
  open(file=filename,unit=myrank+1)           
  div=0.0d0+ic*0.0d0

  d=0.0d0+ic*0.0d0
  d=u
  call deriv(d,1)
  div=d
  d=0.0d0+ic*0.0d0
  d=v
  call deriv(d,2)
  div=div+d
  d=0.0d0+ic*0.0d0
  d=w
  call deriv(d,3)
  div=div+d
  d=0.0d0+ic*0.0d0

  d = w
  call deriv(d, 2)
  vortx = vortx+d
  d = v
  call deriv(d, 3)
  vortx = vortx-d

  d = u
  call deriv(d, 3)
  vorty = vorty+d
  d = w
  call deriv(d, 1)
  vorty = vorty-d

  d = v
  call deriv(d, 1)
  vortz = vortz+d
  d = u
  call deriv(d, 2)
  vortz = vortz-d

  ! vort is real matrix
  vort_matrix = (vortx)**2 + (vorty)**2 + (vortz)**2 

  write(myrank+1,*)  'variables= "x","y","z","u","v","w","ro","p","T","e","et","div"'
  write(myrank+1,*)  'zone I=',lz-fz+1,' J=',ly-fy+1,' K=',lx-fx+1,' F=POINT'
  do i=fx,lx
     do j=fy,ly
        do k=fz,lz
           write(myrank+1,'(1x,13E15.6)')x(i),y(j),z(k),real(u(i,j,k)),real(v(i,j,k)),real(w(i,j,k)),real(ro(i,j,k)),real(p(i,j,k)),real(T(i,j,k)),real(e(i,j,k)),real(et(i,j,k)),real(div(i,j,k)),real(vort_matrix(i,j,k))
        end do
     end do
  end do
  close(10)

  deallocate(in)
  deallocate(out)
  call dfftw_destroy_plan(plan)
  call dfftw_destroy_plan(plan1)

contains

  subroutine touprime()

    use mpidef
    use gdef
    implicit none

    txxprime = 2.0d0*dudpdx-dvdpdy-dwdpdz
    tyyprime = 2.0d0*dvdpdy-dwdpdz-dudpdx
    tzzprime = 2.0d0*dwdpdz-dudpdx-dvdpdy

    txyprime = dvdpdx+dudpdy
    tyzprime = dwdpdy+dvdpdz
    tzxprime = dudpdz+dwdpdx

    tyxprime = txyprime
    tzyprime = tyzprime
    txzprime = tzxprime

  end subroutine touprime

  subroutine convectioncalc()

    implicit none

    !Calculating convection term

    tke_matrix = 0.50d0*(udp**2 + vdp**2 + wdp**2)
    call favreavgfunc(tke_matrix, ro, favreavg)
    tke_favreavg = favreavg
    convec_1 = -1.0d0*robar*tke_favreavg*ubar
    convec_2 = -1.0d0*robar*tke_favreavg*vbar
    convec_3 = -1.0d0*robar*tke_favreavg*wbar

    tke_temp_real = 0.0d0

    d = convec_1
    call deriv(d,1)
    tke_temp_real = tke_temp_real + real(d)

    d = convec_2
    call deriv(d,2)
    tke_temp_real = tke_temp_real + real(d)

    d = convec_3
    call deriv(d,3)
    tke_temp_real = tke_temp_real + real(d)

    convection = real(tke_temp_real)

  end subroutine convectioncalc

  subroutine productioncalc()

    implicit none

    !Calculating production
    !Calculating reynolds stress tensor (complex matrix)
    rst = 0.0d0

    tke_temp_real = ro*udp*udp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,u,favreavg)
    d = favreavg
    call deriv(d,1)
    rst = rst + d*rst_temp

    tke_temp_real = ro*udp*vdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,u,favreavg)
    d = favreavg
    call deriv(d,2)
    rst = rst + d*rst_temp

    tke_temp_real = ro*udp*wdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,u,favreavg)
    d = favreavg
    call deriv(d,3)
    rst = rst + d*rst_temp

    tke_temp_real = ro*vdp*udp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,v,favreavg)
    d = favreavg
    call deriv(d,1)
    rst = rst + d*rst_temp

    tke_temp_real = ro*vdp*vdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,v,favreavg)
    d = favreavg
    call deriv(d,2)
    rst = rst + d*rst_temp

    tke_temp_real = ro*vdp*wdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,v,favreavg)
    d = favreavg
    call deriv(d,3)
    rst = rst + d*rst_temp

    tke_temp_real = ro*wdp*udp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,w,favreavg)
    d = favreavg
    call deriv(d,1)
    rst = rst + d*rst_temp

    tke_temp_real = ro*wdp*vdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,w,favreavg)
    d = favreavg
    call deriv(d,2)
    rst = rst + d*rst_temp

    tke_temp_real = ro*wdp*wdp
    call getavg(tke_temp_real)
    rst_temp = val1
    call favreavgfunc(ro,w,favreavg)
    d = favreavg
    call deriv(d,3)
    rst = rst + d*rst_temp

    production = -1.0d0*real(rst)

  end subroutine productioncalc

  subroutine dissipationcalc()

    implicit none

    !Calculating dissipation
    tke_temp_real = 0.0d0

    d = udp
    call deriv(d,1)
    tke_temp_real = tke_temp_real + real(txx)*real(d)
    d = udp
    call deriv(d,2)
    tke_temp_real = tke_temp_real + real(txy)*real(d)
    d = udp
    call deriv(d,3)
    tke_temp_real = tke_temp_real + real(txz)*real(d)

    d = vdp
    call deriv(d,1)
    tke_temp_real = tke_temp_real + real(tyx)*real(d)
    d = vdp
    call deriv(d,2)
    tke_temp_real = tke_temp_real + real(tyy)*real(d)
    d = vdp
    call deriv(d,3)
    tke_temp_real = tke_temp_real + real(tyz)*real(d)

    d = wdp
    call deriv(d,1)
    tke_temp_real = tke_temp_real + real(tzx)*real(d)
    d = wdp
    call deriv(d,2)
    tke_temp_real = tke_temp_real + real(tzy)*real(d)
    d = wdp
    call deriv(d,3)
    tke_temp_real = tke_temp_real + real(tzz)*real(d)

    call getavg(tke_temp_real)
    dissipation = -1.0d0*val1

  end subroutine dissipationcalc

  subroutine diffusioncalc()

    implicit none

    !Calculating diffusion
    diffusion = 0.0d0
    diffusion_1 = 0.0d0
    diffusion_2 = 0.0d0
    diffusion_3 = 0.0d0

    tke_temp_real = udp*real(txx) + vdp*real(tyx) + wdp*real(tzx)
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,1)
    diffusion_1 = diffusion_1 + d

    tke_temp_real = udp*real(txy) + vdp*real(tyy) + wdp*real(tzy)
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,2)
    diffusion_1 = diffusion_1 + d

    tke_temp_real = udp*real(txz) + vdp*real(tyz) + wdp*real(tzz)
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,3)
    diffusion_1 = diffusion_1 + d

    tke_temp_real = tke_matrix*ro*udp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,1)
    diffusion_2 = diffusion_2 - d

    tke_temp_real = tke_matrix*ro*vdp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,2)
    diffusion_2 = diffusion_2 - d

    tke_temp_real = tke_matrix*ro*wdp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,3)
    diffusion_2 = diffusion_2 - d

    tke_temp_real = (p-pbar(nt))*udp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,1)
    diffusion_3 = diffusion_3 - d

    tke_temp_real = (p-pbar(nt))*vdp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,2)
    diffusion_3 = diffusion_3 - d

    tke_temp_real = (p-pbar(nt))*wdp
    call getavg(tke_temp_real)
    d = val1
    call deriv(d,3)
    diffusion_3 = diffusion_3 - d

    diffusion = real(diffusion_1) + real(diffusion_2) + real(diffusion_3)

  end subroutine diffusioncalc

  subroutine pressureworkcalc()

    implicit none

    !Calculating pressure work
    pressure_work = 0.0d0

    d = pbar(nt)
    call deriv(d,1)
    pressure_work = pressure_work - real(d)*udpbar

    d = pbar(nt)
    call deriv(d,2)
    pressure_work = pressure_work - real(d)*vdpbar

    d = pbar(nt)
    call deriv(d,3)
    pressure_work = pressure_work - real(d)*wdpbar

  end subroutine pressureworkcalc

  subroutine pressuredilatationcalc()

    implicit none

    !Calculating pressure dilatation
    tke_temp_real = 0.0d0

    d = udp
    call deriv(d,1)
    tke_temp_real = tke_temp_real + real(d)

    d = vdp
    call deriv(d,2)
    tke_temp_real = tke_temp_real + real(d)

    d = wdp
    call deriv(d,3)
    tke_temp_real = tke_temp_real + real(d)

    tke_temp_real = tke_temp_real*(real(p)-pbar(nt))
    call getavg(tke_temp_real)
    pressure_dilatation = val1

  end subroutine pressuredilatationcalc

  subroutine totaltkeequationcalc()

    implicit none

    tke_lhs = ( (robar*tke_favreavg) - (robar_prev*tke_favreavg_prev) )/dt
    tke_rhs = maxval(convection) + maxval(production) + maxval(dissipation) + maxval(diffusion) + maxval(pressure_work) + maxval(pressure_dilatation)

  end subroutine totaltkeequationcalc

  subroutine enstrophyrhscalc1()

    implicit none

    realtemp = real(vortdp1)**2*(2.0d0*real(dudpdx)) + real(vortdp1)*real(vortdp2)*(real(dudpdy) + real(dvdpdx)) + real(vortdp1)*real(vortdp3)*(real(dudpdz) + real(dwdpdx)) + real(vortdp2)*real(vortdp1)*(real(dvdpdx) + real(dudpdy)) + real(vortdp2)**2*(2.0d0*real(dvdpdy)) + real(vortdp2)*real(vortdp3)*(real(dvdpdz) + real(dwdpdy)) + real(vortdp3)*real(vortdp1)*(real(dwdpdx) + real(dudpdz)) + real(vortdp3)*real(vortdp2)*(real(dwdpdy) + real(dvdpdz)) + real(vortdp3)**2*(2.0d0*real(dwdpdz))

    call getavg(realtemp)
    enstrophy_rhs1 = val1

  end subroutine enstrophyrhscalc1

  subroutine enstrophyrhscalc2()

    implicit none

    realtemp = ( real(vortdp1)**2 + real(vortdp2)**2 + real(vortdp3)**2 ) * ( real(dudpdx) + real(dvdpdy) + real(dwdpdz) )

    call getavg(realtemp)
    enstrophy_rhs2 = -1.0d0*val1

  end subroutine enstrophyrhscalc2

  subroutine enstrophyrhscalc3()

    implicit none

    realtemp = 0.0d0
    d  = mu*( 2.0d0*dudpdx - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( dudpdy + dvdpdx )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( dudpdz + dwdpdx )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = realtemp*real(dvortdp3dy)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

    realtemp = 0.0d0
    d  = mu*( dvdpdx + dudpdy )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( 2.0d0*dudpdy - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( dvdpdz + dwdpdy )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = realtemp*real(dvortdp1dz)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

    realtemp = 0.0d0
    d  = mu*( dwdpdx + dudpdz )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( dwdpdy + dvdpdz )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( 2.0d0*dwdpdz - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = realtemp*real(dvortdp2dx)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

    realtemp = 0.0d0
    d  = mu*( 2.0d0*dudpdx - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( dudpdy + dvdpdx )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( dudpdz + dwdpdx )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = -1.0d0*realtemp*real(dvortdp2dz)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

    realtemp = 0.0d0
    d  = mu*( dvdpdx + dudpdy )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( 2.0d0*dudpdy - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( dvdpdz + dwdpdy )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = -1.0d0*realtemp*real(dvortdp3dx)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

    realtemp = 0.0d0
    d  = mu*( dwdpdx + dudpdz )
    call deriv(d,1)
    realtemp = realtemp + real(d)
    d = mu*( dwdpdy + dvdpdz )
    call deriv(d,2)
    realtemp = realtemp + real(d)
    d = mu*( 2.0d0*dwdpdz - (2.0d0/3.0d0)*(dudpdx + dvdpdy + dwdpdz) )
    call deriv(d,3)
    realtemp = realtemp + real(d)
    realtemp = -1.0d0*realtemp*real(dvortdp1dy)*(2.0d0/Re)
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/real(ro(i,j,k))
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs3 = enstrophy_rhs3 + val1

  end subroutine enstrophyrhscalc3

  subroutine enstrophyrhscalc4()

    implicit none

    realtemp = 2.0d0*( real(vortdp1)*(real(drody)*real(dpdz) - real(drodz)*real(dpdy)) + real(vortdp2)*(real(drodz)*real(dpdx) - real(drodx)*real(dpdz)) + real(vortdp3)*(real(drodx)*real(dpdy) - real(drody)*real(dpdx)) )
    do k = fz, lz
       do j = fy, ly
          do i = fx, lx
             realtemp(i,j,k) = realtemp(i,j,k)/(real(ro(i,j,k))**2)
          end do
       end do
    end do
    call getavg(realtemp)
    enstrophy_rhs4 = val1

  end subroutine enstrophyrhscalc4

  subroutine enstrophyrhscalc5()

    implicit none

    realtemp = 2.0d0*( real(vortdp1)**2 * real(A11) + real(vortdp2)**2 * real(A22) + real(vortdp3)**2 * real(A33) + real(vortdp1)*real(vortdp2)*(real(A12) + real(A21)) + real(vortdp1)*real(vortdp3)*(real(A13) + real(A31)) + real(vortdp2)*real(vortdp3)*(real(A23) + real(A32)) )
    call getavg(realtemp)
    enstrophy_rhs5 = val1

  end subroutine enstrophyrhscalc5

  subroutine enstrophyrhscalc6()

    implicit none

    realtemp = -2.0d0*( real(vortdp1)**2 + real(vortdp2)**2 + real(vortdp3)**2  )*( real(A11) + real(A22) + real(A33) )
    call getavg(realtemp)
    enstrophy_rhs6 = val1

  end subroutine enstrophyrhscalc6

  subroutine enstrophyrhscalc7()

    implicit none

    realtemp = (real(A32) - real(A23))*( 2.0d0*real(vortdp1)*real(dudpdx) + real(vortdp2)*(real(dudpdy) + real(dvdpdx)) + real(vortdp3)*(real(dudpdz) + real(dwdpdx)) ) + (real(A13) - real(A31))*( 2.0d0*real(vortdp2)*real(dvdpdy) + real(vortdp1)*(real(dvdpdx) + real(dudpdy)) + real(vortdp3)*(real(dvdpdz) + real(dwdpdy)) ) + (real(A21) - real(A12))*( 2.0d0*real(vortdp3)*real(dwdpdz) + real(vortdp1)*(real(dwdpdx) + real(dudpdz)) + real(vortdp2)*(real(dwdpdy) + real(dvdpdz)) )
    call getavg(realtemp)
    enstrophy_rhs7 = val1

  end subroutine enstrophyrhscalc7

  subroutine enstrophyrhscalc8()

    implicit none

  end subroutine enstrophyrhscalc8

  subroutine enstrophyrhscalc9()

    implicit none

  end subroutine enstrophyrhscalc9

  subroutine totalenstrophyrhscalc()

    implicit none

    enstrophy_rhs = enstrophy_rhs1 + enstrophy_rhs2 + enstrophy_rhs3 + enstrophy_rhs4 + enstrophy_rhs5 + enstrophy_rhs6 + enstrophy_rhs7 + enstrophy_rhs8 + enstrophy_rhs9

  end subroutine totalenstrophyrhscalc

  subroutine favrevelocitygradients()

    implicit none

    d = udp
    call deriv(d,1)
    dudpdx = d

    d = udp
    call deriv(d,2)
    dudpdy = d

    d = udp
    call deriv(d,3)
    dudpdz = d

    d = vdp
    call deriv(d,1)
    dvdpdx = d

    d = vdp
    call deriv(d,2)
    dvdpdy = d

    d = vdp
    call deriv(d,3)
    dvdpdz = d

    d = wdp
    call deriv(d,1)
    dwdpdx = d

    d = wdp
    call deriv(d,2)
    dwdpdy = d

    d = wdp
    call deriv(d,3)
    dwdpdz = d

  end subroutine favrevelocitygradients

  subroutine favrevelocityvorticity()

    implicit none

    vortdp1 = 0.0d0
    vortdp2 = 0.0d0
    vortdp3 = 0.0d0

    d = wdp
    call deriv(d, 2)
    vortdp1 = vortdp1+d
    d = vdp
    call deriv(d, 3)
    vortdp1 = vortdp1-d

    d = udp
    call deriv(d, 3)
    vortdp2 = vortdp2+d
    d = wdp
    call deriv(d, 1)
    vortdp2 = vortdp2-d

    d = vdp
    call deriv(d, 1)
    vortdp3 = vortdp3+d
    d = udp
    call deriv(d, 2)
    vortdp3 = vortdp3-d

    ! vort is real matrix
    vort = real(vortx)**2 + real(vorty)**2 + real(vortz)**2

  end subroutine favrevelocityvorticity

  subroutine Amatrix()

    implicit none

    call favreavgfunc(u, ro, favreavg)
    d = favreavg
    call deriv(d,1)
    A11 = d

    call favreavgfunc(u, ro, favreavg)
    d = favreavg
    call deriv(d,2)
    A12 = d

    call favreavgfunc(u, ro, favreavg)
    d = favreavg
    call deriv(d,3)
    A13 = d

    call favreavgfunc(v, ro, favreavg)
    d = favreavg
    call deriv(d,1)
    A21 = d

    call favreavgfunc(v, ro, favreavg)
    d = favreavg
    call deriv(d,2)
    A22 = d

    call favreavgfunc(v, ro, favreavg)
    d = favreavg
    call deriv(d,3)
    A23 = d

    call favreavgfunc(w, ro, favreavg)
    d = favreavg
    call deriv(d,1)
    A31 = d

    call favreavgfunc(w, ro, favreavg)
    d = favreavg
    call deriv(d,2)
    A32 = d

    call favreavgfunc(w, ro, favreavg)
    d = favreavg
    call deriv(d,3)
    A33 = d

  end subroutine Amatrix

  subroutine favrevorticitygradients()

    d = vortdp1
    call deriv(d,1)
    dvortdp1dx = d

    d = vortdp1
    call deriv(d,2)
    dvortdp1dy = d

    d = vortdp1
    call deriv(d,3)
    dvortdp1dz = d

    d = vortdp2
    call deriv(d,1)
    dvortdp2dx = d

    d = vortdp2
    call deriv(d,2)
    dvortdp2dy = d

    d = vortdp2
    call deriv(d,3)
    dvortdp2dz = d

    d = vortdp3
    call deriv(d,1)
    dvortdp3dx = d

    d = vortdp3
    call deriv(d,2)
    dvortdp3dy = d

    d = vortdp3
    call deriv(d,3)
    dvortdp3dz = d

  end subroutine favrevorticitygradients

  subroutine densitypressuregradients()

    d = ro
    call deriv(d,1)
    drodx = d

    d = ro
    call deriv(d,2)
    drody = d

    d = ro
    call deriv(d,3)
    drodz = d

    d = p
    call deriv(d,1)
    dpdx = d

    d = p
    call deriv(d,2)
    dpdy = d

    d = p
    call deriv(d,3)
    dpdz = d

  end subroutine densitypressuregradients

  subroutine favreavgfunc(f, density, favreavg)

    use mpidef
    use gdef
    implicit none

    double precision, dimension(fx:lx,fy:ly,fz:lz) :: temp_matrix
    double complex, dimension(fx:lx,fy:ly,fz:lz) :: density, f
    double precision favreavg

    temp_matrix = real(density)*real(f)
    call getavg(temp_matrix)
    favreavg = val1

    temp_matrix = real(density)
    call getavg(temp_matrix)
    favreavg = favreavg/val1

  end subroutine favreavgfunc

  !******************************************
  ! COMPUTES AVERAGE
  !******************************************


  subroutine getavg(f)

    use mpidef
    use gdef
    implicit none

    double precision f(fx:lx,fy:ly,fz:lz),fbar,val2
    integer i,j,k

    fbar=0.0d0;val1=0.0d0;val2=0.0d0
    do k=fz,lz
       do j=fy,ly
          do i=fx,lx
             fbar=fbar+f(i,j,k)
          end do
       end do
    end do

    val2=fbar/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
    call mpi_allreduce(val2,val1,1,mpi_double_precision,mpi_sum,comm2d,ierr)
    val1=val1/npy/npz

  end subroutine getavg

  !******************************************
  ! COMPUTES ROOT MEAN SQUARE
  !******************************************

  subroutine getrms(f)

    use mpidef
    use gdef
    implicit none

    double complex f(fx:lx,fy:ly,fz:lz)
    integer i,j,k
    double precision fbar,fsbar,val3,val4,val5,val6,fp,fps

    fbar=0.0d0;val=0.0d0;val3=0.0d0;val4=0.0d0;val5=0.0d0;val6=0.0d0;fp=0.0d0;fps=0.0d0

    do i=fx,lx
       do j=fy,ly
          do k=fz,lz
             fbar=fbar+real(f(i,j,k))
          end do
       end do
    end do

    fbar=fbar/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
    val3=fbar

    call mpi_allreduce(val3,val4,1,mpi_double_precision,mpi_sum,comm2d,ierr) 
    val4=val4/npy/npz

    do i=fx,lx
       do j=fy,ly
          do k=fz,lz
             fp=f(i,j,k)-val4
             fps=fps+fp**2
          end do
       end do
    end do

    fps=fps/(lx-fx+1)/(ly-fy+1)/(lz-fz+1)
    val5=fps

    call mpi_allreduce(val5,val6,1,mpi_double_precision,mpi_sum,comm2d,ierr)
    val6=val6/npy/npz

    val=dsqrt(val6)

  end subroutine getrms

  subroutine commons()

    use mpidef
    use gdef
    implicit none

    droudx = 0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=rou
    call deriv(d,1)
    droudx=droudx + d
    d=u
    call deriv(d,1)
    droudx=droudx + ro*d
    d=ro
    call deriv(d,1)
    droudx=droudx + u*d
    droudx=0.50d0*droudx

    drovdy = 0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=rov
    call deriv(d,2)
    drovdy=drovdy + d
    d=v
    call deriv(d,2)
    drovdy=drovdy + ro*d
    d=ro
    call deriv(d,2)
    drovdy=drovdy + v*d
    drovdy=0.50d0*drovdy

    drowdz=0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=row
    call deriv(d,3)
    drowdz=drowdz + d
    d=w
    call deriv(d,3)
    drowdz=drowdz + ro*d
    d=ro
    call deriv(d,3)
    drowdz=drowdz + w*d
    drowdz=0.50d0*drowdz

    d=0.0d0+ic*0.0d0
    d=u
    call deriv(d,1)
    dudx=d

    d=0.0d0+ic*0.0d0
    d=u
    call deriv(d,2)
    dudy=d

    d=0.0d0+ic*0.0d0
    d=u
    call deriv(d,3)
    dudz=d

    d=0.0d0+ic*0.0d0
    d=v
    call deriv(d,1)
    dvdx=d

    d=0.0d0+ic*0.0d0
    d=v
    call deriv(d,2)
    dvdy=d

    d=0.0d0+ic*0.0d0
    d=v
    call deriv(d,3)
    dvdz=d

    d=0.0d0+ic*0.0d0
    d=w
    call deriv(d,1)
    dwdx=d

    d=0.0d0+ic*0.0d0
    d=w
    call deriv(d,2)
    dwdy=d

    d=0.0d0+ic*0.0d0
    d=w
    call deriv(d,3)
    dwdz=d

    d=0.0d0+ic*0.0d0

  end subroutine commons

  !***********************************************
  !   COMPUTES THE RHS OF THE CONTINUITY EQUATION
  !***********************************************

  subroutine rhsct()

    use mpidef
    use gdef
    implicit none

    rhsc = -1.0d0*(droudx+drovdy+drowdz)

  end subroutine rhsct


  !******************************
  !  MOMENTUM EQUATION
  !******************************************************
  ! Computes the RHS of the momentum equation
  !******************************************************

  subroutine rhsmom()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    call mconv()
    call pgrad()
    call tou()
    call mvisc()

    rmx = -px-cx+mvx
    rmy = -py-cy+mvy
    rmz = -pz-cz+mvz    

  end subroutine rhsmom

  !******************************************************
  ! Computes nine components of the Stress Tensor (Tou)
  !******************************************************

  subroutine tou()

    use mpidef
    use gdef
    implicit none

    txx = 2.*dudx-dvdy-dwdz
    tyy = 2.*dvdy-dwdz-dudx
    tzz = 2.*dwdz-dudx-dvdy

    txy = dvdx+dudy
    tyz = dwdy+dvdz
    tzx = dudz+dwdx

    tyx = txy
    tzy = tyz
    txz = tzx

  end subroutine tou

  !**************************************************************
  ! Computes gradients of tous
  !**************************************************************
  subroutine mvisc()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    mvx=0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=txx
    call deriv(d,1)
    mvx=mvx+0.50d0*(2.0d0/3.0d0)*(mu/Re)*d
    d=(2.0d0/3.0d0)*(mu/Re)+ic*0.0d0
    call deriv(d,1)
    mvx=mvx+0.50d0*txx*d
    txx=(2.0d0/3.0d0)*(mu/Re)*txx
    d=txx
    call deriv(d,1)
    mvx = mvx + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=txy
    call deriv(d,2)
    mvx=mvx+0.50d0*(mu/Re)*d
    d=(mu/Re)+ic*0.0d0
    call deriv(d,2)
    mvx=mvx+0.50d0*txy*d
    txy=(mu/Re)*txy
    d=txy
    call deriv(d,2)
    mvx = mvx + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=txz
    call deriv(d,3)
    mvx=mvx+0.50d0*(mu/Re)*d
    d=(mu/Re)+ic*0.0d0
    call deriv(d,3)
    mvx=mvx+0.50d0*txz*d
    txz=(mu/Re)*txz
    d=txz
    call deriv(d,3)
    mvx = mvx + 0.50d0*d

    mvy = 0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=tyx
    call deriv(d,1)
    mvy = mvy + 0.50d0*(mu/Re)*d
    d=(mu/Re)+ic*0.0d0
    call deriv(d,1)
    mvy = mvy + 0.50d0*tyx*d
    tyx=(mu/Re)*tyx
    d=tyx
    call deriv(d,1)
    mvy = mvy + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=tyy
    call deriv(d,2) 
    mvy = mvy + 0.50d0*(2.0d0/3.0d0)*(mu/Re)*d
    d=(2.0d0/3.0d0)*(mu/Re)+ic*0.0d0
    call deriv(d,2)
    mvy = mvy + 0.50d0*tyy*d
    tyy=(2.0d0/3.0d0)*(mu/Re)*tyy
    d=tyy
    call deriv(d,2)
    mvy = mvy + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=tyz
    call deriv(d,3)
    mvy = mvy + 0.50d0*(mu/Re)*d
    d=(mu/Re)+ic*0.0d0
    call deriv(d,3)
    mvy = mvy + 0.50d0*tyz*d
    tyz=(mu/Re)*tyz
    d=tyz
    call deriv(d,3)
    mvy = mvy + 0.50d0*d

    mvz=0.0d0+ic*0.0d0
    d=0.0d0+ic*0.0d0
    d=tzx
    call deriv(d,1)
    mvz = mvz + 0.50d0*(mu/Re)*d
    d = (mu/Re)+ic*0.0d0
    call deriv(d,1)
    mvz = mvz + 0.50d0*tzx*d
    tzx=(mu/Re)*tzx
    d=tzx
    call deriv(d,1)
    mvz = mvz + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=tzy
    call deriv(d,2)
    mvz = mvz + 0.50d0*(mu/Re)*d
    d=(mu/Re)+ic*0.0d0
    call deriv(d,2)
    mvz = mvz + 0.50d0*tzy*d
    tzy=(mu/Re)*tzy
    d=tzy
    call deriv(d,2)
    mvz = mvz + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d = tzz
    call deriv(d,3)
    mvz = mvz + 0.50d0*(2.0d0/3.0d0)*(mu/Re)*d
    d = (2.0d0/3.0d0)*(mu/Re)+ic*0.0d0
    call deriv(d,3)
    mvz = mvz + 0.50d0*tzz*d
    tzz=(2.0d0/3.0d0)*(mu/Re)*tzz
    d=tzz
    call deriv(d,3)
    mvz = mvz + 0.50d0*d

    d=0.0d0+ic*0.0d0

  end subroutine mvisc

  !******************************************************
  ! Computes Pressure Gradients
  !******************************************************

  subroutine pgrad()

    use mpidef
    use gdef
    implicit none

    px=p
    call deriv(px,1)

    py=p
    call deriv(py,2)

    pz=p
    call deriv(pz,3)

  end subroutine pgrad

  !***********************************************************
  ! Computes the Convection terms for the Momentum Equation
  !***********************************************************

  subroutine mconv()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    cx=0.0d0+ic*0.0d0;cy=0.0d0+ic*0.0d0;cz=0.0d0+ic*0.0d0

    call con1(rou,u,1,1)

    d=0.0d0+ic*0.0d0
    d=u*droudx
    cx=cx+0.50d0*d

    d=0.0d0+ic*0.0d0
    d=rou*dudx
    cx=cx+0.50d0*d

    call con1(rou,v,2,1)
    call con2(v,rou,2,1)

    d=0.0d0+ic*0.0d0
    d=rou*dvdy
    cx=cx+0.50d0*d

    call con1(rou,w,3,1)
    call con2(w,rou,3,1)

    d=0.0d0+ic*0.0d0
    d=rou*dwdz
    cx=cx+0.50d0*d

    call con1(rov,u,1,2)
    call con2(u,rov,1,2)

    d=0.0d0+ic*0.0d0
    d=rov*dudx
    cy=cy+0.50d0*d

    call con1(rov,v,2,2)

    d=0.0d0+ic*0.0d0
    d=v*drovdy
    cy=cy+0.50d0*d

    d=0.0d0+ic*0.0d0
    d=rov*dvdy
    cy=cy+0.50d0*d

    call con1(rov,w,3,2)
    call con2(w,rov,3,2)

    d=0.0d0+ic*0.0d0
    d=rov*dwdz
    cy=cy+0.50d0*d

    call con1(row,u,1,3)
    call con2(u,row,1,3)

    d=0.0d0+ic*0.0d0
    d=row*dudx
    cz=cz+0.50d0*d

    call con1(row,v,2,3)
    call con2(v,row,2,3)

    d=0.0d0+ic*0.0d0
    d=row*dvdy
    cz=cz+0.50d0*d

    call con1(row,w,3,3)

    d=0.0d0+ic*0.0d0
    d=w*drowdz
    cz=cz+0.50d0*d

    d=0.0d0+ic*0.0d0
    d=row*dwdz
    cz=cz+0.50d0*d

  end subroutine mconv

  !************************************************************
  ! Conservative Convection Terms(Sub-part of conv())
  !************************************************************

  subroutine con1(a,b,fl,fl1)

    use mpidef
    use gdef
    implicit none

    double complex a(fx:lx,fy:ly,fz:lz),b(fx:lx,fy:ly,fz:lz)
    double complex c(fx:lx,fy:ly,fz:lz)
    integer i,j,k,fl,fl1

    d=0.0d0+ic*0.0d0;c=0.0d0+ic*0.0d0
    d=a*b
    call deriv(d,fl)
    call incconv(fl1)

  end subroutine con1

  !************************************************************
  ! Non-conservative Convection Terms(Sub-part of conv())
  !************************************************************

  subroutine con2(a,b,fl,fl1)

    use mpidef
    use gdef
    implicit none

    integer fl,fl1
    double complex a(fx:lx,fy:ly,fz:lz),b(fx:lx,fy:ly,fz:lz)
    double complex c(fx:lx,fy:ly,fz:lz)

    d=0.0d0+ic*0.0d0;c=0.0d0+ic*0.0d0

    d=b
    call deriv(d,fl)
    c=a*d
    d=0.0d0+ic*0.0d0
    d=c
    call incconv(fl1)

  end subroutine con2

  !*************************************************************
  ! Increments Convection Terms(Sub-Part of Conv)
  !*************************************************************
  subroutine incconv(fl1)

    use mpidef
    use gdef
    implicit none	

    integer fl1

    if(fl1 == 1)then
       cx = cx + 0.50d0*d
    end if

    if(fl1 == 2)then
       cy = cy + 0.50d0*d
    end if

    if(fl1 == 3)then
       cz = cz + 0.50d0*d
    end if

    if(fl1 == 4)then
       ec = ec + 0.50d0*d
    end if

    if(fl1 == 5)then
       ev = ev + d
    end if

  end subroutine incconv


  !*********************************
  !        ENERGY EQUATION
  !******************************************************
  ! Computes RHS of the Energy Equation
  !******************************************************

  subroutine rhseng()

    use mpidef
    use gdef
    implicit none

    call epgrad()
    call qflux()
    call divq()
    call econv()
    call evisc()

    reng = -ep-ec+ev+eq

  end subroutine rhseng

  !******************************************************
  ! Computes Pressure gradient terms in the energy eqn
  !******************************************************		  

  subroutine epgrad()

    use mpidef
    use gdef
    implicit none

    ep = 0.0d0+ic*0.0d0

    d=0.0d0+ic*0.0d0
    d=p
    call deriv(d,1)
    ep = ep + 0.50d0*u*d
    d=u
    call deriv(d,1)
    ep = ep + 0.50d0*p*d
    d=p*u
    call deriv(d,1)
    ep = ep + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=p
    call deriv(d,2)
    ep = ep + 0.50d0*v*d
    d=v
    call deriv(d,2)
    ep = ep + 0.50d0*p*d
    d=p*v
    call deriv(d,2)
    ep = ep + 0.50d0*d

    d=0.0d0+ic*0.0d0
    d=p
    call deriv(d,3)
    ep = ep + 0.50d0*w*d
    d=w
    call deriv(d,3)
    ep = ep + 0.50d0*p*d
    d=p*w
    call deriv(d,3)
    ep = ep + 0.50d0*d

    d=0.0d0+ic*0.0d0

  end subroutine epgrad

  !******************************************************
  ! Computes the Heat fluxes
  !******************************************************

  subroutine qflux()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    qx=T
    call deriv(qx,1)
    qx=-qx

    qy=T
    call deriv(qy,2)
    qy=-qy

    qz=T
    call deriv(qz,3)
    qz=-qz


  end subroutine qflux

  !******************************************************
  ! Computes Divergence of Heat Fluxes
  !******************************************************

  subroutine divq()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    eq=0.0d0+ic*0.0d0

    d=qx
    call deriv(d,1)
    eq = eq - 0.50d0*Kt*d
    d=Kt
    call deriv(d,1)
    eq = eq - 0.50d0*qx*d
    qx=Kt*qx 
    call deriv(qx,1)
    eq=eq-0.50d0*qx

    d=qy
    call deriv(d,2)
    eq = eq - 0.50d0*Kt*d
    d=Kt
    call deriv(d,2)
    eq = eq - 0.50d0*qy*d
    qy=Kt*qy
    call deriv(qy,2)
    eq=eq-0.50d0*qy

    d=qz
    call deriv(d,3)
    eq = eq - 0.50d0*Kt*d
    d=Kt
    call deriv(d,3)
    eq = eq - 0.50d0*qz*d
    qz=Kt*qz
    call deriv(qz,3)
    eq=eq-0.50d0*qz

    d = 0.0d0 + ic*0.0d0

  end subroutine divq

  !******************************************************
  ! Computes convection of energy terms(Skew Symm Form)
  !******************************************************

  subroutine econv()

    use mpidef
    use gdef
    implicit none

    integer i,j,k

    ec=0.0d0+ic*0.0d0

    call con1(roet,u,1,4)
    call con2(u,roet,1,4)

    d=0.0d0+ic*0.0d0
    d=roet*dudx
    ec=ec+0.50d0*d

    call con1(roet,v,2,4)
    call con2(v,roet,2,4)

    d=0.0d0+ic*0.0d0
    d=roet*dvdy
    ec=ec+0.50d0*d

    call con1(roet,w,3,4)
    call con2(w,roet,3,4)

    d=0.0d0+ic*0.0d0
    d=roet*dwdz
    ec=ec+0.50d0*d

  end subroutine econv

  !******************************************************
  ! Computes Viscous Terms in the Energy Equation (SSF)
  !******************************************************

  subroutine evisc()

    use mpidef
    use gdef
    implicit none

    ev=0.0d0+ic*0.0d0

    call con1(txx,u,1,5) 
    call con1(txy,v,1,5) 
    call con1(txz,w,1,5)

    call con1(tyx,u,2,5)
    call con1(tyy,v,2,5)
    call con1(tyz,w,2,5)

    call con1(tzx,u,3,5)
    call con1(tzy,v,3,5)
    call con1(tzz,w,3,5)

  end subroutine evisc


  !*****************************
  ! FOURTH ORDER RUNGA KUTTA
  !*****************************************************
  ! Integrates the RHS of 5 equations Simultaneously
  !*****************************************************

  subroutine rungakutta()

    use mpidef
    use gdef
    implicit none

    integer i,j,k
    real*8 dt

    dt=0.0050d0

    call commons()
    call RHS()

    ! Fourth Order Runga Kutta
    ! Step 1

    co=0
    co=co+1
    !write(*,*)nt,co
    kc1=rhsc*dt
    rod=ro
    ro=ro+0.50d0*kc1

    kmx1=rmx*dt
    kmy1=rmy*dt
    kmz1=rmz*dt

    roud=rou
    rovd=rov
    rowd=row

    rou=rou+0.50d0*kmx1
    rov=rov+0.5d0*kmy1
    row=row+0.5d0*kmz1

    ke1=reng*dt
    roetd=roet
    roet=roet+0.50d0*ke1

    call primitive()

    call commons()
    call RHS()

    ! Step 2

    co=co+1
    !write(*,*)nt,co

    kc2=rhsc*dt
    ro=rod+0.5d0*kc2

    kmx2=rmx*dt
    kmy2=rmy*dt
    kmz2=rmz*dt

    rou=roud+0.50d0*kmx2
    rov=rovd+0.50d0*kmy2
    row=rowd+0.50d0*kmz2

    ke2=reng*dt
    roet=roetd+0.50d0*ke2

    call primitive()

    call commons()
    call RHS()

    ! Step 3 

    co=co+1
    !write(*,*)nt,co

    kc3=rhsc*dt
    ro=rod+kc3

    kmx3=rmx*dt
    kmy3=rmy*dt
    kmz3=rmz*dt

    rou=roud+kmx3
    rov=rovd+kmy3
    row=rowd+kmz3

    ke3=reng*dt
    roet=roetd+ke3

    call primitive()

    call commons()
    call RHS()

    ! Step 4

    co=co+1
    !write(*,*)nt,co

    kc4=rhsc*dt
    ro=rod+(1.0d0/6.0d0)*(kc1+2*kc2+2*kc3+kc4)

    kmx4=rmx*dt
    kmy4=rmy*dt
    kmz4=rmz*dt

    rou=roud+(1.0d0/6.0d0)*(kmx1+2*kmx2+2*kmx3+kmx4)
    rov=rovd+(1.0d0/6.0d0)*(kmy1+2*kmy2+2*kmy3+kmy4)
    row=rowd+(1.0d0/6.0d0)*(kmz1+2*kmz2+2*kmz3+kmz4)

    ke4=reng*dt
    roet=roetd+(1.0d0/6.0d0)*(ke1+2*ke2+2*ke3+ke4)

    call primitive()

  end subroutine rungakutta

  !***********************************************************
  ! Computes the RHS for all equations
  !***********************************************************

  subroutine RHS()

    use mpidef
    use gdef
    implicit none

    call rhsct()
    call rhsmom()
    call rhseng()

  end subroutine RHS

  !***********************************************************
  ! Computes Primitive Variables ro,u,v,w,T and p in physical
  ! space from integral variables ro,rou,rov,row and roet in 
  ! physical space.
  !***********************************************************


  subroutine primitive()

    use mpidef
    use gdef
    implicit none

    double complex eps

    ! Solve for primitive variables

    u=rou/ro
    v=rov/ro
    w=row/ro

    roe=roet-0.50d0*ro*(u**2+v**2+w**2)
    e=roe/ro
    et=roet/ro

    T=gamma*(gamma-1)*e
    p=(ro*T)/gamma

    mu=T**no
    Kt=(mu)/(Re*Pr*(gamma-1))

  end subroutine primitive

  !******************************************************
  ! COMPUTES DERIVATIVE OF A 3D FUNCTION
  ! F=3D FUCTION WHOSE DERIVATIVE IS TO BE COMPUTED
  ! FLAG = 1,2,3 FOR X,Y,Z DERIVATIVES 
  !******************************************************

  subroutine deriv(f,flag)

    integer flag,i,j,k
    !integer*8 plan
    !double complex in(1:NX),out(1:NX)
    double complex f(fx:lx,fy:ly,fz:lz)

    ! fft along x

    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan) 
          f(1:NX,j,k)=out(1:NX)/NX
       end do
    end do

    ! Transpose Array(x,y,z --> y,x,z)

    dx=dy
    nax=npy
    call transposexy(f)

    ! fft along y

    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan) 
          f(1:NX,j,k)=out(1:NX)/NX
       end do
    end do

    ! Transpose Array(y,x,z --> z,x,y) 

    dx=dz
    nax=npz
    call transposexz(f)

    ! fft along z

    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan) 
          f(1:NX,j,k)=out(1:NX)/NX
       end do
    end do

    ! Derivative along x

    if(flag==1)then
       do k=fz,lz
          do i=fx,lx
             do j=fy,ly
                f(i,j,k)=ic*kx(j)*f(i,j,k);
             end do
             if(((NX/2+1)>=fy) .and. ((NX/2+1)<=ly))then
                f(i,NX/2+1,k)=0.0d0+ic*0.0d0
             end if
          end do
       end do
    end if

    ! Derivative along y

    if(flag==2)then
       do i=fx,lx
          do j=fy,ly
             do k=fz,lz
                f(i,j,k)=ic*ky(k)*f(i,j,k);
             end do
             if(((NY/2+1)>=fz) .and. ((NY/2+1)<=lz))then 
                f(i,j,NY/2+1)=0.0d0+ic*0.0d0
             end if
          end do
       end do
    end if

    ! Derivative along z

    if(flag==3)then
       do j=fy,ly
          do k=fz,lz
             do i=fx,lx
                f(i,j,k)=ic*kz(i)*f(i,j,k);
             end do
             f(NZ/2+1,j,k)=0.0d0+ic*0.0d0
          end do
       end do
    end if

    ! Inverse fft along z

    !tl=NX
    !vl=1
    !vs=tl
    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan1) 
          f(1:NX,j,k)=out(1:NX)
       end do
    end do

    ! Transpose Array (z,x,y --> y,x,z)

    dx=dz
    nax=npz
    call transposexz(f)

    ! Inverse fft along y

    !tl=NX
    !vl=1
    !vs=tl
    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan1) 
          f(1:NX,j,k)=out(1:NX)
       end do
    end do

    ! Transpose Array (y,x,z --> x,y,z)

    dx=dy
    nax=npy
    call transposexy(f)

    ! Inverse fft along x

    !tl=NX
    !vl=1
    !vs=tl
    do k=fz,lz
       do j=fy,ly
          in(1:NX)=f(1:NX,j,k)
          !call dfftw_plan_dft_1d(plan,NX,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
          call dfftw_execute(plan1) 
          f(1:NX,j,k)=out(1:NX)
       end do
    end do

    !call dfftw_destroy_plan(plan)

  end subroutine deriv
  !  (x,y,z <--> y,x,z)

  !*************************************************
  ! COMPUTES TRANSPOSE OF A 3D COMPLEX FUNCTION f
  ! ALONG X AND Y
  !*************************************************

  subroutine transposexy(f)

    use mpidef
    use gdef
    implicit none

    double complex f(fx:lx,fy:ly,fz:lz),send(dx,dy,dz),recv(dx,dy,dz)
    integer i,j,k,lfx,stack(nax),p,sendc,recvc,dest,source,count,index
    integer status(MPI_STATUS_SIZE)

    sendc=dx*dy*dz
    recvc=dx*dy*dz

    index=myrank/npz
    count=npy-index

    do i=1,nax
       stack(i)=count
       count=count+1
       if(count>nax)then
          count=1
       end if
    end do

    do p=1,nax
       ix=stack(p)
       lfx=dx*(ix-1)+1

       do k=1,dz
          do j=1,dy
             do i=1,dx
                send(j,i,k)=f(lfx+i-1,fy+j-1,fz+k-1)
                recv(i,j,k)=0.0d0
             end do
          end do
       end do

       source=npz*(ix-1)+mod(myrank,npz)
       dest=npz*(ix-1)+mod(myrank,npz)

       call mpi_sendrecv(send,sendc,mpi_double_complex,dest,1,recv,&
            recvc,mpi_double_complex,source,1,comm2d,status,ierr)

       do k=1,dz
          do j=1,dy
             do i=1,dx
                f(lfx+i-1,fy+j-1,fz+k-1)=recv(i,j,k)
             end do
          end do
       end do

    end do

  end subroutine transposexy


  ! (y,x,z <--> z,x,y>

  !***********************************************
  !   COMPUTES TRANSPOSE OF 3D COMPLEX FUNCTION f
  !   ALONG X AND Z
  !***********************************************

  subroutine transposexz(f)

    use mpidef
    use gdef
    implicit none

    double complex f(fx:lx,fy:ly,fz:lz),send(dx,dy,dz),recv(dx,dy,dz)
    integer i,j,k,lfx,stack(nax),p,sendc,recvc,dest,source,count,index
    integer status(MPI_STATUS_SIZE)

    sendc=dx*dy*dz
    recvc=dx*dy*dz

    index=mod(myrank,npz)
    count=npz-index

    do i=1,nax
       stack(i)=count
       count=count+1
       if(count>nax)then
          count=1
       end if
    end do

    do p=1,nax
       ix=stack(p)
       lfx=dx*(ix-1)+1

       do k=1,dz
          do j=1,dy
             do i=1,dx
                send(k,j,i)=f(lfx+i-1,fy+j-1,fz+k-1)
                recv(i,j,k)=0.0d0
             end do
          end do
       end do

       source=(ix-1)+npz*(myrank/npz)
       dest=(ix-1)+npz*(myrank/npz)

       call mpi_sendrecv(send,sendc,mpi_double_complex,dest,1,recv,&
            recvc,mpi_double_complex,source,1,comm2d,status,ierr)

       do k=1,dz
          do j=1,dy
             do i=1,dx
                f(lfx+i-1,fy+j-1,fz+k-1)=recv(i,j,k)
             end do
          end do
       end do

    end do

  end subroutine transposexz

end subroutine isoturb





























