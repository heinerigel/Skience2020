!---------------------------------------------------------------
module module_sem1d
!---------------------------------------------------------------
  implicit none

  public :: modname,load_sem_dat,read_model,get_model,model &
            ,mesh_model,mesh1d,simu,newmark,ACC,velocity_stress,init_correctors,correctors
  private
  type mesh1d
!ngll: element degree
!nelem: numbre fo elements in the mesh
!coord(i,j): coordinates of gll point number i in element j
!jacob: jacobian of a transformation of an lement
     integer :: ngll,nelem
     doubleprecision :: cfl
     integer, dimension(:)  ,  pointer :: elem2part 
     doubleprecision, dimension(:)  ,  pointer  :: xgll,wgll,inver
     doubleprecision, dimension(:,:),  pointer :: deriv
     doubleprecision, dimension(:,:), pointer :: coord,mu,rho
     doubleprecision, dimension(:), pointer :: jacob     
  end type mesh1d
  type model
!nbp: number of discontinuities+1 (number of parts)
!position_interface(nbp+1): discontiniuities position
!np: number of discretisation points in a part
!rho, v, mu: density, velocity, elastic parameter
     integer :: nbp
     doubleprecision :: k2 !for homogeneized model
     integer, dimension(:), pointer :: np
     doubleprecision :: xmin,xmax
     doubleprecision, dimension(:), pointer ::  position_interface
     doubleprecision, dimension(:,:), pointer :: xp,rho,v,mu     
     doubleprecision, dimension(:,:,:), pointer :: y2
  end type model
  type correctors
     logical :: allcorrectors
     integer :: ngll,nelem,horder,nloc
     doubleprecision, dimension(:,:), pointer :: coord,x1,x2,x2p
     doubleprecision, dimension(:), pointer :: coordloc,x1loc
  end type correctors
  type simu
     type(mesh1d) :: msh    
     doubleprecision :: xsource,fmax_mesh,dt,timel,f0,l0,t0,tsnap,abc1,abc2
     doubleprecision, dimension(:,:), pointer :: Fext,iM,rec_interp,traces
     doubleprecision, dimension(:), pointer :: g,xrec
     integer, dimension(:), pointer :: rec_elem
     integer :: Nt,isnap,Nrec,scheme,horder
     logical :: t0shift
     character(len=20) :: receivers_file,correctorfile
     type(correctors) :: corr,corrmsh
    logical :: abc
  end type simu
!
  character(len=30) :: modname
  doubleprecision, parameter :: beta=0.0d0,gamma=0.5d0,alpha=0.5d0,courant=0.81
  integer, parameter :: ACC=1,VSTRESS=2  
!--------------------------------------------------------------------
contains
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine newmark(sim)
!--------------------------------------------------------------------
    implicit none
!
    type(simu), intent(inout) :: sim
    doubleprecision :: coef1,coef2,coef3,coef4,dt,t1,t2,v1,v2,eps
    doubleprecision, dimension(sim%msh%ngll,sim%msh%nelem) :: U,Fint,V,A,Up,Vp,Uc
    integer :: i,ie,j,k,i_interne,ng,ne
    character(len=13) :: name
!
    doubleprecision, dimension(sim%msh%NGLL,sim%msh%nelem) :: U_X
    logical :: save_deformation=.true.

    
!
    call init_simu(sim)
    call f_ext(sim)
!
    ng=sim%msh%ngll; ne=sim%msh%nelem
    dt=sim%dt
    coef1    =dt*(1.d0 - gamma)
    coef2    =dt*dt*(0.5d0 - beta)
    coef3    =dt*gamma
    coef4    =dt*dt*beta
    if (sim%abc.and.coef4>1.d-8) STOP 'beta must be zero to use abc'
!
    U=0.d0;V=0.d0; A=0.d0;k=0
    open(12,file="snapshot_table")
    write(12,'("snatpshot number, corresponding simulation time")') 
   call cpu_time(t1)
   do i=1,sim%Nt
     call save_traces(sim,i,U)
     if (mod(i,sim%isnap)==0) then
        if (sim%horder>0.and.sim%corr%allcorrectors)  then
           call apply_correctors_snapshot(sim,U,Uc) 
        else
           Uc(:,:)=U(:,:)
        endif
        k=k+1
        write(12,*) k,(i-1)*sim%dt
        write(name,'("snapshot",i3.3)')k
        open(11,file=name)
        do ie=1,sim%msh%nelem
           do j=1,sim%msh%ngll-1
              write(11,*) sim%msh%coord(j,ie),Uc(j,ie)
           enddo
        enddo
        close(11)
        
        if (save_deformation) then
           write(name,'("snapshot_e",i3.3)')k
           call derive_part(Uc,U_X,sim%msh)
           open(11,file=name)
           do ie=1,sim%msh%nelem
              do j=1,sim%msh%ngll-1
                 write(11,*) sim%msh%coord(j,ie),U_X(j,ie)
              enddo
           enddo
           close(11)
           write(name,'("snapshot_s",i3.3)')k
           open(11,file=name)
           do ie=1,sim%msh%nelem
              do j=1,sim%msh%ngll-1
                 write(11,*) sim%msh%coord(j,ie),U_X(j,ie)*sim%msh%mu(j,ie)
              enddo
           enddo
           close(11)
        endif
     endif
     if (mod(i,sim%Nt/10)==0) then
        call cpu_time(t2)
        print*,'Time step n: ',i,sngl(t2-t1),sngl(maxval(abs(sim%traces)))
        t1=t2 
     endif
!
! prediction 
! 
     Up(:,:)=U(:,:) + dt*V(:,:) + coef2*A(:,:)
     Vp(:,:)=V(:,:) + coef1*A(:,:)
!
! billan des forces
!
     call force_interne(sim,Up,Fint)
     call assemble(sim%msh,Fint)
     i_interne=0
     eps=1
     do while (eps>1.d-4)
        if (sim%abc) then
           if (i_interne==0) &
               call apply_clay_abc_newmark(sim,Fint,U,V,A,coef3,coef4,v1,v2,eps,0)
        else
           eps=0.d0
        endif
        if (i_interne==0) then
           A(:,:)=(sim%Fext(:,:)*sim%g(i)-Fint(:,:))*sim%iM(:,:)        
!
! correction
!
           V(:,:)=Vp(:,:)+coef3*A(:,:)
           U(:,:)=Up(:,:)+coef4*A(:,:)
        else
           if (.not.sim%abc) STOP 'I should not be here ...'
           call apply_clay_abc_newmark(sim,Fint,U,V,A,coef3,coef4,v1,v2,eps,1)
        endif
         i_interne= i_interne+1
     enddo
!
  enddo
  close(12)
  print*,'Dumping traces'
  call dump_traces(sim)
  print*,'Newmark finished' 
!
!--------------------------------------------------------------------
  end subroutine newmark
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine velocity_stress(sim)
!--------------------------------------------------------------------
    implicit none
!
    type(simu), intent(inout) :: sim
    doubleprecision :: coef1,coef2,coef3,coef4,dt,t1,t2
!these field are not all necessary. Could use  less memory!
    doubleprecision, dimension(sim%msh%ngll,sim%msh%nelem) :: Un,Unp1s2,Vn,Fint,dV
    integer :: i,ie,j,k
    character(len=13) :: name
!
    call init_simu(sim)
    call f_ext(sim)
!
    dt=sim%dt
!
    Unp1s2=0.d0;Vn=0.d0;k=0
    open(12,file="snapshot_table")
    write(12,'("snatpshot number, corresponding simulation time")') 
    call cpu_time(t1)
    if (sim%horder>0) STOP 'the homogeneization correctors are not implemented yet for the velocity stress time scheme'
   do i=1,sim%Nt
!saving:
     Un(:,:)=Unp1s2(:,:)
!prediction Un+1/2=Un-1/2+dt*Vn
     Unp1s2(:,:)=Unp1s2(:,:)+dt*Vn(:,:)
!resolution
     call force_interne(sim,Unp1s2,Fint)
     call assemble(sim%msh,Fint)
     dV(:,:)=(sim%Fext(:,:)*sim%g(i)-Fint(:,:))*sim%iM(:,:)*dt        
!correction Vn+1=Vn+dV
     Vn(:,:)=Vn(:,:)+dV(:,:)
!for outpout:Un=(Un-1/2+Un+1/2):
     Un(:,:)=(Un(:,:)+Unp1s2(:,:))/2.d0
!
     call save_traces(sim,i,Un)
!
     if (mod(i,sim%isnap)==0) then
        k=k+1
        write(12,*) k,(i-1)*sim%dt
        write(name,'("snapshot",i3.3)')k
        open(11,file=name)
        do ie=1,sim%msh%nelem
           do j=1,sim%msh%ngll-1
              write(11,*) sim%msh%coord(j,ie),Un(j,ie)
           enddo
        enddo
        close(11)
     endif
!
     if (mod(i,sim%Nt/10)==0) then
        call cpu_time(t2)
        print*,'Time step n: ',i,sngl(t2-t1),sngl(maxval(abs(sim%traces)))
        t1=t2 
     endif
  enddo
  call dump_traces(sim)
  close(12)
!--------------------------------------------------------------------
  end subroutine velocity_stress
!-------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine init_receiver(sim)
!--------------------------------------------------------------------
    use funaro
    implicit none
    type(simu), intent(inout) :: sim
!
    integer :: i,n
    doubleprecision :: x,x_loc
    doubleprecision, dimension(sim%msh%NGLL) :: wx
!
    open(11,file=sim%receivers_file,status='old')
    read(11,*) sim%Nrec
    allocate(sim%xrec(sim%Nrec))
    do i=1,sim%Nrec
       read(11,*) sim%xrec(i)
    enddo
    close(11)
    allocate(sim%traces(sim%Nt,sim%Nrec),sim%rec_interp(sim%msh%ngll,sim%Nrec) &
             ,sim%rec_elem(sim%Nrec))
    sim%rec_elem(:)=-1
    do i=1,sim%Nrec
       x=sim%xrec(i)
       do n=1,sim%msh%nelem
          if (x>=sim%msh%coord(1,n).and.x<sim%msh%coord(sim%msh%NGLL,n)) sim%rec_elem(i)=n
       enddo       
       if (sim%rec_elem(i)<0) then
          print*,'receiver:',i,' xrec=',x
          STOP 'init_receiver can''t locate a receiver in the mesh'
       endif
       x_loc=(x-sim%msh%coord(1,sim%rec_elem(i)))/sim%msh%jacob(sim%rec_elem(i))-1.0d0
       call def_heta(x_loc,sim%msh%NGLL,wx)
       sim%rec_interp(:,i)=wx(:)
    enddo
!!$!test
!!$do i=1,300
!!$x_loc=(i-1)*2./299.-1
!!$call def_heta(x_loc,9,wx)
!!$do n=1,9
!!$write(300+n,*) x_loc,wx(n)
!!$enddo
!!$enddo
!!$!fin test
!--------------------------------------------------------------------
  end subroutine init_receiver
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine save_traces(sim,it,U)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim
    integer:: it
    doubleprecision, dimension(:,:), intent(in) :: U
!    
    integer :: i
!
    do i=1,sim%Nrec
       sim%traces(it,i)=SUM(sim%rec_interp(:,i)*U(:,sim%rec_elem(i)))
    enddo
!--------------------------------------------------------------------
  end subroutine save_traces
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine dump_traces(sim)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim    
!
    character(len=20) :: name
    integer :: it,i
!
    do i=1,sim%Nrec
       print*,i
       write(name,'("trace",i4.4)') i
       open(11,file=name)
       if (sim%t0shift) then
          do it=1,sim%Nt
             write(11,*) (it-1)*sim%dt-sim%t0,sim%traces(it,i)
          enddo
       else
          do it=1,sim%Nt
             write(11,*) (it-1)*sim%dt,sim%traces(it,i)
          enddo
       endif
       close(11)
    enddo
!--------------------------------------------------------------------
  end subroutine dump_traces
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine allocate_simu(sim)
!--------------------------------------------------------------------
    implicit none
    type(simu)   , intent(out):: sim
!
    integer :: ng,ne
    ne=sim%msh%nelem
    ng=sim%msh%ngll
    allocate(sim%iM(ng,ne),sim%Fext(ng,ne))
!--------------------------------------------------------------------
  end subroutine allocate_simu
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine load_sem_dat(sim)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim
    open(11,file="sem1d.dat",status='old')
    read(11,*)
    read(11,*) modname
    read(11,*) 
    read(11,*) sim%msh%ngll
    read(11,*) 
    read(11,*) sim%fmax_mesh
    read(11,*) 
    read(11,*) sim%xsource
    read(11,*) 
    read(11,*) sim%f0
    read(11,*) 
    read(11,*) sim%t0shift
    read(11,*) 
    read(11,*) sim%timel
    read(11,*) 
    read(11,*) sim%dt
    read(11,*) 
    read(11,*) sim%tsnap
    read(11,*) 
    read(11,*) sim%receivers_file
    read(11,*) 
    read(11,*) sim%scheme
    read(11,*) 
    read(11,*) sim%abc
    read(11,*) 
    read(11,*) sim%horder
    sim%corr%horder=sim%horder
    if (sim%horder>0) then 
       read(11,*) 
       read(11,*) sim%correctorfile
    else
       sim%corr%horder=-1
    endif
    close(11)
!    sim%l0=2.5d0/(4.d0*atan(1.d0)*sim%f0) !2.5*pi/f0
    sim%l0=1.d0/(4.d0*atan(1.d0)*sim%f0) !1./pi/f0
    sim%t0=5.*sim%l0    
    print*,'********>>>> source t0=',sim%t0
    if (sim%scheme==ACC) then
       print*,'The acceleration time evolution scheme has been selected'
    else if (sim%scheme==VSTRESS) then
       print*,'The velocity-stress time evolution scheme has been selected'
    else
       STOP 'The  time evolution scheme must be either acceleration (1) or velocity-stress (2)'
    endif
!--------------------------------------------------------------------
  end subroutine load_sem_dat
!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine mesh_model(mod,sim)
!loading model file "model_name" into mod
!--------------------------------------------------------------------
  use funaro
  type(model), intent(in) :: mod
  type(simu), intent(out) :: sim
!
  integer :: i,k,nelem,j,kelem,ngll
  doubleprecision :: lambda_min,dxmin,vmax,fmax
  doubleprecision, parameter :: cst=6.d0
  fmax=sim%fmax_mesh
  ngll=sim%msh%ngll
  allocate(sim%msh%xgll(ngll),sim%msh%wgll(ngll),sim%msh%deriv(ngll,ngll))
!initialisation des polynomes de GLL
  call def_xgll (sim%msh%xgll ,NGLL)
  call def_wgll (sim%msh%wgll ,NGLL)
  call def_deriv(sim%msh%deriv,NGLL)
!
  sim%msh%cfl=1.d30
  nelem=0
  do i=1,mod%nbp
     if (sim%horder>=0) then
        lambda_min=1./(fmax/minval(mod%v(1:mod%np(i),i))+mod%k2)
     else
        lambda_min=minval(mod%v(1:mod%np(i),i))/fmax
     endif
     nelem=nelem+max(int((mod%position_interface(i+1)-mod%position_interface(i))/lambda_min*cst/real(NGLL-1)),1) 
  enddo
  print*,'lambda_min estimate:',lambda_min
  sim%msh%nelem=nelem
  allocate(sim%msh%coord(ngll,nelem),sim%msh%jacob(nelem),sim%msh%rho(ngll,nelem),sim%msh%mu(ngll,nelem))
  allocate(sim%msh%elem2part(nelem),sim%msh%inver(nelem))
  kelem=0
  do i=1,mod%nbp
     if (sim%horder>=0) then
        lambda_min=1./(fmax/minval(mod%v(1:mod%np(i),i))+mod%k2)
     else
        lambda_min=minval(mod%v(1:mod%np(i),i))/fmax
     endif
     nelem=max(int((mod%position_interface(i+1)-mod%position_interface(i))/lambda_min*cst/real(NGLL-1)),1)
     do k=1,nelem
        kelem=kelem+1
        sim%msh%elem2part(kelem)=i
        sim%msh%jacob(kelem)=(mod%position_interface(i+1)-mod%position_interface(i))/real(nelem)/2.d0
        sim%msh%inver(kelem)=1.d0/sim%msh%jacob(kelem)
        do j=1,NGLL
           if (kelem/=1) then
              sim%msh%coord(j,kelem)=sim%msh%jacob(kelem)*(sim%msh%xgll(j)+1.d0)+sim%msh%coord(NGLL,kelem-1)
           else
              sim%msh%coord(j,kelem)=sim%msh%jacob(kelem)*(sim%msh%xgll(j)+1.d0)+mod%position_interface(i)
           endif
        enddo
     enddo 
  enddo
  open(11,file='density')
  open(12,file='velocity')
  do j=1,sim%msh%nelem
     do k=1,ngll
         call get_model(mod,sim%msh%coord(k,j),sim%msh%rho(k,j),sim%msh%mu(k,j),sim%msh%elem2part(j))
         if (sim%msh%mu(k,j)<0.d0) STOP 'negative mu in mesh_model !!!'
        write(11,*) sim%msh%coord(k,j),sim%msh%rho(k,j)
        write(12,*) sim%msh%coord(k,j),sqrt(sim%msh%mu(k,j)/sim%msh%rho(k,j))
      enddo
      dxmin=sim%msh%coord(2,j)-sim%msh%coord(1,j)
      vmax=maxval(sqrt(sim%msh%mu(:,j)/sim%msh%rho(:,j)))
      if (dxmin<0) then
         print*,j,k,sim%msh%coord(2,j),sim%msh%coord(1,j)
         STOP 'dx<0!!!'
      endif
      sim%msh%cfl=min(sim%msh%cfl,dxmin/vmax)
   enddo
   close(11); close(12)
   sim%abc1=sim%msh%rho(1,1)*sqrt(sim%msh%mu(1,1)/sim%msh%rho(1,1))
   sim%abc2=sim%msh%rho(ngll,sim%msh%nelem)*sqrt(sim%msh%mu(ngll,sim%msh%nelem)/sim%msh%rho(ngll,sim%msh%nelem))

   print*,'-----Mesh information:'
   print*,'Number of elements:',sim%msh%nelem
   print*,'cfl=',sim%msh%cfl
   print*,'dxmin=',dxmin
   print*,'-----'
   
!--------------------------------------------------------------------
end subroutine mesh_model
!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine read_model(model_name,mod,homogenized)
!loading model file "model_name" into mod
!--------------------------------------------------------------------
  use module_spline
  implicit none
  character(len=*), intent(in) :: model_name
  type(model), intent(out) :: mod
  logical, optional, intent(in) :: homogenized
!
  doubleprecision :: x1,x
  integer :: i,j,n,npmax
!
  open(11,file=model_name,status='old')  
  if (present(homogenized)) then
     if (homogenized) then
        read(11,*) mod%nbp,mod%k2
     else
        read(11,*) mod%nbp
     endif
  else
     read(11,*) mod%nbp
  endif
  allocate(mod%np(mod%nbp),mod%position_interface(mod%nbp+1))
  do i=1,mod%nbp
     read(11,*) mod%position_interface(i),mod%np(i)
     if (mod%np(i)<2) STOP 'read_model: nbp should be >=2 !!!'
     do j=1,mod%np(i)
        read(11,*)
     enddo
  enddo
  npmax=maxval(mod%np)
  allocate(mod%xp(npmax,mod%nbp),mod%rho(npmax,mod%nbp),mod%v(npmax,mod%nbp),mod%mu(npmax,mod%nbp))
  rewind(11)
  read(11,*)  
  do i=1,mod%nbp
     read(11,*)
     do j=1,mod%np(i)
        read(11,*) mod%xp(j,i),mod%rho(j,i),mod%v(j,i)
        mod%mu(j,i)=mod%rho(j,i)*mod%v(j,i)**2
     enddo
  enddo
  mod%xmin=mod%xp(1,1)
  mod%xmax=mod%xp(mod%np(mod%nbp),mod%nbp) 
  mod%position_interface(mod%nbp+1)=mod%xmax
  close(11)
  allocate(mod%y2(npmax,mod%nbp,3))
  do i=1,mod%nbp 
     n=mod%np(i)
     call spline(mod%xp(1:n,i),mod%rho(1:n,i),1.d30,1.d30,mod%y2(1:n,i,1))
     call spline(mod%xp(1:n,i),mod%v  (1:n,i),1.d30,1.d30,mod%y2(1:n,i,2))
     call spline(mod%xp(1:n,i),mod%mu (1:n,i),1.d30,1.d30,mod%y2(1:n,i,3))
  enddo
!
  
!--------------------------------------------------------------------
end subroutine read_model
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine get_model(mod,x,rho,mu,npart)
!--------------------------------------------------------------------
  use module_spline
  implicit none
  type(model ), intent(in) :: mod
  doubleprecision, intent(in ) :: x
  doubleprecision, intent(out) :: rho,mu
  integer, optional, intent(in) :: npart
!
  integer :: n,np
!
  if (x<=mod%xmin) then
     rho=mod%rho(1,1)
     mu=mod%mu(1,1)
  else if( x>=mod%xmax) then
     rho=mod%rho(mod%np(mod%nbp),mod%nbp)
     mu =mod%mu (mod%np(mod%nbp),mod%nbp)
  else
     if (present(npart)) then
        np=npart
     else
        np=locate(mod%position_interface,mod%nbp+1,x)
     endif
     n=mod%np(np)
     rho=splint(mod%xp(1:n,np),mod%rho(1:n,np),mod%y2(1:n,np,1),x)
     mu =splint(mod%xp(1:n,np),mod%mu (1:n,np),mod%y2(1:n,np,3),x)
  endif
!--------------------------------------------------------------------
end subroutine get_model
!--------------------------------------------------------------------
!-----------------------------------------------------------------
integer function locate(inter,n,x)
!-----------------------------------------------------------------
  implicit none
  integer, intent(in) :: n
  doubleprecision, dimension(n), intent(in) :: inter
  doubleprecision, intent(in) :: x
!
  integer :: i
  i=1
  if (abs (x-inter(1))/maxval(inter(:)) <1.d-10) then
     i=1
  else if ((x-inter(n))/maxval(inter(:))  >0.0) then
     i=n
  else
     do while (  x <= inter(i) .or. x > inter(i+1) )
        i=i+1
        if (i > n -1) then
           print*,'i=',i,'n=',n
           print*,'x=',x
           print*,'inter=',inter
           stop 'locate: failed to locate x in inter!'
        endif
     enddo
  endif
  locate=i
!-----------------------------------------------------------------
end function locate
!-----------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine init_simu(sim)
!construction de l'inverse de la matrice de masse
!-------------------------------------------------------------------------
      implicit none
      type(simu)   , intent(inout):: sim
!
      integer :: i
      doubleprecision :: t
!
      call allocate_simu(sim)
      call init_masse(sim)
      if (sim%dt<0.d0) then
         sim%dt=courant*sim%msh%cfl
!it seems the clayton and Egquist are debcreasing the courant number ...? strange.
         if (sim%abc)  sim%dt=sim%dt*0.68 
         print*,'Computed time step:',sim%dt
      else
!         sim%t0shift=.false.
         print*,'Maximum advised time step:',courant*sim%msh%cfl
         print*,'Input time step:',sim%dt
      endif
      sim%Nt=sim%timel/sim%dt+1
      sim%isnap=int(sim%tsnap/sim%dt)
      print*,'Number of time step:',sim%Nt
      print*,'Snapshot every :',sim%isnap
      allocate(sim%g(sim%Nt))
      open(11,file='source.gnu')
      do i=1,sim%Nt
         if (sim%scheme==ACC) then
!g(i)=time step i+1
            t=(i-1+1)*sim%dt
         else
!!g(i)=time step i+1/2
            t=(i-1+0.5)*sim%dt
         endif
         sim%g(i)=ggauss(t,sim%t0,sim%l0)
         write(11,*) t,sim%g(i)
      enddo
      close(11)   
      call init_receiver(sim)
!-------------------------------------------------------------------------
    end subroutine init_simu
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine init_masse(sim)
!construction de l'inverse de la matrice de masse
!-------------------------------------------------------------------------
      implicit none
      type(simu)   , intent(inout):: sim
!
      integer :: n,i,nelem,ngll
!
      nelem=sim%msh%nelem
      ngll =sim%msh%ngll
      do n=1,nelem
         do i=1,NGLL
            sim%iM(i,n)=sim%msh%wgll(i)*sim%msh%rho(i,n)*sim%msh%jacob(n)
         enddo
      enddo
      call assemble(sim%msh,sim%iM)
      sim%iM(:,:)=1.d0/sim%iM(:,:)
!-------------------------------------------------------------------------
    end subroutine init_masse
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine assemble(msh,U)
!-------------------------------------------------------------------------
      implicit none
!
      type(mesh1d), intent(in) :: msh
      doubleprecision, dimension(:,:), intent(inout)  :: U
!
      integer :: i,nelem,ngll
!
      nelem=msh%nelem
      ngll =msh%ngll
      do i=2,nelem
         U(NGLL,i-1)=U(NGLL,i-1)+U(1,i)
         U(1   ,i  )=U(NGLL,i-1)
      enddo
!-------------------------------------------------------------------------
    end subroutine assemble
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine f_ext(sim)
!-------------------------------------------------------------------------
      use funaro
      implicit none
    type(simu), intent(inout) :: sim
!
      integer :: nelem,ngll,xelem,n
      doubleprecision :: x,x_loc 
      doubleprecision, dimension(sim%msh%ngll) :: wx,d,dd,h
!
      nelem=sim%msh%nelem
      ngll =sim%msh%ngll
!
! localisation de x
!
      x=sim%xsource
      xelem=0
      do n=1,nelem
         if (x>=sim%msh%coord(1,n).and.x<sim%msh%coord(NGLL,n)) xelem=n
      enddo
      print*,'The source is in element number: ',xelem
      if (xelem==0) stop 'f_ext, can''t locate source element (out of domain??)!'
!
      sim%Fext(:,:)=0.0d0
      x_loc=(x-sim%msh%coord(1,xelem))/sim%msh%jacob(xelem)-1.0d0
      call def_heta(x_loc,NGLL,wx)
      sim%Fext(:,xelem) =wx(:)! pas de *jacob(xelem), c'est in dirac: \int h \delta dx = h
      if (sim%corr%horder>0) then

         call def_derivx(x_loc,h,d,dd,ngll)
          sim%Fext(:,xelem) = sim%Fext(:,xelem)+sim%corr%x1loc(1)*d(:)*sim%msh%inver(xelem)
      endif
      call assemble(sim%msh,sim%Fext)
!
!-------------------------------------------------------------------------
    end subroutine f_ext
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine derive_part(U,U_X,msh)
!-------------------------------------------------------------------------
      implicit none
!
      type(mesh1d), intent(in) :: msh
      doubleprecision, dimension(:,:), intent(in)  :: U
      doubleprecision, dimension(:,:), intent(out) :: U_X
!
      integer :: n,i
!
      do n=1,msh%nelem
         do i=1,msh%NGLL
            U_X(i,n)=SUM(U(:,n)*msh%deriv(:,i))*msh%inver(n)
         enddo
      enddo
!
!-------------------------------------------------------------------------
    end subroutine derive_part
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine force_interne(sim,U,F)
!-------------------------------------------------------------------------
      implicit none
      type(simu), intent(inout) :: sim
      doubleprecision, dimension(:,:), intent(in)  :: U
      doubleprecision, dimension(:,:), intent(out) :: F
!

      integer :: i,j,k,l,n
      doubleprecision, dimension(sim%msh%NGLL,sim%msh%nelem) :: U_X
!      
      call derive_part(U,U_X,sim%msh)
      do n=1,sim%msh%nelem
         do i=1,sim%msh%NGLL
            F(i,n)=SUM(sim%msh%mu(:,n)*U_X(:,n)*sim%msh%deriv(i,:)*sim%msh%wgll(:))*sim%msh%inver(n)*sim%msh%jacob(n)
         enddo
      enddo
!-------------------------------------------------------------------------
    end subroutine force_interne
!-------------------------------------------------------------------------
!----------------------------------------------------------------------
  doubleprecision function ggauss(t,t0,L0)
!----------------------------------------------------------------------
    doubleprecision :: t,t0,L0
    ggauss=(1.d0-2.d0*((t-t0)/L0)**2)*exp(-((t-t0)/L0)**2)
!
!----------------------------------------------------------------------
  end function ggauss
!----------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine read_correctors(sim)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim
    integer :: i,j
    open(11,file=sim%correctorfile,status='old')
    read(11,*) sim%corr%allcorrectors,sim%corr%horder,sim%corr%nloc,sim%corr%ngll,sim%corr%nelem
    select case(sim%corr%horder)
    case(1)
       allocate(sim%corr%coordloc(sim%corr%nloc),sim%corr%x1loc(sim%corr%nloc))
       do i=1,sim%corr%nloc
          read(11,*) sim%corr%coordloc(i),sim%corr%x1loc(i)
       enddo
       if (sim%corr%allcorrectors) then
          allocate(sim%corr%coord(sim%corr%ngll,sim%corr%nelem),sim%corr%x1(sim%corr%ngll,sim%corr%nelem))
          do j=1,sim%corr%nelem
             do i=1,sim%corr%ngll
                read(11,*) sim%corr%coord(i,j),sim%corr%x1(i,j)
             enddo
          enddo
       endif
    case(2)
       STOP 'read_corrector: can''t be there!'
    case default
       STOP 'read_corrector: can''t be there!'
    end select
    close(11)

!--------------------------------------------------------------------
  end subroutine read_correctors
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine init_correctors(sim)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim
!
    call read_correctors(sim)
    sim%corrmsh%ngll=sim%msh%ngll
    sim%corrmsh%nelem=sim%msh%nelem
    if (sim%corr%allcorrectors) then
       allocate(sim%corrmsh%coord(sim%msh%ngll,sim%msh%nelem),sim%corrmsh%x1(sim%msh%ngll,sim%msh%nelem))
       sim%corrmsh%coord(:,:)=sim%msh%coord(:,:)
       call interp_correctors(sim%corr,sim%corrmsh)
    endif
!--------------------------------------------------------------------
  end subroutine init_correctors
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine interp_correctors(corr1,corr2)
!--------------------------------------------------------------------
    use funaro
    implicit none
    type(correctors), intent(in   ) :: corr1
    type(correctors), intent(inout) :: corr2
!
    integer :: ie1,ig,ie2,ig2
    doubleprecision, dimension(corr1%nelem+1) :: interfaces
    doubleprecision, dimension(corr1%ngll) :: wx
    doubleprecision, parameter :: eps=1.d-2
    doubleprecision :: x_loc
    
!
    interfaces(1)=corr1%coord(1,1)-eps
    interfaces(2:corr1%nelem+1)=corr1%coord(corr1%ngll,:)
    interfaces(corr1%nelem+1)=interfaces(corr1%nelem+1)+eps
!
    do ie2=1,corr2%nelem
       do ig2=1,corr2%ngll
          ie1=locate(interfaces,corr1%nelem+1,corr2%coord(ig2,ie2))
          x_loc=2.d0*(corr2%coord(ig2,ie2)-interfaces(ie1))/(interfaces(ie1+1)-interfaces(ie1))-1.d0
          call def_heta(x_loc,corr1%ngll,wx)
          select case(corr1%horder)
          case(1)
            corr2%x1(ig2,ie2)=SUM(wx(:)*corr1%x1(:,ie1)) 
          end select
       enddo
    enddo 
!--------------------------------------------------------------------
  end subroutine interp_correctors
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine apply_correctors_snapshot(sim,U,Uc)
!works for samemsh
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(in) :: sim
    doubleprecision, dimension(:,:), intent(in ) :: U
    doubleprecision, dimension(:,:), intent(out) :: Uc
!
    doubleprecision, dimension(sim%msh%ngll,sim%msh%nelem) :: Ux
!
    call derive_part(U,Ux,sim%msh)
    Uc(:,:)=U(:,:)+sim%corrmsh%x1(:,:)*Ux(:,:)
!--------------------------------------------------------------------
  end subroutine apply_correctors_snapshot
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine apply_clay_abc_newmark(sim,F,U,V,A,coef3,coef4,v1,v2,eps,flag)
!--------------------------------------------------------------------
    implicit none
    type(simu), intent(inout) :: sim
    doubleprecision, dimension(:,:), intent(inout) :: V,F,A,U
    doubleprecision, intent(in) :: coef3,coef4
    doubleprecision, intent(inout) :: v1,v2,eps
    integer, intent(in) :: flag
!
    integer :: ne,ng
    doubleprecision :: dv1,dv2,da1,da2,df1,df2
!
    ng=sim%msh%ngll; ne=sim%msh%nelem
    if (flag==0) then
       v1=V( 1, 1)
       v2=V(ng,ne)
       F( 1, 1)=F( 1, 1)+sim%abc1*v1
       F(ng,ne)=F(ng,ne)+sim%abc2*v2
!
       df1=+sim%abc1*v1
       df2=+sim%abc1*v2
       dA1=-df1*sim%iM( 1, 1)        
       dA2=-df2*sim%iM(ng,ne)  
        eps=max(abs(da1/(A( 1, 1)+1.d-12)),abs(da2/(A(ng,ne)+1.d-12)))
    else
       dv1=V( 1, 1)-v1
       dv2=V(ng,ne)-v2
!
       df1=+sim%abc1*dv1
       df2=+sim%abc2*dv2
       dA1=-df1*sim%iM( 1, 1)        
       dA2=-df2*sim%iM(ng,ne)        
       v1=V( 1, 1)
       v2=V(ng,ne)
!
! abc correction
!
       V( 1, 1)=V( 1, 1)+coef3*dA1
       U( 1, 1)=U( 1, 1)+coef4*dA1
       A( 1, 1)=A( 1, 1)+      dA1
       V(ng,ne)=V(ng,ne)+coef3*dA2
       U(ng,ne)=U(ng,ne)+coef4*dA2
       A(ng,ne)=A(ng,ne)+      dA2
       eps=max(abs(da1/(A( 1, 1)+1.d-12)),abs(da2/(A(ng,ne)+1.d-12)))
    endif
!--------------------------------------------------------------------
  end subroutine apply_clay_abc_newmark
!--------------------------------------------------------------------
!---------------------------------------------------------------
end module module_sem1d
!---------------------------------------------------------------
