 


 subroutine wacomm(udims,vdims,wdims,                               &
                  uwacom,vwacom,wwacom,akt,mask_u,mask_v,mask_rho, &
                  conc,zeta,s_w,s_rho,lon_u,lat_v,                 &
                  dti,deltat,ninterp,ncalc)

!--------------------------------------------------------------------------!
 use omp_lib
 use  comblk

 implicit none

 integer :: ninterp, ncalc

! velocity
 integer :: udims(4),vdims(4),wdims(4)

 real :: uwacom(  udims(1)  ,0:udims(2)-1,-udims(3)+1:0)
  

 real :: vwacom(0:vdims(1)-1,  vdims(2)  ,-vdims(3)+1:0)
 real :: wwacom(0:wdims(1)-1,0:wdims(2)-1,-wdims(3)+1:0)
 real :: akt  (0:wdims(1)-1,0:wdims(2)-1,-wdims(3)+1:0)

! mask
 real*8 :: mask_u  (  udims(1),0:vdims(2))
 real*8 :: mask_v  (0:udims(1),  vdims(2))
 real*8 :: mask_rho(0:udims(1),0:vdims(2))

! time step
 real :: dti, deltat

! conc
 real :: conc(0:udims(1),0:vdims(2),-wdims(3)+2:0)

! vertical stretch and depth
 real :: zeta(0:udims(1),0:vdims(2))
 real*8 :: s_w(-wdims(3)+1:0), s_rho(-udims(3)+1:0)
 real*8 :: depth(-wdims(3)+2:0)

! lat/lon coord
 real*8 :: lon_u(1:udims(1),0:vdims(2))
 real*8 :: lat_v(0:udims(1),1:vdims(2))

! miscellaneous
 real*8 :: rn
 integer :: time(8), seed(8)

 integer :: ixx, iyy, izz
 integer :: i, j, k, iint, ll
 real :: vs, fxx, fyy, fzz
 real :: uu, vv, ww, w
 real :: u1, u2, u3, u4
 real :: r, xleap, yleap, zleap, sigmaprof, zdet, g, rxleap, &
         ryleap, rzleap, crid, d1, d2, dist
 real :: xold, yold, zold, tau, t90, xdet, ydet

! u,v,w interpolated on interior rho points
 real :: ucomp(udims(1), vdims(2),-udims(3)+1:0)
 real :: vcomp(udims(1), vdims(2),-vdims(3)+1:0)
 real :: wcomp(udims(1), vdims(2),-wdims(3)+1:0)
 real :: uw1(-udims(3)+1:0), uw2(-udims(3)+1:0)

! akt interpolated on interior rho points
 real :: aktcomp(udims(1), vdims(2),-wdims(3)+1:0)
 
 integer :: res_u=0, res_v=0, res_nsources=0, nsources_loc=0, id_thread=0, num_threads=0, udims_loc=0, vdims_loc=0, res_i=0
 integer :: npart_loc=0, res_npart=0, res_iint=0, iint_loc=0, N=0
 real :: start1=0, end1=0, time_tot=0;
 
 start1=omp_get_wtime() 

 call OMP_SET_NUM_THREADS(ninterp)
 

!$OMP PARALLEL PRIVATE (num_threads, i, j, id_thread, k, u1, u2, u3, uw1, uw2) SHARED( mask_rho, wdims, vdims, udims) 

id_thread=OMP_GET_THREAD_NUM()
num_threads=OMP_GET_NUM_THREADS()
!vdims_loc=vdims(2)/num_threads
!udims_loc=udims(1)/num_threads
!res_v=MOD(vdims(2),num_threads)
!res_u=MOD(udims(1),num_threads)
!write(*,*) "vdims(2): ", vdims(2)
!write(*,*) "udims(1): ", udims(1)

!write(*,*) "id_thread FUORI: ", id_thread
!write(*,*) "num_threads: ", num_threads

if(id_thread == 0) then
  
!write(*,*) "ESEGUO id_thread : ", id_thread 
!write(*,*) "lavoro su u e sono il thread: ", id_thread

!write(*,*) "Thread 0 vdims_loc: ",vdims_loc
!write(*,*) "Thread 0 res_v: ", res_v
!write(*,*) "Thread 0 udims(1): ", udims(1)

 
do j=1, vdims(2)
   do i=1, udims(1)
   
      if ( mask_rho(i,j) > 0.0 ) then
      if ( mask_u(i,j) > 0.0 ) then
      
      uw1=uwacom(i,j,:)
      
      else
      
      uw1=0.0
      
      endif
      if ( mask_u(i,j-1) > 0.0 ) then
    
      uw2=uwacom(i,j-1,:)
      
      else
      
      uw2=0.0
      
      endif
      
       ucomp(i,j,:)=0.5*(uw1+uw2)
       
      else      
        
       ucomp(i,j,:)=0.0
       
     endif
    ! write(*,*) "Thread 0: ucomp= ", ucomp(i,j,:)
    enddo
  enddo
  

endif 

  
 ! write(*,*) "Thread 1 vdims_loc: ",vdims_loc
 ! write(*,*) "Thread 1 res_v: ", res_v
 ! write(*,*) "Thread 1 udims(1): ", udims(1)
 
 !do j=(vdims_loc*id_thread)+res_v+1, (vdims_loc*(id_thread+1))+res_v
  !do i=1, udims(1)
  
  !if ( mask_rho(i,j) > 0.0 ) then
    !if ( mask_u(i,j) > 0.0 ) then
      
          
         ! uw1=uwacom(i,j,:)
          
          !else
          
         ! uw1=0.0
          
          !endif
          !if ( mask_u(i,j-1) > 0.0 ) then
          
         ! uw2=uwacom(i,j-1,:)
          
          !else
          
          !uw2=0.0
        
          !endif
          
          !ucomp(i,j,:)=0.5*(uw1+uw2)
          
         !else      
           
          !ucomp(i,j,:)=0.0
          
        !endif
        !write(*,*) "Thread 1: ucomp= ", ucomp(i,j,:)
       !enddo
      !enddo
       
!endif
  
  

 ! write(*,*) "lavoro su v e sono il thread: ", id_thread
  
if (id_thread == 1) then
!write(*,*) "ESEGUO id_thread : ", id_thread
!interpolate v on rho points
 do j=1, vdims(2)
   do i=1, udims(1)
      if ( mask_rho(i,j) > 0.0 ) then
       if ( mask_v(i,j) > 0.0 ) then
        uw1=vwacom(i,j,:)
      else
      uw1=0.0
      endif
      if ( mask_v(i-1,j) > 0.0 ) then
      uw2=vwacom(i-1,j,:)
       else
       uw2=0.0
     endif
      vcomp(i,j,:)=0.5*(uw1+uw2)
      else
      vcomp(i,j,:)=0.0
      endif
   enddo
  enddo

endif

!else
  !do j=(id_thread*vdims_loc)+res_v+1, (vdims_loc*(id_thread+1))+res_v
  ! do i=(id_thread*udims_loc)+res_u+1, (udims_loc*(id_thread+1))+res_u
    !if(mask_rho(i,j) > 0.0) then
    !if(mask_v(i,j) > 0.0) then
    !uw1=vwacom(i,j,:);
    !else
    !uw1=0;
    !endif
    !if(mask_v(i-1,j) > 0.0) then
    !uw2=vwacom(i-1,j,:)
    !else
    !uw2=0.0
    !endif
    !vcomp(i,j,:)=0.5*(uw1+uw2)
    !else
    !vcomp(i,j,:)=0.0
   !endif
  !enddo
 !enddo

 

 if (id_thread == 2) then
!interpolate w on rho points
!write(*,*) "ESEGUO id_thread : ", id_thread
 

 
 	do k=-wdims(3)+1,0
 		do j=1, vdims(2)
 			do i=1, udims(1)
 				if ( mask_rho(i,j) > 0.0 ) then
 					if ( mask_rho(i-1,j) > 0.0 ) then
 						u1=wwacom(i-1,j,k)
 					else
 						u1=0.0
 					endif
 					if ( mask_rho(i,j-1) > 0.0 ) then
 						u2=wwacom(i,j-1,k)
 					else
 						u2=0.0
 					endif
 					if ( mask_rho(i-1,j-1) > 0.0 ) then
 						u3=wwacom(i-1,j-1,k)
 					else
 						u3=0.0
 					endif
 					wcomp(i,j,k)=0.25*(u1+u2+u3+wwacom(i,j,k))
 				else
     					wcomp(i,j,k)=0.0
    				endif
   			enddo
  		enddo
 	enddo
 
      endif


      if (id_thread == 3) then
        
! interpolate akt on rho points
 !write(*,*) "sono in akt e sono il thread: ", id_thread 
 
 
 	do k=-wdims(3)+1,0
 		do j=1, vdims(2)
 			do i=1, udims(1)
 				if ( mask_rho(i,j) > 0.0 ) then
 					if ( mask_rho(i-1,j) > 0.0 ) then
 						u1=akt(i-1,j,k)
 					else
 						u1=0.0
 					endif
 					if ( mask_rho(i,j-1) > 0.0 ) then
 						u2=akt(i,j-1,k)
 					else
 						u2=0.0
 					endif
 					if ( mask_rho(i-1,j-1) > 0.0 ) then
 						u3=akt(i-1,j-1,k)
 					else
 						u3=0.0
 					endif
 					aktcomp(i,j,k)=0.25*(u1+u2+u3+akt(i,j,k))
 				else
 					aktcomp(i,j,k)=0.0
 				endif
 			enddo
 		enddo
 	enddo
      
        endif
  !$OMP END PARALLEL
 
        



! initialize level depths
 do i=-wdims(3)+2, 0
 	depth(i)=s_w(i)-s_w(i-1)
 enddo

! initialite sed. velocity
 vs=0.0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! seed initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call date_and_time(values=time)     ! get the current time
 seed(1)=time(4)*(360000*time(5)+6000*time(6)+100*time(7)+time(8))
 call random_seed(put=seed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          dispersion                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(*,*) "Dispersion" 
 iint=deltat/dti ! time steps per hour

 call OMP_SET_NUM_THREADS(ncalc)


 !$OMP PARALLEL PRIVATE (i, ip, id_thread, num_threads, nsources_loc, res_nsources) SHARED( nsources, npart, npartsperhour) 
 
 num_threads=OMP_GET_NUM_THREADS()
 id_thread=OMP_GET_THREAD_NUM()
 nsources_loc=nsources/num_threads
 res_nsources=MOD(nsources, num_threads)
 
  
 
 if(id_thread==0) then
 
 do i=1, nsources_loc+res_nsources
      do ip=npart+1, npart+npartsperhour(i)
 	      health(ip)=health0
 		if ( mode(i) == 1 ) then
 		      call random_number(rn)
 			xpart(ip)=i_source(i) + rn*0.5 - 0.25
                          !xpart(ip)=i_source(i) + 0.5*0.5 - 0.25
 				call random_number(rn)
 				ypart(ip)=j_source(i) + rn*0.5 - 0.25
                                !ypart(ip)=j_source(i) + 0.5*0.5 - 0.25
 				call random_number(harvest = rn)
 				zpart(ip)=k_source(i) + rn*0.5 - 0.25
                                !zpart(ip)=k_source(i) + 0.5*0.5 - 0.25
 			else if ( mode(i) <= 0 ) then
 				health(ip)=0.0
 			else
 				print*, "release mode not defined"
 				print*, "mode: ", mode
 				stop 77777
 			endif
 		enddo
 		do ll=0, iint-1
 			do ip=npart+1+ll*npartsperhour(i)/iint, npart+(ll+1)*npartsperhour(i)/iint
 				tpart(ip)=ll*dti
 			enddo
 		enddo
 		npart=npart+npartsperhour(i)
                enddo
              
  
              else


  do i=(nsources_loc*id_thread)+res_nsources+1, (nsources_loc*(id_thread+1))+res_nsources
      do ip=npart+1, npart+npartsperhour(i)
          health(ip)=health0
          if ( mode(i) == 1 ) then
            call random_number(rn)
            xpart(ip)=i_source(i) + rn*0.5 - 0.25
            !xpart(ip)=i_source(i) + 0.5*0.5 - 0.25  
            call random_number(rn)
            ypart(ip)=j_source(i) + rn*0.5 - 0.25
            !ypart(ip)=j_source(i) + 0.5*0.5 - 0.25  
            call random_number(rn)
            zpart(ip)=k_source(i) + rn*0.5 - 0.25
            !zpart(ip)=k_source(i) + 0.5*0.5 - 0.25  
          else if ( mode(i) <= 0 ) then  
              health(ip)=0.0
            else
              print*, "release mode not defined"
              print*, "mode: ", mode
              stop 77777
            endif
        enddo
       do ll=0, iint-1
        do ip=npart+1+ll*npartsperhour(i)/iint, npart+(ll+1)*npartsperhour(i)/iint
         tpart(ip)=ll*dti
         enddo
        enddo
        npart=npart+npartsperhour(i)
        enddo
              
              
    endif
  !$OMP END PARALLEL

  
  
  !call OMP_SET_NUM_THREADS(N)
  
   
! Time steps

 !OMP PARALLEL PRIVATE(ll, res_iint, iint_loc, ip, id_thread, num_threads), &
 !OMP& SHARED(wdims, vdims, udims), &
 !OMP& SHARED(iint, npart, ixx, iyy, zdet, sigmaprof), &
 !OMP& SHARED(health, pstatus, xpart, ypart, zpart, mask_rho), &
 !OMP& SHARED(ucomp, vcomp, wcomp, aktcomp, lat_v, lon_u) &
 !OMP& SHARED(depth, zeta, tpart)
 !OMP& SHARED(g, health0, tau0, u4), &
 !OMP& SHARED(xold, yold, zold, ixx, iyy, izz, fxx, fyy, fzz, u1, u2, u3), &
 !OMP& SHARED( uu, vv, ww, xpart, zpart, ypart, tpart), &
 !OMP& SHARED(d1, d2, dist, xdet, ydet, zdet), &
 !OMP& SHARED(wcomp, aktcomp, lat_v, lon_u, depth, zeta), &
 !OMP& SHAREDxleap, yleap, zleap, rxleap, ryleap, rzleap), &
 !OMP& SHARED(survprob, idead, iouter, dti, sigmaprof, mask_rho, ucomp)

 num_threads=OMP_GET_NUM_THREADS()
 id_thread=OMP_GET_THREAD_NUM()
 iint_loc=iint/num_threads
 res_iint=MOD(iint,num_threads)

 !write(*,*) "numero threads: ", num_threads
 !write(*,*) "id_thread: ", id_thread
 !write(*,*) "iint: ", iint
 !write(*,*) "iint_loc: ", iint_loc
 !write(*,*) "res_iint: ", res_iint

 !if(id_thread == 0) then
 
   
 do ll=1, iint

!##################### outer cycle over all particles ###################
  
!write(*,*) "outer cycle over all particles:", ll
 
 

 	do ip=1,npart
 		if (health(ip) <= survprob) then  ! dead particle
 			if ( pstatus(ip) > 0 ) idead=idead+1
 			pstatus(ip)=0
 			cycle
 		endif
    
 		xold=xpart(ip)
 		yold=ypart(ip)
 		zold=zpart(ip)

 		ixx=xold     ! integer part
 		fxx=xold-ixx ! remainder

 		iyy=yold
 		fyy=yold-iyy

 		izz=zold
 		fzz=abs(zold-izz)
                
 		if ( (ixx < 1) .or. (ixx >= udims(1)) ) then
 			health(ip)=-1.0
 			iouter=iouter+1
 			pstatus(ip)=0
 			cycle
 		endif
                
 		if ( (iyy < 1) .or. (iyy >= vdims(2)) ) then
 			health(ip)=-1.0
 			iouter=iouter+1
 			pstatus(ip)=0
 			cycle
 		endif
                
	 	if ( mask_rho(ixx,iyy) <= 0.0 ) then
 			health(ip)=-1.0
 			iouter=iouter+1
 			pstatus(ip)=0
 			cycle
 		endif
 !OMP CRITICAL
 		u1=ucomp(ixx  ,iyy  ,izz)*(1.0-fxx)*(1.0-fyy)
 		u2=ucomp(ixx  ,iyy+1,izz)*(1.0-fxx)*     fyy
 		u3=ucomp(ixx+1,iyy+1,izz)*     fxx *     fyy
 		u4=ucomp(ixx+1,iyy  ,izz)*     fxx *(1.0-fyy)
 		uu=u1+u2+u3+u4

 		u1=vcomp(ixx  ,iyy  ,izz)*(1.0-fxx)*(1.0-fyy)
 		u2=vcomp(ixx  ,iyy+1,izz)*(1.0-fxx)*     fyy
 		u3=vcomp(ixx+1,iyy+1,izz)*     fxx *     fyy
 		u4=vcomp(ixx+1,iyy  ,izz)*     fxx *(1.0-fyy)
 		vv=u1+u2+u3+u4

 		u1=wcomp(ixx  ,iyy  ,izz  )*(1.0-fxx)*(1.0-fyy)*(1.0-fzz)
 		u2=wcomp(ixx  ,iyy+1,izz  )*(1.0-fxx)*     fyy *(1.0-fzz)
 		u3=wcomp(ixx+1,iyy+1,izz  )*     fxx *     fyy *(1.0-fzz)
 		u4=wcomp(ixx+1,iyy  ,izz  )*     fxx *(1.0-fyy)*(1.0-fzz)
 		ww=u1+u2+u3+u4

 		u1=wcomp(ixx  ,iyy  ,izz-1)*(1.0-fxx)*(1.0-fyy)*     fzz
 		u2=wcomp(ixx  ,iyy+1,izz-1)*(1.0-fxx)*     fyy *     fzz
 		u3=wcomp(ixx+1,iyy+1,izz-1)*     fxx *     fyy *     fzz
 		u4=wcomp(ixx+1,iyy  ,izz-1)*     fxx *(1.0-fyy)*     fzz
 		ww=ww+u1+u2+u3+u4
 
 		if ( abs(uu) > 10.0 .or. abs(vv) > 10.0 .or. abs(ww) > 10.0 ) then
 			write(*,*) " warning! speeds: ", uu, vv, ww
 		endif

!-----particle leap (deterministic component)

 		xleap=uu*dti
 		yleap=vv*dti
 		zleap=(vs+ww)*dti

! stochastic counterpart

! sigma profile
                
		sigmaprof=3.46*(1+zdet/wdims(3))
  
!OMP END CRITICAL


! pseudo-random numbers
!write(*,*) "pseudo-random numbers"
 		
 		g=0.
 		do i=1,12
 			call random_number(r)
 			g=g+r-0.5
                        !g=g+0.5-0.5
 		enddo
 		rxleap=g*sigmaprof
 		g=0.
 		do i=1,12
 			call random_number(r)
 			g=g+r-0.5
                        !g=g+0.5-0.5
 		enddo
 		ryleap=g*sigmaprof
 		g=0.
 		do i=1,12
 			call random_number(r)
 			g=g+r-0.5
                        !g=g+0.5-0.5
 		enddo
 !OMP CRITICAL
 		u1=aktcomp(ixx  ,iyy  ,izz  )*(1.0-fxx)*(1.0-fyy)*(1.0-fzz)
 		u2=aktcomp(ixx  ,iyy+1,izz  )*(1.0-fxx)*     fyy *(1.0-fzz)
 		u3=aktcomp(ixx+1,iyy+1,izz  )*     fxx *     fyy *(1.0-fzz)
 		u4=aktcomp(ixx+1,iyy  ,izz  )*     fxx *(1.0-fyy)*(1.0-fzz)
 		ww=u1+u2+u3+u4

 		u1=aktcomp(ixx  ,iyy  ,izz-1)*(1.0-fxx)*(1.0-fyy)*     fzz
 		u2=aktcomp(ixx  ,iyy+1,izz-1)*(1.0-fxx)*     fyy *     fzz
 		u3=aktcomp(ixx+1,iyy+1,izz-1)*     fxx *     fyy *     fzz
 		u4=aktcomp(ixx+1,iyy  ,izz-1)*     fxx *(1.0-fyy)*     fzz
 		ww=ww+u1+u2+u3+u4
 		rzleap=g*ww*crid ! reduction coefficient

!-----move particle
!write(*,*) "move particle"
 		xleap=xleap+rxleap ! add random displacement
 		yleap=yleap+ryleap
 		zleap=zleap+rzleap

 		d1=(lat_v(ixx,iyy+1)-lat_v(ixx,iyy))
 		d2=(lon_u(ixx+1,iyy)-lon_u(ixx,iyy))
 		d1=sin(0.5*d1)**2 + (sin(0.5*d2)**2)* &
 		cos(lat_v(ixx,iyy+1))*cos(lat_v(ixx,iyy))
 		dist=2.0*atan2(sqrt(d1),sqrt(1.0-d1))
 		dist=6371.0*dist
 		xdet=xold+0.001*xleap/dist

 		d1=(lat_v(ixx,iyy+1)-lat_v(ixx,iyy))
 		d2=(lon_u(ixx+1,iyy)-lon_u(ixx,iyy))
 		d1=sin(0.5*d1)**2 + (sin(0.5*d2)**2)* &
 		cos(lat_v(ixx,iyy+1))*cos(lat_v(ixx,iyy))
 		dist=2.0*atan2(sqrt(d1),sqrt(1.0-d1))
 		dist=6371.0*dist
 		ydet=yold+0.001*yleap/dist

 		dist=depth(izz)*zeta(ixx,iyy)
 		if ( abs(zleap) .GT. abs(dist) ) zleap=sign(dist, zleap)
 		zdet=zold+zleap/dist

! reflect if out-of-column
!write(*,*) "reflect if out-of-column"
        
              
 		if ( zdet < -wdims(3)+2 ) zdet=2.0*(-wdims(3)+2)-zdet
 		if ( zdet > 0. ) zdet=-zdet
                

! reflect if crossed the coastline
!write(*,*) "reflect if crossed the coastline"
 		if ( mask_rho(int(xdet),int(ydet)) <= 0.0 ) then
 			if ( int(xdet) < ixx ) then
 				xdet=real(ixx) + abs(xold-xdet)
 			elseif ( int(xdet) > ixx ) then
 				xdet=real(int(xdet)) - mod(xdet,1.0)
 			endif
 			if ( int(ydet) < iyy ) then
 				ydet=real(iyy) + abs(yold-ydet)
 			elseif ( int(ydet) > iyy ) then
 				ydet=real(int(ydet)) - mod(ydet,1.0)
 			endif
 		endif
! update particle position
!write(*,*) "update particle position"
 		xpart(ip)=xdet
 		ypart(ip)=ydet
 		zpart(ip)=zdet

 		tpart(ip)=tpart(ip)+dti
 		health(ip)=health0*exp(-tpart(ip)/tau0)
                
    !OMP END CRITICAL

 	enddo ! cycle over particles
!write(*,*) "end cycle over particles"


 enddo ! cycle over time steps
 
 
 !else
  !OMP CRITICAL
   !do ll=(id_thread*iint_loc)+res_iint+1, (iint_loc*(id_thread+1))+res_iint
   
  ! do ip=1,npart
   !if (health(ip) <= survprob) then
    ! if ( pstatus(ip) > 0 ) idead=idead+1
     !pstatus(ip)=0
     !cycle
   !endif

   !xold=xpart(ip)
   !yold=ypart(ip)
   !zold=zpart(ip)
   !ixx=xold 
   !fxx=xold-ixx

   !iyy=yold
   !fyy=yold-iyy

   !izz=zold
   !fzz=abs(zold-izz)

   
   !if ( (ixx < 1) .or. (ixx >= udims(1)) ) then
    ! health(ip)=-1.0
     !iouter=iouter+1
     !pstatus(ip)=0  
     !cycle
  ! endif
   

  !if ( (iyy < 1) .or. (iyy >= vdims(2)) ) then
   ! health(ip)=-1.0
    !iouter=iouter+1
    !pstatus(ip)=0
     ! cycle
    !endif

   !if ( mask_rho(ixx,iyy) <= 0.0 ) then
    ! health(ip)=-1.0
     !iouter=iouter+1
     !pstatus(ip)=0
     !cycle
   !endif

   !u1=ucomp(ixx  ,iyy  ,izz)*(1.0-fxx)*(1.0-fyy)
   !u2=ucomp(ixx  ,iyy+1,izz)*(1.0-fxx)*     fyy
   !u3=ucomp(ixx+1,iyy+1,izz)*     fxx *     fyy
   !u4=ucomp(ixx+1,iyy  ,izz)*     fxx *(1.0-fyy)
   !uu=u1+u2+u3+u4

   !u1=vcomp(ixx  ,iyy  ,izz)*(1.0-fxx)*(1.0-fyy)
   !u2=vcomp(ixx  ,iyy+1,izz)*(1.0-fxx)*     fyy
   !u3=vcomp(ixx+1,iyy+1,izz)*     fxx *     fyy
   !u4=vcomp(ixx+1,iyy  ,izz)*     fxx *(1.0-fyy)
   !vv=u1+u2+u3+u4

   !u1=wcomp(ixx  ,iyy  ,izz  )*(1.0-fxx)*(1.0-fyy)*(1.0-fzz)
   !u2=wcomp(ixx  ,iyy+1,izz  )*(1.0-fxx)*     fyy *(1.0-fzz)
   !u3=wcomp(ixx+1,iyy+1,izz  )*     fxx *     fyy *(1.0-fzz)
   !u4=wcomp(ixx+1,iyy  ,izz  )*     fxx *(1.0-fyy)*(1.0-fzz)
   !ww=u1+u2+u3+u4

   !u1=wcomp(ixx  ,iyy  ,izz-1)*(1.0-fxx)*(1.0-fyy)*     fzz
   !u2=wcomp(ixx  ,iyy+1,izz-1)*(1.0-fxx)*     fyy *     fzz
   !u3=wcomp(ixx+1,iyy+1,izz-1)*     fxx *     fyy *     fzz
   !u4=wcomp(ixx+1,iyy  ,izz-1)*     fxx *(1.0-fyy)*     fzz
   !ww=ww+u1+u2+u3+u4

   !if ( abs(uu) > 10.0 .or. abs(vv) > 10.0 .or. abs(ww) > 10.0 ) then
    ! write(*,*) " warning! speeds: ", uu, vv, ww
   !endif
   
   !-----particle leap (deterministic component)

   !xleap=uu*dti
   !yleap=vv*dti
   !zleap=(vs+ww)*dti

   ! stochastic counterpart

   ! sigma profile
   
   
   !sigmaprof=3.46*(1+zdet/wdims(3))
   

   ! pseudo-random numbers
   !write(*,*) "pseudo-random numbers"

   !g=0.
   !do i=1,12
   !call random_number(r)
   !g=g+r-0.5
   !g=g+0.5-0.5
   !enddo
   !rxleap=g*sigmaprof

   !g=0.
   !do i=1,12
   !call random_number(r)
   !g=g+r-0.5
   !g=g+0.5-0.5
   !enddo
   !ryleap=g*sigmaprof

   !g=0.
   !do i=1,12
   !call random_number(r)
   !g=g+r-0.5
   !g=g+0.5-0.5
   !enddo

   !u1=aktcomp(ixx  ,iyy  ,izz  )*(1.0-fxx)*(1.0-fyy)*(1.0-fzz)
   !u2=aktcomp(ixx  ,iyy+1,izz  )*(1.0-fxx)*     fyy *(1.0-fzz)
   !u3=aktcomp(ixx+1,iyy+1,izz  )*     fxx *     fyy *(1.0-fzz)
   !u4=aktcomp(ixx+1,iyy  ,izz  )*     fxx *(1.0-fyy)*(1.0-fzz)
   !ww=u1+u2+u3+u4

   !u1=aktcomp(ixx  ,iyy  ,izz-1)*(1.0-fxx)*(1.0-fyy)*     fzz
   !u2=aktcomp(ixx  ,iyy+1,izz-1)*(1.0-fxx)*     fyy *     fzz
   !u3=aktcomp(ixx+1,iyy+1,izz-1)*     fxx *     fyy *     fzz
   !u4=aktcomp(ixx+1,iyy  ,izz-1)*     fxx *(1.0-fyy)*     fzz
   !ww=ww+u1+u2+u3+u4

   !rzleap=g*ww*crid ! reduction coefficient

   !-----move particle
   !write(*,*) "move particle"
   !xleap=xleap+rxleap ! add random displacement
   !yleap=yleap+ryleap
   !zleap=zleap+rzleap

   !d1=(lat_v(ixx,iyy+1)-lat_v(ixx,iyy))
   !d2=(lon_u(ixx+1,iyy)-lon_u(ixx,iyy))
   !d1=sin(0.5*d1)**2 + (sin(0.5*d2)**2)* &
   !cos(lat_v(ixx,iyy+1))*cos(lat_v(ixx,iyy))
   !dist=2.0*atan2(sqrt(d1),sqrt(1.0-d1))
   !dist=6371.0*dist
   !xdet=xold+0.001*xleap/dist

   !d1=(lat_v(ixx,iyy+1)-lat_v(ixx,iyy))
   !d2=(lon_u(ixx+1,iyy)-lon_u(ixx,iyy))
   !d1=sin(0.5*d1)**2 + (sin(0.5*d2)**2)* &
   !cos(lat_v(ixx,iyy+1))*cos(lat_v(ixx,iyy))
   !dist=2.0*atan2(sqrt(d1),sqrt(1.0-d1))
   !dist=6371.0*dist
   !ydet=yold+0.001*yleap/dist
   !dist=depth(izz)*zeta(ixx,iyy)
   !if ( abs(zleap) .GT. abs(dist) ) zleap=sign(dist, zleap)
   !zdet=zold+zleap/dist

   ! reflect if out-of-column
   !write(*,*) "reflect if out-of-column"

    
    !if ( zdet < -wdims(3)+2 ) zdet=2.0*(-wdims(3)+2)-zdet
    !if ( zdet > 0. ) zdet=-zdet
   

   !if ( mask_rho(int(xdet),int(ydet)) <= 0.0 ) then
    ! if ( int(xdet) < ixx ) then
     !  xdet=real(ixx) + abs(xold-xdet)
     !elseif ( int(xdet) > ixx ) then
     !xdet=real(int(xdet)) - mod(xdet,1.0)
    !endif
    !if ( int(ydet) < iyy ) then
     ! ydet=real(iyy) + abs(yold-ydet)
    !elseif ( int(ydet) > iyy ) then
    !ydet=real(int(ydet)) - mod(ydet,1.0)
  !endif
!endif
! update particle position
!write(*,*) "update particle position"

!xpart(ip)=xdet
!ypart(ip)=ydet
!zpart(ip)=zdet

!tpart(ip)=tpart(ip)+dti
!health(ip)=health0*exp(-tpart(ip)/tau0)

!enddo ! cycle over particles
!write(*,*) "end cycle over particles"

 
!enddo ! cycle over time steps
!OMP END CRITICAL

! endif!************************************************************
 

 !OMP END PARALLEL
 
!write(*,*) "end cycle over time steps"

!-----zero-th the concentration matrix
 do k=-wdims(3)+2,0
 	do j=1,vdims(2)
 		do i=1,udims(1)
 			conc(i,j,k)=0.
 		enddo
 	enddo
 enddo
!-----update the concentration matrix
 idead=0
 iouter=0

 call OMP_SET_NUM_THREADS(ncalc)

 !$OMP PARALLEL PRIVATE(id_thread, num_threads, ip, npart_loc, res_npart, ixx, iyy, izz) SHARED(npart)
 
 num_threads=OMP_GET_NUM_THREADS()
 id_thread=OMP_GET_THREAD_NUM()
 npart_loc=npart/num_threads
 res_npart=MOD(npart,num_threads)
 
 if(id_thread == 0) then

 do ip=1,npart_loc+res_npart
 	if(health(ip) > survprob) then
 		ixx=xpart(ip)
 		iyy=ypart(ip)
 		izz=zpart(ip)
 		if ( mask_rho(ixx,iyy) > 0.0 ) then
 			conc(ixx,iyy,izz)=conc(ixx,iyy,izz)+1.0
 		else
 			health(ip)=-2.0
 			idead=idead+1
 			iouter=iouter+1
 		endif
 	else
 		health(ip)=-2.0
 		idead=idead+1
 	endif
 enddo
 
 else
 do ip=(npart_loc*id_thread)+res_npart+1 , (npart_loc*(id_thread+1))+res_npart
    if(health(ip) > survprob) then
      ixx=xpart(ip)
      iyy=ypart(ip)
      izz=zpart(ip)
      if ( mask_rho(ixx,iyy) > 0.0 ) then
        conc(ixx,iyy,izz)=conc(ixx,iyy,izz)+1.0
      else
        health(ip)=-2.0
        idead=idead+1
        iouter=iouter+1
      endif
    else
      health(ip)=-2.0
      idead=idead+1
    endif
    enddo

 endif 


 !$OMP END PARALLEL

 end1=omp_get_wtime()
 time_tot=time_tot+(end1-start1)
 write(*,*) "Tempo trascorso: ", time_tot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           print section                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 write(*,*)
 write(*,'(a24,i6)')'total # of parts ->',npart
 write(*,'(a24,i6)')'dead parts       ->',idead
 write(*,'(a24,i6)')'outer parts      ->',iouter

! write(*,*)'conc             ->',conc(98,337,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 return 
 end subroutine wacomm
