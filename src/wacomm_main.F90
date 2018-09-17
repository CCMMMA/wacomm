      PROGRAM WACOM_main

      use netcdf
      use comblk
      use omp_lib


      implicit none

      integer clock
      
      INTEGER TID
      
    

      integer seed
      allocatable seed(:)

       
! namelist vars
      integer ninterp, ncalc
      character*200 nc_input, nc_output_root, restartfile,              &
                    historyfile, inithsttime, rstfiledate, namelistfile
      integer io_form_wacomm

      namelist /io/ nc_input, nc_output_root, io_form_wacomm
      namelist /ems/ nsources, i_source, j_source, k_source,            &
                     id_source, npartsperhour, mode,                    &
                     source_start, source_end
      namelist /chm/ survprob, tau0
      namelist /rst/ restart, restartfile, interval
      namelist /hst/ history, historyfile, outfreq, inithsttime
      namelist /omp/ ninterp, ncalc

      integer dims
      allocatable dims(:)

      integer ncid, dimid, varid, xtype, nshp
      integer  ndims, nvars, natts, id_unlim_dim
      character*132 dim_name, var_name, att_name
      integer diml(100)

      real u, v, w, akt
      allocatable u(:,:,:,:)
      allocatable v(:,:,:,:)
      allocatable w(:,:,:,:)
      allocatable akt(:,:,:,:)

      real*8 mask_u, mask_v, mask_rho
      allocatable mask_u(:,:), mask_v(:,:), mask_rho(:,:)

      real*8 lat_rho, lon_rho
      allocatable lat_rho(:,:)
      allocatable lon_rho(:,:)
      real*8 lat_v, lon_u
      allocatable lat_v(:,:)
      allocatable lon_u(:,:)

      real*8 ocean_time
      allocatable ocean_time(:)

      real conc
      allocatable conc(:,:,:)

      real*8 s_rho, s_w
      allocatable s_rho(:)
      allocatable s_w(:)

      real zeta
      real*8 bath
      allocatable zeta(:,:,:), bath(:,:)

      integer udims(4), vdims(4), wdims(4), aktdims(4)

      real dti, deltat

      logical restart, history
      integer nrstpart
      real outfreq, interval
      integer syear, smon, sday, shour, rstcounter
      double precision jdref, jd

      integer ncmode

      real pi
      integer i, n, k


      call getarg(1, namelistfile)
      print*, "Open namelist ", trim(namelistfile)
      open(8,file=trim(namelistfile))
      read(8,nml=io)
      read(8,nml=ems)
      read(8,nml=chm)
      read(8,nml=rst)
      read(8,nml=hst)
      close(8)

      if ( restart ) then
       write(*,*)
       write(*,*) "reading restart from file ", trim(restartfile)
       open(unit=20,file=trim(restartfile),status='old',action='read')
       read(20,*) nrstpart
      else
       nrstpart=0
      endif

      read(inithsttime,60) syear, smon, sday, shour
      call greg2jd(syear, smon, sday, i)
      jdref=real(i) + shour/24.0
60    format(i4,i2,i2,1x,i2)

      print*, "open file: ",nc_input
      call check( nf90_open(trim(nc_input), nf90_nowrite, ncid) )
      print*, "inquire about dimension of ", trim(nc_input),            &
        " file..."
      call check( nf90_inquire(ncid, ndims, nvars, natts,               &
                               id_unlim_dim) )
      write(*,10) "ndims=",ndims, "nvars=",nvars, "ngatts=",natts,      &
                  "id_unlimited_dimension=",id_unlim_dim
 10   format(4(a10,i3))

      allocate(dims(ndims))

      print*, "*** dimensions:"
      do dimid = 1, ndims
         call check( nf90_inquire_dimension(ncid, dimid, dim_name,      &
                     diml(dimid)) )
         write(*,20) 'dimid', dimid, 'dim_name=', dim_name, 'dimlen=',  &
                     diml(dimid)
      end do
20    format(a10,i3,a10,a30,a10,i3)

      print*, "*** variables: "
      do varid = 1, nvars
       call check ( nf90_inquire_variable(ncid, varid, var_name,        &
                                          xtype, nshp, dims) )

       call tolower(var_name)

       if ( trim(var_name) == 'u' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 4 ) allocate(u(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 3 ) allocate(u(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 2 ) allocate(u(diml(dims(1)),diml(dims(2)),1,1))
        if ( nshp == 1 ) allocate(u(diml(dims(1)),1,1,1))
        if ( nshp == 0 ) allocate(u(1,1,1,1))
        write(*,*) "dims:", diml(dims(1:4))
        udims=diml(dims(1:4))
       endif
       if ( var_name == 'v' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 4 ) allocate(v(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 3 ) allocate(v(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 2 ) allocate(v(diml(dims(1)),diml(dims(2)),1,1))
        if ( nshp == 1 ) allocate(v(diml(dims(1)),1,1,1))
        if ( nshp == 0 ) allocate(v(1,1,1,1))
        write(*,*) "dims:", diml(dims(1:4))
        vdims=diml(dims(1:4))
       endif
       if ( var_name == 'w' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 4 ) allocate(w(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 3 ) allocate(w(diml(dims(1)),diml(dims(2)),        &
                                    diml(dims(3)),1))
        if ( nshp == 2 ) allocate(w(diml(dims(1)),diml(dims(2)),1,1))
        if ( nshp == 1 ) allocate(w(diml(dims(1)),1,1,1))
        if ( nshp == 0 ) allocate(w(1,1,1,1))
        write(*,*) "dims:", diml(dims(1:4))
        wdims=diml(dims(1:4))
       endif

       if ( var_name == 'akt' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 4 ) allocate(akt(diml(dims(1)),diml(dims(2)),      &
                                      diml(dims(3)),1))
        if ( nshp == 3 ) allocate(akt(diml(dims(1)),diml(dims(2)),      &
                                      diml(dims(3)),1))
        if ( nshp == 2 ) allocate(akt(diml(dims(1)),diml(dims(2)),1,1))
        if ( nshp == 1 ) allocate(akt(diml(dims(1)),1,1,1))
        if ( nshp == 0 ) allocate(akt(1,1,1,1))
        write(*,*) "dims:", diml(dims(1:4))
        aktdims=diml(dims(1:4))
       endif

       if ( var_name == 'mask_u' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(mask_u(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, mask_u) )
        write(*,*) "dims:", diml(dims(1:2))
       endif
       if ( var_name == 'mask_v' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(mask_v(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, mask_v) )
        write(*,*) "dims:", diml(dims(1:2))
       endif
       if ( var_name == 'mask_rho' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 )                                                &
         allocate(mask_rho(0:diml(dims(1))-1,0:diml(dims(2))-1))
        call check( nf90_get_var(ncid, varid, mask_rho) )
        write(*,*) "dims:", diml(dims(1:2))
       endif

       if ( var_name == 'lat_v' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(lat_v(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, lat_v) )
        write(*,*) "dims:", diml(dims(1:2))
       endif

       if ( var_name == 'lon_u' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(lon_u(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, lon_u) )
        write(*,*) "dims:", diml(dims(1:2))
       endif

       if ( var_name == 'lat_rho' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(lat_rho(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, lat_rho) )
        write(*,*) "dims:", diml(dims(1:2))
       endif
       if ( var_name == 'lon_rho' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(lon_rho(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, lon_rho) )
        write(*,*) "dims:", diml(dims(1:2))
       endif

       if ( var_name == 'zeta' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 3 )                                                &
         allocate(zeta(0:diml(dims(1))-1,0:diml(dims(2))-1,1))
        write(*,*) "dims:", diml(dims(1:3))
       endif
       if ( var_name == 'h' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 2 ) allocate(bath(diml(dims(1)),diml(dims(2))))
        call check( nf90_get_var(ncid, varid, bath) )
        write(*,*) "dims:", diml(dims(1:2))
       endif

       if ( var_name == 'ocean_time' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 1 ) allocate(ocean_time(diml(dims(1))))
        call check( nf90_get_var(ncid, varid, ocean_time) )
        write(*,*) "dims:", diml(dims(1:1))
       endif

       if ( var_name == 's_rho' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 1 ) allocate(s_rho(diml(dims(1))))
        call check( nf90_get_var(ncid, varid, s_rho) )
        write(*,*) "dims:", diml(dims(1:1))
       endif
       if ( var_name == 's_w' ) then
        write(*,40) 'variable', var_name, ' is', nshp, '-d shaped'
        if ( nshp == 1 ) allocate(s_w(diml(dims(1))))
        call check( nf90_get_var(ncid, varid, s_w) )
        write(*,*) "dims:", diml(dims(1:1))
       endif
      enddo
30    format(a10,i3,3a10,i3,a10,10i3)
40    format(a11,1x,a12,a3,i2,a9)


      allocate(conc(0:udims(1),0:vdims(2),-wdims(3)+2:0))

      totpart=size(ocean_time,1)*sum(npartsperhour(1:nsources)) +       &
        nrstpart
      allocate(xpart(totpart))
      allocate(ypart(totpart))
      allocate(zpart(totpart))
      allocate(health(totpart))
      allocate(tpart(totpart))
      allocate(pstatus(totpart))

      write(*,*) "end variable allocation"

! initialize random seed
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
      deallocate(seed)

      dti=30.0

      pi=4.0*atan(1.0)
      lat_rho=lat_rho*pi/180.0
      lon_rho=lon_rho*pi/180.0
      lat_v=lat_v*pi/180.0
      lon_u=lon_u*pi/180.0

      health=0.0
      tpart=0.0
      pstatus=1
      health0=1.0
      npart=0

      deltat=ocean_time(2)-ocean_time(1)
      write(*,*) "using a time step of:", deltat
      do i=1, nsources
       npartsperhour=npartsperhour*nint(3600.0/deltat)
       write(*,*) "using npartsperhour", npartsperhour(i),              &
                  "for source", i
      enddo

      if ( restart ) then
       do i=1, nrstpart
        read(20,*) xpart(i), ypart(i), zpart(i), health(i), tpart(i)
       enddo
       npart=nrstpart
       close(20)
      endif

      ncmode = 1
      call fgennc(ncmode,trim(nc_output_root)//'001.nc',                &
          vdims(1), udims(2), udims(3), udims(4), nsources, 0,          &
          lat_rho, lon_rho, mask_rho, zeta, bath, s_rho, ocean_time,    &
          conc)
      ncmode = 2

      iouter=0
      idead=0
      rstcounter=0

! call wacom
      do i=1, udims(4)
       write(*,*) "reading time level:", i
       call check( nf90_inq_varid(ncid, "u", varid) )
       call check( nf90_get_var(ncid, varid, u, start=(/1,1,1,i/),      &
         count = (/ udims(1), udims(2), udims(3), 1 /)) )
       call check( nf90_inq_varid(ncid, "v", varid) )
       call check( nf90_get_var(ncid, varid, v, start=(/1,1,1,i/),      &
         count = (/ vdims(1), vdims(2), vdims(3), 1 /)) )
       call check( nf90_inq_varid(ncid, "w", varid) )
       call check( nf90_get_var(ncid, varid, w, start=(/1,1,1,i/),      &
         count = (/ wdims(1), wdims(2), wdims(3), 1 /)) )
       call check( nf90_inq_varid(ncid, "zeta", varid) )
       call check( nf90_get_var(ncid, varid, zeta, start=(/1,1,i/),     &
         count = (/ wdims(1), wdims(2), 1 /)) )
       call check( nf90_inq_varid(ncid, "AKt", varid) )
       call check( nf90_get_var(ncid, varid, akt, start=(/1,1,1,i/),    &
         count = (/ wdims(1), wdims(2), wdims(3), 1 /)) )

       if ( i == 1 ) zeta(:,:,1)=zeta(:,:,1)+bath
       call wacomm(udims,vdims,wdims,u,v,w,akt,                          &
                  mask_u,mask_v,mask_rho,conc,                          &
                  zeta,s_w,s_rho,lon_u,lat_v,                           &
                  dti,deltat,ninterp,ncalc)
       if ( i == 1 ) zeta(:,:,1)=zeta(:,:,1)-bath

       call fgennc(ncmode,trim(nc_output_root)//'001.nc',               &
           vdims(1), udims(2), udims(3), udims(4), nsources, i,         &
           lat_rho, lon_rho, mask_rho, zeta(1,1,1), bath, s_rho,        &
           ocean_time(i), conc)

       if ( i < udims(4) .AND. history ) then
        if ( int((ocean_time(i+1)-ocean_time(1))/outfreq) >             &
              rstcounter ) then
         jd=jdref + (ocean_time(i+1)-ocean_time(1))/(3600.0*24.0)
         call jd2greg(int(jd*10.0)/10, syear, smon, sday)
         shour=nint(24.0*(jd-int(jd)))
         rstcounter=rstcounter+1
         write(rstfiledate,70) syear, smon, sday, 'Z', shour
         rstfiledate=trim(historyfile)//trim(rstfiledate)//'.txt'
         call writerestartfile(npart, health, survprob, mask_rho,       &
           xpart, ypart, zpart, tpart, shape(mask_rho), rstfiledate)
        endif
       endif
 70    format(i4,i2.2,i2.2,a1,i2.2)

       
      enddo


      ncmode = 3
      call fgennc(ncmode,trim(nc_output_root)//'001.nc',                &
          vdims(1), udims(2), udims(3), udims(4), nsources, 0,          &
          lat_rho, lon_rho, mask_rho, zeta(1,1,i), bath, s_rho,         &
          ocean_time, conc)


#if DEBUG
       open(unit=66,file="partpos2"//".txt",action="write")
       k=0
       do i=udims(4), udims(4)
        do Ip=1, totpart
         if ( health(Ip) .GT. survprob ) k=k+1
        enddo
       enddo

       write(66,'(I6)') k

       do i=udims(4), udims(4)
        do Ip=1, totpart
         if ( health(Ip) .GT. survprob ) then
          write(66,'(3(1x,f14.8))')                                     &
             xpart(Ip)*pi/180.0, ypart(Ip)*pi/180.0, zpart(Ip)
         endif
        enddo
       enddo
       close(66)
#endif


      end program

      subroutine check(status)

      include "netcdf.inc"
      integer status

      if(status .ne. NF_NOERR) then
       print*, trim(NF_STRERROR(status))
       stop ":-("
      end if
      end subroutine check


      subroutine greg2jd(i, j, k, jd)

      implicit none

      integer i, j, k
      integer jd, jd1, jd2, jd3

      jd1 = k-32075+1461*(i+4800+(j-14)/12)/4
      jd2 = 367*(j-2-(j-14)/12*12)/12
      jd3 = 3*((i+4900+(j-14)/12)/100)/4
      jd = jd1+jd2-jd3

      end subroutine greg2jd

      subroutine jd2greg(jd, i, j, k)

      implicit none

      integer jd, i, j, k
      integer l, n

      l = jd + 68569
      n = 4*l/146097
      l = l - (146097*n + 3)/4
      i = 4000*(l + 1)/1461001
      l = l - 1461*i/4 + 31
      j = 80*l/2447
      k = l - 2447*j/80
      l = j/11
      j = j + 2 - 12*l
      i = 100*(n - 49) + i + l

      end subroutine



      subroutine writerestartfile(npart, health, survprob, mask_rho,    &
         xpart, ypart, zpart, tpart, dimsrho, restartfile)

      implicit none

      integer npart, nrstpart, ixx, iyy, ip, dimsrho(2)
      real health(*), xpart(*), ypart(*), zpart(*), tpart(*), survprob
      real*8 mask_rho(0:dimsrho(1)-1, 0:dimsrho(2)-1)
      character*(*) restartfile

      write(*,*) "writing restart file ", trim(restartfile)
      open(unit=20,file=trim(restartfile),action='write')
      nrstpart=0
      do ip=1, npart
       if ( health(ip) .gt. survprob ) then
        ixx=int(xpart(ip))
        iyy=int(ypart(ip))
        if ( mask_rho(ixx,iyy) .gt. 0.0 ) then
         nrstpart=nrstpart+1
        endif
       endif
      enddo

      write(20,*) nrstpart

      do ip=1, npart
       if ( health(ip) .gt. survprob ) then
        ixx=int(xpart(ip))
        iyy=int(ypart(ip))
        if ( mask_rho(ixx,iyy) .gt. 0.0 ) then
         write(20,50) xpart(ip),ypart(ip),zpart(ip),health(ip),tpart(ip)
        endif
       endif
      enddo
      close(20)
50    format(5(e12.6,1x))

      end subroutine


      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. nf_noerr) then
      print *, nf_strerror(iret)
      stop
      endif
      end

      subroutine fgennc(mode,ncfile,                                    &
       xi_rho, eta_rho, s_rho, ntimes, nsources, slice,                 &
       lat_rho, lon_rho, mask_rho, zeta, h, z, time, conc)

      implicit none

      include 'netcdf.inc'

      integer mode
      character*(*) ncfile
      integer xi_rho, eta_rho, s_rho, ntimes, nsources

      integer xi_rho_dim
      integer eta_rho_dim
      integer s_rho_dim
      integer source_dim
      integer ocean_time_dim

      integer  conc_rank
      integer  lon_rho_rank
      integer  lat_rho_rank
      integer  mask_rho_rank
      integer  zeta_rank
      integer  h_rank, h_id
      integer  s_rho_rank
      integer  ocean_time_rank

      parameter (conc_rank = 4)
      parameter (lon_rho_rank = 2)
      parameter (lat_rho_rank = 2)
      parameter (mask_rho_rank = 2)
      parameter (zeta_rank = 3)
      parameter (h_rank = 2)
      parameter (s_rho_rank = 1)
      parameter (ocean_time_rank = 1)

      integer  dims(conc_rank)
      integer  zeta_dims(zeta_rank)
      integer  s_rho_dims(s_rho_rank)
      integer  ocean_time_dims(ocean_time_rank)
! error status return
      integer  iret

! id netcdf file
      integer ncid

! id vars
      integer conc_id, lon_rho_id, lat_rho_id, mask_rho_id, zeta_id, &
      & h_i, d, s_rho_id, ocean_time_id
      integer i,j,k

      real conc(xi_rho, eta_rho, s_rho)
      real*8 lat_rho(xi_rho, eta_rho)
      real*8 lon_rho(xi_rho, eta_rho)
      real*8 mask_rho(xi_rho, eta_rho)
      real*8 h(xi_rho, eta_rho)
      real*8 z(s_rho)
      real*8 time
      real zeta(xi_rho, eta_rho)

      integer start(4), count(4), slice
      real pi, realval(1)

save ncid, conc_id, zeta_id, ocean_time_id

      goto ( 10,20,30 ) mode

10    continue

! enter define mode
      iret = nf_create(trim(ncfile), nf_clobber, ncid)
      call check_err(iret)

! define dimensions
      iret = nf_def_dim(ncid, 'xi_rho', xi_rho, xi_rho_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'eta_rho', eta_rho, eta_rho_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 's_rho', s_rho, s_rho_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'ocean_time', NF_UNLIMITED, ocean_time_dim)
      call check_err(iret)

! define vars
      dims(4) = ocean_time_dim
      dims(3) = s_rho_dim
      dims(2) = eta_rho_dim
      dims(1) = xi_rho_dim

      iret = nf_def_var(ncid, 'conc', nf_real, conc_rank, dims, conc_id)
      iret = nf_def_var(ncid, 'lon_rho', NF_DOUBLE, lon_rho_rank,       &
      (/xi_rho_dim, eta_rho_dim /), lon_rho_id)
      call check_err(iret)

      iret = nf_def_var(ncid, 'lat_rho', NF_DOUBLE, lat_rho_rank,       &
      (/xi_rho_dim, eta_rho_dim /), lat_rho_id)

      call check_err(iret)
      iret = nf_def_var(ncid, 'mask_rho', NF_DOUBLE, mask_rho_rank,     &
      (/xi_rho_dim, eta_rho_dim /), mask_rho_id)
      call check_err(iret)

      zeta_dims(3) = ocean_time_dim
      zeta_dims(2) = eta_rho_dim
      zeta_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta', NF_REAL, zeta_rank, zeta_dims,    &
      zeta_id)
      iret = nf_def_var(ncid, 'h', NF_DOUBLE, h_rank, dims, h_id)
      call check_err(iret)

      s_rho_dims(1) = s_rho_dim
      iret = nf_def_var(ncid, 's_rho', NF_DOUBLE, s_rho_rank, s_rho_dims&
      , s_rho_id)
      call check_err(iret)

      ocean_time_dims(1) = ocean_time_dim
      iret = nf_def_var(ncid, 'ocean_time', NF_DOUBLE, ocean_time_rank, &
      ocean_time_dims, ocean_time_id)
      call check_err(iret)

! assign attributes
      iret = nf_put_att_text(ncid, conc_id, 'long_name', 22, 'particle-c&
      oncentration')
      call check_err(iret)
      iret = nf_put_att_text(ncid, conc_id, 'units', 1, '1')
      call check_err(iret)
      iret = nf_put_att_text(ncid, conc_id, 'long_name', 46, 'concentrat&
      ion_of_suspended_matter_in_sea_water')
      call check_err(iret)
      iret = nf_put_att_text(ncid, conc_id, 'time', 10, 'ocean_time'    &
      )
      call check_err(iret)
      iret = nf_put_att_text(ncid, conc_id, 'coordinates', 21, 'lon_rho &
      lat_rho s_rho')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_rho_id, 'long_name', 23, 'longitu&
      de of rho-points')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_rho_id, 'units', 11, 'degree_east&
      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_rho_id, 'standard_name', 9, 'long&
      itude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_rho_id, 'field', 15, 'lon_rho, sc&
      alar')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_rho_id, '_coordinateaxistype', 3,&
       'lon')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_rho_id, 'long_name', 22, 'latitud&
      e of rho-points')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_rho_id, 'units', 11, 'degree_east&
      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_rho_id, 'standard_name', 8, 'lati&
      tude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_rho_id, 'field', 15, 'lat_rho, sc&
      alar')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_rho_id, '_coordinateaxistype', 3,&
       'lat')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, 'long_name', 26, 'S-coordin&
      ate at RHO-points')
      call check_err(iret)
      !iret = nf_put_att_text(ncid, s_rho_id, 'valid_min', 4, '-1.0')
      !call check_err(iret)
      !iret = nf_put_att_text(ncid, s_rho_id, 'valid_max', 3, '0.0')
      !call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, 'positive', 2, 'up')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, 'standard_name', 22, 'ocean&
      _sigma_coordinate')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, 'field', 13, 's_rho, scalar&
      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, '_CoordinateTransformType',&
       8, 'Vertical')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, '_CoordinateAxes', 5, 's_rh&
      o')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, '_CoordinateAxisType', 4, '&
      GeoZ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, s_rho_id, '_CoordinateZisPositive', 2&
      , 'up')
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_rho_id, 'long_name', 18, 'mask o&
      n RHO-points')
      call check_err(iret)
!     iret = nf_put_att_text(ncid, mask_rho_id, 'flag_values',  9, '0.0 &
!     , 1.0')
!     call check_err(iret)
!     iret = nf_put_att_text(ncid, mask_rho_id, 'flag_meanings', 10, 'la&
!     nd water')
!     call check_err(iret)
      iret = nf_put_att_text(ncid, mask_rho_id, 'coordinates', 15, 'lon_&
      rho lat_rho')
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_rho_id, 'units', 1, '1')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'long_name', 12, 'free-surfa&
      ce')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'time', 10, 'ocean_time')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'coordinates', 26, 'lon_rho &
      lat_rho ocean_time')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'field', 28, 'free-surface, &
      scalar, series')
      call check_err(iret)
      realval(1) = 9.9999999e+36
      iret = nf_put_att_real(ncid, zeta_id, '_FillValue', NF_REAL, 1,   &
      realval)
      call check_err(iret)
      iret = nf_put_att_text(ncid, h_id, 'long_name', 24, 'bathymetry at&
       RHO-points')
      call check_err(iret)
      iret = nf_put_att_text(ncid, h_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid, h_id, 'coordinates', 15, 'lon_rho lat&
      _rho')
      call check_err(iret)
      iret = nf_put_att_text(ncid, h_id, 'field', 12, 'bath, scalar')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ocean_time_id, 'long_name', 25, 'time&
       since initialization')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ocean_time_id, 'units', 37, 'seconds &
      since 1968-05-23 00:00:00 GMT')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ocean_time_id, 'calendar', 9, 'gregor&
      ian')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ocean_time_id, 'field', 20, 'time, sc&
      alar, series')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ocean_time_id, '_CoordinateAxisType',&
       4, 'Time')
      call check_err(iret)


      iret = nf_put_att_text(ncid, nf_global, 'format', 26, 'netcdf-3 64&
      bit offset file')
      call check_err(iret)
      iret = nf_put_att_text(ncid, nf_global, 'Conventions', 6, 'CF-1.4'&
      )
      call check_err(iret)
      iret = nf_put_att_text(ncid, nf_global, 'type', 22, 'roms/toms his&
      tory file')
      call check_err(iret)

! leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)

! write time-independent fields
      pi=4.0*atan(1.0)
      lat_rho=lat_rho*180./pi
      lon_rho=lon_rho*180./pi

      start(1)=1 ; start(2)=1
      count(1)=xi_rho ; count(2)=eta_rho
      iret=nf_put_vara_double(ncid,lat_rho_id,start,count,lat_rho)
      call check_err(iret)

      iret=nf_put_vara_double(ncid,lon_rho_id,start,count,lon_rho)
      call check_err(iret)

      iret=nf_put_vara_double(ncid,h_id,start,count,h)
      call check_err(iret)

      lat_rho=lat_rho*pi/180.
      lon_rho=lon_rho*pi/180.

      count(1)=s_rho
      iret=nf_put_vara_double(ncid,s_rho_id,start,count,z)
      call check_err(iret)

      return

 20   continue

! write time-dependent fields
      start(1)=1 ; start(2)=1 ; start(3)=1 ; start(4)=slice
      count(1)=xi_rho ; count(2)=eta_rho
      count(3)=s_rho ; count(4)=1
      !do k=1, count(3)
      !  do j=1, count(2)
      !    do i=1, count(1)
      !      if ( conc(i,j,k) > 0.0 ) print*, i, j, k, conc(i,j,k)
      !    enddo
      !  enddo
      !enddo
      iret=nf_put_vara_real(ncid,conc_id,start,count,conc)
      call check_err(iret)

      start(1)=1 ; start(2)=1 ; start(3)=slice
      count(1)=xi_rho ; count(2)=eta_rho ; count(3)=1
      iret=nf_put_vara_real(ncid,zeta_id,start,count,zeta)
      call check_err(iret)

      start(1)=slice
      count(1)=1
      iret=nf_put_vara_double(ncid,ocean_time_id,start,count,time)
      call check_err(iret)

      return

 30   continue

      iret = nf_close(ncid)
      call check_err(iret)


      end subroutine

      subroutine tolower(line)

      implicit none

      character*(*) line

      integer del
      integer i
!
!     Convert input file lines to lowercase, and save to output file.
!
      DEL = IACHAR('a') - IACHAR('A')                                   !    ! fid ASCII position diff between 'A' and 'a'

      DO I = 1, LEN_TRIM(LINE)                                          !    scan each character in line
       IF (LGE(LINE(I:I),'A') .AND. LLE(LINE(I:I),'Z'))  &              !    distance between 'A' and 'Z'..
        line(I:I) = ACHAR(IACHAR(LINE(I:I)) + DEL)                      !    then convert to lowercase
      END DO

      end subroutine
