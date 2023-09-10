      module TestCases

      use precis_param

      use io_units

      use Atmosphere,    only : &
         R_GAS_DRYAIR

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public set_TestCase_globvars,set_TestCase_windfield,DistSource,&
             Testcase_CalcErrors

        ! Publicly available variables

      integer :: TestCase
      logical :: fullASCIIOutput
      integer :: SubCase

      real(kind=ip) :: TC4_conc_1 = 1.0_ip  ! in kg/km3
      real(kind=ip) :: TC4_conc_2 = 0.0_ip  ! in kg/km3
      real(kind=ip) :: TC4_k_1 = 100.0_ip   ! in m2/s
      real(kind=ip) :: TC4_k_2 = 50.0_ip   ! in m2/s

      !real(kind=ip) :: MMS_eta0      = 1.72e-5_ip   ! Ref visc (Pa s)
      !real(kind=ip) :: MMS_suthcons  = 117.0_ip     ! Sutherland Constant (K)
      !real(kind=ip) :: MMS_suthtref  = 273.0_ip     ! Sutherland Ref
      !temperature (K)
      !real(kind=ip) :: MMS_eta0      = 1.8325e-5_ip   ! Ref visc (Pa s)
      !real(kind=ip) :: MMS_suthcons  = 120.0_ip     ! Sutherland Constant (K)
      !real(kind=ip) :: MMS_suthtref  = 296.16_ip     ! Sutherland Ref temperature (K)

      real(kind=ip) :: MMS_ztrop     = 10000.0_ip   ! Height of troposphere (m)
      !real(kind=ip) :: MMS_Ffac1     = 1.39896599613964_ip  ! WilsonHuang fac1
      !real(kind=ip) :: MMS_Ffac2     = 0.635137780328017_ip ! WilsonHuang fac2
      !real(kind=ip) :: MMS_Rs        = R_GAS_DRYAIR
      !real(kind=ip) :: MMS_temper0   = 300.0_ip     ! Surf temp (K)
      !real(kind=ip) :: MMS_pres0     = 100000.0_ip     ! Surf pres (Pa)
      !real(kind=ip) :: MMS_skinz     = 7000.0_ip    ! pres skin depth (m)
      !real(kind=ip) :: MMS_lpsr      = -0.007_ip    ! lapse rate (K/m)
      !real(kind=ip) :: MMS_diam      = 0.0001_ip    ! grain size (m)
      !real(kind=ip) :: MMS_rhom      = 2000.0_ip    ! density (kg/m3)
      real(kind=ip) :: MMS_U0        = -10.0_ip      ! ref x vel (m/s)
      real(kind=ip) :: MMS_V0        = -10.0_ip      ! ref y vel (m/s)
      real(kind=ip) :: MMS_W0        = 1.0_ip       ! ref z vel (m/s)
      real(kind=ip) :: MMS_Q0        = 1.0e10_ip    ! ref concen (kg/m3)
      real(kind=ip) :: MMS_zeta0     = 200000.0_ip    ! ref zeta (m)
      real(kind=ip) :: MMS_qxylen    = 200000.0_ip     ! horiz length of solution (m)
      !real(kind=ip) :: MMS_vzlen     = 100000.0_ip     ! horiz length of vz (m)

      logical :: MMS_USE_X    = .true.
      logical :: MMS_USE_Y    = .true.
      logical :: MMS_USE_Z    = .true.
      logical :: MMS_USE_T    = .true.
      logical :: MMS_USE_DifF = .true.

      !real(kind=ip) :: MMS_U0        = 0.0_ip      ! ref x vel (m/s)
      !real(kind=ip) :: MMS_V0        = 0.0_ip      ! ref y vel (m/s)
      !real(kind=ip) :: MMS_W0        = 1.0_ip       ! ref z vel (m/s)

      real(kind=ip) :: total_mass
      real(kind=ip) :: L1_toterror
      real(kind=ip) :: L2_toterror

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine set_TestCase_globvars

      use global_param,  only : &
         useVertAdvect,useHorzAdvect,&
         !DT_MIN,DT_MAX,&
         KM2_2_M2,HR_2_S

      use mesh,          only : &
        ZPADDING

      implicit none

#ifndef TESTCASE_1
#ifndef TESTCASE_2
#ifndef TESTCASE_3
#ifndef TESTCASE_4
#ifndef TESTCASE_5
#ifndef TESTCASE_6
      TestCase = 0
      useVertAdvect = .true.
      useHorzAdvect = .true.
      fullASCIIOutput = .false.
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef TESTCASE_1
      ! Horizontal advection:
      !  subcase 1 : x+y0
      !  subcase 2 : x-y0
      !  subcase 3 : x0y+
      !  subcase 4 : x0y-
      !  subcase 5 : x+y+
      !  subcase 6 : x-y-
      !  subcase 7 : x-y+
      !  subcase 8 : x+y-
      TestCase = 1
      useVertAdvect = .false.
      useHorzAdvect = .true.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
#endif
#ifdef TESTCASE_2
      ! Vertical advection:
      !  subcase 1 : Wind blows up (no fall velocity)
      !  subcase 2 : Wind blows down (no fall velocity)
      !  subcase 3 : No z wind (fall velocity +)
      !  subcase 4 : No z wind (fall velocity -)
      TestCase = 2
      useVertAdvect = .true.
      useHorzAdvect = .false.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
#endif
#ifdef TESTCASE_3
      ! Horizontal rotation of block and cone
      TestCase = 3
      useVertAdvect = .false.
      useHorzAdvect = .true.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
#endif
#ifdef TESTCASE_4
      ! Diffusion
      TestCase = 4
      useVertAdvect = .false.
      useHorzAdvect = .false.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
      TC4_k_1 = TC4_k_1 * HR_2_S/KM2_2_M2 ! in km^2/hr
      TC4_k_2 = TC4_k_2 * HR_2_S/KM2_2_M2 ! in km^2/hr
      ! more setting below after subcase setting
#endif
#ifdef TESTCASE_5
      ! Horizontal rotational shear
      TestCase = 5
      useVertAdvect = .false.
      useHorzAdvect = .true.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
#endif
#ifdef TESTCASE_6
      ! Method of manufactured solutions
      TestCase = 6
      useVertAdvect = .true.
      useHorzAdvect = .true.
      fullASCIIOutput = .true.
      ZPADDING  = 1.0_ip
#endif

#ifndef SUBCASE_1
#ifndef SUBCASE_2
#ifndef SUBCASE_3
#ifndef SUBCASE_4
#ifndef SUBCASE_5
#ifndef SUBCASE_6
#ifndef SUBCASE_7
#ifndef SUBCASE_8
      SubCase = 0
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#ifdef SUBCASE_1
      SubCase = 1
#endif
#ifdef SUBCASE_2
      SubCase = 2
#endif
#ifdef SUBCASE_3
      SubCase = 3
#endif
#ifdef SUBCASE_4
      SubCase = 4
#endif
#ifdef SUBCASE_5
      SubCase = 5
#endif
#ifdef SUBCASE_6
      SubCase = 6
#endif
#ifdef SUBCASE_7
      SubCase = 7
#endif
#ifdef SUBCASE_8
      SubCase = 8
#endif

!#ifdef TESTCASE_4
      ! Diffusion
      !if (SubCase.gt.3)then
      !  ! For CN integrations with no advection, we need to specify the
      !  ! time step explicitly
      !  DT_MIN    = 1.0e-5_ip
      !  DT_MAX    = 1.0e-5_ip
      !endif
!#endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"       TestCase = ",TestCase
        write(outlog(io),*)"        SubCase = ",SubCase
        write(outlog(io),*)"  useVertAdvect = ",useVertAdvect
        write(outlog(io),*)"  useHorzAdvect = ",useHorzAdvect
        write(outlog(io),*)"fullASCIIOutput = ",fullASCIIOutput
        write(outlog(io),*)"       ZPADDING = ",ZPADDING
      endif;enddo

      end subroutine set_TestCase_globvars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine set_TestCase_windfield

      use global_param,  only : &
         KM_2_M,MPS_2_KMPHR,PI

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd,z_cc_pd,&
         x_cc_pd,y_cc_pd,ZPADDING

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd

      use time_data,     only : &
         time,dt

      implicit none

      integer :: i,j,k

      real(kind=ip) :: lon_pole,lat_pole  ! pole of rotation
      real(kind=ip) :: omega              ! rotational velocity (rev/hour)
      real(kind=ip) :: ve_out,vn_out      ! output easterly and northerly velocities
      real(kind=ip) :: fac1, fac2
      real(kind=ip) :: dist

      real(kind=ip) :: xcc,ycc,zcc

      if(TestCase.eq.1)then
        omega = 0.1_ip ! concentration peak should travel 36 degrees/hour
        ! Test case 1 has a wind field SubCase
        if(SubCase.eq.1)then
          if(IsLatLon)then
            ! wind blows to the east
            ! Set pole of rotation to north pole
            lon_pole = 0.0_ip
            lat_pole = 90.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) =  1.0_ip
            vy_pd(:,:,:) =  0.0_ip
          endif
        elseif(SubCase.eq.2)then
          if(IsLatLon)then
            ! wind blows to the east
            ! Set pole of rotation to south pole
            lon_pole =  0.0_ip
            lat_pole = -90.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) = -1.0_ip
            vy_pd(:,:,:) =  0.0_ip
          endif
        elseif(SubCase.eq.3)then
          if(IsLatLon)then
            ! wind blows to the north
            ! Set pole of rotation to the Galapagos
            lon_pole = -90.0_ip
            lat_pole =  0.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) =  0.0_ip
            vy_pd(:,:,:) =  1.0_ip
          endif
        elseif(SubCase.eq.4)then
          if(IsLatLon)then
            ! wind blows to the south
            ! Set pole of rotation to the eastern Indian Ocean
            lon_pole = 90.0_ip
            lat_pole =  0.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) =  0.0_ip
            vy_pd(:,:,:) = -1.0_ip
          endif
        elseif(SubCase.eq.5)then
          if(IsLatLon)then
            ! wind blows to the northeast
            ! Set pole of rotation to the Athen, Wisconsin
            lon_pole = -90.0_ip
            lat_pole =  45.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                  !vx(i,j,k) = 0.0
                  !vy(i,j,k) = 0.0
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) =  1.0_ip
            vy_pd(:,:,:) =  1.0_ip
          endif
        elseif(SubCase.eq.6)then
          if(IsLatLon)then
            ! wind blows to the southwest
            ! Set pole of rotation to the southeast Austrailia/Antarctic Ocean
            lon_pole = 90.0_ip
            lat_pole = -45.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) = -1.0_ip
            vy_pd(:,:,:) = -1.0_ip
          endif
        elseif(SubCase.eq.7)then
          if(IsLatLon)then
            ! wind blows to the northwest
            ! Set pole of rotation to west of Chile
            lon_pole = -90.0_ip
            lat_pole = -45.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) = -1.0_ip
            vy_pd(:,:,:) =  1.0_ip
          endif
        elseif(SubCase.eq.8)then
          if(IsLatLon)then
            ! wind blows to the southeast
            ! Set pole of rotation to the Gobi desert
            lon_pole =  90.0_ip
            lat_pole =  45.0_ip
            do i=-1,nxmax+2
              do j=-1,nymax+2
                do k=-1,nzmax+2
                  call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                   lon_pole,lat_pole,omega,        &
                                 ve_out,vn_out)
                  vx_pd(i,j,k) = ve_out
                  vy_pd(i,j,k) = vn_out
                enddo
              enddo
            enddo
          else
            vx_pd(:,:,:) =  1.0_ip
            vy_pd(:,:,:) = -1.0_ip
          endif
        else
          write(global_info,*)"Invalid SubCase"
          stop 1
        endif
      endif

      if(TestCase.eq.2)then
        if(ZPADDING.gt.1.0_ip)then
          write(global_info,*)"ERROR: ZPADDING != 1.0"
          stop 1
        endif
        ! Test case 2 has a wind field SubCase
        if(SubCase.eq.1)then
            ! Wind blows up (no fall velocity)
          vx_pd(:,:,:) =  0.0_ip
          vy_pd(:,:,:) =  0.0_ip
          vz_pd(:,:,:) =  1.0_ip
          vf_pd(:,:,:,:) =  0.0_ip
        elseif(SubCase.eq.2)then
            ! Wind blows down (no fall velocity)
          vx_pd(:,:,:) =  0.0_ip
          vy_pd(:,:,:) =  0.0_ip
          vz_pd(:,:,:) = -1.0_ip
          vf_pd(:,:,:,:) =  0.0_ip
        elseif(SubCase.eq.3)then
            ! No z wind (fall velocity +)
          vx_pd(:,:,:) =  0.0_ip
          vy_pd(:,:,:) =  0.0_ip
          vz_pd(:,:,:) =  0.0_ip
          vf_pd(:,:,:,:) =  1.0_ip
        elseif(SubCase.eq.4)then
            ! No z wind (fall velocity -)
          vx_pd(:,:,:) =  0.0_ip
          vy_pd(:,:,:) =  0.0_ip
          vz_pd(:,:,:) =  0.0_ip
          vf_pd(:,:,:,:) =  -1.0_ip
        endif
      endif

      if(TestCase.eq.3)then
        ! Test case 3 has circular wind field
        if(IsLatLon)then
          ! Set angular velocity in rev/hour
          omega = -1.0_ip/PI         ! 1 revolution in pi hours
              ! wind rotates clockwise
             ! Set pole of rotation to middle of domain 
          lon_pole = 360.0_ip
          lat_pole = 0.0_ip
          do i=-1,nxmax+2
            do j=-1,nymax+2
              do k=-1,nzmax+2
                call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                 lon_pole,lat_pole,omega,        &
                               ve_out,vn_out)
                vx_pd(i,j,k) = ve_out
                vy_pd(i,j,k) = vn_out
              enddo
            enddo
          enddo
        else
          do i=-1,nxmax+2
            do j=-1,nymax+2
                ! cartesian grid
              vx_pd(i,j,:) = 40.0_ip * ( y_cc_pd(j))/500.0_ip
              vy_pd(i,j,:) = 40.0_ip * (-x_cc_pd(i))/500.0_ip
              vx_pd(i,j,:) =  2.0_ip *   y_cc_pd(j)
              vy_pd(i,j,:) = -2.0_ip *   x_cc_pd(i)
            enddo
          enddo
        endif
      endif

      if(TestCase.eq.4)then
        ! Test case 4 has zero wind field
        vx_pd = 0.0_ip
        vy_pd = 0.0_ip
        vz_pd = 0.0_ip
      endif

      if(TestCase.eq.5)then
        ! Test case 5 has transient circular wind field
        ! See Numerical Methods for Wave Equations in Geophysical Fluid
        ! Dynamics; Durran p284; 1998
        if(IsLatLon)then
             ! Set pole of rotation to middle of domain 
          lon_pole = 360.0_ip
          lat_pole = 0.0_ip
          do i=-1,nxmax+2
            do j=-1,nymax+2
              do k=-1,nzmax+2
                ! Omega is a function of the distance from the pole of rotation
                dist = sqrt((lon_cc_pd(i)-lon_pole)**2.0_ip + &
                            (lat_cc_pd(j)-lat_pole)**2.0_ip)
                  ! This gives a maximum of 1rev/2hours
                omega = -sin(dist/36.0_ip * 2.0_ip*PI) * &
                         cos(PI*time*0.25_ip) / 2.0_ip
                call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                 lon_pole,lat_pole,omega,        &
                               ve_out,vn_out)
                vx_pd(i,j,k) = ve_out
                vy_pd(i,j,k) = vn_out

                ! Calculate dv_dt by finite differencing
                omega = -sin(dist/36.0_ip * 2.0_ip*PI) * &
                         cos(PI*(time+0.25_ip*dt)*0.25_ip) / 2.0_ip
                call set_LL_wind(lon_cc_pd(i),lat_cc_pd(j),z_cc_pd(k),  &
                                 lon_pole,lat_pole,omega,        &
                               ve_out,vn_out)
              enddo
            enddo
          enddo
        else
          do i=-1,nxmax+2
            do j=-1,nymax+2
                ! cartesian grid
              fac1 = (sin(PI*x_cc_pd(i))**2.0_ip) * &
                     (sin(2.0_ip*PI*y_cc_pd(j)))
              fac2 = -(sin(PI*y_cc_pd(j))**2.0_ip) * &
                      (sin(2.0_ip*PI*x_cc_pd(i))) 
              vx_pd(i,j,:) = fac1 * cos(PI*time*0.2_ip)
              vy_pd(i,j,:) = fac2 * cos(PI*time*0.2_ip)
            enddo
          enddo
        endif
      endif

      if(TestCase.eq.6)then
        ! MMS
        if(IsLatLon)then
          write(global_info,*)"MMS solution only in Cartesian."
          stop 1
        else
          ! Vx is constant
          if(MMS_USE_X)then
            vx_pd = MMS_U0 * MPS_2_KMPHR
          else
            vx_pd = 0.0_ip
          endif
          ! Vy only varies in z
          do k=-1,nzmax+2
            zcc = z_cc_pd(k) * KM_2_M
            if(MMS_USE_Y)then
              vy_pd(:,:,k) = MMS_V0 * 0.5_ip * (1.0_ip + tanh(zcc-MMS_ztrop))
            else
              vy_pd = 0.0_ip
            endif
          enddo
          do i=-1,nxmax+2
            xcc = x_cc_pd(i) * KM_2_M
            do j=-1,nymax+2
              ycc = y_cc_pd(j) * KM_2_M
              ! z-velocity in m/s
              if(MMS_USE_Z)then
                vz_pd(i,j,:) = -MMS_W0*cos(xcc*PI/MMS_qxylen)*cos(ycc*PI/MMS_qxylen)
              else
                vz_pd = 0.0_ip
              endif
            enddo
          enddo
        endif
      endif


      end subroutine set_TestCase_windfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine set_LL_wind(lon_pt,lat_pt,z_pt,             &
                             lon_pole,lat_pole,omega,        &
                             ve_out,vn_out)

      use global_param,  only : &
         PI,DEG2RAD,EPS_TINY

      use projection

      implicit none

      real(kind=ip),intent(in)  :: lon_pt,lat_pt,z_pt ! coordinates of input point
      real(kind=ip),intent(in)  :: lon_pole,lat_pole  ! pole of rotation
      real(kind=ip),intent(in)  :: omega              ! rotational velocity (rev/hour)
      real(kind=ip),intent(out) :: ve_out,vn_out      ! output easterly and northerly velocities

      real(kind=ip) :: phi_in,lam_in,rad_in
      !real(kind=ip) :: phir
      real(kind=ip) :: lamr,cophir
      real(kind=ip) :: slmr,clmr,scphr,ccphr
      real(kind=ip) :: slm,clm,sph,cph
      real(kind=ip),dimension(3,3) :: R1,R2,R1b,R2b
      real(kind=ip),dimension(3) :: xvec,xnew
      real(kind=ip),dimension(3) :: v,vnew
      real(kind=ip),dimension(3) :: vr,vtan,ve_norm,vn_norm,xvec_norm
      real(kind=ip),dimension(3) :: ve_vec,vn_vec
      real(kind=ip) :: sxlm,cxlm,cxph
      real(kind=ip) :: xcolat,xlm,xph
      real(kind=ip) :: eqvel                          ! equitorial velocity of rotation (km/hr)

      rad_in = PJ_Re + z_pt
      eqvel = 2.0_ip*PI*rad_in*omega
      phi_in  = lat_pt   * DEG2RAD;
      lam_in  = lon_pt   * DEG2RAD;
      !phir    = lat_pole * DEG2RAD;
      lamr    = lon_pole * DEG2RAD;
      cophir = (90.0_ip-lat_pole)*DEG2RAD;
      slmr = sin(lamr);     clmr = cos(lamr);
      scphr = sin(cophir); ccphr = cos(cophir);

      slm = sin(lam_in)
      clm = cos(lam_in)
      sph = sin(phi_in);
      cph = cos(phi_in);

      ! Rotate phi_r about -y (from pole towards x-axis)
      !R1 = [[ ccphr 0.0_ip  scphr];
      !      [  0.0_ip   1.0_ip    0.0_ip  ];
      !      [-scphr 0.0_ip  ccphr]];
      R1(1,1) = ccphr;  R1(1,2) = 0.0_ip; R1(1,3) = scphr;
      R1(2,1) = 0.0_ip; R1(2,2) = 1.0_ip; R1(2,3) = 0.0_ip;
      R1(3,1) = -scphr; R1(3,2) = 0.0_ip; R1(3,3) = ccphr;

      ! Rotate lambda_r about z'
      !R2 = [[ clmr -slmr 0.0_ip  ];
      !      [ slmr  clmr 0.0_ip  ];
      !      [  0.0_ip     0.0_ip   1.0_ip  ]];
      R2(1,1) = clmr;   R2(1,2) = -slmr;  R2(1,3) = 0.0_ip;
      R2(2,1) = slmr;   R2(2,2) =  clmr;  R2(2,3) = 0.0_ip;
      R2(3,1) = 0.0_ip; R2(3,2) = 0.0_ip; R2(3,3) = 1.0_ip;

      !Backwards rotations
      ! Rotate - lambda_r about z'
      !R2b = [[ clmr  slmr 0.0_ip  ];
      !       [-slmr  clmr 0.0_ip  ];
      !       [  0.0_ip     0.0_ip   1.0_ip  ]];
      R2b(1,1) = clmr;   R2b(1,2) =  slmr;  R2b(1,3) = 0.0_ip;
      R2b(2,1) =-slmr;   R2b(2,2) =  clmr;  R2b(2,3) = 0.0_ip;
      R2b(3,1) = 0.0_ip; R2b(3,2) = 0.0_ip; R2b(3,3) = 1.0_ip;

      ! Rotate -phi_r about -y 
      !R1b = [[ ccphr 0.0_ip -scphr];
      !       [  0.0_ip   1.0_ip    0.0_ip  ];
      !       [ scphr 0.0_ip  ccphr]];
      R1b(1,1) = ccphr;  R1b(1,2) = 0.0_ip; R1b(1,3) = -scphr;
      R1b(2,1) = 0.0_ip; R1b(2,2) = 1.0_ip; R1b(2,3) = 0.0_ip;
      R1b(3,1) = scphr;  R1b(3,2) = 0.0_ip; R1b(3,3) = ccphr;

        ! Get coordinates in standard orientation
      xvec(1) = rad_in*clm*cph
      xvec(2) = rad_in*slm*cph
      xvec(3) = rad_in*sph
        ! Normalize position vector
      xvec_norm = xvec/sqrt(dot_product(xvec,xvec))
        ! back-rotate position vector
      xnew = matmul(R1b,matmul(R2b,xvec))
        ! Get lon/lat for this new position
      xcolat = acos(xnew(3)/rad_in);
      xph    = (90.0_ip*DEG2RAD)-xcolat;
      xlm    = atan2(xnew(2),xnew(1));
      sxlm = sin(xlm); cxlm = cos(xlm);
      cxph = cos(xph);

        ! Get velocity (fraction of equitorial velocity)
      v(1) = -sxlm*cxph*eqvel
      v(2) =  cxlm*cxph*eqvel
      v(3) = 0.0_ip
        ! Rotate velocity to get the correct x,y,z velocity
      vnew = matmul(R2,matmul(R1,v))
        ! Now convert this 3D velocity to vn,ve
        ! Get radial velocity (should be zero)
      vr   = dot_product(vnew,xvec_norm)*xvec_norm
      vtan = vnew - vr

        ! To get unit vector ve, cross z with r
      ve_norm(1) = -xvec_norm(2)
      ve_norm(2) =  xvec_norm(1)
      ve_norm(3) = 0.0_ip
      ve_norm = ve_norm/(sqrt(dot_product(ve_norm,ve_norm)))
        ! To get unit vector vn, cross r with ve
      vn_norm(1) = (xvec_norm(2)*ve_norm(3) - xvec_norm(3)*ve_norm(2))
      vn_norm(2) =-(xvec_norm(1)*ve_norm(3) - xvec_norm(3)*ve_norm(1))
      vn_norm(3) = (xvec_norm(1)*ve_norm(2) - xvec_norm(2)*ve_norm(1))
      vn_norm = vn_norm/(sqrt(dot_product(vn_norm,vn_norm)))

        ! To get ve, project tangent v onto ve_norm
      ve_vec = dot_product(vtan,ve_norm)*ve_norm
        ! To get vn, subtract vn from tangent v
      vn_vec = vtan - ve_vec

        ! All that needs to be returned are the signed lengths
      if(dot_product(ve_vec,ve_norm).gt.EPS_TINY)then
        ve_out = sqrt(dot_product(ve_vec,ve_vec))
      elseif(dot_product(ve_vec,ve_norm).lt.-EPS_TINY)then
        ve_out = -sqrt(dot_product(ve_vec,ve_vec))
      else
        ve_out = 0.0_ip
      endif
      if(dot_product(vn_vec,vn_norm).gt.EPS_TINY)then
        vn_out = sqrt(dot_product(vn_vec,vn_vec))
      elseif(dot_product(vn_vec,vn_norm).lt.-EPS_TINY)then
        vn_out = -sqrt(dot_product(vn_vec,vn_vec))
      else
        vn_out = 0.0_ip
      endif

      end subroutine set_LL_wind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DistSource

      use global_param,  only : &
         PI

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,lat_cc_pd,lon_cc_pd,kappa_pd,z_cc_pd,&
         x_cc_pd,y_cc_pd,ts0,IsLatLon !,dx,dy

      use solution,      only : &
         concen_pd

      use Diffusion,     only : &
         kx,ky,kz

      implicit none

      integer :: i,j,k,n
      real(kind=ip) :: r

      real(kind=ip) :: xpeak,ypeak,zpeak
      real(kind=ip) :: char_len,Dist
      real(kind=ip) :: fac_1         = 0.0_ip
      real(kind=ip) :: temp1         = 0.0_ip
      real(kind=ip) :: Const1_1d     = 1.0_ip
      real(kind=ip) :: Const1_2d     = 0.6820926132509800_ip
      !real(kind=ip) :: Const1_3d     = 0.477464829275686_ip
      real(kind=ip) :: lon_shifted

      !real(kind=ip) ::  Vs,Vs_p1,Vs_m1,!dws_dz
      !real(kind=ip) ::  kap
      !real(kind=ip) ::  zeta
      !real(kind=ip) ::  xcc,ycc,zcc
      !real(kind=ip) ::  x_qxylen,y_qxylen
      !real(kind=ip) ::  tanhx,tanhy,tanhz
      !real(kind=ip) ::  sechx,sechy
      !real(kind=ip) ::  x_vzlen,y_vzlen
      !real(kind=ip) ::  expzcc,expzccm
      !real(kind=ip) ::  MMS_TrueSol,MMS_Source
      !real(kind=ip) ::  mms_src
      !real(kind=ip) ::  init_sol
      !real(kind=ip) ::  Dq_dt,Dq_dx,Dq_dy,Dq_dz
      !real(kind=ip) ::  D2q_dx2,D2q_dy2,D2q_dz2

      total_mass = 0.0_ip

      if(TestCase.eq.1)then
        if(IsLatLon)then
          ! Setting up concentration for 2D pulse to be advected
          ! horizontally
          xpeak = 360.0_ip
          ypeak = 0.0_ip
          char_len = 6.0_ip
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                  ! This is a 2D cubic spline concentration centered at
                  ! 0,0 with a characteristic width of h=0.1 and
                  ! normalized to unit volume

                  if ( (abs(lon_cc_pd(i)-xpeak).lt.2.0_ip*char_len) .or.  &
                       (abs(lat_cc_pd(j)-ypeak).lt.2.0_ip*char_len) ) then

       Dist = sqrt(  (lon_cc_pd(i)-xpeak)**2.0_ip + (lat_cc_pd(j)-ypeak)**2.0_ip )
       r = Dist/char_len
       fac_1 = Const1_2d/(char_len*char_len)
       if (r>2.0_ip) then
         concen_pd(i,j,k,n,:) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         concen_pd(i,j,k,n,:) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         concen_pd(i,j,k,n,:) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif

                  else
                    concen_pd(i,j,k,n,:) = 0.0_ip
                  endif
                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        else
          ! Setting up concentration for 2D pulse to be advected
          ! horizontally
          xpeak = 0.0_ip
          ypeak = 0.0_ip
          char_len = 0.1_ip
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                    ! This is a 2D cubic spline concentration centered at
                    ! 0,0 with a characteristic width of h=0.1 and
                    ! normalized to unit volume
  
                  if ( (abs(x_cc_pd(i)-xpeak).lt.2.0_ip*char_len) .or.  &
                       (abs(y_cc_pd(j)-ypeak).lt.2.0_ip*char_len) ) then

       Dist = sqrt(  (x_cc_pd(i)-xpeak)**2.0_ip + (y_cc_pd(j)-ypeak)**2.0_ip )
       r = Dist/char_len
       fac_1 = Const1_2d/(char_len*char_len)
       if (r>2.0_ip) then
         concen_pd(i,j,k,n,:) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         concen_pd(i,j,k,n,:) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         concen_pd(i,j,k,n,:) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif

                  else
                    concen_pd(i,j,k,n,:) = 0.0_ip
                  endif
  
                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        endif !IsLatLon
      endif !TestCase.eq.1

      if(TestCase.eq.2)then
        ! Setting up concentration for 1D pulse to be advected
        ! vertically
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Setting up concentration for TestCase 2"
        endif;enddo
        zpeak = 1.0_ip
        char_len = 0.1
        do n=1,nsmax
          do k=-1,nzmax+2
                  ! This is a 1D cubic spline concentration centered at
                  ! 0,0 with a characteristic width of h=0.1 and
                  ! normalized to unit volume
                 Dist = abs(z_cc_pd(k)-zpeak)
                if ( Dist.lt.2.0_ip*char_len) then

       r = Dist/char_len
       fac_1 = Const1_1d/(char_len)
       if (r>2.0_ip) then
         concen_pd(:,:,k,n,:) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         concen_pd(:,:,k,n,:) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         concen_pd(:,:,k,n,:) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif

                else
                  concen_pd(:,:,k,n,:) = 0.0_ip
                endif
          enddo ! loop over k
        enddo ! loop over n
      endif

      if(TestCase.eq.3)then
        ! Setting up concentration for 2D cone/box to be advected
        ! horizontally
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Setting up concentration for TestCase 3"
        endif;enddo
        if(IsLatLon)then
          ! Lat/lon grid
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                    ! A similar concentration is described in Example 20.1 of
                    ! LeVeque's Finite Volume book (p.460)
                      ! Set up "square" 16deg x 16deg at lon = 16
                  lon_shifted = lon_cc_pd(i) - 360.0_ip
                  if ((lon_shifted.lt.24.0_ip).and.(lon_shifted.gt.8.0_ip).and.&
                      (lat_cc_pd(j).gt.-8.0_ip).and.(lat_cc_pd(j).lt.8.0_ip))then
                    concen_pd(i,j,k,n,:) = 1.0_ip
                  else
                    concen_pd(i,j,k,n,:) = 0.0_ip
                  endif
                      ! Set up "cone" at lon = -16 with radius of 8deg
                  r = sqrt((lon_shifted+16.0_ip)**2 + lat_cc_pd(j)**2)
                  if (r .lt. 8.0_ip) then
                    concen_pd(i,j,k,n,:) = 1.0_ip - r/8.0_ip
                  endif
                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        else
          ! This is the projected grid
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                    ! This concentration is described in Example 20.1 of
                    ! LeVeque's Finite Volume book (p.460)
                  if ((x_cc_pd(i).lt.0.6_ip).and.(x_cc_pd(i).gt.0.1_ip).and.      &
                      (y_cc_pd(j).gt.-0.25_ip).and.(y_cc_pd(j).lt.0.25_ip)) then
                    concen_pd(i,j,k,n,:) = 1.0_ip
                  else
                    concen_pd(i,j,k,n,:) = 0.0_ip
                  endif
                  r = sqrt((x_cc_pd(i)+0.45_ip)**2 + y_cc_pd(j)**2)
                  if (r .lt. 0.35_ip) then
                    concen_pd(i,j,k,n,:) = 1.0_ip - r/0.35_ip
                  endif
                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        endif
      endif ! TestCase.eq.3

      if(TestCase.eq.4)then
        ! Setting up halfspace concentration for diffusion testing
                !  With diffusivity_horz=12000 and dx=dy=dz=0.5, use dt =
                !  5.0e-4_ip
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Setting up concentration for TestCase 4"
        endif;enddo
        if(SubCase.eq.1.or.SubCase.eq.4)then
          concen_pd = 0.0_ip
          ky = 0.0_ip
          kz = 0.0_ip
          do i=1,nxmax
            if(x_cc_pd(i).lt.1.0_ip) then
              concen_pd(i,:,:,:,:) = TC4_conc_1
              kx(i,:,:) = TC4_k_1
            else
              concen_pd(i,:,:,:,:) = TC4_conc_2
              kx(i,:,:) = TC4_k_2
            endif
          enddo
        elseif(SubCase.eq.2.or.SubCase.eq.5)then
          concen_pd = 0.0_ip
          kx = 0.0_ip
          kz = 0.0_ip
          do j=1,nymax
            if(y_cc_pd(j).lt.1.0_ip) then
              concen_pd(:,j,:,:,:) = TC4_conc_1
              ky(:,j,:) = TC4_k_1
            else
              concen_pd(:,j,:,:,:) = TC4_conc_2
              ky(:,j,:) = TC4_k_2
            endif
          enddo
        elseif(SubCase.eq.3.or.SubCase.eq.6)then
          concen_pd = 0.0_ip
          kx = 0.0_ip
          ky = 0.0_ip
          do k=1,nzmax
            if(z_cc_pd(k).lt.1.0_ip) then
              concen_pd(:,:,k,:,:) = TC4_conc_1
              kz(:,:,k) = TC4_k_1
            else
              concen_pd(:,:,k,:,:) = TC4_conc_2
              kz(:,:,k) = TC4_k_2
            endif
          enddo
        endif
      endif

      if(TestCase.eq.5)then
        ! Setting up concentration pulse for circular shear advection
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Setting up concentration for TestCase 5"
        endif;enddo
        if(IsLatLon)then
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                    ! This concentration is similar to Example 5.4.4 of
                    ! Durran's Finite Volume book (p.284)
                  lon_shifted = lon_cc_pd(i) - 360.0_ip
                  r = sqrt((lon_shifted+16.0_ip)**2.0_ip + &
                           (lat_cc_pd(j))**2.0_ip) / 16.0_ip
                  r = min(1.0_ip,r)
                  concen_pd(i,j,k,n,:) = 0.5_ip*(1.0_ip+cos(PI*r))

                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        else
          ! Setting up concentration pulse for circular shear advection
          do n=1,nsmax
            do k=-1,nzmax+2
              do j=-1,nymax+2
                do i=-1,nxmax+2
                    ! This concentration is described in Example 5.4.4 of
                    ! Durran's Finite Volume book (p.284)
                  r = 4.0_ip*sqrt((x_cc_pd(i)-0.25_ip)**2.0_ip + &
                                  (y_cc_pd(j)-0.25_ip)**2.0_ip)
                  r = min(1.0_ip,r)
                  concen_pd(i,j,k,n,:) = 0.5_ip*(1.0_ip+cos(PI*r))
                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
        endif
      endif ! TestCase.eq.5


!      if(TestCase.eq.6)then
!!        ! Setting up MMS
!        if(IsLatLon)then
!          write(outlog(io),*)"MMS not set up for Lat/Lon grids."
!          stop 1
!        else
!        do n=1,nsmax
!          do k=1,nzmax
!            ! Fall velocity is only a function of z
!            Vs_m1 = vf(1,1,k-1,n)
!            Vs    = vf(1,1,k  ,n)
!            Vs_p1 = vf(1,1,k+1,n)
!            do j=1,nymax
!             do i=1,nxmax
!               ! Set up initial distribution
!               if(itime.eq.0)then
!                 init_sol = MMS_TrueSol(x_cc(i),y_cc(j),z_cc(k),time)
!                 concen_pd(i,j,k,n,ts0) = init_sol
!               endif
!
!            ! Get approximate z derivative of vertical velocities
!            dws_dz = ((Vs_p1+vz(i,j,k+1))-(Vs_m1+vz(i,j,k-1)))/(2.0_ip*dz_vec(k))
!
!            mms_src = MMS_Source(x_cc(i),y_cc(j),z_cc(k),time, &
!                                 vx(i,j,k),vy(i,j,k),vz(i,j,k),Vs, &
!                                 dws_dz)
!
!            concen_pd(i,j,k,n,ts1) =  concen_pd(i,j,k,n,ts0) + dt*mms_src 
!
!              enddo ! loop over i
!            enddo ! loop over j
!          enddo ! loop over k
!        enddo ! loop over n
!         ! copy q-star slice back to q
!        concen_pd(1:nxmax,1:nymax,1:nzmax,:,ts0) = &
!          concen_pd(1:nxmax,1:nymax,1:nzmax,:,ts1)
!
!        endif
!      endif

      do n=1,nsmax
        do k=1,nzmax
          do j=1,nymax
            do i=1,nxmax
              total_mass = total_mass + concen_pd(i,j,k,n,ts0)*kappa_pd(i,j,k)
            enddo ! loop over i
          enddo ! loop over j
        enddo ! loop over k
      enddo ! loop over n

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Total mass = ",total_mass
      endif;enddo

      end subroutine DistSource


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Testcase_CalcErrors

      use global_param,  only : &
         PI,KM3_2_M3

      use time_data,     only : &
         time

      use solution,      only : &
         outflow_xy1_pd,outflow_xy2_pd,outflow_xz1_pd,outflow_xz2_pd, &
         outflow_yz1_pd,outflow_yz2_pd,concen_pd

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dz_vec_pd,kappa_pd,x_cc_pd,y_cc_pd,z_cc_pd,lat_cc_pd,&
         lon_cc_pd,ts1,IsLatLon,dx,dy

      implicit none

      integer :: i,j,k,n
      integer :: lx,ly,lz
      real(kind=ip) :: r

      real(kind=ip) :: xpeak,ypeak,zpeak
      real(kind=ip) :: char_len,Dist
      real(kind=ip) :: fac_1         = 0.0_ip
      real(kind=ip) :: temp1         = 0.0_ip
      real(kind=ip) :: Const1_1d     = 1.0_ip
      real(kind=ip) :: Const1_2d     = 0.6820926132509800_ip
      !real(kind=ip) :: Const1_3d     = 0.477464829275686_ip
      real(kind=ip) :: tsol
      real(kind=ip),dimension(nxmax)       :: tsolx
      real(kind=ip),dimension(nymax)       :: tsoly
      real(kind=ip),dimension(nzmax)       :: tsolz
      real(kind=ip),dimension(nxmax,nymax) :: err,truesol
      real(kind=ip),dimension(nzmax)       :: errz,truesolz
      real(kind=ip),dimension(nxmax,nymax,nzmax) :: err3D
      !real(kind=ip),dimension(nxmax,nymax,nzmax) :: truesol3D

      real(kind=ip) :: lon_pole,lat_pole
      real(kind=ip) :: deg_rot
      real(kind=ip) :: lon_br,lat_br

      character(len=30) :: ofile1,ofile2
      real(kind=ip) :: lon_shifted

      real(kind=ip) :: MassConsError
      real(kind=ip) :: TotalVol
      real(kind=ip) :: eta1,eta2
      real(kind=ip) :: xm,ym,zm,Tc

      !real(kind=ip) :: MMS_TrueSol

      !Note: concentrations should just be set for concen(:,:,:,:,ts0),
      !but this seems to cause problems with the semi-lagrange routine

      L1_toterror = 0.0_ip
      L2_toterror = 0.0_ip
      MassConsError = 0.0_ip
      TotalVol    = 0.0_ip

      do k=1,nzmax
        do j=1,nymax
          do i=1,nxmax
            TotalVol = TotalVol + kappa_pd(i,j,k)
          enddo ! loop over i
        enddo ! loop over j
      enddo ! loop over k

      if(TestCase.eq.1)then
        if(IsLatLon)then
          ! To compare with the true solution, we'll back-rotate each
          ! cell-center point and to the original position and calculate
          ! the correct concentration from the original profile.
          ! Each subcase has a different pole of rotations.
          if(SubCase.eq.1)then
            lon_pole = 0.0_ip
            lat_pole = 90.0_ip
          elseif(SubCase.eq.2)then
            lon_pole =  0.0_ip
            lat_pole = -90.0_ip
          elseif(SubCase.eq.3)then
            lon_pole = -90.0_ip
            lat_pole =  0.0_ip
          elseif(SubCase.eq.4)then
             lon_pole = 90.0_ip
             lat_pole =  0.0_ip
          elseif(SubCase.eq.5)then
            lon_pole = -90.0_ip
            lat_pole =  45.0_ip
          elseif(SubCase.eq.6)then
            lon_pole = 90.0_ip
            lat_pole = -45.0_ip
          elseif(SubCase.eq.7)then
            lon_pole = -90.0_ip
            lat_pole = -45.0_ip
          elseif(SubCase.eq.8)then
            lon_pole =  90.0_ip
            lat_pole =  45.0_ip
          endif

          if(time.lt.0.6_ip)then
             deg_rot = -18.0_ip
          else
             deg_rot = -36.0_ip
          endif

          xpeak = 360.0_ip
          ypeak = 0.0_ip
          char_len = 6.0_ip
        do n=1,nsmax
          do k=1,nzmax
            do j=1,nymax
              do i=1,nxmax
                  ! This is a 2D cubic spline concentration centered at
                  ! 0,0 with a characteristic width of h=0.1 and
                  ! normalized to unit volume
                call back_rotate(lon_cc_pd(i),lat_cc_pd(j),    &
                                 lon_pole,lat_pole,deg_rot, &
                                 lon_br,lat_br)
                if (lon_br.lt.180.0_ip)lon_br=lon_br+360.0_ip
                if ( (abs(lon_br-xpeak).lt.2.0_ip*char_len) .or.  &
                     (abs(lat_br-ypeak).lt.2.0_ip*char_len) ) then

       Dist = sqrt(  (lon_br-xpeak)**2.0_ip + (lat_br-ypeak)**2.0_ip)
       r = Dist/char_len
       fac_1 = Const1_2d/(char_len*char_len)
       if (r>2.0_ip) then
         truesol(i,j) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         truesol(i,j) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         truesol(i,j) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif

                else
                  truesol(i,j) = 0.0_ip
                endif

                ! Zero nearly zero concentration values since these cause problems
                ! when written out.
                if (abs(concen_pd(i,j,k,n,ts1)).lt.1.0e-50_ip)then
                  concen_pd(i,j,k,n,ts1) = 0.0_ip
                endif
                err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)

                L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

              enddo ! loop over i
            enddo ! loop over j
          enddo ! loop over k
        enddo ! loop over n
        ! Now account for mass conservation at boundaries
        do n=1,nsmax
          do i=1,nxmax
            ! Adding bottom and top
            do j=1,nymax
              MassConsError = MassConsError + outflow_xy1_pd(i,j,n)*kappa_pd(i,j,      0) &
                                            + outflow_xy2_pd(i,j,n)*kappa_pd(i,j,nzmax+1)
            enddo
            do k=1,nzmax
              MassConsError = MassConsError + outflow_xz1_pd(i,k,n)*kappa_pd(i,     0,k) &
                                            + outflow_xz2_pd(i,k,n)*kappa_pd(i,nymax+1,k)
            enddo
          enddo

          do j=1,nymax
            do k=1,nzmax
              MassConsError = MassConsError + outflow_yz1_pd(j,k,n)*kappa_pd(      0,j,k) &
                                            + outflow_yz2_pd(j,k,n)*kappa_pd(nxmax+1,j,k)
            enddo
          enddo
        enddo

        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror /  nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        if(time.lt.0.6)then
          write(ofile1,'(a10,i1,a12)') &
                'TC1_LL_Sub',SubCase,'_St1_err.dat'
          write(ofile2,'(a10,i1,a12)') &
                'TC1_LL_Sub',SubCase,'_St1_sol.dat'
        else
          write(ofile1,'(a10,i1,a12)') &
                'TC1_LL_Sub',SubCase,'_St2_err.dat'
          write(ofile2,'(a10,i1,a12)') &
                'TC1_LL_Sub',SubCase,'_St2_sol.dat'
        endif
        open(200,file=ofile1,status='replace')
        open(201,file=ofile2,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError
        do i=1,nxmax
          do j=1,nymax
            write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
          enddo
        enddo
        close(200)
        close(201)

        else
        ! Setting up concentration for 2D pulse to be advected
        ! horizontally
        if(SubCase.eq.1)then
          xpeak = 1.0_ip*time
          ypeak = 0.0_ip*time
        elseif(SubCase.eq.2)then
          xpeak =-1.0_ip*time
          ypeak = 0.0_ip*time
        elseif(SubCase.eq.3)then
          xpeak = 0.0_ip*time
          ypeak = 1.0_ip*time
        elseif(SubCase.eq.4)then
          xpeak = 0.0_ip*time
          ypeak =-1.0_ip*time
        elseif(SubCase.eq.5)then
          xpeak = 1.0_ip*time
          ypeak = 1.0_ip*time
        elseif(SubCase.eq.6)then
          xpeak =-1.0_ip*time
          ypeak =-1.0_ip*time
        elseif(SubCase.eq.7)then
          xpeak =-1.0_ip*time
          ypeak = 1.0_ip*time
        elseif(SubCase.eq.8)then
          xpeak = 1.0_ip*time
          ypeak =-1.0_ip*time
        endif
        char_len = 0.1
        do n=1,nsmax
          do k=1,nzmax
            do j=1,nymax
              do i=1,nxmax
                  ! This is a 2D cubic spline concentration centered at
                  ! 0,0 with a characteristic width of h=0.1 and
                  ! normalized to unit volume

                if ( (abs(x_cc_pd(i)-xpeak).lt.2.0_ip*char_len) .or.  &
                     (abs(y_cc_pd(j)-ypeak).lt.2.0_ip*char_len) ) then

       Dist = sqrt(  (x_cc_pd(i)-xpeak)**2.0_ip + (y_cc_pd(j)-ypeak)**2.0_ip )
       r = Dist/char_len
       fac_1 = Const1_2d/(char_len*char_len)
       if (r>2.0_ip) then
         truesol(i,j) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         truesol(i,j) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         truesol(i,j) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif

                else
                  truesol(i,j) = 0.0_ip
                endif

                ! Zero nearly zero concentration values since these cause problems
                ! when written out.
                if (abs(concen_pd(i,j,k,n,ts1)).lt.1.0e-50_ip)then
                  concen_pd(i,j,k,n,ts1) = 0.0_ip
                endif
                err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)

                L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

              enddo ! loop over i
            enddo ! loop over j
          enddo ! loop over k
        enddo ! loop over n
        ! Now account for mass conservation at boundaries
        ! HFS Fix this to be in terms of kappa_pd(i,j,k)
        do n=1,nsmax
          do i=1,nxmax
            do j=1,nymax
              MassConsError = MassConsError + outflow_xy1_pd(i,j,n)*dx*dy*dz_vec_pd(0) &
                                            + outflow_xy2_pd(i,j,n)*dx*dy*dz_vec_pd(nzmax+1)
            enddo
            do k=1,nzmax
              MassConsError = MassConsError + outflow_xz1_pd(i,k,n)*dx*dy*dz_vec_pd(k) &
                                            + outflow_xz2_pd(i,k,n)*dx*dy*dz_vec_pd(k)
            enddo
          enddo

          do j=1,nymax
            do k=1,nzmax
              MassConsError = MassConsError + outflow_yz1_pd(j,k,n)*dx*dy*dz_vec_pd(k) &
                                            + outflow_yz2_pd(j,k,n)*dx*dy*dz_vec_pd(k)
            enddo
          enddo
        enddo

        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror /  nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        if(time.lt.0.6)then
          write(ofile1,'(a7,i1,a12)')'TC1_Sub',SubCase,'_St1_err.dat'
          write(ofile2,'(a7,i1,a12)')'TC1_Sub',SubCase,'_St1_sol.dat'
        else
          write(ofile1,'(a7,i1,a12)')'TC1_Sub',SubCase,'_St2_err.dat'
          write(ofile2,'(a7,i1,a12)')'TC1_Sub',SubCase,'_St2_sol.dat'
        endif
        open(200,file=ofile1,status='replace')
        open(201,file=ofile2,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError
        do i=1,nxmax
          do j=1,nymax
            write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
          enddo
        enddo
        close(200)
        close(201)
        endif
      endif

      if(TestCase.eq.2)then
        ! Setting up concentration for 1D pulse to be advected
        ! vertically
       
        if(SubCase.eq.1.or.SubCase.eq.3)then
          zpeak =1.0_ip+1.0_ip*time
        elseif(SubCase.eq.2.or.SubCase.eq.4)then
          zpeak =1.0_ip-1.0_ip*time
        endif

        char_len = 0.1
        do n=1,nsmax
          do i=1,nxmax
            do j=1,nymax
              do k=1,nzmax
                  ! This is a 1D cubic spline concentration centered at
                  ! 0,0 with a characteristic width of h=0.1 and
                  ! normalized to unit volume

                 Dist = abs(z_cc_pd(k)-zpeak)
                if ( Dist.lt.2.0_ip*char_len) then
       r = Dist/char_len
       fac_1 = Const1_1d/(char_len)
       if (r>2.0_ip) then
         truesolz(k) = 0.0_ip
       elseif (r>1.0_ip) then
         temp1 = 2.0_ip-r
         truesolz(k) = fac_1*((temp1*temp1*temp1)/6.0_ip)
       else
         truesolz(k) = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
       endif
                else
                  truesolz(k) = 0.0_ip
                endif

                ! Zero nearly zero concentration values since these cause problems
                ! when written out.
                if (abs(concen_pd(i,j,k,n,ts1)).lt.1.0e-50_ip)then
                  concen_pd(i,j,k,n,ts1) = 0.0_ip
                endif
                errz(k)=truesolz(k)-concen_pd(10,10,k,1,ts1)

                L1_toterror = L1_toterror + abs(errz(k))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + errz(k)*errz(k)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,1,ts1)*kappa_pd(i,j,k)
              enddo ! loop over k
            enddo
          enddo
        enddo ! loop over n
        ! Now account for mass conservation at boundaries
        do n=1,nsmax
          do i=1,nxmax
            do j=1,nymax
              MassConsError = MassConsError + outflow_xy1_pd(i,j,n)*dx*dy*dz_vec_pd(0) &
                                            + outflow_xy2_pd(i,j,n)*dx*dy*dz_vec_pd(nzmax+1)
            enddo
          enddo
        enddo

        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror /  nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        if(time.lt.0.6)then
          write(ofile1,'(a7,i1,a12)')'TC2_Sub',SubCase,'_St1_err.dat'
          write(ofile2,'(a7,i1,a12)')'TC2_Sub',SubCase,'_St1_sol.dat'
        else
          write(ofile1,'(a7,i1,a12)')'TC2_Sub',SubCase,'_St2_err.dat'
          write(ofile2,'(a7,i1,a12)')'TC2_Sub',SubCase,'_St2_sol.dat'
        endif
        open(200,file=ofile1,status='replace')
        open(201,file=ofile2,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError
        do k=1,nzmax
          write(201,'(3g20.10)')truesolz(k),concen_pd(10,10,k,1,ts1),errz(k)
        enddo
        close(200)
        close(201)
      endif

      if(TestCase.eq.3)then
        ! Setting up concentration for 2D cone/box to be advected
        ! horizontally
        if(IsLatLon)then
          do n=1,nsmax
            do k=1,nzmax
              do j=1,nymax
                do i=1,nxmax
                    ! This concentration is described in Example 20.1 of
                    ! LeVeque's Finite Volume book (p.460)
                  lon_shifted = lon_cc_pd(i) - 360.0_ip
                  if ((lon_shifted.lt.24.0_ip).and.(lon_shifted.gt.8.0_ip).and.&
                      (lat_cc_pd(j).gt.-8.0_ip).and.(lat_cc_pd(j).lt.8.0_ip))then
                    truesol(i,j) = 1.0_ip
                  else
                    truesol(i,j) = 0.0_ip
                  endif
                      ! Set up "cone" at lon = -16 with radius of 8deg
                  r = sqrt((lon_shifted+16.0_ip)**2 + lat_cc_pd(j)**2)
                  if (r .lt. 8.0_ip) then
                    truesol(i,j) = 1.0_ip - r/8.0_ip
                  endif

                  err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)

                  L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                  L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                  MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n

          L1_toterror = L1_toterror / (nsmax * TotalVol)
          L2_toterror = L2_toterror / nsmax
          L2_toterror = sqrt(L2_toterror) / TotalVol
          write(global_info,*)"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
          write(global_info,*)" Original total mass = ",total_mass
          write(global_info,*)" Calculated mass = ",MassConsError
          MassConsError = abs((MassConsError-total_mass)/total_mass)
          write(global_info,*)" Mass ratio = ",MassConsError
          write(global_info,*)"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"

          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
            write(outlog(io),*)" Original total mass = ",total_mass
            write(outlog(io),*)" Calculated mass = ",MassConsError
            write(outlog(io),*)" Mass ratio = ",MassConsError
            write(outlog(io),*)"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
          endif;enddo
          write(ofile1,'(a14)')'TC3_LL_err.dat'
          open(200,file=ofile1,status='replace')
          write(200,*)L1_toterror,L2_toterror,MassConsError
          close(200)

          write(ofile2,'(a14)')'TC3_LL_sol.dat'
          open(201,file=ofile2,status='replace')
          do i=1,nxmax
            do j=1,nymax
              write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
            enddo
          enddo
          close(201)
        else ! IsLatLon :: this is the xy case
          do n=1,nsmax
            do k=1,nzmax
              do j=1,nymax
                do i=1,nxmax
                    ! This concentration is described in Example 20.1 of
                    ! LeVeque's Finite Volume book (p.460)
                  if ((x_cc_pd(i).lt.0.6_ip).and.(x_cc_pd(i).gt.0.1_ip).and.&
                      (y_cc_pd(j).gt.-0.25_ip).and.(y_cc_pd(j).lt.0.25_ip))then
                    truesol(i,j) = 1.0_ip
                  else
                    truesol(i,j) = 0.0_ip
                  endif
                  r = sqrt((x_cc_pd(i)+0.45_ip)**2 + y_cc_pd(j)**2)
                  if (r .lt. 0.35_ip) then
                    truesol(i,j) = 1.0_ip - r/0.35_ip
                  endif

                  err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)

                  L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                  L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                  MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n

          L1_toterror = L1_toterror / (nsmax * TotalVol)
          L2_toterror = L2_toterror /  nsmax
          L2_toterror = sqrt(L2_toterror) / TotalVol
          MassConsError = abs((MassConsError-total_mass)/total_mass)
          write(ofile1,'(a11)')'TC3_err.dat'
          write(ofile2,'(a11)')'TC3_sol.dat'
          open(200,file=ofile1,status='replace')
          open(201,file=ofile2,status='replace')
          write(200,*)L1_toterror,L2_toterror,MassConsError
          do i=1,nxmax
            do j=1,nymax
              write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
            enddo
          enddo
          close(200)
          close(201)
        endif
      endif ! TestCase.eq.3

      if(TestCase.eq.4)then
        ! This solution is from Carslaw and Jaeger 1959, p88        
        xm = 1.0
        ym = 1.0
        zm = 1.0

        Tc = (TC4_conc_2-TC4_conc_1)*sqrt(TC4_k_2) / &
                (sqrt(TC4_k_2) + sqrt(TC4_k_1))
        do i=1,nxmax
          eta1 = 0.5_ip*(x_cc_pd(i)-xm)/(sqrt(TC4_k_1*time))
          eta2 = 0.5_ip*(x_cc_pd(i)-xm)/(sqrt(TC4_k_2*time))
          if (x_cc_pd(i).lt.xm) then
            tsolx(i) = TC4_conc_1 + Tc * &
                        (1.0_ip + erf(eta1))
          else
            tsolx(i) = TC4_conc_1 + Tc * &
                        (1.0_ip + sqrt(TC4_k_1/TC4_k_2)*erf(eta2))
          endif
        enddo

        do j=1,nymax
          eta1 = 0.5_ip*(y_cc_pd(j)-ym)/(sqrt(TC4_k_1*time))
          eta2 = 0.5_ip*(y_cc_pd(j)-ym)/(sqrt(TC4_k_2*time))
          if (y_cc_pd(j).lt.ym) then
            tsoly(j) = TC4_conc_1 + Tc * &
                        (1.0 + erf(eta1))
          else
            tsoly(j) = TC4_conc_1 + Tc * &
                        (1.0 + sqrt(TC4_k_1/TC4_k_2)*erf(eta2))
          endif
        enddo

        do k=1,nzmax
          eta1 = 0.5_ip*(z_cc_pd(k)-zm)/(sqrt(TC4_k_1*time))
          eta2 = 0.5_ip*(z_cc_pd(k)-zm)/(sqrt(TC4_k_2*time))
          if (z_cc_pd(k).lt.zm) then
            tsolz(k) = TC4_conc_1 + Tc * &
                        (1.0 + erf(eta1))
          else
            tsolz(k) = TC4_conc_1 + Tc * &
                        (1.0 + sqrt(TC4_k_1/TC4_k_2)*erf(eta2))
          endif
        enddo

        do n=1,nsmax
          do k=1,nzmax
            do j=1,nymax
              do i=1,nxmax
                if (SubCase.eq.1.or.SubCase.eq.4) then
                  err3D(i,j,k)=tsolx(i)-concen_pd(i,j,k,n,ts1)
                elseif (SubCase.eq.2.or.SubCase.eq.5) then
                  err3D(i,j,k)=tsoly(j)-concen_pd(i,j,k,n,ts1)
                elseif (SubCase.eq.3.or.SubCase.eq.6) then
                  err3D(i,j,k)=tsolz(k)-concen_pd(i,j,k,n,ts1)
                else
                  write(global_error,*)"Subcase not known.  Stopping program."
                  stop 1
                endif

                L1_toterror = L1_toterror + abs(err3D(i,j,k))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + err3D(i,j,k)*err3D(i,j,k)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

              enddo ! loop over i
            enddo ! loop over j
          enddo ! loop over k
        enddo ! loop over n

        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror /  nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        write(ofile1,'(a11)')'TC4_err.dat'
        write(ofile2,'(a11)')'TC4_sol.dat'
        open(200,file=ofile1,status='replace')
        open(201,file=ofile2,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError

        lx=ceiling(nxmax*0.5)
        ly=ceiling(nymax*0.5)
        lz=ceiling(nzmax*0.5)
        if (SubCase.eq.1.or.SubCase.eq.4) then
          do i=1,nxmax
            write(201,'(4f12.7)')x_cc_pd(i),tsolx(i),concen_pd(i,ly,lz,1,ts1),err3D(i,ly,lz)
          enddo
        elseif (SubCase.eq.2.or.SubCase.eq.5) then
          do j=1,nymax
            write(201,'(4f12.7)')y_cc_pd(j),tsoly(j),concen_pd(lx,j,lz,1,ts1),err3D(lx,j,lz)
          enddo
        elseif (SubCase.eq.3.or.SubCase.eq.6) then
          do k=1,nzmax
            write(201,'(4f12.7)')z_cc_pd(k),tsolz(k),concen_pd(lx,ly,k,1,ts1),err3D(lx,ly,k)
          enddo
        endif
        close(200)
        close(201)
        write(global_info,*)"Wrote to files TC4_err.dat and TC4_sol.dat"
      endif

      if(TestCase.eq.5)then
        if(IsLatLon)then
          do n=1,nsmax
            do k=1,nzmax
              do j=1,nymax
                do i=1,nxmax
                 
                    ! This concentration is similar to Example 5.4.4 of
                    ! Durran's Finite Volume book (p.284)
                  lon_shifted = lon_cc_pd(i) - 360.0_ip
                  r = sqrt((lon_shifted+16.0_ip)**2.0_ip + &
                           (lat_cc_pd(j))**2.0_ip) / 16.0_ip
                  r = min(1.0_ip,r)
                  truesol(i,j) = 0.5_ip*(1.0_ip+cos(PI*r))
 
                  err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)
                  
                  L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                  L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                  MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

                enddo ! loop over i
              enddo ! loop over j
            enddo ! loop over k
          enddo ! loop over n
          L1_toterror = L1_toterror / (nsmax * TotalVol)
          L2_toterror = L2_toterror / nsmax
          L2_toterror = sqrt(L2_toterror) / TotalVol
          MassConsError = abs((MassConsError-total_mass)/total_mass)
          if(time.lt.3.0)then
            write(ofile1,'(a18)')'TC5_LL_mid_err.dat'
            write(ofile2,'(a18)')'TC5_LL_mid_sol.dat'
          else
            write(ofile1,'(a14)')'TC5_LL_err.dat'
            write(ofile2,'(a14)')'TC5_LL_sol.dat'
          endif
          open(200,file=ofile1,status='replace')
          open(201,file=ofile2,status='replace')
          write(200,*)L1_toterror,L2_toterror,MassConsError
          do i=1,nxmax
            do j=1,nymax
              write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
            enddo
          enddo
          close(200)
          close(201)

        else
        ! Setting up concentration pulse for circular shear advection
        do n=1,nsmax
          do k=1,nzmax
            do j=1,nymax
              do i=1,nxmax
                  ! This concentration is described in Example 5.4.4 of
                  ! Durran's Finite Volume book (p.284)
                r = 4.0_ip*sqrt((x_cc_pd(i)-0.25_ip)**2.0_ip + &
                                (y_cc_pd(j)-0.25_ip)**2.0_ip)
                r = min(1.0_ip,r)

                truesol(i,j) = 0.5_ip*(1.0_ip+cos(PI*r))

                err(i,j)=truesol(i,j)-concen_pd(i,j,k,n,ts1)

                L1_toterror = L1_toterror + abs(err(i,j))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + err(i,j)*err(i,j)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

              enddo ! loop over i
            enddo ! loop over j
          enddo ! loop over k
        enddo ! loop over n
        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror / nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        write(ofile1,'(a11)')'TC5_err.dat'
        write(ofile2,'(a11)')'TC5_sol.dat'
        open(200,file=ofile1,status='replace')
        open(201,file=ofile2,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError
        do i=1,nxmax
          do j=1,nymax
            write(201,'(3f12.7)')truesol(i,j),concen_pd(i,j,1,1,ts1),err(i,j)
          enddo
        enddo
        close(200)
        close(201)
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Wrote to files TC5_err.dat TC5_sol.dat"
        endif;enddo
        endif

      endif

      if(TestCase.eq.6)then
        do n=1,nsmax
          do k=1,nzmax
            do j=1,nymax
              do i=1,nxmax
                tsol = MMS_TrueSol(x_cc_pd(i),y_cc_pd(j),z_cc_pd(k),time)
                err3D(i,j,k)=tsol-concen_pd(i,j,k,n,ts1)

                L1_toterror = L1_toterror + abs(err3D(i,j,k))*kappa_pd(i,j,k)
                L2_toterror = L2_toterror + err3D(i,j,k)*err3D(i,j,k)*kappa_pd(i,j,k)
                MassConsError = MassConsError + concen_pd(i,j,k,n,ts1)*kappa_pd(i,j,k)

              enddo ! loop over i
            enddo ! loop over j
          enddo ! loop over k
        enddo ! loop over n
        L1_toterror = L1_toterror / (nsmax * TotalVol)
        L2_toterror = L2_toterror / nsmax
        L2_toterror = sqrt(L2_toterror) / TotalVol
        MassConsError = abs((MassConsError-total_mass)/total_mass)
        write(ofile1,'(a11)')'TC6_err.dat'
        open(200,file=ofile1,status='replace')
        write(200,*)L1_toterror,L2_toterror,MassConsError
        close(200)
      endif

      end subroutine Testcase_CalcErrors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function MMS_TrueSol(x,y,z,t)

      use global_param, only : &
        HR_2_S,KM_2_M,KM3_2_M3

      implicit none

      real(kind=ip) :: MMS_TrueSol
      real(kind=ip) :: x
      real(kind=ip) :: y
      real(kind=ip) :: z
      real(kind=ip) :: t

      real(kind=ip) :: zeta
      real(kind=ip) :: xcc,ycc,zcc
      !real(kind=ip) :: expzcc, expzccm
      real(kind=ip) :: sechz
      real(kind=ip) :: y_qxylen, sechy
      real(kind=ip) :: x_qxylen, sechx

            ! Here's the chosen solution:
            !   q = MMS_Q0 * sech(xcc/qxylen) *
            !                sech(ycc/qxylen) * 
            !                zcc * exp(-zcc/zeta);
            !      where zeta = MMS_W0 * time*HR_2_S + 1000.0_ip

      if(MMS_USE_T)then
        zeta = MMS_W0 * t * HR_2_S + MMS_zeta0
      else
        zeta = MMS_zeta0
      endif
      zcc = z * KM_2_M
      ycc = y * KM_2_M
      xcc = x * KM_2_M

      !expzcc = exp(zcc/zeta)
      !expzccm= 1.0_ip/expzcc
      y_qxylen = ycc/MMS_qxylen
      sechy = 1.0_ip/cosh(y_qxylen)
      x_qxylen = xcc/MMS_qxylen
      sechx = 1.0_ip/cosh(x_qxylen)

!      MMS_TrueSol = MMS_Q0
!      if(MMS_USE_X) MMS_TrueSol = MMS_TrueSol * sechx
!      if(MMS_USE_Y) MMS_TrueSol = MMS_TrueSol * sechy
!      if(MMS_USE_Z) MMS_TrueSol = MMS_TrueSol * zcc/MMS_ztrop * expzccm
!
!      MMS_TrueSol = MMS_TrueSol / KM3_2_M3


            ! Here's the chosen solution:
            !   q = MMS_Q0 * sech(xcc/qxylen) *
            !                sech(ycc/qxylen) * 
            !                sech(zcc/zeta);
            !      where zeta = MMS_W0 * time*HR_2_S + 1000.0_ip
      sechz = 1.0_ip/cosh(zcc/zeta)

      MMS_TrueSol = MMS_Q0
      if(MMS_USE_X) MMS_TrueSol = MMS_TrueSol * sechx
      if(MMS_USE_Y) MMS_TrueSol = MMS_TrueSol * sechy
      if(MMS_USE_Z) MMS_TrueSol = MMS_TrueSol * sechz

      MMS_TrueSol = MMS_TrueSol / KM3_2_M3

      end function MMS_TrueSol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function MMS_Source(x,y,z,t,velx,vely,velz,vels,dws_dz)

      use global_param, only : &
        MPS_2_KMPHR,MPS_2_KMPHR,KM_2_M,HR_2_S,KM3_2_M3

      use Diffusion,    only : &
        diffusivity_horz

      implicit none

      real(kind=ip) :: MMS_Source
      real(kind=ip) :: x
      real(kind=ip) :: y
      real(kind=ip) :: z
      real(kind=ip) :: t
      real(kind=ip) :: velx
      real(kind=ip) :: vely
      real(kind=ip) :: velz
      real(kind=ip) :: vels
      real(kind=ip) :: dws_dz

      real(kind=ip) :: true_sol
      real(kind=ip) :: kap
      real(kind=ip) :: zeta
      real(kind=ip) :: xcc,ycc,zcc
      real(kind=ip) :: x_qxylen,y_qxylen
      real(kind=ip) :: tanhx,tanhy,tanhz
      real(kind=ip) :: sechx,sechy,sechz
      !real(kind=ip) :: x_vzlen,y_vzlen
      real(kind=ip) :: expzcc,expzccm
      real(kind=ip) :: Dq_dt,Dq_dx,Dq_dy,Dq_dz
      real(kind=ip) :: D2q_dx2,D2q_dy2,D2q_dz2
      real(kind=ip) :: Lap_q

            ! Source for the MMS solution is obtained by plugging the
            ! chosen solution into the governing equation
            ! Here's the chosen solution:
            !   q = MMS_Q0 * sech(xcc/qxylen) *
            !                sech(ycc/qxylen) * 
            !                zcc * exp(-zcc/zeta);
            !      where zeta = MMS_W0 * time*HR_2_S + 1000.0_ip
            !
            ! And the governing equation:
            ! Dq_dt + Vx Dq_dx + Vy Dq_dy + (Vz+Vs) Dq_dz + q dws_dz
            !      -  Div.(kappa Grad q)  = Q(x,y,z,t)

      ! First get true solution at this point and time
      true_sol = MMS_TrueSol(x,y,z,t)
      true_sol = true_sol * KM3_2_M3

      dws_dz = dws_dz/MPS_2_KMPHR/KM_2_M

      if(.not.MMS_USE_DifF)diffusivity_horz=0.0_ip
        ! Convert diffusivity_horz back to m2/s
      kap  = diffusivity_horz * KM_2_M * KM_2_M / HR_2_S
      if(MMS_USE_T)then
        zeta = MMS_W0 * t             * HR_2_S + MMS_zeta0
        !zeta = MMS_W0 * (t-0.5_ip*dt) * HR_2_S + MMS_zeta0
      else
        zeta = MMS_zeta0
      endif
      ! Setting up MMS
          ! NOTE: when the source is transient, use RK-2 (Eq 17.40 of
          ! LeV) to integrate instead of FE
      zcc = z * KM_2_M
      expzcc = exp(zcc/zeta)
      expzccm= 1.0_ip/expzcc
      sechz = 1.0_ip/cosh(zcc/zeta)
      tanhz =        tanh(zcc/zeta)
      ycc =  y * KM_2_M
      y_qxylen = ycc/MMS_qxylen
      !y_vzlen  = ycc/MMS_vzlen
      sechy = 1.0_ip/cosh(y_qxylen)
      tanhy =        tanh(y_qxylen)
      xcc =  x * KM_2_M
      x_qxylen = xcc/MMS_qxylen
      !x_vzlen  = xcc/MMS_vzlen
      sechx = 1.0_ip/cosh(x_qxylen)
      tanhx =        tanh(x_qxylen)

      ! Get Derivatives for solution #1
      Dq_dt   = true_sol * (zcc*MMS_W0/(zeta*zeta))
      Dq_dx   = true_sol * (-tanhx/MMS_qxylen)
      D2q_dx2 = true_sol * (tanhx*tanhx - sechx*sechx) / &
                               (MMS_qxylen  * MMS_qxylen)
      Dq_dy   = true_sol * (-tanhy/MMS_qxylen)
      D2q_dy2 = true_sol * (tanhy*tanhy - sechy*sechy) / &
                               (MMS_qxylen  * MMS_qxylen)
      Dq_dz   = true_sol * (1.0_ip/zcc - 1.0_ip/zeta)
      D2q_dz2 = true_sol * (1.0_ip/(zeta*zeta) - 1.0_ip/(zeta*zcc))

      ! Get Derivatives for solution #2
      Dq_dt   = true_sol * (tanhz*zcc*MMS_W0/(zeta*zeta))
      Dq_dx   = true_sol * (-tanhx/MMS_qxylen)
      D2q_dx2 = true_sol * (tanhx*tanhx - sechx*sechx) / &
                               (MMS_qxylen  * MMS_qxylen)
      Dq_dy   = true_sol * (-tanhy/MMS_qxylen)
      D2q_dy2 = true_sol * (tanhy*tanhy - sechy*sechy) / &
                               (MMS_qxylen  * MMS_qxylen)
      Dq_dz   = true_sol * (-tanhz/zeta)
      D2q_dz2 = true_sol * ((tanhz*tanhz - sechz*sechz)/(zeta*zeta))

      Lap_q = D2q_dx2 + D2q_dy2 + D2q_dz2

      if(.not.MMS_USE_X)Dq_dx = 0.0_ip
      if(.not.MMS_USE_Y)Dq_dy = 0.0_ip
      if(.not.MMS_USE_Z)Dq_dz = 0.0_ip
      if(.not.MMS_USE_Z)Dq_dt = 0.0_ip
      if(.not.MMS_USE_T)Dq_dt = 0.0_ip
      !if(.not.MMS_USE_DifF)then
      !  D2q_dx2 = 0.0_ip
      !  D2q_dy2 = 0.0_ip
      !  D2q_dz2 = 0.0_ip
      !endif

      ! Assemble source term in kg/m3/s
      MMS_Source = 0.0_ip
      if(MMS_USE_X) MMS_Source = MMS_Source + &
                              (velx/MPS_2_KMPHR) * Dq_dx
      if(MMS_USE_Y) MMS_Source = MMS_Source + &
                              (vely/MPS_2_KMPHR) * Dq_dy
      if(MMS_USE_Z) MMS_Source = MMS_Source + &
                              ((vels+velz)/MPS_2_KMPHR) * Dq_dz    + &
                              dws_dz                       * true_sol + &
                              Dq_dt

      if(MMS_USE_DifF) MMS_Source = MMS_Source - kap * Lap_q

           ! Convert from kg/m3/s to kg/km3/hr
      MMS_Source = MMS_Source * HR_2_S / KM3_2_M3

      end function MMS_Source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine back_rotate(lon,lat,    &
                             lon_pole,lat_pole,deg_rot, &
                             lon_br,lat_br)

      use global_param, only : &
        DEG2RAD

      implicit none

      real(kind=ip) :: lon,lat
      real(kind=ip) :: lon_pole,lat_pole
      real(kind=ip) :: deg_rot
      real(kind=ip) :: lon_br,lat_br

      real(kind=ip),dimension(3,3) :: Rot
      real(kind=ip),dimension(3)   :: xvec,xnew

      real(kind=ip) :: clm,clmr,comg,cph,cphr
      real(kind=ip) :: lam,lamr,omega
      real(kind=ip) :: phi,phir
      real(kind=ip) :: slm,slmr,somg,sph,sphr
      real(kind=ip) :: urx,ury,urz
      real(kind=ip) :: x_usph,y_usph,z_usph
      real(kind=ip) :: xcolat,xlm,xph

      ! back_rotate rotates the point lon,lat by an angle deg_rot about
      ! the axis lon_pole,lat_pole
      ! returns the lon,lat coordinates of the rotated point

      ! Get details of position vector
      phi  = lat*DEG2RAD
      lam  = lon*DEG2RAD
      slm = sin(lam)
      clm = cos(lam)
      sph = sin(phi)
      cph = cos(phi)

      ! Get details of rotation axis
      phir = lat_pole*DEG2RAD
      lamr = lon_pole*DEG2RAD
      slmr = sin(lamr)
      clmr = cos(lamr)
      sphr = sin(phir)
      cphr = cos(phir)
      urx = clmr*cphr
      ury = slmr*cphr
      urz = sphr

      omega = deg_rot * DEG2RAD
      somg = sin(omega)
      comg = cos(omega)

      !% build rotation matrix
      !% (see http://en.wikipedia.org/wiki/Rotation_matrix : Rotation matrix
      !given an axis and an angle)
      Rot(1,1) = urx*urx + (1.0_ip-urx*urx)*comg
      Rot(1,2) = urx*ury * (1.0_ip-comg) - urz*somg
      Rot(1,3) = urx*urz * (1.0_ip-comg) + ury*somg
      Rot(2,1) = urx*ury * (1.0_ip-comg) + urz*somg
      Rot(2,2) = ury*ury + (1.0_ip-ury*ury)*comg
      Rot(2,3) = ury*urz * (1.0_ip-comg) - urx*somg
      Rot(3,1) = urx*urz * (1.0_ip-comg) - ury*somg
      Rot(3,2) = ury*urz * (1.0_ip-comg) + urx*somg
      Rot(3,3) = urz*urz + (1.0_ip-urz*urz)*comg

      !  % Get unit sphere coordinates in standard orientation
      x_usph = clm*cph
      y_usph = slm*cph
      z_usph = sph
      !  % Back-rotate coordinates
      xvec(1) = x_usph
      xvec(2) = y_usph
      xvec(3) = z_usph
      xnew = matmul(Rot,xvec)
      !  % Get lon/lat for this new position
      xcolat = acos(xnew(3))
      xph    = (90.0_ip*DEG2RAD)-xcolat
      xlm    = atan2(xnew(2),xnew(1))

      lon_br = xlm/DEG2RAD
      lat_br = xph/DEG2RAD


      end subroutine back_rotate

      end module TestCases

