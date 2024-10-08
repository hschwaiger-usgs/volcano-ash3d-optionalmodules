!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Ash3d is a program for modeling volcanic ash transport and dispersion.
!
!  This software is written in Frotran 2003 and is designed for use on a Linux
!  operating system.
!  
!  This software, along with auxillary USGS libraries and related repositories,
!  can be found at https://code.usgs.gov/vsc/ash3d
!
!  Installation instructions are given in the README.md file of this repository.
!  Basic usage instructions are given in doc/UsersGuide.md.
!
!  The program description and numerical methodology employed is described in:
!   Schwaiger, H.F., R.P. Denlinger, and L.G. Mastin, 2012, Ash3d, a finite-
!     volume, conservative numerical model for ash transport and tephra
!     deposition, Journal of Geophysical Research, 117, B04204,
!     doi:10.1029/2011JB008968
!
!  A complete user's guide and reference manual is available at
!
!  The USGS provides a web-interface to this software at:
!    https://vsc-ash.wr.usgs.gov
!  with instructions on web-interface usage provided by:
!   Mastin, L.G., M.J. Randall, H.F. Schwaiger and R.P. Denlinger, 2021, User's
!     Guide and Reference to Ash3d-- A Three-Dimensional Model for Eulerian
!     Atmospheric Tephra Transport and Deposition, USGS Open-File Report
!     2013-1122, doi:10.3133/ofr20131122.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program Ash3d

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useDiffusion,useHorzAdvect,useVertAdvect,&
         HR_2_S,useTemperature,DT_MIN,EPS_TINY,EPS_SMALL,&
         nmods,OPTMOD_names,StopConditions,CheckConditions      

      use mesh,          only : &
         ivent,jvent,nxmax,nymax,nzmax,nsmax,ts0,ts1,ZPADDING,dz_vec_pd,z_cc_pd

      use solution,      only : &
         concen_pd,DepositGranularity,StopValue_FracAshDep,aloft_percent_remaining, &
         SourceCumulativeVol,dep_vol,aloft_vol,outflow_vol,tot_vol,vf_pd

      use Output_Vars,   only : &
         DepositAreaCovered,DepositThickness,LoadVal,CloudLoadArea,&
         Calculated_Cloud_Load,Calculated_AshThickness,Calc_vprofile, &
           Allocate_Output_UserVars, &
           Allocate_NTime,   &
           Allocate_Profile, &
           Gen_Output_Vars,&
           FirstAsh

      use io_data,       only : &
         Called_Gen_Output_Vars,isFinal_TS,LoadConcen,log_step,&
         Output_at_logsteps,Output_at_WriteTimes,Output_every_TS,&
         NextWriteTime,iTimeNext,nvprofiles,nWriteTimes,&
         Write_PT_Data,Write_PR_Data

      use time_data,     only : &
         time,dt,Simtime_in_hours,t0,t1,t2,ntmax,tcount1,tcount2,&
         tcount_rate,tcount_max,tw_tot

      use Ash3d_Program_Control, only : &
           Parse_Command_Line, &
           Set_OS_Env, &
           Read_Control_File

      use Source,        only : &
         SourceNodeFlux,e_EndTime_final,e_Volume,&
         SourceType,Source_in_dt,IsCustom_SourceType, &
           Calc_Normalized_SourceCol,&
           EruptivePulse_MassFluxRate,&
           CheckEruptivePulses,&
           TephraSourceNodes,&
           SourceVolInc

      use Source_Umbrella, only : &
         ibase,itop,SourceNodeFlux_Umbrella, &
           Allocate_Source_Umbrella,&
           TephraSourceNodes_Umbrella,&
           SourceVolInc_Umbrella,&
           AvgCon_Umbrella

      use Tephra,        only : &
         n_gs_max,n_gs_aloft,Tephra_gsdiam,&
           Allocate_Tephra,&
           Allocate_Tephra_Met,&
           Prune_GS
   
      use Atmosphere,    only : &
           Allocate_Atmosphere_Met

      use AdvectionHorz, only : &
           AdvectHorz

      use AdvectionVert_DCU, only : &
           advect_z

      use Diffusion,     only : &
           DiffuseHorz,&
           DiffuseVert

      use Airports,      only : &
         Airport_thickness_TS,Airport_thickness,nairports,&
           ReadAirports

      use Ash3d_ASCII_IO,  only : &
           vprofilewriter

#ifdef USENETCDF
      use Ash3d_Netcdf_IO,only : &
           NC_RestartFile_LoadConcen
#endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert 'use' statements here
!
#ifdef TESTCASES
      use TestCases
#endif
!------------------------------------------------------------------------------

      implicit none

      integer               :: itime
      integer               :: i,k,isize
      real(kind=dp)         :: Interval_Frac
      real(kind=ip)         :: falltime
      logical               :: Load_MesoSteps
      logical               :: StopTimeLoop   = .false.
      real(kind=ip)         :: MassConsErr
      real(kind=dp)         :: dt_TC

      INTERFACE
!        subroutine input_data_ResetParams
!        end subroutine input_data_ResetParams
        subroutine alloc_arrays
        end subroutine alloc_arrays
        subroutine calc_mesh_params
        end subroutine calc_mesh_params
        subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac)
          integer,parameter  :: dp         = 8 ! Double precision
          real(kind=dp),intent(in)    :: TimeNow
          real(kind=dp),intent(out)   :: Interval_Frac
          logical      ,intent(inout) :: Load_MesoSteps
        end subroutine MesoInterpolater
         ! We do need to call Adjust_DT from this top-level file
        subroutine Adjust_DT(mesostep)
          logical, intent(in), optional :: mesostep
        end subroutine
        subroutine output_results
        end subroutine output_results
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
        subroutine TimeStepTotals(itime)
          integer, intent(in) :: itime
        end subroutine TimeStepTotals
        subroutine dealloc_arrays
        end subroutine dealloc_arrays
      END INTERFACE

      ! Start time logging
      call cpu_time(t0) !time is a scaler real
      call system_clock(tcount1,tcount_rate,tcount_max)

      ! First, parse the command line
      call Parse_Command_Line

      ! Before we do anything, get the state of the executable, system, environment and run
      call Set_OS_Env

      aloft_percent_remaining = 1.0_ip
      SourceCumulativeVol     = 0.0_ip
      MassConsErr             = 0.0_ip

#ifdef TESTCASES
      call set_TestCase_globvars
#endif

        ! input data for ash transport
      call Read_Control_File

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to custom input blocks here
!
!  Loop through all the optional modules compiled and test against the list
!  from the input file (e.g. OPTMOD=TOPO), then call the special input reader
!  for that block
!  Do a sanity check on optional module requested in the input file v.s. those
!  compiled in this executable and for consistency among modules:
!  e.g. SRC_RESUSP will require the VARDIFF and LC be set
!
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now looping through optional modules found in input file"
      endif;enddo
      do i=1,nmods
        do io=1,2;if(VB(io).le.verbosity_essential)then
          write(outlog(io),*)"Testing for ",OPTMOD_names(i),i
        endif;enddo
        !if(OPTMOD_names(i).eq.'RESETPARAMS')then
        !  do io=1,2;if(VB(io).le.verbosity_essential)then
        !    write(outlog(io),*)"  Reading input block for RESETPARAMS"
        !  endif;enddo
        !  call input_data_ResetParams
        !endif
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then    
        write(outlog(io),*)"Finished reading all specialized input blocks"
      endif;enddo
!
!------------------------------------------------------------------------------

        ! Read airports/POI and allocate/initilize arrays
        ! We only need to do this if an output variable demands it since this is
        ! a burden every time step
      if(Output_every_TS) &
        call ReadAirports

      call alloc_arrays
        ! Set up grids for solution and Met data
      call calc_mesh_params

      if(((SourceType.eq.'umbrella').or.(SourceType.eq.'umbrella_air')))then
        call Allocate_Source_Umbrella(nxmax,nymax,nzmax)
      endif
      if(.not.IsCustom_SourceType)then
        call Calc_Normalized_SourceCol
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize concen and any special source terms here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(LoadConcen)then
        ! We are initializing the concentration and time from an output file
        ! Currently, Ash3d assumes the concentration file is compatible with
        ! the computational grid and grainsize distribution
#ifdef USENETCDF
        call NC_RestartFile_LoadConcen
#else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Loading concentration files requires previous netcdf"
          write(errlog(io),*)"       output.  This Ash3d executable was not compiled with"
          write(errlog(io),*)"       netcdf support.  Please recompile Ash3d with"
          write(errlog(io),*)"       USENETCDF=T, or select another source."
        endif;enddo
        stop 1
#endif
      else
        ! Initialize arrays if we haven't already loaded the concentration from
        ! a previous run
        concen_pd = 0.0_ip
        DepositGranularity = 0.0_ip
      endif
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert special source terms here
!
!------------------------------------------------------------------------------

      if(useTemperature)then
        call Allocate_Atmosphere_Met
#ifdef TESTCASES
        if(TestCase.eq.6)then
          ! Setup up full atmosphere here
          call Set_MMS_Atmos
          CheckConditions(1) = .false. ! less than 1% of ash aloft
          CheckConditions(2) = .true.  ! time.ge.Simtime_in_hours
          CheckConditions(3) = .false. ! n_gs_aloft.eq.0
          CheckConditions(4) = .false. ! MassConsErr.gt.1.0e-3_ip
          CheckConditions(5) = .false. ! if any volume measure is negative
        endif
#endif
      endif
      ! Now we can allocate the variables that live on the met grid
      if(useCalcFallVel)then
        call Allocate_Tephra_Met
      endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional variable allocation subroutines here
!
!#ifdef TESTCASES
!#endif
!------------------------------------------------------------------------------
      ! Allocate all the output variables
      call Allocate_Output_UserVars(nxmax,nymax,nzmax,nsmax)

      ! Now that we have the Met grids initialized, get the state variables
      ! interpolated on the start time
      time           = 0.0_ip
      Load_MesoSteps = .true.
      Interval_Frac  = 0.0_8
!****
!****
!      call MesoInterpolater(time , Load_MesoSteps , Interval_Frac)
!****

      ! Calculate the fall time of each grain size
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5020)
      endif;enddo
      do isize = 1,n_gs_max
        falltime = 0.0_ip
        do k = nzmax,1,-1
          if(z_cc_pd(k)+0.5_ip*dz_vec_pd(k).lt.z_cc_pd(nzmax)/ZPADDING)then
            if(abs(vf_pd(ivent,jvent,k,isize)).gt.EPS_SMALL)then
              falltime = falltime - dz_vec_pd(k)/vf_pd(ivent,jvent,k,isize)
            else
              falltime = 0.0_ip
            endif
          endif
        enddo
        do io=1,2;if(VB(io).le.verbosity_info)then
          if(falltime.lt.EPS_SMALL)then
            write(outlog(io),*)isize," Tracer particle; no appreciable fall velocity."
          else
            write(outlog(io),5021)isize,Tephra_gsdiam(isize)*1000.0_ip,falltime
          endif
        endif;enddo
      enddo

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
#ifdef TESTCASES
        call set_TestCase_windfield
        call Adjust_DT(.false.)
        ! Log the initial dt for test cases.  We need this for TC5 where the
        ! wind speed drops to zero, then reverses
        dt_TC = dt
#endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to prep user-specified output
!
!------------------------------------------------------------------------------

        ! Call output_results before time loop to create output files
      call output_results

      ntmax = max(1,4*int(Simtime_in_hours/dt))
      call Allocate_NTime(ntmax)
      if (Write_PR_Data)then
        call Allocate_Profile(nzmax,ntmax,nvprofiles)
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5007)
      endif;enddo

      ! Calculate mass flux and end times of each eruptive pulse
      call EruptivePulse_MassFluxRate
#ifdef TESTCASES
        call DistSource
#endif

      ! Write out starting volume, max time steps, and headers for the table that follows
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5001) tot_vol,ntmax
      endif;enddo

      ! Get the cpu time for the start of the time loop
      call cpu_time(t1) !time is a scaler real

      ! ************************************************************************
      ! ****** begin time simulation *******************************************
      ! ************************************************************************
      itime = 0
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Starting time loop."
      endif;enddo

      do while (StopTimeLoop.eqv..false.)
        ! Note: stop conditions are evaluated at the end of the time loop

        ! Copy previous t=n+1 to t=n
        concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)
          ! re-initialize slice (1)
        concen_pd(:,:,:,1:nsmax,ts1) = 0.0_ip

        Called_Gen_Output_Vars  = .false.
        Calculated_Cloud_Load   = .false.
        Calculated_AshThickness = .false.

        itime = itime + 1
        if(itime.gt.ntmax)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: The number of time steps attempted exceeds 3x that anticipated."
            write(outlog(io),*)"         Check that the winds are stable"
            write(outlog(io),*)"        Simtime_in_hours = ",Simtime_in_hours
            write(outlog(io),*)"                   ntmax = ",ntmax
            write(outlog(io),*)"            current step = ",itime
          endif;enddo
        endif

          ! find the wind field at the current time
!****
!****
!        call MesoInterpolater(time , Load_MesoSteps , Interval_Frac)
!****
!****
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
#ifdef TESTCASES
        call set_TestCase_windfield
        call Adjust_DT(.false.)
        dt = min(dt_TC,dt)
        if(TestCase.eq.6)then
          call DistSource
        endif
#endif
!------------------------------------------------------------------------------

          ! Determine if (and which) eruptive pulses are active in the current dt
        call CheckEruptivePulses

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to specialized MassFluxRate calculations here
!
!------------------------------------------------------------------------------

        ! Add source term
        if(Source_in_dt) then
          ! Check if the source type is one of the standard types with a 1-node column
          if ((SourceType.eq.'point')  .or. &
              (SourceType.eq.'line')   .or. &
              (SourceType.eq.'profile').or. &
              (SourceType.eq.'suzuki'))then

            ! Calculating the flux into the vent column
            call TephraSourceNodes

            ! Most standard source types (point, line, profile, suzuki) are
            ! integrated as follows.
            concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) =  &
                concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) +  &
                  real(dt,kind=ip) * &
                  SourceNodeFlux(1:nzmax+1,1:n_gs_max)

            ! Keep track of the accumulated source inserted for mass conservation error-checking
            SourceCumulativeVol = SourceCumulativeVol + SourceVolInc(dt)

          elseif (SourceType.eq.'umbrella'.or. &
                 (SourceType.eq.'umbrella_air')) then
            ! Umbrella clouds have a special source insertion with a 3x3 column

            ! Umbrella sources still need the Suzuki distribution of mass above the vent
            ! Calculating the flux into the vent column
            call TephraSourceNodes

            ! Now modify the above result to distribute the total SourceNodeFlux
            ! into the umbrella form
            call TephraSourceNodes_Umbrella

            ! Before source insertion, we smooth over the concentration over the
            ! 3x3 patch within the umbrella zone using the weighted averaging stencil
            do k=ibase,itop
              do isize=1,n_gs_max
                concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,k,isize,ts0) = &
                  AvgCon_Umbrella(concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,k,isize,ts0),k)
              enddo
            enddo
            ! Here the integration is the same as for standard sources (point, line, profile, suzuki),
            ! except we use a 3x3 column around ivent,jvent instead of just a single-node column.
            concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,1:nzmax+1,1:n_gs_max,ts0) =  &
                concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,1:nzmax+1,1:n_gs_max,ts0)  +  &
                  real(dt,kind=ip) * &
                  SourceNodeFlux_Umbrella(1:3,1:3,1:nzmax+1,1:n_gs_max)

            ! Keep track of the accumulated source inserted for mass conservation error-checking
            SourceCumulativeVol = SourceCumulativeVol + SourceVolInc_Umbrella(dt)

          else
            ! This is not a standard source.
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: source type is non-standard"
            endif;enddo
            stop 1
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional sources here
!         These subroutines need to calculate the mass of a species inserted
!         into specific cells and update concen accordingly
!
!          elseif () then
!
!------------------------------------------------------------------------------
          endif
        endif ! MassFluxRate_dt1.gt.0.0_ip

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set Boundary Conditions
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional boundary conditions here
!
!------------------------------------------------------------------------------
        call Set_BC(1)
#ifdef TESTCASES
        if(TestCase.eq.6)then
          call Set_MMS_BC
        endif
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Advection / Diffusion / Deposition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional advection/diffusion routines here
!
!------------------------------------------------------------------------------

        if(useHorzAdvect) call AdvectHorz(itime)

        if(useVertAdvect) call advect_z

        if(useDiffusion)then
          call Set_BC(2)
#ifdef TESTCASES
        if(TestCase.eq.6)then
          call Set_MMS_BC
        endif
#endif
          call DiffuseVert
          call DiffuseHorz(itime)
        endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional deposition routines here
!
!------------------------------------------------------------------------------

        ! Advance time
        time = time + dt

        ! If there is any time-series output requiring evaluation at each
        ! time-step, then extract output variables from concen here
        if(Output_every_TS)then
          call Gen_Output_Vars

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every timestep) here
!
!------------------------------------------------------------------------------

            ! See whether the ash has hit any airports/POI
          call FirstAsh

            ! Track ash on vertical profiles
          if (Write_PR_Data)then
            call Calc_vprofile(itime)
            call vprofilewriter(itime)     !write out vertical profiles
          endif
        endif

        ! Go to output results if we're at the next output stage
        ! Note that dt was set in Adjust_DT so that it is no larger than
        ! DT_MIN, but may be adjusted down so as to land on the next
        ! output time.  time has already been integrated forward so
        ! NextWriteTime-time should be near zero for output steps.
        if(Output_at_WriteTimes.and.(abs(NextWriteTime-time).lt.DT_MIN))then
            ! Generate output variables if we haven't already
          if(.not.Called_Gen_Output_Vars)then
            call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every output-step) here
!
!------------------------------------------------------------------------------
          endif
          call output_results
#ifdef TESTCASES
          call Testcase_CalcErrors
#endif
          !if ((WriteAirportFile_ASCII.or.WriteAirportFile_KML).and. &
          if (Write_PT_Data.and. &
              (iTimeNext.lt.nWriteTimes)) then
            do i=iTimeNext,nWriteTimes
              Airport_Thickness_TS(1:nairports,i) = Airport_Thickness(1:nairports)
            enddo
          endif
        endif

        if(Output_at_logsteps)then  
          ! Write summary information on mass conservation every log_step time steps
          if(mod(itime,log_step).eq.0) then
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every log-step) here
!
!------------------------------------------------------------------------------
            endif
            call TimeStepTotals(itime)
          endif

            ! Only consider reducing GS bins at logsteps
          if(time.gt.e_EndTime_final)then
            ! Doesn't make sense to flag bins as flushed out while the eruption is on-going
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
            endif
            call Prune_GS
          endif
        else
            ! If we are not monitoring deposits through logsteps, then set
            ! tot_vol to 0 and rely on the simulation ending through input
            ! duration values
          tot_vol = 0.0_ip
        endif

        if(tot_vol.gt.EPS_SMALL)then
          aloft_percent_remaining = aloft_vol/tot_vol
        else
          aloft_percent_remaining = 1.0_ip
        endif

        ! Check stop conditions
        !  If any of these is true, then the time loop will stop
           ! Stops if there is less than 1% of ash aloft in the domain
        StopConditions(1) = (aloft_percent_remaining.lt.(1.0_ip-StopValue_FracAshDep))
           ! Normal stop condition if simulation exceeds alloted time
        StopConditions(2) = (time.ge.Simtime_in_hours)
           ! Normal stop condition when nothing is left to advect
        StopConditions(3) = (n_gs_aloft.eq.0)
        if(SourceCumulativeVol.gt.EPS_TINY)then
          MassConsErr = abs(SourceCumulativeVol-tot_vol)/SourceCumulativeVol
        endif
           ! Error stop condition if the concen and outflow do not match the source,
           ! but only trigger this condition if not a restart case (until outflow is tracked)
        if(.not.LoadConcen) &
          StopConditions(4) = (MassConsErr.gt.1.0e-3_ip)
           ! Error stop condition if any volume measure is negative
        StopConditions(5) = (dep_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (aloft_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (outflow_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (SourceCumulativeVol.lt.-1.0_ip*EPS_SMALL)

        if((CheckConditions(1).eqv..true.).and.&
           (StopConditions(1).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(2).eqv..true.).and.&
               (StopConditions(2).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(3).eqv..true.).and.&
               (StopConditions(3).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(4).eqv..true.).and.&
               (StopConditions(4).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(5).eqv..true.).and.&
               (StopConditions(5).eqv..true.))then
          StopTimeLoop = .true.
        else
          StopTimeLoop = .false.
        endif
      enddo  !loop over itime
              !  ((dep_percent_accumulated.le.StopValue_FracAshDep).and. &
              !    (time.lt.Simtime_in_hours)        .and. &
              !    (n_gs_aloft.gt.0))

      ! Reset ntmax to the actual number of time steps
      ntmax = itime

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Time integration completed for the following reason:"
      endif;enddo
      if((CheckConditions(1).eqv..true.).and.&
         (StopConditions(1).eqv..true.))then
        ! Normal stop condition set by user tracking the deposit
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a35,f8.3)')"Percent accumulated/exited exceeds ",StopValue_FracAshDep
        endif;enddo
      endif
      if((CheckConditions(2).eqv..true.).and.&
         (StopConditions(2).eqv..true.))then
        ! Normal stop condition if simulation exceeds alloted time
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a24)')"time.ge.Simtime_in_hours"
          write(outlog(io),'(a21,f15.3)')"              Time = ",time
          write(outlog(io),'(a21,f15.3)')"  Simtime_in_hours = ",Simtime_in_hours
        endif;enddo
      endif
      if((CheckConditions(3).eqv..true.).and.&
         (StopConditions(3).eqv..true.))then
        ! Normal stop condition when nothing is left to advect
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"No ash species remain aloft."
        endif;enddo
      endif
      if((CheckConditions(4).eqv..true.).and.&
         (StopConditions(4).eqv..true.))then
        ! Error stop condition if the concen and outflow do not match the source
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Cummulative source volume does not match aloft + outflow"
          write(errlog(io),'(a11,f15.5)')" tot_vol = ",tot_vol
          write(errlog(io),'(a11,f15.5)')" SourceCumulativeVol = ",SourceCumulativeVol
          write(errlog(io),'(a14,f15.5)')" Abs. Error = ",&
                               abs((tot_vol-SourceCumulativeVol)/SourceCumulativeVol)
          write(errlog(io),'(a12,f15.5)')" e_Volume = ",e_Volume
        endif;enddo
        stop 1
      endif
      if((CheckConditions(5).eqv..true.).and.&
         (StopConditions(5).eqv..true.))then
        ! Error stop condition if any volume measure is negative
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"One of the volume measures is negative."
          write(errlog(io),'(a30,f13.5)')"                    dep_vol = ",dep_vol
          write(errlog(io),'(a30,f13.5)')"                  aloft_vol = ",aloft_vol
          write(errlog(io),'(a30,f13.5)')"                outflow_vol = ",outflow_vol
          write(errlog(io),'(a30,f13.5)')"        SourceCumulativeVol = ",SourceCumulativeVol
        endif;enddo
        stop 1
      endif

      ! ************************************************************************
      ! ****** end time simulation *********************************************
      ! ************************************************************************

      isFinal_TS = .true.
      Called_Gen_Output_Vars  = .false.
      Calculated_Cloud_Load   = .false.
      Calculated_AshThickness = .false.

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5012)   ! put footnotes below output table
        write(outlog(io),'(a5,f10.3,a5,f10.3)')'time=',time,', dt=',dt
        write(outlog(io),'(a26,g15.5)')"Mass Conservation Error = ",MassConsErr
      endif;enddo

        ! Make sure we have the latest output variables and go to write routines
      call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines here
!
!------------------------------------------------------------------------------

      call output_results
      ! Write results to log and standard output
      call cpu_time(t2) ! time is a scalar real
      call system_clock(tcount2,tcount_rate,tcount_max)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5003) t1-t0,tw_tot,t2-t1,&
                               real(tcount2-tcount1,kind=dp)/real(tcount_rate,kind=dp)
      endif;enddo
      call TimeStepTotals(itime)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5005) dep_vol
        write(outlog(io),5006) tot_vol
        write(outlog(io),5009) maxval(DepositThickness), DepositAreaCovered
        write(outlog(io),5034)       ! write out area of cloud at different thresholds
        do i=1,5
          write(outlog(io),5035) LoadVal(i), CloudLoadArea(i)
        enddo
        write(outlog(io),5033)       ! write "normal completion"
      endif;enddo

      ! Format statements
      ! Starting with 5000
5001  format(/,5x,'Starting volume (km3 DRE)    = ',f11.4,       &
             /,5x,'maximum number of time steps = ',i8,          &
             //,21x,'Time',19x,                                  &
                '|--------------------Volume (km3 DRE)-------------------|', &
                3x,'Cloud Area',                                 &
              /,7x,'step',8x,'(hrs)',2x,'yyyymmddhh:mm',         &
                5x,'Source',6x,'Deposit',7x,'Aloft',5x,'Outflow',&
                7x,'Total',10x,'km2')

5003  format(/,5x,'Set-up time (cpu)         = ',f15.4,' seconds',/&
               5x,'MetReader time (cpu)      = ',f15.4,' seconds',/&
               5x,'Total solver time (cpu)   = ',f15.4,' seconds',/&
              ,5x,'Wall clock time           = ',f15.4,' seconds') 
!5003  format(/,5x,'Set-up time              = ',f15.4,' seconds',/&
!               5x,'Execution time           = ',f15.4,' seconds',/&
!              ,5x,'Simulation time          = ',f15.4,' seconds')      
!5004  format(  5x,'Execution time/CPU time  = ',f15.4)
5005  format(  5x,'Ending deposit volume    = ',f15.4,' km3 DRE')       
5006  format(  5x,'Ending total volume      = ',f15.4,' km3 DRE')       
5007  format(  5x,'Building time array of plume height & eruption rate')
5009  format(/,5x,'Maximum deposit thickness (mm)   = ',f10.4, &
             /,5x,'Area covered by >0.01 mm (km2)   = ',f10.1,/)
5012  format(4x,'*=files written out')

5020  format('Calculating fall time from plume top',/,&
              5x,'GS index',5x,'diam (mm)',5x,'fall time (hours)')
5021  format(5x,i4,7x,f8.3,10x,f15.1)

5033  format(/,5x,'Normal completion')
5034  format(/,'  Ash load   cloud area',/, &
               '      T/km2         km2')
5035  format(2f11.1)

      ! clean up memory
      call dealloc_arrays

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls deallocation routines here
!
!------------------------------------------------------------------------------

      close(fid_logfile)       !close log file 

      end program Ash3d
