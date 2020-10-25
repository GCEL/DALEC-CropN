! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! Code used for inferring crop leaf nitrogen per area (LNA) based on the calibration of a semi-empirical critical N dilution function
! derived from the ATEC project experimental field trials. Code adapted from DALEC_GSI_DFOL_FR_CROP.f90 - key changes are the use of 
! N dilution slope and intercept parameters applied in lines 376 to 389. From a calibration of parameters 1-34 from a prior calibraiton 
! step the CARDAMOM MDF framework optimises the slope and intercept of the N dilusion (parameters 35 and 36) in order to find an N dilution
! fit consistent with an LAI timeseries.
! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! POOLS:  1.labile                                                              PARAMETERS:  1.Decomposition rate
!         2.foliar 																			 2.Fraction of GPP to respired
!         3.root                                                                             3.Maximum developmental rate: DS (0->1) 
!         4.stem                                                                             4.Maximum developmental rate: DS (1->2)
!         5.litter                                                                           5.Turnover rate of foliage pool
!		  6.SOM                  													         6.Turnover rate of stem pool
!		  7.autotrophic respiration 													     7.Maximum rate of foliar turnover due to self-shading
!         8.storage                                                                          8.Effective vernalisation days when plants are 50% vernalised	
! ---------------------------------------------------------------------------                9.Mineralisation rate of soil organic matter (SOM)
! FLUXES: 1.GPP (gC.m-2.day-1)																10.Mineralisation rate of litter
!         2.temprate (i.e. temperature modified rate of metabolic activity)                 11.Photosynthetic nitrogen use efficiency
!         3.autotrophic respiration (gC.m-2.day-1)                                          12.Sowing day
!         4.leaf production rate (gC.m-2.day-1)                                             13.Respiratory cost of labile transfer
!         5.labile production (gC.m-2.day-1)                                                14.phenological heat units required for emergence
!         6.root production (gC.m-2.day-1)                                                  15.Harvest day
!         7.stem production (gC.m-2.day-1)                                                  16.Plough day
!         8.labile consumption (gC.m-2.day-1)                                               17.Leaf C per area (LCA)
!         9.allocation to storage organ (gC.m-2.day-1)										18.Initial labile pool
!        10.total leaf litter production (gC.m-2.day-1)									    19.Initial foliar pool
!        11.stem litter production (gC.m-2.day-1)  									        20.Initial root pool
!        12.root litter production (gC.m-2.day-1)											21.Initial stem pool
!        13.respiration heterotrophic litter (gC.m-2.day-1)								    22.Initial litter pool
!        14.respiration heterotrophic SOM (gC.m-2.day-1)									23.Initial SOM pool
!        15.litter to SOM rate (gC.m-2.day-1)											    24.autotrophic pool
!        16.allocation to autotrophic pool (gC.m-2.day-1)									25.Initial storage pool  																					 
!                                                                                           26.Minimum temperature for development
! 																							27.Maximum temperature for development
!                                                                                           28.Optimum temperature for development
!																							29.Minimum temperature for vernalisation
!																						    30.Maximum temperature for vernalisation															
!																						    31.Optimum temperature for vernalisation
!																							32.Critical photoperiod for development
!																							33.Photoperiod sensitivity
!																							34.turnover rate of labile 
!																							35.Slope of nitrogen dilution 
!																							36.Intercept of nitrogen dilution
! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

module CARBON_MODEL_CROP_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL_CROP    &
         ,acm, management_dates,harvest &
         ,plough, interpolate, vernalization &
         ,development_stage, photoperiod_impact, temperature_impact &
         ,calc_pools_crops, carbon_alloc_fractions & 
         ,resp_rate_temp_coeff & 
         ,ts_length            &
         ,sec_in_day           &
         ,sec_in_hour

! declare module level variables
double precision, parameter :: sec_in_day  = 86400.0
double precision, parameter :: sec_in_hour = 3600.0
double precision, parameter :: pi          = 3.14159265
double precision, parameter :: deg_to_rad  = pi/180d0

! variables local to this module
integer ::   plough_day, & ! day-of-year when field is ploughed  
                sow_day, & ! day-of-year when field is sown      
            harvest_day, & ! day-of-year when field is harvested 
                stmob=0, & ! remoblise stem C to labile (1 = on)
 turnover_labile_switch    ! begin turnover of labile C

logical :: vernal_calcs, &
               ploughed, & ! was the soil ploughed ?
        use_seed_labile, & ! whether to use seed labile for growth 
                   sown, & ! has farmer sown the crop ?
                emerged    ! has the crop emerged  ?

double precision ::           ts_length, & ! time step length in hours
                            step_of_day, & ! current step of the day (default = 1)
                           steps_in_day, & ! number of steps in a day (default = 1)
                                    doy, & ! decimal doy of year 
                                gpp_acm, & ! gross primary productivity (gC.m-2.day-1)
                    stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m--2)
                     stock_dead_foliage, & ! dead but still standing foliage (gC.m--2)
                        stock_resp_auto, & ! autotrophic respiration pool (gC.m--2)
                           stock_labile, & ! labile C pool (gC.m--2)
                          stock_foliage, & ! foliage C pool (gC.m--2)
                             stock_stem, & ! stem C pool (gC.m--2)
                            stock_roots, & ! roots C pool (gC.m--2)
                           stock_litter, & ! litter C pool (gC.m--2)
                    stock_soilOrgMatter, & ! SOM C pool (gC.m--2)
                              resp_auto, & ! autotrophic respiration (gC.m-2.t-1)
                          resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.t-1)
                   resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2)
                                    npp, & ! net primary productivity (gC.m-2.t-1)
                              nee_dalec, & ! net ecosystem exchange (gC.m-2.t-1)
                              lai_local, & !
                                     DS, & ! Developmental state and initial condition
                                    LCA, & ! leaf carbon area (gC._leafm-2)
            mean_alloc_to_storage_organ, & ! rolling average allocation of GPP to storage organ (gC.m-2)
        mean_alloc_to_storage_organ_old, & ! ...same but previous value...
                     decomposition_rate, & ! decomposition rate (frac / hr)
                     frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                  turnover_rate_foliage, & ! turnover rate of foliage (frac/hr)
                     turnover_rate_stem, & ! same for stem
                   turnover_rate_labile, & ! same for labile 
                turnover_rate_resp_auto, & ! same for autotrophic C pool
                 resp_cost_labile_trans, & ! labile lost to respiration per gC labile to GPP
             mineralisation_rate_litter, & ! mineralisation rate of litter
      mineralisation_rate_soilOrgMatter, & ! mineralisation rate of SOM
                                  PHUem, & ! emergance value for phenological heat units
                                    PHU, & ! phenological heat units
                                 DR_pre, & ! development rate coefficient DS 0->1
                                DR_post, & ! development rate coefficient DS 1->2
                                   tmin, & ! min temperature for development
                                   tmax, & ! max temperature for development
                                   topt, & ! optimum temperature for development
                                 tmin_v, & ! min temperature for vernalisation
                                 tmax_v, & ! max temperature for vernalisation
                                 topt_v, & ! optimim temperature for vernalisation
                                    VDh, & ! effective vernalisation days when plants are 50 % vernalised 
                                     VD, & ! count of vernalisation days
                               RDRSHMAX, & ! maximum rate of self shading turnover
                                   PHCR, & ! critical value of photoperiod for development
                                   PHSC, & ! photoperiod sensitivity
                             raso = 0.0, & ! rolling average for alloc to storage organ
                         max_raso = 0.0, & ! maximum value for rolling average alloc to storage organ
                                  BM_EX, & ! 
                                     HI, & ! Harvest Index
                                  yield, & ! crop yield (gC.m-2)
                     alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                    turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                gso_max, & !
                     max_raso_old = 0.0, & !
                        raso_old  = 0.0, & !
            resp_cost_labile_to_foliage, & ! respiratory cost of moving carbon..from labile to foliage pools
            resp_cost_foliage_to_labile, & ! ..from foliage to labile pools
                              resp_rate, & ! rate of respiration at given temperature
                                 Cshoot, & !
                                     DR, & !
                        fol_frac_intpol, & ! Fraction of NPP C that goes to foliage pool
                       stem_frac_intpol, & ! Fraction of NPP C that goes to stem pool
                       root_frac_intpol, & ! Fraction of NPP C that goes to root pool
                      grain_frac_intpol, & ! Fraction of NPP C that goes to grain/storage_organ pool
                               fP,fT,fV, & !                      
                                  remob, & !
                                 avtemp, & !
                 alloc_to_storage_organ, & ! Ammount of C that is allocated to grain/storage_organ pool
                     litterfall_foliage, & !
                        litterfall_stem, & !
                       litterfall_roots, & !
                          decomposition, & !
                              npp_shoot, & ! Fraction of NPP C that goes to shoot (stem + leaves + grain)
                      alloc_from_labile, & ! Amount of C that is allocated from labile pool
                        alloc_to_labile, & ! Amount of C that is allocated to labile pool
                         alloc_to_roots, & ! Ammount of C that is allocated to roots pool
                       alloc_to_foliage, & ! Ammount of C that is allocated to leaves pool
                          alloc_to_stem, & ! Ammount of C that is allocated to stem pool
                                raremob, & !
                                  RDRSH, & !
                                  RDRDV, & !
                                    RDR, & !
                              daylength 
                                  

  ! defines Q10 = 2 in exponential temperature response for heterotrophic respiration
  double precision, parameter :: resp_rate_temp_coeff = 0.0693
  ! residue fraction of leaves left post harvest
  double precision, parameter :: lv_res = 0.1
  ! residue fraction of stem left post harvest
  double precision, parameter :: st_res = 0.1
  ! LAI above which self shading turnover occurs
  double precision, parameter :: LAICR = 4.0
  ! allocation to storage organ relative to GPP 
  double precision, parameter :: rel_gso_max = 0.35
  ! initial C available AR: hard code seed labile
  double precision, parameter :: stock_seed_labile = 9.0

save 

contains

!
!-------------------------------------------------------------------------------------------------------------------
!

  subroutine CARBON_MODEL_CROP(start,finish,met,pars,deltat,nodays,lat, &
                               LAI,NEE,FLUXES,POOLS,GPP, LNA,           &
                               nopars,nomet,nopools,nofluxes)

    ! The Data Assimilation Linked Ecosystem Carbon (DALEC) model. 
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP 
    ! and partitions between various ecosystem carbon pools. 
    ! These pools are subject to turnovers/decompostion resulting in 
    ! ecosystem phenology and fluxes of CO2
    ! use DALEC_CROP_DEV_VARIABLES ! Now need to define the variables here

    implicit none

    integer, intent(in) :: start    & 
                          ,finish   &
                          ,nopars   & ! number of paremeters in vector
						  ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation
                          
                         

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)              & ! time step in decimal days
                         ,pars(nopars)                & ! number of parameters
                         ,lat                           ! site latitude (degrees)

    double precision,  intent(out) :: lai(nodays) &   ! LAI
                                     ,GPP(nodays) &   ! GPP
                                     ,NEE(nodays) &   ! NEE
                                     ,LNA(nodays)     ! Leaf nitrogen per area (LNA)

    double precision, intent(out)  :: POOLS((nodays+1),nopools)! vector of ecosystem pools
 
    double precision, intent(out)  :: FLUXES(nodays,nofluxes)! vector of ecosystem fluxes
    
    ! Local only
    double precision :: gpppars(11)       & ! ACM inputs (LAI+met)
                       ,constants(10)     & ! parameters for ACM
                       ,slope_n           & ! slope of N dilution function
                       ,intercept_n         ! intercept of N dilution function   

    integer :: n
    
    
    double precision, allocatable, dimension(:) ::  DS_shoot, & !
                                                     DS_root, & !
                                                    fol_frac, & !
                                                   stem_frac, & !
                                                   root_frac
                                                             
    save

    gpppars(7)  = lat
    gpppars(9)  = -2.0 !-2.060814 ! leafWP-soilWP
    gpppars(10) = 1.0 ! 0.2 !1.0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters -- cropland
    constants(1)  = pars(11)
    constants(2)  = 0.04
    constants(3)  = 2.70e-4  
    constants(4)  = 83.18
    constants(5)  = 0.03     
    constants(6)  = 4.54
    constants(7)  = 3.86     
    constants(8)  = 4.10e-3
    constants(9)  = 0.38     
    constants(10) = 2.72e-8

    ! create DS arrays ---------------------------------------------------

    ! Crop C allocation per DS -- arrays   
    if (.not. allocated(DS_shoot)) allocate( DS_shoot(10) , fol_frac(10) , stem_frac(10)  ) ! Shoot 
    if (.not. allocated(DS_root)) allocate( DS_root(5) , root_frac(5) ) ! Root 
    
    DS_shoot  = (/ 0.0,0.33,0.43,0.53,0.62,0.77,0.95,1.14,1.38,2.1 /) ! DS 
    fol_frac  = (/ 0.9,0.85,0.83,0.75,0.56,0.20,0.10,0.05,0.0,0.0 /) ! Leaves 
    stem_frac = (/ 0.1,0.15,0.17,0.25,0.44,0.80,0.60,0.55,0.0,0.0 /) ! Stem

    DS_root   = (/ 0.0,0.33,0.53,1.0,2.1 /) ! DS 
    root_frac = (/ 0.5,0.5,0.25,0.0,0.0 /) ! Roots    
    ! ------------------------------------------------------------------  

    ! length of time step in hours
    ts_length = ((sum(deltat)/nodays) * sec_in_day) / sec_in_hour
    
    ! temporal resolution of ins/outs
    steps_in_day = 1 / (sum(deltat)/nodays) ! number of steps in a day 
    step_of_day  = steps_in_day             ! current step of the day 

    ! DALEC parameters 
    ! parameters from file
    decomposition_rate                = pars(1)  ! decomposition rate (frac / hr)
    frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    turnover_rate_foliage             = pars(5)  !pars(5)  ! turnover_rate of foliage (frac/hr)
    turnover_rate_stem                = pars(6)  ! turnover rate of stem (frac/hr)
    RDRSHMAX                          = pars(7)  ! maximum rate of foliar turnover due to self shading
    VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised 
    mineralisation_rate_soilOrgMatter = pars(9)  ! mineralisation rate som
    mineralisation_rate_litter        = pars(10) ! mineralisation rate litter ! PARS(11) = NUE, defined earler for constant(1)
    sow_day                           = nint(pars(12)) ! sow day (doy)
    resp_cost_labile_trans            = pars(13) ! labile lost to respiration per gC labile to GPP
    PHUem                             = pars(14) ! phenological heat units required for emergence
    harvest_day                       = nint(pars(15)) ! nint(mod(pars(15),365.25)) ! harvest day (doy)
    plough_day                        = nint(pars(16)) ! nint(mod(pars(16),365.25)) ! plough day (doy)
    LCA                               = pars(17) ! leaf mass area (gC.m-2)
    tmin                              = pars(26) ! min temperature for development
    tmax                              = pars(27) ! max temperature for development
    topt                              = pars(28) ! optimum temperature for development
    tmin_v                            = pars(29) ! min temperature for vernalisation
    tmax_v                            = pars(30) ! max temperature for vernalisation
    topt_v                            = pars(31) ! optimim temperature for vernalisation
    PHCR                              = pars(32) ! critical value of photoperiod for development
    PHSC                              = pars(33) ! photoperiod sensitivity
    turnover_rate_labile              = pars(34) ! turnover rate labile C
    turnover_rate_resp_auto           = pars(35) ! turnover rate of autotrophic carbon for respiration
    slope_n	                      	  = pars(36) ! Estimated slope for N dilution function
    intercept_n 	                  = pars(37) ! LNA estimated intecept for N dilution function


    if (start == 1) then

        yield = 0.0
        DS = -1.0 
        mean_alloc_to_storage_organ = 0.0
        mean_alloc_to_storage_organ_old = 0.0
        PHU = 0.0
        VD = 0.0
        BM_EX = 0.0
        HI = 0.0
        ploughed = .false.
        vernal_calcs = .true. 
        sown = .false.
        use_seed_labile = .true.  
        emerged = .false.
        stock_dead_foliage=0.0
        DR = 0.0
        alloc_to_labile = 0.0
        stmob = 0.0
        max_raso = 0.0
        DR = 0.0
        raso = 0.0
        RDRDV = 0.0 
        
        ! assigning initial conditions
        POOLS(1,1) = nint(pars(18)) ! labile C
        POOLS(1,2) = nint(pars(19)) ! foliar C
        POOLS(1,3) = nint(pars(20)) ! root C
        POOLS(1,4) = nint(pars(21)) ! stem C
        POOLS(1,5) = nint(pars(22)) ! litter C
        POOLS(1,6) = nint(pars(23)) ! som C
        POOLS(1,7) = pars(24)       ! autotrophic resp pool
        POOLS(1,8) = pars(25)       ! storage organ

        stock_labile           = POOLS(1,1) 
        stock_foliage          = POOLS(1,2) 
        stock_roots            = POOLS(1,3) 
        stock_stem             = POOLS(1,4) 
        stock_litter           = POOLS(1,5) 
        stock_soilOrgMatter    = POOLS(1,6) 
        stock_resp_auto        = POOLS(1,7) 
        stock_storage_organ    = POOLS(1,8) 


    endif 

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish 
       
      ! Calculate LAI 
      lai(n)    = POOLS(n,2) / LCA
      lai_local = lai(n)

      ! ACM's parameters 
      gpppars(1)   = lai_local                         ! LAI 
      gpppars(2)   = met(3,n)                          ! max T
      gpppars(3)   = met(2,n)                          ! min T

	  ! ------------------------------------------------------
	  ! CALCULATION OF LNA (gpppars(4)) 
	  ! The calculation of LNA (used in ACM) is based on the destructively sampled winter wheat C and N data, which was acquired during key growth stages in accordance with
	  ! the UK AHDB winter wheat growth guide (https: //projectblue.blob.core.windows.net/media/Default/Imported%20Publication%20Docs/Wheat%20growth%20guide.pdf): from early tillering (GS 20) to 
	  ! post-anthesis/reprodictive (GS 75).
	  
	  if (DS .lt. 0.15) then                           ! DS < 0.15 corresponds to the growth stage at beginning of the UK reccommended period of      
         gpppars(4) = intercept_n                      ! N fertiliser application for winter wheat (Zodocks growth stage 20) - the early tillering stage (typically mid-march to April)
      else if (DS .ge. 0.15) then                      
         gpppars(4) = (slope_n*POOLS(n,2))+intercept_n ! NOTE: The slope_n parameter can be included in the MDF optimisation. The value for this parameter has also been observed to be 
      end if                                           !       similar across different N application rates and so a mean value of -0.02448 could also be used.

      if (DS .gt. 1.15) then                           ! Set LNA to 0.5 after anthesis (Zodocks growth stage 75)  
         gpppars(4) = 0.5
      end if
      ! ------------------------------------------------------
	  
      gpppars(5)   = met(5,n)                          ! co2
      gpppars(6)   = ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)   = met(4,n)                          ! radiation     
      
      ! Day of Year
      doy = met(6,n)

      ! Calculate GPP (gC.m-2.day-1)
      if (lai(n) > 1e-5) then
          gpp_acm = acm(gpppars,constants)
      else
          gpp_acm = 0.0
      end if

      ! Calculate mean Temperature
      avtemp = 0.5 * ( met(3,n) + met(2,n) )

      ! Calculate weighted air temperature value based on daily minimum, maximum and means
      resp_rate = 0.5 * exp( resp_rate_temp_coeff * avtemp )

      ! daily avg of allocation to storage organ (needed to determine max storage organ growth rate)
      mean_alloc_to_storage_organ_old = mean_alloc_to_storage_organ
      mean_alloc_to_storage_organ     = 0.0

      ! Turnover rate of roots        
      if ( DS .lt. 1.1 ) then 
        ! turnover_rate_foliage = 0.0
        turnover_rate_roots   = 0.0
      else
        ! turnover_rate_foliage = pars(5)  
        turnover_rate_roots   = 0.00462
      endif

      ! determine development stage (DS) 
      call development_stage(deltat(n))
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)
      ! begin carbon allocation for crops
      call calc_pools_crops()
      ! conduct management updates (at the end of timestep)
      call management_dates(stock_seed_labile,deltat(n))

      NEE(n) = nee_dalec
      GPP(n) = gpp_acm
      LNA(n) = gpppars(4)

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = GPP(n)
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = resp_rate
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = resp_auto * steps_in_day
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = alloc_to_foliage * steps_in_day
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (alloc_to_labile + remob) * steps_in_day
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = alloc_to_roots * steps_in_day
      ! wood production 
      FLUXES(n,7) = alloc_to_stem * steps_in_day
      ! alloc to storage organ
      FLUXES(n,9) = alloc_to_storage_organ * steps_in_day
      ! alloc to autotrophic pool
      FLUXES(n,16) = frac_GPP_resp_auto * gpp_acm
      ! labile production
      FLUXES(n,8) = (alloc_from_labile + resp_cost_labile_to_foliage) * steps_in_day
      ! total leaf litter production
      FLUXES(n,10) = litterfall_foliage * steps_in_day
      ! total wood production
      FLUXES(n,11) = litterfall_stem * steps_in_day
      ! total root litter production
      FLUXES(n,12) = litterfall_roots * steps_in_day
      ! respiration heterotrophic litter
      FLUXES(n,13) = resp_h_litter * steps_in_day
      ! respiration heterotrophic som
      FLUXES(n,14) = resp_h_soilOrgMatter * steps_in_day
      ! litter to som 
      FLUXES(n,15) = decomposition * steps_in_day
      ! alloc to autotrophic pool
      FLUXES(n,16) = frac_GPP_resp_auto * GPP(n)

      ! labile pool
      POOLS(n+1,1) = stock_labile
      ! foliar pool
      POOLS(n+1,2) = stock_foliage
      ! wood pool
      POOLS(n+1,4) = stock_stem
      ! root pool
      POOLS(n+1,3) = stock_roots
      ! litter pool
      POOLS(n+1,5) = stock_litter
      ! som pool
      POOLS(n+1,6) = stock_soilOrgMatter
      ! autotrophic pool 
      POOLS(n+1,7) = stock_resp_auto
      ! storage organ pool
      POOLS(n+1,8) = stock_storage_organ

      ! harvest index 
      !HarvestIndex = HI
      
    end do 

  end subroutine CARBON_MODEL_CROP

!
!-------------------------------------------------------------------------------------------------------------------
!
  
  double precision function acm(drivers,constants)

    ! The Aggregated Canopy Model is a GPP emulator which operates at a daily time step.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(11) & ! drivers
                                   ,constants(10) ! parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, ci, e0, cps, dec, nit &
             ,trange, sinld, cosld,aob,pi &
             ,mint,maxt,radiation,co2,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! initial values
    gc=0.0 ; pp=0.0 ; qq=0.0 ; ci=0.0 ; e0=0.0 ; cps=0.0 ; dec=0.0 ; nit=1.

    ! load driver values 
    maxt = drivers(2)
    mint = drivers(3)
    nit = drivers(4)    ! Linked to N dilution function
    co2 = drivers(5)
    radiation = drivers(8)
    lat = drivers(7)

    ! load parameters
    pi = drivers(11)
    deltaWP = drivers(9)
    Rtot = drivers(10)
    NUE = constants(1)  
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

    ! determine temperature range 
    trange=0.5*(maxt-mint)
    ! daily canopy conductance 
    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange)) ! default
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn=lai_local*nit*NUE*exp(temp_exponent*maxt) ! default
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci=0.5*(co2+qq-pp+((co2+qq-pp)**2-4.0*(co2*qq-pp*co2_comp_point))**0.5)
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0=lai_coef*lai_local**2/(lai_local**2+lai_const)
    ! calculate day length (in hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    sinld = sin( lat*deg_to_rad ) * sin( dec )
    cosld = cos( lat*deg_to_rad ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    daylength = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )
    ! calculate CO2 limited rate of photosynthesis
    pd=gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps=e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm=cps*(dayl_coef*daylength+dayl_const)

    return

  end function acm
  
!
!-------------------------------------------------------------------------------------------------------------------
!
  
  subroutine calc_pools_crops()

    ! Allocation of GPP to NPP and various carbon pools.   !
    ! Based on physiological responses to temperature      !
    ! vernalisation, and photoperiod.                      !

    implicit none

    ! if sown turn on labile / seed turnover for growth
    if ( sown ) then
      ! turnover on
      turnover_labile_switch = 1
    else
      ! turnover off
      turnover_labile_switch = 0
    endif

    ! Initialise..
    resp_cost_foliage_to_labile = 0.0 

    ! respiratory cost of C transfer from labile pool to short-term pool (NPP)
    resp_cost_labile_to_foliage = turnover_rate_labile * stock_labile * resp_cost_labile_trans * resp_rate &
                                   * ts_length * real(turnover_labile_switch)
    ! resp_cost_labile_to_foliage = stock_labile * min(1.0,resp_cost_labile_to_foliage)

    ! allocation flux from labile C pool to NPP
    alloc_from_labile = turnover_rate_labile * stock_labile * ( 1.0 - resp_cost_labile_trans ) * resp_rate &
                                  * ts_length * real(turnover_labile_switch)
    ! alloc_from_labile = stock_labile * min(1.0,alloc_from_labile)

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    if ( ( gpp_acm .gt. alloc_from_labile ) .and. ( use_seed_labile ) ) then
      stock_litter = stock_litter + stock_labile
      stock_labile = 0.0
      use_seed_labile = .false.
    endif

    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation..
    npp = ( 1.0 - frac_GPP_resp_auto ) * gpp_acm + alloc_from_labile

    root_frac_intpol  = max(0d0,min(1d0,root_frac_intpol))
    alloc_to_roots    = root_frac_intpol * npp         !
    npp_shoot         = npp - alloc_to_roots           ! NPP remaining after root growth = SHOOT frac
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot   !
    alloc_to_stem     = stem_frac_intpol * npp_shoot   !
    alloc_to_storage_organ = max(0d0,npp_shoot - alloc_to_foliage - alloc_to_stem)
    
    if ( alloc_to_storage_organ > 0.0 ) then  ! allocation flux to storage organ limited by maximum growth rate
        gso_max  = ( stock_storage_organ + 0.5 ) * rel_gso_max / steps_in_day
        alloc_to_storage_organ = min( alloc_to_storage_organ , gso_max )
        
        if ( sown ) then
          alloc_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                              * ( 1.0 - resp_cost_labile_trans )
          resp_cost_foliage_to_labile =  ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                              * resp_cost_labile_trans
        else
          alloc_to_labile             = 0.0
          resp_cost_foliage_to_labile = 0.0
        endif
    
    endif
    
    mean_alloc_to_storage_organ = mean_alloc_to_storage_organ + alloc_to_storage_organ

    ! set switches to (de)activate leaf, root and stem remobliization
    if ( step_of_day .eq. steps_in_day ) then
      mean_alloc_to_storage_organ = mean_alloc_to_storage_organ / steps_in_day
      raso_old = raso
      ! running average of growth rate of storage organ..
      raso = ( mean_alloc_to_storage_organ + mean_alloc_to_storage_organ_old ) * 0.5
      max_raso_old = max_raso
      max_raso = max( raso , max_raso_old )
      ! Stem remobilisation triggered once running average of storage organ growth declines
      ! Second part prevents premature remobilisation
      if ( ( raso .lt. raso_old ) .and. &
            ( mean_alloc_to_storage_organ .gt. ( mean_alloc_to_storage_organ_old + 0.5 ) / steps_in_day ) ) then
          stmob = 1.
      else 
          stmob = 0.
      endif
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a 
    ! function of shading (RDRSH) or developmental stage (RDRT)

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    RDRSH = min( RDRSHMAX , max( 0d0 , RDRSHMAX * ( lai_local - LAICR ) / LAICR ) )
    if ( DS < 1.0 ) then
       RDRDV = 0.0
    else
       RDRDV = turnover_rate_foliage * ( 1.0 / ( ( max( 2.0 - DS , 0.1 ) ) * 8.0 ) ) ** 2
    ENDIF

    ! relative leaf death rate is the maximum value of the arguments RDRSH and RDRDV
    RDR = max( RDRSH , RDRDV )

    litterfall_foliage = stock_foliage * ts_length * RDR
    litterfall_stem    = stock_stem    * ts_length * DR * turnover_rate_stem * real(stmob) 
    litterfall_roots   = stock_roots   * ts_length * turnover_rate_roots

    ! remobilized C to NPP (from both leaves and stems)
    remob   = ( litterfall_foliage * 0.5 + litterfall_stem ) * ( 1.0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates)
    Raremob = ( litterfall_foliage * 0.5 + litterfall_stem ) * resp_cost_labile_trans

    ! heterotrophic respiration component 1: mineralisation of litter C pool
    resp_h_litter = mineralisation_rate_litter * stock_litter * resp_rate * ts_length
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool
    resp_h_soilOrgMatter = mineralisation_rate_soilOrgMatter * stock_soilOrgMatter * resp_rate * ts_length
    ! decomposition of litter to soil organic matter
    decomposition = decomposition_rate * stock_litter * resp_rate * ts_length

    ! Recalculate Carbon Pools...
    stock_foliage       = max(0d0, stock_foliage + alloc_to_foliage - litterfall_foliage)
    stock_stem          = max(0d0, stock_stem + alloc_to_stem - litterfall_stem)
    stock_storage_organ = max(0d0, stock_storage_organ + alloc_to_storage_organ)
    stock_roots         = max(0d0, stock_roots + alloc_to_roots   - litterfall_roots)
    stock_litter        = max(0d0, stock_litter + litterfall_roots - resp_h_litter - decomposition)
    stock_soilOrgMatter = max(0d0, stock_soilOrgMatter + decomposition - resp_h_soilOrgMatter)
    stock_dead_foliage  = max(0d0, stock_dead_foliage  + litterfall_foliage * 0.5)
    stock_labile        = max(0d0, stock_labile + alloc_to_labile - alloc_from_labile - resp_cost_labile_to_foliage + remob)

    ! respiratory pool: new photosynthates are added 
    stock_resp_auto = stock_resp_auto + frac_GPP_resp_auto * gpp_acm
    ! autotrophic respiration; Ra (7% of respiratory pool)
    resp_auto = stock_resp_auto * turnover_rate_resp_auto * ts_length
    ! respiratory pool reduced by Ra (amount of C respired by plant)
    stock_resp_auto = max(0d0, stock_resp_auto - resp_auto)
    ! respiratory cost of C transfer from labile pool to short-term pool added
    ! to yield total autotrophic respiration
    resp_auto = resp_auto + resp_cost_labile_to_foliage + resp_cost_foliage_to_labile + Raremob
    ! nee (gC.m-2.t-1)
    nee_dalec = (resp_auto + resp_h_litter + resp_h_soilOrgMatter) - gpp_acm

  end subroutine calc_pools_crops
  
!
!-------------------------------------------------------------------------------------------------------------------
!

  subroutine carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

    double precision, dimension(:), intent(in) ::   DS_shoot, &
                                                     DS_root, & 
                                                    fol_frac, & 
                                                   stem_frac, & 
                                                   root_frac  
						

    double precision, dimension(:), allocatable :: frac_shoot, frac_root

    if ( sown ) then ! after sowing

       ! leaf development stages and corresponding fractions
       frac_shoot = fol_frac
       ! interpolate between PdV allocation values with reference to developmental stage 
       fol_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )
       ! stem DS and corresponding fractions
       frac_shoot = stem_frac
       stem_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )
       ! root DS and corresponding fractions
       frac_root = root_frac
       root_frac_intpol = interpolate( DS , DS_root , frac_root , size(DS_root) )
	 

    endif

  end subroutine carbon_alloc_fractions

!
!-------------------------------------------------------------------------------------------------------------------
!
  
  subroutine development_stage(days_in_step)

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

    ! agruments
    double precision :: days_in_step    
    
    ! local variables..
    double precision ::  doptmin, & ! Difference between optimum and minimum temperature
                         dmaxmin, & ! Difference between maximum and minimum temperature
                          dttmin, & ! Difference between daiy average and minimum temperatures
                       doptmin_v, & ! Difference between optimum and minimum vernalization temperatures
                       dmaxmin_v, & ! Difference between maximum and minimum vernalization temperatures
                       dttmin_v     ! Difference between daily average and minimum vernalization temperatures

    doptmin   = topt   - tmin   ! difference between optimal and minimum cardinal temperatures
    dmaxmin   = tmax   - tmin   ! difference between maximum and minimum cardinal temperatures
    dttmin    = avtemp - tmin   ! difference between daily average and minimum cardinal temperatures
    doptmin_v = topt_v - tmin_v ! same as above,
    dmaxmin_v = tmax_v - tmin_v !       but for vernalization 
    dttmin_v  = avtemp - tmin_v ! cardinal temperatures
        
    ! Calculation of developmental function values: vernalization (fV), temperature (fT) and
    ! photoperiod (fP) these values are multiplicative factors of DRmax (maximum developmental
    ! rate), each ranging between 0 (no development) and 1 (unrestricted development).

    ! Summation of vernalization days (VD), not before sowing and only if
    ! average temperature is within min and max cardinal temperatures..
    if ( ( avtemp .gt. tmin_v ) .and. ( avtemp .lt. tmax_v ) .and. sown .and. .not. emerged ) then
       fV = vernalization( doptmin_v , dmaxmin_v , dttmin_v, days_in_step)
    endif     
    
    ! Only calculate temperature coefficient if avtemp lies within (tmin,tmax)
    ! range.
    if ( (avtemp .gt. tmin ) .and. ( avtemp .lt. tmax ) ) then
      fT = temperature_impact( doptmin , dmaxmin , dttmin )
    else
      fT = 0.0
    endif

    
    fP = photoperiod_impact( PHCR , PHSC ) ! calculation of photoperiod coefficient    

    if ( emerged .and. ( DS .lt. 2.0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( DS .lt. 1.0 ) then  ! in the vegetative phase (before flowering):

          DR = DR_pre * fT * fP   ! DR is affected by temperature, photoperiod...

          if ( vernal_calcs ) DR = DR * fV ! ...and vernalization (for winter cereals)
          DS = DS + (DR * days_in_step)  ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          DR = DR_post * fT   ! DR is affected only by temperature
          DS = DS + (DR * days_in_step)

       endif

    endif
  
  end subroutine development_stage
  
!
!-------------------------------------------------------------------------------------------------------------------
!
  
  subroutine management_dates (stock_seed_labile,days_in_step)

    ! use CARBON_MODEL_CROP_MOD

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

    ! arguments
    double precision, intent(in) :: stock_seed_labile, days_in_step
    ! local variables
    double precision :: tmp
    logical :: plough_sanity, sow_sanity, harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.

    ! spring crops
    ! if (sow_day < harvest_day .and. nint(doy) < harvest_day) sow_sanity = .true.
    ! if (plough_day < harvest_day .and. nint(doy) < harvest_day) plough_sanity = .true.
    ! if (harvest_day > sow_day) harvest_sanity = .true.
    ! winter crops
    if (sow_day > harvest_day) sow_sanity = .true. 
    if (plough_day > harvest_day) plough_sanity = .true.
    if (harvest_day < plough_day .and. nint(doy) < plough_day) harvest_sanity = .true.   
    
    if ( .not. sown ) then ! fresh field...

      if ( plough_sanity .and. .not.ploughed .and. nint(doy) >= plough_day ) then
        ! the field needs ploughing..
        call plough 

      elseif ( sow_sanity .and. nint(doy) >= sow_day ) then
                 
        ! ensure that the field has indeed been ploughed
        if (.not.ploughed) then
         call plough
        end if 

        ! the field needs sowing..
        sown = .true.
                
        ! this switch controls whether the labile C within the seed is used for growth
        use_seed_labile = .true.
        stock_labile = stock_seed_labile

      endif ! plough or sow?

    else
      
      ! crop in field... calculate when crop emerges
      if ( .not. emerged ) then

         ! estimate emergence date based on the accumulated phenological heat units (PHU)
         ! where PHU is the (positive) heat over tmin..
         tmp = max( avtemp - tmin_v , 0.0 ) * days_in_step
         PHU = PHU + tmp         
         
         ! set the development stage and emergence..
         if ( PHU .ge. PHUem ) then
           emerged = .true.
           DS = 0.0
         else
           emerged = .false.
           DS = -1.0
         endif
                  
      endif ! emerged or not

      ! note that in this case harvest day has been fixed relative to the sow day
      if ( harvest_sanity .and. nint(doy) >= harvest_day) then ! the field needs harvesting
         call harvest
      endif

   endif ! sown or not

  end subroutine management_dates

!
!-------------------------------------------------------------------------------------------------------------------
!

  subroutine harvest()

    implicit none

    ! shoot biomass..
    Cshoot = stock_foliage + stock_stem + stock_storage_organ + stock_labile

    ! determine harvest index..
    HI = stock_storage_organ / Cshoot

    ! the stuff we actually want from the harvest...
    yield  = stock_storage_organ

    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = stock_foliage * ( 1.0 - lv_res )          &
              + stock_stem * ( 1.0 - st_res )          &
               + stock_dead_foliage * ( 1.0 - lv_res ) &
                + stock_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter                     &
                    + stock_resp_auto                &
                     + stock_foliage * lv_res        &
                      + stock_stem * st_res          &
                       + stock_dead_foliage * lv_res

    ! empty the plant stocks..
    stock_storage_organ = 0.0
    stock_foliage       = 0.0
    stock_stem          = 0.0
    stock_dead_foliage  = 0.0
    stock_labile        = 0.0
    stock_resp_auto     = 0.0

    ! reset logical variables..
    sown     = .false.
    emerged  = .false.
    ploughed = .false.
    DS = -9999 ; fV = 0.0 ; fT = 0.0 ; fP = 0.0 ; VD = 0.0

  end subroutine harvest

!
!-------------------------------------------------------------------------------------------------------------------
!

  double precision function photoperiod_impact( PH_crit , PH_sens )

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity
    
    photoperiod_impact = max(0.0, 1.0 - exp ( - PH_Sens * ( daylength - PH_crit ) ))

  end function photoperiod_impact

!
!-------------------------------------------------------------------------------------------------------------------
!
  
  subroutine plough

    implicit none

    stock_litter        = stock_litter + stock_dead_foliage &
                          + stock_foliage + stock_labile    &
                           + stock_roots + stock_stem       &
                            + stock_storage_organ

    stock_dead_foliage  = 0.0
    stock_foliage       = 0.0
    stock_labile        = 0.0
    stock_roots         = 0.0
    stock_stem          = 0.0
    stock_storage_organ = 0.0

    ploughed = .true. ; DS = -1.0 ; PHU = 0.0
    max_raso = 0.0 ; raso = 0.0 ; max_raso_old = 0.0 ; raso_old = 0.0
    mean_alloc_to_storage_organ_old = 0.0 ; mean_alloc_to_storage_organ = 0.0

  end subroutine plough

!
!-------------------------------------------------------------------------------------------------------------------
!
  
  double precision function temperature_impact( doptmin , dmaxmin , dttmin )

    ! Function to determine the coefficent for temperature impact on developmental rate. !
    ! From Streck et al., 2003.                                                          !

    implicit none

    ! arguments..
    double precision,intent(in) :: doptmin , dmaxmin , dttmin   ! temperature differences
    ! local variables..
    double precision :: a , nmr , dnr

    a   = log( 2.0 ) / ( log( ( dmaxmin ) / doptmin ) )
    nmr = 2.0 * ( ( dttmin ) ** a ) * ( doptmin ** a ) - ( ( dttmin ) ** ( 2.0 * a ) )
    dnr = doptmin ** ( 2.0 * a )
    temperature_impact = nmr / dnr

  return

  end function temperature_impact
  
!
!-------------------------------------------------------------------------------------------------------------------
!
  
  double precision function vernalization( doptmin_v , dmaxmin_v , dttmin_v, days_in_step )

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

    double precision,intent(in) :: dmaxmin_v , doptmin_v , dttmin_v, days_in_step ! temperature differences
    
    double precision :: a , dnr , fvn , nmr

    a   = log( 2.0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2.0 * ( ( dttmin_v ) ** a ) * ( doptmin_v ** a ) - ( ( dttmin_v ) ** (2.0 * a ) )
    dnr = doptmin_v ** ( 2.0 * a )
    fvn = nmr / dnr
    VD = VD + (fvn * days_in_step)
    ! final output value : 
    vernalization = max( 0d0 , min( 1d0 , ( VD ** 5.0 ) / ( ( VDh ** 5.0 ) + (VD ** 5.0 ) ) ) )

  return
  
  end function vernalization

  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  
  double precision function interpolate( x , reference_x , reference_y , row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    integer,intent(in)                         :: row
    double precision,intent(in)                :: x
    double precision,dimension(row),intent(in) :: reference_x , reference_y

    integer::i

    do i = 1 , row

       if ( x .le. reference_x(1) ) then
          interpolate = reference_y(1)
          exit
       endif

       ! cycling means growth rate remains constant between DS levels
       if ( ( x .gt. reference_x(i) ) .and. ( i .lt. row ) ) cycle

       if ( x .eq. reference_x(i) ) then
          interpolate = reference_y(i)
          exit
       endif

       if ( x .lt. reference_x(i) ) then
          interpolate = reference_y(i-1) + ( x - reference_x(i-1) ) &
                       * ( reference_y(i) - reference_y(i-1) )      &
                       / ( reference_x(i) - reference_x(i-1) )
          exit
       else
          interpolate = reference_y(row)
       endif

    enddo
    
  return

  end function interpolate
  
!
!-------------------------------------------------------------------------------------------------------------------
!

end module CARBON_MODEL_CROP_MOD
