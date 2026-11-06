!DEC$ FREEFORM
program main
    use DDMParameters
    use strainFatigue_
    use stressFatigue_ 
    use thermalFatigue_
    use goldenGIc_
    use bisectionGIcPrime_
    use Data_
    implicit none
    ! (c) 2026 Ever Barbero, ORCID 0000-0001-9121-6277
    !
    integer, allocatable :: LSS(:)          !degrees
    real(8), allocatable :: Orientations(:) !lamina stacking sequence in radians 
    real(8), allocatable :: Thicknesses(:) 
    real(8), allocatable :: trueStateV(:,:)
    real(8), allocatable :: Damages(:, :)
    real(8), allocatable :: StressInLaminaCS(:,:)
    real(8), allocatable :: PROPS(:), STATEV(:)
    !
    real(8), allocatable :: EnergyState(:,:)!activation function g, ERR GI, and GII, in this order, for all laminas
    ! Activation function, GI and GII for all laminas. If commented, not used at this time. 
    ! real(8) GIenergy(NL), GIIenergy(NL), gFunction(NL) 
    ! Next, to calculate laminate CTE. If commented, not used at this time. 
    ! real(8) alphalaminate(2,3), Genergy(2)
    ! 
    ! UGENS/UMAT variables. Additional variables declared in downstream procedures.
    integer :: NPROPS   ! = nLaminateProps + NL * nLaminaProps   !...Props from DDMParameters.mod
    integer :: NSTATV   ! = nStateVarPerPly*NL
    integer, parameter :: NTENS = 3  !passed to UMAT for explicit declaration of arrays
    integer, parameter :: NSECV = 6  !passed to UGENS for explicit declaration of arrays
    real(8) :: TEMP=0.D0, DTEMP=0.D0 !initialized in main.f
    real(8) :: STRAN(NSECV),DSTRAN(NSECV),TIME(2),DTIME
    integer :: BC(NSECV)                ! BOUNDARY CONDITIONS KEYS
    real(8) :: FC(NSECV), deltaFC(NSECV)! BOUNDARY CONDITIONS VALUES (forces and strains) 
    ! material properties & state variables
    real(8) :: PlyThickness
    real(8) :: TempDepProps(5,3)! formerly Laminae Properties(5,3)   *
    real(8) :: TempDepCTE(2,4)  ! formerly Laminae CTE(2,4)          *
    real(8) :: GIc(5), GIIc(5), TempRange(5), ThermFatigG(5)
    real(8) :: F1t, F1c, F2c, mWeibull, SFT
    ! local
    integer :: ios, NL !, analysisType moved to Data_
    real(8) :: targetS(NSECV), targetTemp    ! 3 variables read from input file, the remaining 3 not used
    real(8) :: Nx, Ny, Nxy, LT  ! stress resultants and laminate thickness
    integer, parameter :: nSteps = 200
    real(8), parameter :: pi=4d0*atan(1d0)
    ! Key = -2
    real(8) :: dummy, lambda_fit
    integer, parameter :: ndp=30  ! number of data points
    !
    ! Executable
    write(*,*) 'DDM6-gfortran-v1.0.0 (C) 2016 Ever Barbero'
    open(5,FILE='FILE5.INP',STATUS='OLD',iostat=ios)
    if (ios /= 0) stop '*** ERROR OPENING the INP FILE ***'
    read(5,*) NL
    !
    ! allocate requires NL, so it must be below read(5,*) NL
    NPROPS = nLaminateProps + NL * nLaminaProps ! nLaminateProps & nLaminaProps from DDMParameters
    NSTATV = nStateVarPerPly*NL                 ! nStateVarPerPly from DDMParameters
    allocate (LSS(NL))                          ! integer, degrees
    allocate (PROPS(NPROPS))
    allocate (STATEV(NSTATV))
    allocate (Orientations(NL))                 ! real, radians
    allocate (Thicknesses(NL)) 
    allocate (trueStateV(NL,nTrueStateVar))    
    allocate (Damages(NL, nDamages))
    allocate (StressInLaminaCS(NL, nStressComp))
    allocate (EnergyState(NL, nEnergies))
    !
    read(5,*) PlyThickness
    read(5,*) LSS   ! Orientations later
    read(5,*) Thicknesses
    read(5,*) GIc(1), GIc(2), GIc(3), GIc(4), GIc(5)  !GIc(A,B,C,TGImax,TGImin) temp-depend. properties
    read(5,*) GIIc(1), GIIc(2), GIIc(3), GIIc(4), GIIc(5) !GIIc(A,B,C,TGIImax,TGIImin) temp-depend. properties
    read(5,*) TempRange(1), TempRange(2), TempRange(3), TempRange(4), TempRange(5) !TempRange(SFT,Temax,Temin,Ttmax,Ttmin)
    read(5,*) ThermFatigG(1), ThermFatigG(2), ThermFatigG(3), ThermFatigG(4), ThermFatigG(5)    ! betaI, betaII, Ncycles, aI, aII
    SFT = TempRange(1)
    read(5,*) F1t
    read(5,*) F1c
    read(5,*) F2c
    read(5,*) mWeibull
    ! E1,E2,G12,v12,v23, same for all laminas, 2nd order polynomial
    ! * Do not use TempDepProps & TempDepCTE outside main and PackProperties, use LaminaeProperties(5) and laminaeCTE(3) instead *
    read(5,*) TempDepProps(1,1), TempDepProps(1,2) , TempDepProps(1,3) !E1=1 --> A,B,C temp-depend. properties
    read(5,*) TempDepProps(2,1), TempDepProps(2,2) , TempDepProps(2,3) !E2=2 --> A,B,C temp-depend. properties  
    read(5,*) TempDepProps(3,1), TempDepProps(3,2) , TempDepProps(3,3) !G12=3 --> A,B,C temp-depend. properties  
    read(5,*) TempDepProps(4,1), TempDepProps(4,2) , TempDepProps(4,3) !nu12=4 --> A,B,C temp-depend. properties 
    read(5,*) TempDepProps(5,1), TempDepProps(5,2) , TempDepProps(5,3) !nu23=5 --> A,B,C temp-depend. properties 
    ! CTE1, CTE2:  3rd order polynomial
    read(5,*) TempDepCTE(1,1), TempDepCTE(1,2), TempDepCTE(1,3) , TempDepCTE(1,4) ! tangentCTE Longitudinal=1 --> A,B,C,D temp-depend. properties
    read(5,*) TempDepCTE(2,1), TempDepCTE(2,2), TempDepCTE(2,3) , TempDepCTE(2,4) ! tangentCTE Transverse=2 --> A,B,C,D temp-depend. properties
    !
    ! read target epsx & final temperature
    read(5,*) targetS(1:3), targetTemp, analysisType  !analysisType=1: targetS=Strains, analysisType=2: targetS=Stresses,
    write(*,"(A,F5.2,A,L2)") " Seed crack density =",InitialLambda," useJacobian =",useJacobian
    ! write(*,*) "If analysisType = 0 > thermal-fatigue at Ncycles, = 1 strain-fatigue(*), = 2 stress-fatigue(*), (*)w/fixed deltaT=targetTemp-SFT"
    ! write(*,'(" targetS =",G12.4," targetTemp=",G12.4,"SFT =",G12.4,"Ncycles=",I10)') targetS(1), targetTemp, TempRange(1), int(ThermFatigG(3))
    if (targetTemp > SFT) stop "*** targetTemp must be less than or equal to SFT ***"
    if (all(targetS(1:3) == 0.0d0) .and. analysisType /= 0 ) stop "* non-zero stress or strain must be present unless KEY = 0"
    close(5)
    open(7,FILE='FILE7.out',FORM='FORMATTED',ACTION='WRITE',STATUS='REPLACE',iostat=ios)
    if (ios /= 0) stop '*** ERROR OPENING the OUTPUT FILE ***'
    !
    Orientations = LSS*pi/180d0     !radians
    Thicknesses(:) = Thicknesses(:)*PlyThickness
    LT = sum(Thicknesses(:))   ! laminate thickness (half LSS)
    trueStateV(:,1) = InitialLambda ! seed damages from DDMparamters.mod
    ! Must pack the initial damages into STATEV & properties into PROPS
    call PackStateVariables(trueStateV, Damages, nstatv, NL, STATEV, StressInLaminaCS, EnergyState)!output => STATEV
    call PackProperties(PROPS, NPROPS, NL, GIc, GIIc, F1t, F1c, F2c, mWeibull, &
        Orientations, Thicknesses, TempDepProps, TempDepCTE, TempRange, ThermFatigG)!output => PROPS
    call InitializeStateVariables(STATEV, NSTATV)   !needed, do not remove
    !
    if (analysisType == 0) then
        ! thermal-only static/fatigue analysys w/targets: 0.000 0.0 0.0 Tmin 0
        if ((ThermFatigG(1)<=0.d0 .or. ThermFatigG(2)<=0.d0) .and. ThermFatigG(3)>1.d0) stop "*** betaI and betaII must be > 1.0 for thermal fatigue ***"
        write(*,'(12A12)') '   TEMP        ', ' Strain     ', 'lambda i  '
        write(7,'(12A12)') '   TEMP        ', ' Strain     ', 'lambda i  '
        STRAN = 0.0D0
        DSTRAN = 0.D0 
        TEMP = SFT          !thermal analysis starts at SFT
        DTEMP = -1.0D0      !must be negative
        ! BC(): BOUNDARY CONDITIONS KEYS: 0=FORCE, 1=STRAIN 
        ! FC(): BOUNDARY CONDITIONS VALUES: FORCE:stress-resultants(#) &/or DISPLACEMENT:strains(#). (#) in this work.  
        ! BC(:) = (/0,0,0,0,0,0/)         ! only forces specified
        ! BC(:) = (/1,1,1,1,1,1/)         ! only displacements (strains) specified
        BC(:)   = (/0,0,0,1,1,1/)         ! force, force, force, kappa, kappa, kappa (do not restrict gamma_xy, Javier did.)
        FC(:)   = (/0.,0.,0.,0.,0.,0./)   ! no forces applied Nx,Ny,Nxy==0, all curvatures restrained: 0.,0.,0.
        call thermalFatigue(PROPS,STATEV,STRAN,DSTRAN,TEMP,DTEMP,TIME,DTIME,NSECV,NSTATV,NPROPS,NL,NTENS,FC,BC,targetTemp)
    else if (analysisType == 1) then
        ! strain-driven static/fatigue analysis w/targets: epsx, epsy, gammaxy, Tmin, 1    
        if ((ThermFatigG(1)<=0.d0 .or. ThermFatigG(2)<=0.d0) .and. ThermFatigG(3)>1.d0) stop "*** betaI and betaII must be > 1.0 for strain-driven fatigue ***"
        write(*,'(12A12)') 'Stress_x    ', 'Strain_idx  ', 'lambda i  '
        write(7,'(12A12)') 'Stress_x    ', 'Strain_idx  ', 'lambda i  '
        STRAN(:) = 0.d0
        DSTRAN(:) = targetS(:)/nSteps
        TEMP = targetTemp       !strain/stress analyses happen at Temp=targetTemp=constant
        DTEMP = 0.d0
        call strainFatigue(PROPS,STATEV,STRAN,DSTRAN,TEMP,DTEMP,TIME,DTIME,NSECV,NSTATV,NPROPS,NL,NTENS,targetS)
    else if (analysisType == 2) then
        ! stress-driven fatigue analysis w/targets: sigx, sigy, sigxy, Tmin, 2
        if ((ThermFatigG(1)<=0.d0 .or. ThermFatigG(2)<=0.d0) .and. ThermFatigG(3)>1.d0) stop "*** betaI and betaII must be > 1.0 for stress-driven fatigue ***"
        if ((ThermFatigG(4)<1.d0 .or. ThermFatigG(5)<1.d0) .and. ThermFatigG(3)>1.d0) stop "*** aI and aII must be > 1.0 for stress-driven fatigue ***"
        write(*,'(12A12)') 'Stress      ', 'Strain      ', 'lambda i  '
        write(7,'(12A12)') 'Stress      ', 'Strain      ', 'lambda i  '
        STRAN = 0.0D0
        DSTRAN = 0.D0 
        TEMP = targetTemp       !strain/stress analyses happen at Temp=targetTemp=constant
        DTEMP = 0.d0
        targetS(:)  = targetS(:)*LT !stress to stress-resultants conversion
        Nx = targetS(1); Ny = targetS(2); Nxy = targetS(3)
        BC(:)   = (/0,0,0,1,1,1/)               ! force, force, force, kappa, kappa, kappa 
        FC(:)   = (/Nx,Ny,Nxy,0.d0,0.d0,0.d0/)  ! Nx, Ny, Nxy applied, all curvatures restrained: 0., 0., 0.
        deltaFC(:) = FC(:)/nSteps   ! incremental loading
        FC(:) = deltaFC(:)          ! for the first step, then incremented in stressFatigue()
        call stressFatigue(PROPS,STATEV,STRAN,DSTRAN,TEMP,DTEMP,TIME,DTIME,NSECV,NSTATV,NPROPS,NL,NTENS,FC,BC,deltaFC,targetS,LT)
    else if (analysisType == -1) then
        ! adjust GIc using experimental lambda(epsilon_x) data, doing multiple analysisType=1  
        ! the initial GIc range is (0 to GIc) w/GIc read from .inp line 5 col 1 
        ! provide epsx to cover the range of experimental data 10/23/26
        write(7,'(12A12)') 'Stress      ', 'Strain      ', 'lambda i  '
        PROPS(45) = 1.d0 ! set Ncycles to 1 directly into PROPS in case it is wrong in the .inp file
        STRAN(:) = 0.d0
        DSTRAN(:) = targetS(:)/nSteps
        TEMP = targetTemp       !strain/stress analyses happen at Temp=targetTemp=constant
        DTEMP = 0.d0
        call goldenGIc(PROPS,STATEV,STRAN,DSTRAN,TEMP,DTEMP,TIME,DTIME,NSECV,NSTATV,NPROPS,NL,NTENS,targetS)    ! targetS==strain
        write(*,*) "* Use the new value of GIc to update your .inp file for strain-, stress-, or thermal-analysis *"
    else if (analysisType == -2) then
        ! calculate GIcPrime using experimental lambda(Ncycles) data
        write(7,'(12A12)') 'Stress      ', 'Strain      ', 'lambda i  '
        STRAN(:) = 0.d0     ! If fatigue data from stress tests, calculate targetS as stress/laminateEx0 (line 20 .inp file)
        DSTRAN(:) = targetS(:)/nSteps   ! use targetS = strain from line 20 .inp file 
        TEMP = targetTemp   !strain/stress analyses happen at Temp=targetTemp=constant
        DTEMP = 0.d0
        call bisectionGIcPrime(PROPS,STATEV,STRAN,DSTRAN,TEMP,DTEMP,TIME,DTIME,NSECV,NSTATV,NPROPS,NL,NTENS,targetS)    ! targetS==strain
    else
        stop "*** wrong analysisType ***"
    endif
    write(*,*) 'success, end of pgm ....'
    deallocate (LSS,PROPS,STATEV,Orientations,Thicknesses,trueStatev,Damages,StressInLaminaCS,EnergyState)
end program main
