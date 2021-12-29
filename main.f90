program testRDM_ED_NupNdw
  USE SCIFOR
  !< use COMMON_VARIABLES module to load global variables 
  USE COMMON_VARS
  !
  !< use AUXILIARY FUNCTIONS module to load some auxiliary procedures: binary_decomposition, sector_buildup, etc...
  USE AUX_FUNX
  !
  !< fortran-esque to force explicit declaration of all variables.
  implicit none
  !
  !< define variables to be used in THIS code. Variables are shared with
  !any procedure contained here.
  integer                          :: i,j,k,io,jo
  integer                          :: ibathUp,ibathDw,Ibath,Jbath
  integer                          :: iimpUp,iimpDw,iimp,jimp
  integer                          :: jimpUp,jimpDw
  integer                          :: Ntot,sIup,sIdw,sJup,sJdw
  integer                          :: iup,idw,kup,kdw,bup,bdw,jup,jdw
  integer                          :: Nup,Ndw
  integer                          :: LenInterUp,LenInterDw
  integer,allocatable              :: InterUp(:),InterDw(:)
  integer                          :: Dim
  integer,allocatable              :: ivecUp(:),ivecDw(:)
  integer,allocatable              :: jvecUp(:),jvecDw(:)
  integer                          :: unit
  character(len=32)                :: finput
  !
  !< CUSTOM TYPE: a user defined variable/structure which contains
  !the map between the states of a given sector S and those of the Fock space F
  type(sector_map)                     :: H,Hup,Hdw
  !
  real(8),dimension(:),allocatable     :: state_dvec
  real(8),dimension(:,:),allocatable   :: impurity_density_matrix
  real(8),dimension(:,:),allocatable   :: impurity_density_matrix_
  real(8),dimension(:,:),allocatable   :: reduced_density_matrix
  real(8),dimension(:,:),allocatable   :: big_density_matrix
  real(8),dimension(:,:),allocatable   :: small_density_matrix
  real(8)                              :: peso,delta
  real(8),dimension(:),allocatable     :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable     :: docc

  !  
  !< define some useful parameters:
  call parse_cmd_variable(finput,"FINPUT",default='inputDM')
  call parse_input_variable(Nbath,"Nbath",finput,default=3)
  call parse_input_variable(Nlat,"Nlat",finput,default=2)
  call parse_input_variable(Norb,"Norb",finput,default=1)
  call parse_input_variable(Nup,"Nup",finput,default=4)
  call parse_input_variable(Ndw,"Ndw",finput,default=4)
  call parse_input_variable(verbose,"verbose",finput,default=.false.)
  call parse_input_variable(fast,"fast",finput,default=.true.)
  call parse_input_variable(random,"random",finput,default=.false.)
  !< Number of impurity levels
  Nimp  = Nlat*Norb
  !< Number of spins: keep this low for better visualization
  Ns = Nimp*(1+Nbath)
  !< Total number of bath levels (per spin)
  Nbath_tot = Nimp*Nbath
  if(Nup+Ndw/=Ns)stop "ERROR: The given parameters are inconsistent."
  !
  !< Total number of levels for fermions. 1Fermion=2*spin
  Ntot  = 2*Ns
  !< set default output unit. 6==std.output in fortran-esque
  unit=6
  if(verbose)then
    write(unit,*)""
    write(unit,*)""
    write(unit,*)"Ns      =", Ns
    write(unit,*)"Ntot    =", Ntot
    write(unit,*)"2**Ns   =", 2**Ns
    write(unit,*)"2**Ntot =", 2**Ntot
    write(unit,*)""
    write(unit,*)""
  endif
  !
  !
  allocate(IvecUp(Ns),IvecDw(Ns))
  allocate(JvecUp(Ns),JvecDw(Ns))
  allocate(dens(Nimp),dens_up(Nimp),dens_dw(Nimp))
  allocate(docc(Nimp))





  !> BUILD SECTORS MAP (map + sp)
  write(unit,*)""
  write(unit,"(A,2I3,A,I10)")"Test --- sector: [ Nup Ndw ] = [",Nup,Ndw,"  ] ---- total dimension:",dim_sector_normal(Nup)*dim_sector_normal(Ndw)
  write(unit,*)""
  call build_sector_normal_sigma(Nup,Hup)
  call build_sector_normal_sigma(Ndw,Hdw)

  Dim = dim_sector_normal(Nup)*dim_sector_normal(Ndw)
  !> SECTOR PRINTOUT
  if(verbose)then
    !If dim is small enough list all the states in Human Readable form:
    if(Dim < 65000)then
      write(unit,'(A9,A,A9,A,A1,2A4,A6,2A4,A1)',advance='no')"i"," =>","K"," -- ","(","iup","idw",") => (","kup","kdw",")"
      write(unit,'(A)',advance='no')"|"
      write(unit,"("//str(Ns)//"A1)",advance="no")("b",j=1,Ns)
      write(unit,"(A1)",advance="no")">"
      write(unit,'(A)',advance='no')"|"
      write(unit,"("//str(Ns)//"A1)",advance="no")("b",j=1,Ns)
      write(unit,"(A1)",advance="no")">"
      write(unit,"(A3)",advance="no")" - "
      write(unit,"(A7,A1,"//str(Nimp)//"A1,A1)",advance="no")"kimpUp","|",("b",j=1,Nimp),">"
      write(unit,"(A3)",advance="no")"+"
      write(unit,"(A7,A1,"//str(Nbath_tot)//"A1,A1)",advance="no")"kbathUp","|",("b",j=Nimp+1,Ns),">"
      write(unit,"(A1)",advance="no")" "
      write(unit,"(A7,A1,"//str(Nimp)//"A1,A1)",advance="no")"kimpDw","|",("b",j=1,Nimp),">"
      write(unit,"(A3)",advance="no")"+"
      write(unit,"(A7,A1,"//str(Nbath_tot)//"A1,A1)",advance="no")"kbathDw","|",("b",j=Nimp+1,Ns),">"
      write(unit,"(A1)",advance="yes")""

      do idw = 1,dim_sector_normal(Ndw)
          kdw = Hdw%map(idw)
          ivecDw = bdecomp(kdw,Ns)
          iimpDw = get_imp_state(kdw)
          ibathDw= get_bath_state(kdw)
          !
          do iup = 1,dim_sector_normal(Nup)
            kup    = Hup%map(iup)
            ivecUp = bdecomp(kup,Ns)
            iimpUp = get_imp_state(kup)
            ibathUp= get_bath_state(kup)
            !
            i      = iup + (idw-1)*dim_sector_normal(Nup)
            k      = kup + 2**Ns*kdw
            !
            !
            write(unit,'(i9,A,i9,A)',advance='no')i," =>",k," -- "
            write(unit,'(A1,2i4,A6,2i4,A1)',advance='no')"(",iup,idw,") => (",kup,kdw,")"
            write(unit,'(A)',advance='no')"|"
            write(unit,"("//str(Ns)//"I1)",advance="no")(ivecUp(j),j=1,Ns)
            write(unit,"(A1)",advance="no")">"
            write(unit,'(A)',advance='no')"|"
            write(unit,"("//str(Ns)//"I1)",advance="no")(ivecDw(j),j=1,Ns)
            write(unit,"(A1)",advance="no")">"
            write(unit,"(A3)",advance="no")" - "
            write(unit,"(I7,A1,"//str(Nimp)//"I1,A1)",advance="no")iimpUp,"|",(ivecUp(j),j=1,Nimp),">"
            write(unit,"(A3)",advance="no")"+"
            write(unit,"(I7,A1,"//str(Nbath_tot)//"I1,A1)",advance="no")ibathUp,"|",(ivecUp(j),j=Nimp+1,Ns),">"
            write(unit,"(A1)",advance="no")" "
            write(unit,"(I7,A1,"//str(Nimp)//"I1,A1)",advance="no")iimpDw,"|",(ivecDw(j),j=1,Nimp),">"
            write(unit,"(A3)",advance="no")"+"
            write(unit,"(I7,A1,"//str(Nbath_tot)//"I1,A1)",advance="no")ibathDw,"|",(ivecDw(j),j=Nimp+1,Ns),">"
            write(unit,"(A1)",advance="yes")""
            !
            !
          enddo
      enddo
      write(unit,'(A9,A,A9,A,A1,2A4,A6,2A4,A1)',advance='no')"i"," =>","K"," -- ","(","iup","idw",") => (","Kup","Kdw",")"
      write(unit,'(A)',advance='no')"|"
      write(unit,"("//str(Ns)//"A1)",advance="no")("b",j=1,Ns)
      write(unit,"(A1)",advance="no")">"
      write(unit,'(A)',advance='no')"|"
      write(unit,"("//str(Ns)//"A1)",advance="no")("b",j=1,Ns)
      write(unit,"(A1)",advance="no")">"
      write(unit,"(A3)",advance="no")" - "
      write(unit,"(A7,A1,"//str(Nimp)//"A1,A1)",advance="no")"KimpUp","|",("b",j=1,Nimp),">"
      write(unit,"(A3)",advance="no")"+"
      write(unit,"(A7,A1,"//str(Nbath_tot)//"A1,A1)",advance="no")"KbathUp","|",("b",j=Nimp+1,Ns),">"
      write(unit,"(A1)",advance="no")" "
      write(unit,"(A7,A1,"//str(Nimp)//"A1,A1)",advance="no")"KimpDw","|",("b",j=1,Nimp),">"
      write(unit,"(A3)",advance="no")"+"
      write(unit,"(A7,A1,"//str(Nbath_tot)//"A1,A1)",advance="no")"KbathDw","|",("b",j=Nimp+1,Ns),">"
      write(unit,*)""
    endif
  endif



  
  !> EVALUATE THE INTERSECTION A_up \cap B_up BETWEEN ANY PAIR OF IMP_UP INDICES
  if(verbose)then 
    write(unit,*)""
    do iimpUp=0,Hup%sp%Nimp_state-1
      do jimpUp=0,Hup%sp%Nimp_state-1
          call sp_return_intersection(Hup%sp,iimpUp,jimpUp,InterUp,LenInterUp)
          write(unit,*)"Intersect IimpUP="//str(iimpUp)//", IimpUP="//str(jimpUp)
          write(unit,*)"Size(Intersection)",LenInterUp
          if(LenInterUp==0)cycle
          write(unit,*)"Intersection:"
          write(*,"(10000I5)")(InterUp(i),i=1,size(InterUp))
          write(unit,*)""
      enddo
    enddo
  endif




  !> ALLOCATE AND BUILD A VECTOR OF THE SECTOR: i.e. the linear combination of the basis using random/unit coefficient.
  !> NORMALIZED
  allocate(state_dvec(Dim))
  state_dvec=1d0
  if(random)call mt_random(state_dvec)
  state_dvec = state_dvec/sqrt(dot_product(state_dvec,state_dvec))


  !> EVALUATE OCCUPATION AND DOUBLE OCCUPATIONS FOR THIS STATE
  dens    = 0.d0
  dens_up = 0.d0
  dens_dw = 0.d0
  docc    = 0.d0
  do idw = 1,dim_sector_normal(Ndw)
     kdw = Hdw%map(idw)
     ivecDw = bdecomp(kdw,Ns)
     do iup = 1,dim_sector_normal(Nup)
        kup    = Hup%map(iup)
        ivecUp = bdecomp(kup,Ns)
        !
        i = iup + (idw-1)*dim_sector_normal(Nup)           
        peso=abs(state_dvec(i))**2
        !
        do io=1,Nimp
           dens(io)     = dens(io)      +  (ivecUp(io)+IvecDw(io))*peso
           dens_up(io)  = dens_up(io)   +  IvecUp(io)*peso
           dens_dw(io)  = dens_dw(io)   +  IvecDw(io)*peso
           docc(io)     = docc(io)      +  IvecUp(io)*IvecDw(io)*peso
        enddo
     enddo
  enddo



  !> ALLOCATE IMP DENSITY MATRIX (and a spare copy to test)
  allocate(impurity_density_matrix(4**Nimp,4**Nimp))
  allocate(impurity_density_matrix_(4**Nimp,4**Nimp))


  !>SLOW ALGORITHM:
  print*,"SLOW ALGORITHM"
  if(fast)then
      write(unit,"(A,2I3,A,I10)")"[skipping since <fast> option is enabled]"
      print*, " "
  else
      call start_timer()
      impurity_density_matrix=0d0
      do idw = 1,dim_sector_normal(Ndw)
         ivecDw = bdecomp(Hdw%map(idw),Ns)
         do iup = 1,dim_sector_normal(Nup)
            ivecUp = bdecomp(Hup%map(iup),Ns)
            Iimp = bjoin([IvecUp(1:Nimp),IvecDw(1:Nimp)],2*Nimp) + 1
            Ibath= bjoin([IvecUp(Nimp+1:),IvecDw(Nimp+1:)],2*Nbath_tot) + 1
            !
            do jdw = 1,dim_sector_normal(Ndw)
               jvecDw = bdecomp(Hdw%map(jdw),Ns)
               do jup = 1,dim_sector_normal(Nup)
                  jvecUp = bdecomp(Hup%map(jup),Ns)
                  Jimp = bjoin([JvecUp(1:Nimp),JvecDw(1:Nimp)],2*Nimp) + 1
                  Jbath= bjoin([JvecUp(Nimp+1:),JvecDw(Nimp+1:)],2*Nbath_tot) + 1
                  !
                  if(Jbath/=Ibath)cycle
                  i = iup + (idw-1)*dim_sector_normal(Nup)
                  j = jup + (jdw-1)*dim_sector_normal(Nup)
                  impurity_density_matrix(Iimp,Jimp) = impurity_density_matrix(Iimp,Jimp) + &
                        state_dvec(i)*state_dvec(j)
               enddo
            enddo
         enddo
      enddo
      impurity_density_matrix_ = impurity_density_matrix
      if(Nimp<=2)then
            do i=1,4**Nimp
               !PRINT FULL MATRIX
               write(*,"(1000F7.3)")(impurity_density_matrix(i,j),j=1,4**Nimp)
            enddo
      endif
      call stop_timer()
      print*, " "
  endif


  !>FAST ALGORITHM:
  print*,"FAST ALGORITHM"
  call start_timer()
  impurity_density_matrix=0d0
  do IimpUp=0,Hup%sp%Nimp_state-1
     do JimpUp=0,Hup%sp%Nimp_state-1
        !Finding the unique bath states connecting IimpUp and JimpUp -> InterUp(:)
        call sp_return_intersection(Hup%sp,IimpUp,JimpUp,InterUp,LenInterUp)
        if(LenInterUp==0)cycle  !there are no bath states intersecting IimpUp,JimpUp
        do IimpDw=0,Hdw%sp%Nimp_state-1
           do JimpDw=0,Hdw%sp%Nimp_state-1
              !Finding the unique bath states connecting IimpDw and JimpDw -> InterDw(:)
              call sp_return_intersection(Hdw%sp,IimpDw,JimpDw,InterDw,LenInterDw)
              if(LenInterDw==0)cycle  !there are no bath states intersecting IimpDw,JimpDw
              !------------------------------------------------------------------------------------
              !Arrived here we have finally determined the {IbathUp,IbathDw} states to trace on.
              !We shall use them to build back the Fock space states (IimpSIGMA,IbathSIGMA) and
              !through those search the map to retrieve the formal labels i=(iUP,iDW), j=(jUP,jDW)
              !------------------------------------------------------------------------------------
              do bup=1,LenInterUp !=== >>> THE TRACE OVER THE BATH <<< ============================
                 IbathUp = InterUp(bup)
                 do bdw=1,LenInterDw
                    IbathDw = InterDw(bdw)
                    !-----------------------------------------------------
                    !Allowed spin Fock Indices: 
                    Iup= binary_search(Hup%map,IimpUp + 2**Nimp*IbathUp)
                    Idw= binary_search(Hdw%map,IimpDw + 2**Nimp*IbathDw)
                    !Global sector index:
                    i  = Iup + (Idw-1)*dim_sector_normal(Nup)
                    !-----------------------------------------------------
                    !Allowed spin Fock Jndices: 
                    Jup= binary_search(Hup%map,JimpUp + 2**Nimp*IbathUp)
                    Jdw= binary_search(Hdw%map,JimpDw + 2**Nimp*IbathDw)
                    !Global sector jndex:
                    j  = Jup + (Jdw-1)*dim_sector_normal(Nup)
                    !-----------------------------------------------------
                    !Final Lin-Gubernatis composition for the impurity
                    io = (IimpUp+1) + 2**Nimp*IimpDw
                    jo = (JimpUp+1) + 2**Nimp*JimpDw
                    !------------------------------------------------------------------
                    !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                    impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                         state_dvec(i)*state_dvec(j)
                    !------------------------------------------------------------------
                 enddo
              enddo !==============================================================================
              !
           enddo
        enddo
        !
     enddo
  enddo
  if(Nimp<=2)then
      do i=1,4**Nimp
         !PRINT FULL MATRIX
         write(*,"(1000F7.3)")(impurity_density_matrix(i,j),j=1,4**Nimp)
      enddo
  endif
  call stop_timer()
  print*, " "

  !>CROSSCHECK ALGORITHMS
  if(.not.fast)then
    print*,    "***************************************************"
    if( any(abs(impurity_density_matrix_-impurity_density_matrix)/=0d0) )then
      print*, "ERROR: SLOW and FAST algorithms fail to match"
    else
      print*, "SLOW and FAST algorithms match to machine-precision"
    endif
    print*,    "***************************************************"
  endif

  !>SUBTRACING
  write(*,*)
  write(*,*) "=================="
  write(*,*) "SUBTRACING-ROUTINE"
  write(*,*) "=================="
  write(*,*)
  print*,"ALL REDUCED DENSITY MATRICES"
  big_density_matrix = impurity_density_matrix
  do k=1,Nlat-1
     !<From cluster-dm to all reduced dms [increasing Ntrace]
     call subtrace(reduced_density_matrix,impurity_density_matrix,Nlat,Nlat-k)
     DIM = size(reduced_density_matrix(:,1)) 
     write(*,"(A,2I3,A,I10)") ">rank-"//str(DIM)//" [DIRECT]"
     if(verbose)then 
         do i=1,DIM
            !PRINT REDUCED MATRIX [DIRECT]
            write(*,"(1000F7.3)")(reduced_density_matrix(i,j),j=1,DIM) 
         enddo
     else
         print*, "(no print)"
     endif
     !<From cluster-dm to all reduced dms [site by site]
     call subtrace(small_density_matrix,big_density_matrix,Nlat-k+1,Nlat-k)
     DIM = size(small_density_matrix(:,1)) 
     write(*,"(A,2I3,A,I10)") ">rank-"//str(DIM)//" [ITERATIVE]" 
     if(verbose)then 
         do i=1,DIM
            !PRINT REDUCED MATRIX [ITERATIVE]
            write(*,"(1000F7.3)")(small_density_matrix(i,j),j=1,DIM) 
         enddo
      else
         print*, "(no print)" 
      endif
     !>CROSSCHECK ALGORITHMS
     print*,   "**************************************************************"
     delta = maxval(abs(reduced_density_matrix-small_density_matrix))
     if(delta > 0.00001d0)then
       print*, "ERROR: DIRECT and ITERATIVE subtracings fail to match"
     elseif(delta == 0.d0)then
       print*, "DIRECT and ITERATIVE subtracings match up to machine precision" 
     else
       print*, "DIRECT and ITERATIVE subtracings match up to:",real(delta)
     endif
    print*,    "**************************************************************"
    !
    deallocate(big_density_matrix)
    big_density_matrix = small_density_matrix
    !
  enddo
  !
  !
  if(Norb==1)then
      print*, " "
      print*, "==============================================================="
      print*, "BENCHMARK: LOCAL DENSITY-MATRIX versus SEMI-ANALYTICAL FORMULAE"
      print*, "==============================================================="
      print*, "------ ESTIMATE ------ | -------- ERROR ----------- "
      print*, 1-dens_up(1)-dens_dw(1)+docc(1),abs(1-dens_up(1)-dens_dw(1)+docc(1)-reduced_density_matrix(1,1))
      print*, dens_up(1)-docc(1),abs(dens_up(1)-docc(1)-reduced_density_matrix(2,2))
      print*, dens_dw(1)-docc(1),abs(dens_dw(1)-docc(1)-reduced_density_matrix(3,3))
      print*, docc(1),abs(docc(1)-reduced_density_matrix(4,4))
      print*, "==============================================================="
      print*, " "
  endif

  call map_deallocate(Hup)
  call map_deallocate(Hdw)
  write(unit,*)""
  write(unit,*)""
  write(unit,*)""







end program testRDM_ED_NupNdw








