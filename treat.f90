PROGRAM treat
  !----------------------------------------------------------!
  ! TREAT - tail-regression estimator analysis toolbox       !
  ! =====                                                    !
  ! Tool to evaluate the tail-regression estimator of the    !
  ! mean of a statistical sample that follows a probability  !
  ! distribution with Pareto-type tails of known tail index. !
  !                                                          !
  ! PLR 05.2016                                              !
  !----------------------------------------------------------!
  IMPLICIT NONE

  ! Derived types
  ! =============

  ! * Data (as read) and immutable statistical properties.
  TYPE dataset_type
    CHARACTER(1024) :: description='no data loaded'
    LOGICAL :: have_w=.false.
    INTEGER :: M=0
    DOUBLE PRECISION :: inv_total_w=0.d0
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION,ALLOCATABLE :: A(:),w(:),q(:)
    ! Stats.
    DOUBLE PRECISION :: A_mean=0.d0,A_stderr=0.d0,A_var=0.d0,A_err_var=0.d0,&
       &A_median=0.d0
    ! Parameters for kernel density estimate of P(A).
    ! NB, this is only used for plots, so we use a simple variable-width
    ! kernel which is not guaranteed to yield accurate results.
    DOUBLE PRECISION :: h_KDE_centre=0.d0
    DOUBLE PRECISION :: A_KDE_end(2)=0.d0,A_KDE_onset(2)=0.d0,&
       &h_KDE_end(2)=0.d0
    DOUBLE PRECISION :: max_x=sqrt(-2.d0*log(epsilon(1.d0)))
  END TYPE dataset_type

  ! * Tail form parameters.
  TYPE tail_form_type
    ! Left/right division.
    CHARACTER(32) :: A_centre_def='median'
    DOUBLE PRECISION :: A_centre=0.d0
    INTEGER :: Mhalf(2)=0
    LOGICAL :: self_consistent_centre=.false.
    ! Exponents.
    INTEGER :: nparam(2)=0
    DOUBLE PRECISION :: mu(2)=0.d0, global_delta_mu(2)=1.d0
    ! Fit range.
    CHARACTER :: def_anchor(2)=''
    DOUBLE PRECISION :: anchor(2)=0.d0
    ! Constraints.
    INTEGER :: nconstraint=0
    INTEGER, ALLOCATABLE :: lhs_iparam(:),lhs_itail(:),nrhs(:),&
       &rhs_iparam(:,:),rhs_itail(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: rhs_coeff(:,:)
    ! Whether moments are defined.
    INTEGER :: momdef=0
    ! Result (for self-consistent A_centre).
    LOGICAL :: A_result_set=.false.
    DOUBLE PRECISION :: A_result=0.d0, dA_result=0.d0
    ! Which type of x-dependent fit weight prefactor to use (hill/none).
    CHARACTER(16) :: fit_weight_type='hill'
  END TYPE tail_form_type

  ! Constraint set.
  TYPE constraint_set_type
    INTEGER, ALLOCATABLE :: ndependant(:,:), dependant_itail(:,:,:), &
       &dependant_iparam(:,:,:), ncombine(:,:), combine_itail(:,:,:), &
       &combine_iparam(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: constraint_coeff(:,:,:,:)
  END TYPE constraint_set_type

  ! Plot grid definition.
  TYPE plot_grid_type
    ! Points.
    INTEGER :: npoint=0
    CHARACTER, ALLOCATABLE :: def_point(:)
    DOUBLE PRECISION, ALLOCATABLE :: point(:)
  END TYPE plot_grid_type

  ! Assessment grid definition.
  TYPE assess_grid_type
    ! Per-tail anchors.
    INTEGER :: nanchor=0
    CHARACTER, ALLOCATABLE :: def_anchor(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: anchor(:,:)
    ! Per-tail nparams.
    INTEGER :: nnparam=0,max_nparam=0
    INTEGER, ALLOCATABLE :: nparam(:,:)
  END TYPE assess_grid_type


  ! * Bootstrap parameters for analysis.
  TYPE bs_setup_type
    INTEGER :: nbs=256,nbs_proc=256
  END TYPE bs_setup_type

  ! * Definition of model distribution for tests.
  TYPE model_P_type
    INTEGER :: nterm=0
    DOUBLE PRECISION,ALLOCATABLE :: centrevec(:),muvec(:),avec(:),normvec(:),&
       &lambdavec(:)
  END TYPE model_P_type

  ! Constants
  ! =========
  ! * Labels identifying which moments of the distribution are defined
  !   according to value of mu.
  INTEGER, PARAMETER :: &
     &MOMDEF_NONE=0,&
     &MOMDEF_MOM0=1,&
     &MOMDEF_MOM1_MAYBE=2,&
     &MOMDEF_MOM1=3,&
     &MOMDEF_MOM2=4

  ! Numerical constants.
  DOUBLE PRECISION, PARAMETER :: pi=4.d0*atan(1.d0)
  DOUBLE PRECISION, PARAMETER :: half_sqrt_pi=0.5d0*sqrt(pi)
  DOUBLE PRECISION, PARAMETER :: inv_sqrt_pi=1.d0/sqrt(pi)

  ! Global variables for MPI functionality.
  INTEGER MPI_NPROC,MPI_IPROC
  LOGICAL MPI_MASTER
  ! Load MPI header file.
  INCLUDE "mpif.h"

  ! Initialize MPI.
  call init_mpi(MPI_NPROC,MPI_IPROC,MPI_MASTER)
  ! Initialize random number generator.
  call init_random(1)
  ! Call main driver.
  call treat_main()
  ! Finish MPI.
  call finish_mpi()


CONTAINS


  SUBROUTINE treat_main
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! User interaction variables.
    CHARACTER(8192) command,token
    LOGICAL apply_assessment
    INTEGER ierr,ierror,ifield,it1,ipos,itail
    INTEGER iparam,jparam,jtail,irhs,nrhs,iconstraint
    DOUBLE PRECISION t1,t2
    ! Variables for inspect & load commands.
    INTEGER nline,ncolumn,Mstart,Mend,ipos_A,ipos_w,icol_w,ierr1,ierr2
    INTEGER, ALLOCATABLE :: icol_A(:)
    CHARACTER(256) fname
    ! Model P parameters.
    INTEGER iterm,nsample
    TYPE(model_P_type) model_P
    ! Selected tail form.
    TYPE(tail_form_type) tail_form,tail_form_default
    ! Bootstrap parameters.
    TYPE(bs_setup_type) bs_setup, bs_setup_default
    ! Dataset.
    TYPE(dataset_type) dataset
    ! Variables for plot command.
    TYPE(plot_grid_type) plot_grid
    TYPE(assess_grid_type) assess_grid
    CHARACTER(16) tchar
    CHARACTER(256) plot_fname,plot_P_fname,plot_fit_fname
    INTEGER plot_P_ngrid,plot_fit_ngrid,ngrid
    ! Input echo.
    LOGICAL input_echo

    ! Write header.
    if(MPI_MASTER)then
      write(6,'(a)')'=================================================='
      write(6,'(a)')'TREAT - tail-regression estimator analysis toolbox'
      write(6,'(a)')'=================================================='
      write(6,'()')
      if(MPI_NPROC>1)then
        write(6,'(a)')''//trim(i2s(MPI_NPROC))//' MPI processes available.'
        write(6,'()')
      endif
      write(6,'(a)')'Type "help" for a list of commands.'
      write(6,'()')
    endif

    ! Initialize.
    input_echo=.false.
    call refresh_bs_setup(bs_setup) ! account for MPI_NPROC

    ! Loop over user actions.
    user_loop: do

      ! Read user command.
      if(MPI_MASTER)then
        write(6,'(a)',advance='no')'TREAT> '
        call flush(6)
        read(5,'(a)',iostat=ierr)command
        if(ierr==0)then
          command=adjustl(command)
          if(input_echo.and.command(1:1)/='#')write(6,'(a)')trim(command)
        else
          write(6,'()')
        endif
      endif
      call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
      if(ierr/=0)call quit()
      call mpi_bcast(command,len(command),mpi_character,0,mpi_comm_world,&
         &ierror)
      if(ierror/=0)call quit()
      if(MPI_MASTER)write(6,'()')
      if(command(1:1)=='#')cycle

      ! Execute command.
      select case(trim(field(1,command)))

      case('inspect')
        ! Report number of lines and columns in data file.
        if(MPI_MASTER)then
          fname=trim(field(2,command))
          call parse_file(fname,.true.,nline,ncolumn)
          select case(nline)
          case(-1)
            call master_msg('Could not open file "'//trim(fname)//'".')
          case(-2)
            call master_msg('Problem reading file "'//trim(fname)//'".')
          case(-3)
            call master_msg('Column count problem in file "'//trim(fname)//'".')
          case(0)
            call master_msg('File "'//trim(fname)//'" contains no useful data.')
          case default
            call master_msg('File "'//trim(fname)//'" contains '//&
               &trim(i2s(nline))//' lines with '//trim(i2s(ncolumn))//&
               &' columns.')
          end select
        endif ! MPI_MASTER

      case('load')
        ! Load data from specified columns of data file.
        fname=trim(field(2,command))

        ! Parse subcommands.
        if(allocated(icol_A))deallocate(icol_A)
        allocate(icol_A(1))
        icol_A(1)=1
        icol_w=0
        Mstart=0
        Mend=0
        ifield=2
        do
          ifield=ifield+1
          select case(trim(field(ifield,command)))
          case('type')
            call parse_type_string(field(ifield+1,command),ipos_A,ipos_w,ierr)
            if(ierr/=0)then
              call master_msg('Syntax error in load command: unrecognized &
                 &dataset type.')
              cycle user_loop
            endif
            if(trim(field(ifield+2,command))=='using')then
              ierr1=0
              ierr2=0
              if(allocated(icol_A))deallocate(icol_A)
              icol_w=0
              if(ipos_A>0)call parse_sum(trim(field(ifield+2+ipos_A,command)),&
                 &icol_A)
              if(ipos_w>0)icol_w=int_field(ifield+2+ipos_w,command,ierr2)
              if(.not.allocated(icol_A).or.ierr2/=0)then
                call master_msg('Problem parsing "using" indices.')
                cycle user_loop
              endif
              ifield=ifield+2
              if(ipos_A>0)ifield=ifield+1
              if(ipos_w>0)ifield=ifield+1
            else
              if(allocated(icol_A))deallocate(icol_A)
              allocate(icol_A(1))
              icol_A(1)=ipos_A
              icol_w=ipos_w
              ifield=ifield+1
            endif
          case('trim')
            Mstart=int_field(ifield+1,command,ierr1)
            Mend=int_field(ifield+2,command,ierr2)
            if(ierr1/=0.or.ierr2/=0)then
              call master_msg('Problem parsing "trim" line indices.')
              cycle user_loop
            endif
            if(Mstart<1.or.Mend<Mstart)then
              call master_msg('"trim" line indices out of range or order.')
              cycle user_loop
            endif
            ifield=ifield+2
          case('')
            exit
          case default
            call master_msg('Syntax error in load command: unknown &
               &subcommand "'//trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo ! ifield

        ! Parse file to get number of columns and (if no line range was
        ! provided) number of lines.
        if(MPI_MASTER)call parse_file(fname,Mend==0,nline,ncolumn)
        call mpi_bcast(nline,1,mpi_integer,0,mpi_comm_world,ierror)
        select case(nline)
        case(-1)
          call master_msg('Could not open file "'//trim(fname)//'".')
          cycle user_loop
        case(-2)
          call master_msg('Problem reading file "'//trim(fname)//'".')
          cycle user_loop
        case(-3)
          call master_msg('Column count problem in file "'//trim(fname)//'".')
          cycle user_loop
        case(0)
          if(Mend==0)then
            call master_msg('File "'//trim(fname)//'" contains no useful data.')
            cycle user_loop
          endif
        case default
          if(MPI_MASTER)write(6,'(a)')'File "'//trim(fname)//'" contains '//&
             &trim(i2s(nline))//' lines with '//trim(i2s(ncolumn))//&
             &' columns.'
        end select
        call mpi_bcast(ncolumn,1,mpi_integer,0,mpi_comm_world,ierror)

        ! Set line range.
        if(nline/=0)then
          Mstart=1
          Mend=nline
        endif
        ! Check column indices.
        if(any(icol_A>ncolumn).or.any(icol_A<0).or.icol_w>ncolumn.or.&
           &icol_w<0)then
          call master_msg('Column indices out of range.')
          cycle user_loop
        endif

        ! Load data.
        if(MPI_MASTER)call read_file(fname,icol_A,icol_w,Mstart,Mend,dataset,&
           &ierr)
        call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
        if(ierr/=0)then
          if(MPI_MASTER)then
            if(allocated(dataset%A))deallocate(dataset%A)
            if(allocated(dataset%w))deallocate(dataset%w)
            if(allocated(dataset%indx))deallocate(dataset%indx)
            if(allocated(dataset%q))deallocate(dataset%q)
            call master_msg('Problem reading file - column or line-range &
               &problem?')
          endif
          cycle user_loop
        endif
        if(MPI_MASTER)call master_msg('Loaded '//trim(i2s(dataset%M))//&
           &' data.')

        ! Pass data to slave nodes.
        call mpi_bcast(dataset%description,len(dataset%description),&
           &mpi_character,0,mpi_comm_world,ierror)
        call mpi_bcast(dataset%M,1,mpi_integer,0,mpi_comm_world,ierror)
        call mpi_bcast(dataset%have_w,1,mpi_logical,0,mpi_comm_world,ierror)
        if(.not.MPI_MASTER)then
          if(allocated(dataset%A))deallocate(dataset%A)
          if(allocated(dataset%w))deallocate(dataset%w)
          if(allocated(dataset%indx))deallocate(dataset%indx)
          if(allocated(dataset%q))deallocate(dataset%q)
          allocate(dataset%A(dataset%M))
          if(dataset%have_w)then
            allocate(dataset%w(dataset%M))
            dataset%w=1.d0
          endif
        endif
        call mpi_bcast(dataset%A,dataset%M,mpi_double_precision,0,&
           &mpi_comm_world,ierror)
        if(dataset%have_w)call mpi_bcast(dataset%w,dataset%M,&
           &mpi_double_precision,0,mpi_comm_world,ierror)
        call refresh_dataset(dataset)
        tail_form%A_result_set=.false.
        call refresh_centre(dataset,tail_form)

      case('unload')
        if(dataset%M<1)then
          call master_msg('No data to unload.')
        else
          dataset%description='no data loaded'
          dataset%have_w=.false.
          dataset%M=0
          if(allocated(dataset%A))deallocate(dataset%A)
          if(allocated(dataset%w))deallocate(dataset%w)
          if(allocated(dataset%indx))deallocate(dataset%indx)
          if(allocated(dataset%q))deallocate(dataset%q)
          call master_msg('Data unloaded.')
        endif

      case('generate')
        ! Generate dataset from model distribution.
        ! Get nsample.
        it1=parse_int(field(2,command),ierr)
        if(ierr/=0)then
          call master_msg('Could not parse <nsample> in generate command.')
          cycle user_loop
        endif
        if(it1<1)then
          call master_msg('<nsample> must be a positive integer.')
          cycle user_loop
        endif
        nsample=it1
        ! Get nterm.
        model_P%nterm=nfield(command)-2
        if(model_P%nterm<1)then
          call master_msg('No model specified.')
          cycle user_loop
        endif
        if(allocated(model_P%centrevec))deallocate(model_P%centrevec)
        if(allocated(model_P%muvec))deallocate(model_P%muvec)
        if(allocated(model_P%avec))deallocate(model_P%avec)
        if(allocated(model_P%normvec))deallocate(model_P%normvec)
        if(allocated(model_P%lambdavec))deallocate(model_P%lambdavec)
        allocate(model_P%centrevec(model_P%nterm),&
           &model_P%muvec(model_P%nterm),model_P%avec(model_P%nterm),&
           &model_P%normvec(model_P%nterm),model_P%lambdavec(model_P%nterm))
        ! Parse model parameters.
        do iterm=1,model_P%nterm
          ifield=iterm+2
          call parse_model_P_token(trim(field(ifield,command)),iterm,model_P,&
             &ierr)
          if(ierr/=0)then
            if(ierr==2)then
              call master_msg('<mu> must be >1 in model P generation.')
            else
              call master_msg('Syntax error in model P, term #'//&
                &trim(i2s(iterm)))
            endif
            cycle user_loop
          endif
        enddo ! iterm
        if(model_P%nterm<1)then
          call master_msg('No terms specified for model.')
          cycle user_loop
        endif
        ! Normalize terms.
        call normalize_model_P(model_P)
        ! Set up dataset.
        call generate_model_P(model_P,nsample,dataset)
        call refresh_dataset(dataset)
        tail_form%A_result_set=.false.
        call refresh_centre(dataset,tail_form)

      case('save')
        ! Dump dataset.
        if(MPI_MASTER)then
          if(dataset%M<1)then
            call master_msg('No data to save.')
            cycle user_loop
          endif
          fname=trim(field(2,command))
          if(len_trim(fname)<1)then
            call master_msg('No filename provided.')
            cycle user_loop
          endif
          call dump_dataset(dataset,fname)
          call master_msg('Data dumped to "'//trim(fname)//'".')
        endif

      case('report')
        ! Check there is something to report.
        if(dataset%M<1)then
          call master_msg('No datasets loaded.')
          cycle user_loop
        endif
        ! Select report to display.
        select case(trim(field(2,command)))
        case('stats')
          ! Report 'standard' mean, median, etc, of dataset.
          if(MPI_MASTER)then
            write(6,'(a)')'Basic statistics'
            write(6,'(a)')'================'
            write(6,'("STAT  Mean     ",2(1x,es20.12))')dataset%A_mean,&
               &dataset%A_stderr
            write(6,'("STAT  Variance ",2(1x,es20.12))')dataset%A_var,&
               &dataset%A_err_var
            write(6,'("STAT  Median   ",1x,es20.12)')dataset%A_median
            write(6,'("STAT  Max+     ",1x,es20.12)')maxval(dataset%A)
            write(6,'("STAT  Max-     ",1x,es20.12)')-maxval(-dataset%A)
            write(6,'()')
          endif ! MPI_MASTER
        case('corr')
          ! Report serial correlation info.
          if(MPI_MASTER)call report_serial_correlation(dataset)
        case('')
          call master_msg('No report specified.')
          cycle user_loop
        case default
          call master_msg('Unknown report "'//trim(field(2,command))//'".')
          cycle user_loop
        end select

      case('assess')
        ! Check there are data to assess.
        if(dataset%M<1)then
          call master_msg('No datasets loaded.')
          cycle user_loop
        endif

        ! Check assessment and requirements.
        select case(trim(field(2,command)))
        case('mu')
          continue
        case('fit','nparam,anchor','anchor,nparam')
          if(all(tail_form%mu<=0.d0))then
            call master_msg('mu is not set.')
            cycle user_loop
          endif
        case('')
          call master_msg('No assessemnt specified.')
          cycle user_loop
        case default
          call master_msg('Unknown assessment "'//trim(field(2,command))//'".')
          cycle user_loop
        end select

        ! Parse options.
        call initialize_plot_grid(plot_grid)
        call initialize_assess_grid(assess_grid)
        plot_fname=''
        ifield=3
        apply_assessment=.false.
        do
          select case(trim(field(ifield,command)))
          case('at')
            itail=0
            select case(trim(field(ifield+1,command)))
            case('left')
              itail=1
              ifield=ifield+1
            case('right')
              itail=2
              ifield=ifield+1
            end select
            token=field(ifield+1,command)
            if(len_trim(token)<1)then
              call master_msg('"at" must be followed by a range expression.')
              cycle user_loop
            endif
            select case(trim(field(2,command)))
            case('mu')
              call parse_plot_grid_points(dataset,tail_form,itail,token,&
                 &plot_grid,ierr)
            case('fit','nparam,anchor','anchor,nparam')
              if(itail/=0)then
                call master_msg('No left/right selection for "at" subcommand &
                   &in anchor assessment.')
                cycle user_loop
              endif
              call parse_assess_grid_anchor(dataset,tail_form,token,&
                 &assess_grid,ierr)
            end select
            if(ierr/=0)then
              call master_msg('Could not parse range "'//trim(token)//'".')
              cycle user_loop
            endif
            ifield=ifield+1
          case('for')
            select case(trim(field(2,command)))
            case('mu')
              call master_msg('No "for" subcommand for "mu" assessment.')
              cycle user_loop
            case('fit','nparam,anchor','anchor,nparam')
              token=field(ifield+1,command)
              if(len_trim(token)<1)then
                call master_msg('"for" must be followed by a range &
                   &expression.')
                cycle user_loop
              endif
              call parse_assess_grid_nparam(token,assess_grid,ierr)
              if(ierr/=0)then
                call master_msg('Could not parse range "'//trim(token)//'".')
                cycle user_loop
              endif
            end select
            ifield=ifield+1
          case('to')
            plot_fname=trim(field(ifield+1,command))
            if(len_trim(plot_fname)<1)then
              call master_msg('Must provide a filename after "to".')
              cycle user_loop
            endif
            ifield=ifield+1
          case('apply')
            apply_assessment=.true.
          case('')
            exit
          case default
            call master_msg('Unknown subcommand "'//&
               &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
          ifield=ifield+1
        enddo

        ! Check we have all we need.
        select case(trim(field(2,command)))
        case('mu')
          if(plot_grid%npoint<1)then
            call master_msg('Must set a point range.')
            cycle user_loop
          endif
        case('fit','nparam,anchor','anchor,nparam')
          if(assess_grid%nanchor<1)then
            call master_msg('Must set an anchor range.')
            cycle user_loop
          endif
          if(assess_grid%nnparam<1)then
            call master_msg('Must set an nparam range.')
            cycle user_loop
          endif
        end select

        ! Assess.
        select case(trim(field(2,command)))
        case('mu')
          call assess_mu(dataset,tail_form,bs_setup,plot_grid,&
             &trim(plot_fname),apply_assessment)
        case('fit','nparam,anchor','anchor,nparam')
          call assess_fit(dataset,tail_form,bs_setup,assess_grid,&
             &trim(plot_fname),apply_assessment)
        end select
        if(apply_assessment)call refresh_momdef(tail_form)

      case('plot')
        ! Check there is something to plot.
        if(dataset%M<1)then
          call master_msg('No datasets loaded.')
          cycle user_loop
        endif

        ! Initialize plot_grid.
        call initialize_plot_grid(plot_grid)

        ! Check plot subject.
        plot_fname='plot.dat'
        select case(trim(field(2,command)))
        case('model')
          if(model_P%nterm<1)then
            call master_msg('No model defined.')
            cycle user_loop
          endif
          if(all(tail_form%mu<=0.d0))then
            call master_msg('mu is not set.')
            cycle user_loop
          endif
          plot_fname='TREAT_model_plot.dat'
        case('')
          call master_msg('No plot subject specified.')
          cycle user_loop
        case default
          call master_msg('Unknown plot subject "'//trim(field(2,command))//&
             &'".')
          cycle user_loop
        end select

        ! Parse options.
        ifield=3
        do
          select case(trim(field(ifield,command)))
          case('at')
            itail=0
            select case(trim(field(ifield+1,command)))
            case('left')
              itail=1
              ifield=ifield+1
            case('right')
              itail=2
              ifield=ifield+1
            end select
            token=field(ifield+1,command)
            if(len_trim(token)<1)then
              call master_msg('"at" must be followed by a range expression.')
              cycle user_loop
            endif
            call parse_plot_grid_points(dataset,tail_form,itail,token,&
               &plot_grid,ierr)
            if(ierr/=0)then
              call master_msg('Could not parse range "'//trim(token)//'".')
              cycle user_loop
            endif
            ifield=ifield+1
          case('to')
            plot_fname=trim(field(ifield+1,command))
            if(len_trim(plot_fname)<1)then
              call master_msg('Must provide a filename after "to".')
              cycle user_loop
            endif
            ifield=ifield+1
          case('')
            exit
          case default
            call master_msg('Unknown subcommand "'//&
               &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
          ifield=ifield+1
        enddo

        ! Verify we have a plot_grid.
        if(plot_grid%npoint<1)then
          call master_msg('Must set a plot range.')
          cycle user_loop
        endif

        ! Make the plot.
        select case(trim(field(2,command)))
        case('model')
          call plot_model(model_P,plot_grid,trim(plot_fname))
          call master_msg('Plot written to "'//trim(plot_fname)//'".')
        end select

      case('evaluate')
        ! Check there is anything to evaluate.
        if(dataset%M<1)then
          call master_msg('No datasets loaded.')
          cycle user_loop
        endif

        ! Select what to evaluate.
        select case(trim(field(2,command)))
        case('mean','TRE')
          continue
        case default
          call master_msg('Unknown object "'//trim(field(2,command))//&
             &'" to evaluate.')
          cycle user_loop
        end select

        ! Get additional options.
        plot_P_fname=''
        plot_P_ngrid=0
        plot_fit_fname=''
        plot_fit_ngrid=0
        ifield=2
        do
          select case(trim(field(ifield+1,command)))
          case('plot')
            ifield=ifield+1
            select case(trim(field(ifield+1,command)))
            case('P')
              ifield=ifield+1
              tchar='P'
              ngrid=500
              plot_fname='TREAT_P_plot.dat'
            case('fit')
              ifield=ifield+1
              tchar='fit'
              ngrid=100
              plot_fname='TREAT_fit_plot.dat'
            case default
              call master_msg('Unknown plot subject "'//&
                 &trim(field(ifield+1,command))//'" to evaluate.')
              cycle user_loop
            end select
            do
              select case(trim(field(ifield+1,command)))
              case('grid')
                ifield=ifield+2
                ngrid=int_field(ifield,command,ierr1)
                if(ierr1/=0)then
                  call master_msg('Must provide an integer after "grid".')
                  cycle user_loop
                endif
                if(ngrid<2)then
                  call master_msg('Number of grid points must be > 1.')
                  cycle user_loop
                endif
              case('to')
                ifield=ifield+2
                plot_fname=trim(field(ifield,command))
                if(len_trim(plot_fname)<1)then
                  call master_msg('Must provide a filename after "to".')
                  cycle user_loop
                endif
              case default
                exit
              end select
            enddo
            select case(trim(tchar))
            case('P')
              plot_P_fname=plot_fname
              plot_P_ngrid=ngrid
            case('fit')
              plot_fit_fname=plot_fname
              plot_fit_ngrid=ngrid
            end select
          case('')
            exit
          case default
           call master_msg('Unrecognized subcommand "'//&
              &trim(field(ifield+1,command))//'" in evaluate command.')
           cycle user_loop
          end select
        enddo

        ! Evaluate.
        call evaluate_TRE(dataset,tail_form,bs_setup,trim(plot_P_fname),&
           &plot_P_ngrid,trim(plot_fit_fname),plot_fit_ngrid)

      case('set')

        ! Check if there is a tail selection subcommand.
        itail=0
        ifield=2
        select case(trim(field(ifield,command)))
        case('left')
          itail=1
          ifield=ifield+1
        case('right')
          itail=2
          ifield=ifield+1
        end select

        ! Set variables.
        select case(trim(field(ifield,command)))

        case('Ac','centre')
          ! Centre of distribution.
          if(itail/=0)then
            call master_msg('No left/right selection for variable Ac.')
            cycle user_loop
          endif
          select case(trim(field(ifield+1,command)))
          case('mean','median')
            continue
          case('')
            call master_msg('No value of Ac provided.')
            cycle user_loop
          case default
            t1=dble_field(ifield+1,command,ierr)
            if(ierr/=0)then
              call master_msg('Problem parsing value of Ac - must be &
                 &mean/median or real number.')
              cycle user_loop
            endif
          end select
          tail_form%A_centre_def=field(ifield+1,command)
          tail_form%A_result_set=.false.
          call refresh_centre(dataset,tail_form)
          call master_msg('Set Ac to '//trim(tail_form%A_centre_def)//'.')

        case('self-consistent','sc')
          ! Whether centre is self consistent.
          if(itail/=0)then
            call master_msg('No left/right selection for variable &
               &self-consistent.')
            cycle user_loop
          endif
          tail_form%self_consistent_centre=.true.
          call refresh_centre(dataset,tail_form)
          call master_msg('Enabled self-consistent centre.')

        case('weight','weights')
          ! Which x-dependent fit weight prefactor to use.
          if(itail/=0)then
            call master_msg('No left/right selection for variable weight.')
            cycle user_loop
          endif
          select case(trim(field(3,command)))
          case('hill','none')
            tail_form%fit_weight_type=trim(field(3,command))
          case default
            call master_msg('Unknown value for variable weight.')
            cycle user_loop
          end select
          call master_msg('Set fit weights to '//&
             &trim(tail_form%fit_weight_type)//'.')

        case('mu')
          ! Leading-order exponent.
          t1=parse_dble(field(ifield+1,command),ierr)
          if(ierr/=0)then
            call master_msg('Problem parsing value of mu.')
            cycle user_loop
          endif
          if(t1<=1.d0)then
            call master_msg('mu must be greater than 1.')
            cycle user_loop
          endif
          if(itail==0)then
            tail_form%mu(:)=t1
            call master_msg('Set mu for both tails to '//&
               &trim(field(ifield+1,command))//'.')
          else
            tail_form%mu(itail)=t1
            call master_msg('Set mu for tail #'//trim(i2s(itail))//' to '//&
               &trim(field(ifield+1,command))//'.')
          endif
          call refresh_momdef(tail_form)

        case('dmu')
          ! Global delta mu.
          t1=parse_dble(field(ifield+1,command),ierr)
          if(ierr/=0)then
            call master_msg('Problem parsing value of dmu.')
            cycle user_loop
          endif
          if(t1<=0.d0)then
            call master_msg('dmu must be greater than 0.')
            cycle user_loop
          endif
          if(itail==0)then
            tail_form%global_delta_mu(:)=t1
            call master_msg('Set dmu for both tails to '//&
               &trim(field(ifield+1,command))//'.')
          else
            tail_form%global_delta_mu(itail)=t1
            call master_msg('Set dmu for tail #'//trim(i2s(itail))//' to '//&
               &trim(field(ifield+1,command))//'.')
          endif

        case('nparam')
          ! Number of parameters in fit.
          it1=parse_int(field(ifield+1,command),ierr)
          if(ierr/=0)then
            call master_msg('Problem parsing value of nparam.')
            cycle user_loop
          endif
          if(it1<1)then
            call master_msg('nparam must be equal to or greater than 1.')
            cycle user_loop
          endif
          if(itail==0)then
            tail_form%nparam(:)=it1
            call master_msg('Set nparam for both tails to '//&
               &trim(i2s(it1))//'.')
          else
            tail_form%nparam(itail)=it1
            call master_msg('Set nparam for tail #'//trim(i2s(itail))//&
               &' to '//trim(i2s(it1))//'.')
          endif
          call refresh_momdef(tail_form)

        case('constraint')
          if(itail/=0)then
            call master_msg('No left/right selection for variable constraint.')
            cycle user_loop
          endif
          ! Basic check of left-hand side.
          call parse_constraint_token(field(3,command),iparam,itail,t1)
          if(iparam<1)then
            call master_msg('Syntax error parsing left-hand side of &
               &constraint.')
            cycle user_loop
          endif
          if(are_equal(t1,0.d0))then
            call master_msg('Zero coefficient in left-hand side of &
               &constraint.')
            cycle user_loop
          endif
          if(trim(field(4,command))/='=')then
            call master_msg('Syntax error in constraint: second field is not &
               &"=".')
            cycle user_loop
          endif
          do iconstraint=1,tail_form%nconstraint
            do irhs=1,tail_form%nrhs(iconstraint)
              if(tail_form%rhs_iparam(irhs,iconstraint)==iparam.and.&
                 &tail_form%rhs_itail(irhs,iconstraint)==itail)then
                call master_msg('Coefficient on left-hand side cannot &
                   &appear on right-hand side of other constraints.')
                cycle user_loop
              endif
            enddo ! irhs
          enddo ! iconstraint
          ! Basic check of right-hand side.
          nrhs=nfield(command)-4
          do irhs=1,nrhs
            call parse_constraint_token(field(4+irhs,command),jparam,jtail,t2)
            if(jparam<1)then
              call master_msg('Syntax error parsing token #'//&
                 &trim(i2s(irhs))//' on right-hand side of constraint.')
              cycle user_loop
            endif
            if(itail==jtail.and.iparam==jparam)then
              call master_msg('Coefficient on left-hand side cannot appear on &
                 &right-hand side of constraint.')
              cycle user_loop
            endif
            do iconstraint=1,tail_form%nconstraint
              if(tail_form%lhs_iparam(iconstraint)==jparam.and.&
                 &tail_form%lhs_itail(iconstraint)==jtail)then
                call master_msg('Coefficients on right-hand side cannot &
                   &appear on left-hand side of other constraints.')
                cycle user_loop
              endif
            enddo ! iconstraint
          enddo ! irhs
          ! Store constraint as provided.
          call increment_constraint_storage(nrhs,tail_form)
          iconstraint=tail_form%nconstraint+1
          tail_form%nconstraint=iconstraint
          tail_form%lhs_iparam(iconstraint)=iparam
          tail_form%lhs_itail(iconstraint)=itail
          tail_form%nrhs(iconstraint)=nrhs
          do irhs=1,nrhs
            call parse_constraint_token(field(4+irhs,command),jparam,jtail,t2)
            tail_form%rhs_iparam(irhs,iconstraint)=jparam
            tail_form%rhs_itail(irhs,iconstraint)=jtail
            tail_form%rhs_coeff(irhs,iconstraint)=t2/t1
          enddo ! irhs
          call refresh_momdef(tail_form)
          call master_msg('Defined constraint #'//trim(i2s(iconstraint))//'.')

        case('anchor')
          ! Fit anchor.
          if(dataset%M<1)then
            call master_msg('Cannot set anchor until dataset is loaded.')
            cycle user_loop
          endif
          ! Parse "variable=value" string.
          token=field(ifield+1,command)
          ipos=scan(token,'=')
          if(ipos<1)then
            call master_msg('Problem parsing <variable>=<value> clause.')
            cycle user_loop
          endif
          t1=parse_dble(token(ipos+1:),ierr)
          if(ierr/=0)then
            call master_msg('Problem parsing value of anchor.')
            cycle user_loop
          endif
          select case(trim(token(:ipos-1)))
          case('q','M','mlogq')
            select case(trim(token(:ipos-1)))
            case('M')
              t1=(t1-0.5d0)/dble(dataset%M)
            case('mlogq')
              t1=exp(-t1)
            end select
            if(t1>1.d0-0.5d0/dble(dataset%M).or.t1<0.5d0/dble(dataset%M))then
              call master_msg('anchor q is out of range.')
              cycle user_loop
            endif
            if(itail==0)then
              tail_form%def_anchor(1:2)='q'
              tail_form%anchor(1:2)=t1
            else
              tail_form%def_anchor(itail)='q'
              tail_form%anchor(itail)=t1
            endif
          case('A','1/A','logA')
            select case(trim(token(:ipos-1)))
            case('A')
              continue
            case('1/A')
              if(t1<=0.d0)then
                call master_msg('anchor 1/A must be greater than 0.')
                cycle user_loop
              endif
              t1=1.d0/t1
            case('logA')
              t1=exp(t1)
            end select
            if(itail==0)then
              tail_form%def_anchor(1:2)='A'
              tail_form%anchor(1:2)=t1
            else
              tail_form%def_anchor(itail)='A'
              tail_form%anchor(itail)=t1
            endif
          case('x')
            if(t1<=0.d0)then
              call master_msg('anchor x must be greater than 0.')
              cycle user_loop
            endif
            if(itail==0)then
              tail_form%def_anchor(1:2)='A'
              tail_form%anchor(1:2)=1.d0/t1**tail_form%global_delta_mu(1:2)
            else
              tail_form%def_anchor(itail)='A'
              tail_form%anchor(itail)=1.d0/t1**tail_form%global_delta_mu(itail)
            endif
          end select
          if(itail==0)then
            call master_msg('Set anchor for both tails to '//&
               &trim(token)//'.')
          else
            call master_msg('Set anchor for tail #'//trim(i2s(itail))//&
               &' to '//trim(token)//'.')
          endif

        case('nbs')
          ! Number of bootstrap samples.
          if(itail/=0)then
            call master_msg('No left/right selection for variable nbs.')
            cycle user_loop
          endif
          it1=parse_int(field(ifield+1,command),ierr)
          if(ierr/=0)then
            call master_msg('Problem parsing value of nbs.')
            cycle user_loop
          endif
          if(it1<2)then
            call master_msg('nbs must be equal to or greater than 2.')
            cycle user_loop
          endif
          bs_setup%nbs=it1
          call master_msg('Set nbs to '//trim(i2s(it1))//'.')
          call refresh_bs_setup(bs_setup)
          if(bs_setup%nbs/=it1)call master_msg('nbs re-adjusted up to '//&
             &trim(i2s(bs_setup%nbs)))

        case('random')
          ! Random number seed.
          if(itail/=0)then
            call master_msg('No left/right selection for variable random.')
            cycle user_loop
          endif
          select case(trim(field(3,command)))
          case('timer')
            call init_random(0)
          case default
            it1=parse_int(field(ifield+1,command),ierr)
            if(ierr/=0)then
              call master_msg('Problem parsing value of random.  Must be &
                 &timer/<integer>.')
              cycle user_loop
            endif
            if(it1<1)then
              call master_msg('Value of random must be greater than zero.')
              cycle user_loop
            endif
            call init_random(it1)
          end select
          call master_msg('Set random seed to '//trim(field(3,command))//'.')

        case('echo')
          ! Echo mode.
          if(itail/=0)then
            call master_msg('No left/right selection for variable echo.')
            cycle user_loop
          endif
          input_echo=.true.
          call master_msg('Enabled input echo.')

        case default
          call master_msg('Unknown variable "'//&
             &trim(field(ifield,command))//'".')
          cycle user_loop

        end select

      case('unset')
        ! Unset variables.
        ! Check if there is a tail selection subcommand.
        itail=0
        ifield=2
        select case(trim(field(ifield,command)))
        case('left')
          itail=1
          ifield=ifield+1
        case('right')
          itail=2
          ifield=ifield+1
        end select

        select case(trim(field(ifield,command)))
        case('Ac','centre')
          ! Centre of distribution.
          if(itail/=0)then
            call master_msg('No left/right selection for variable Ac.')
            cycle user_loop
          endif
          tail_form%A_centre_def=tail_form_default%A_centre_def
          tail_form%A_result_set=.false.
          call refresh_centre(dataset,tail_form)
          call master_msg('Reset Ac to '//trim(tail_form%A_centre_def)//'.')

        case('self-consistent','sc')
          ! Whether centre is self consistent.
          if(itail/=0)then
            call master_msg('No left/right selection for variable &
               &self-consistent.')
            cycle user_loop
          endif
          tail_form%self_consistent_centre=.false.
          call refresh_centre(dataset,tail_form)
          call master_msg('Disabled self-consistent centre.')

        case('weight')
          ! Which fit weight to use.
          if(itail/=0)then
            call master_msg('No left/right selection for variable weight.')
            cycle user_loop
          endif
          tail_form%fit_weight_type=tail_form_default%fit_weight_type
          call master_msg('Reset fit weights to '//&
             &trim(tail_form%fit_weight_type)//'.')

        case('mu')
          ! Leading-order exponent.
          if(itail==0)then
            tail_form%mu(:)=0.d0
            call master_msg('Unset mu for both tails.')
          else
            tail_form%mu(itail)=0.d0
            call master_msg('Unset mu for tail #'//trim(i2s(itail))//'.')
          endif
          call refresh_momdef(tail_form)

        case('dmu')
          ! Global delta mu.
          if(itail==0)then
            tail_form%global_delta_mu(:)=tail_form_default%global_delta_mu(:)
            call master_msg('Reset dmu for both tails to 1.')
          else
            tail_form%global_delta_mu(itail)=&
               &tail_form_default%global_delta_mu(itail)
            call master_msg('Reset dmu for tail #'//trim(i2s(itail))//' to 1.')
          endif

        case('nparam')
          ! Number of parameters in fit.
          if(itail==0)then
            tail_form%nparam(:)=0
            call master_msg('Unset nparam for both tails.')
          else
            tail_form%nparam(itail)=0
            call master_msg('Unset nparam for tail #'//trim(i2s(itail))//'.')
          endif
          call refresh_momdef(tail_form)

        case('constraint')
          if(itail/=0)then
            call master_msg('No left/right selection for variable constraint.')
            cycle user_loop
          endif
          ! Reset constraints.
          if(tail_form%nconstraint==1)then
            call master_msg('Removed 1 active constraint.')
          else
            call master_msg('Removed '//trim(i2s(tail_form%nconstraint))//&
               &' active constraints.')
          endif
          if(tail_form%nconstraint>0)deallocate(tail_form%lhs_iparam,&
             &tail_form%lhs_itail,tail_form%nrhs,tail_form%rhs_iparam,&
             &tail_form%rhs_itail,tail_form%rhs_coeff)
          tail_form%nconstraint=0
          call refresh_momdef(tail_form)

        case('anchor')
          ! Fit anchor.
          if(itail==0)then
            tail_form%def_anchor(1:2)=''
            tail_form%anchor(1:2)=0.d0
            call master_msg('Unset anchor for both tails.')
          else
            tail_form%def_anchor(itail)=''
            tail_form%anchor(itail)=0.d0
            call master_msg('Unset anchor for tail #'//trim(i2s(itail))//'.')
          endif

        case('nbs')
          if(itail/=0)then
            call master_msg('No left/right selection for variable nbs.')
            cycle user_loop
          endif
          bs_setup%nbs=bs_setup_default%nbs
          call master_msg('Reset nbs to '//trim(i2s(bs_setup%nbs))//'.')
          call refresh_bs_setup(bs_setup)
          if(bs_setup%nbs/=bs_setup_default%nbs)call master_msg('nbs &
             &re-adjusted up to '//trim(i2s(bs_setup%nbs)))

        case('random')
          if(itail/=0)then
            call master_msg('No left/right selection for variable random.')
            cycle user_loop
          endif
          call master_msg('Reset random seed to 1.')
          call init_random(1)

        case('echo')
          ! Echo mode.
          if(itail/=0)then
            call master_msg('No left/right selection for variable echo.')
            cycle user_loop
          endif
          input_echo=.false.
          call master_msg('Disabled input echo.')

        case default
          call master_msg('Unknown variable "'//&
             &trim(field(ifield,command))//'".')
          cycle user_loop

        end select

      case('status')

        if(MPI_MASTER)then

          ! In dataset.
          write(6,'(a)')'Data:'
          write(6,'(2x,a)')trim(dataset%description)
          if(dataset%have_w)then
            write(6,'(2x,a)')'* Data are weighted'
          else
            write(6,'(2x,a)')'* Data are unweighted'
          endif
          write(6,'()')

          ! In tail_form.
          write(6,'(a)')'Tail form in fits:'
          select case(trim(tail_form%A_centre_def))
          case('mean','median')
            write(6,'(2x,a,es11.4,a)')'Ac = ',tail_form%A_centre,&
               &' ('//trim(tail_form%A_centre_def)//')'
          case default
            write(6,'(2x,a,es11.4,a)')'Ac = ',tail_form%A_centre,&
               &' (user-defined value)'
          end select
          write(6,'(2x,a)')'* Fit weight prefactor in use: '//&
             &trim(tail_form%fit_weight_type)
          do itail=1,2
            write(6,'(2x,a)')'Tail '//trim(i2s(itail))//':'
            if(tail_form%mu(itail)>0.d0)then
              write(6,'(4x,a,es11.4)')'mu = ',tail_form%mu(itail)
            else
              write(6,'(4x,a)')'mu = (unset)'
            endif
            if(tail_form%nparam(itail)>0)then
              write(6,'(4x,a)')'nparam = '//trim(i2s(tail_form%nparam(itail)))
            else
              write(6,'(4x,a)')'nparam = (unset)'
            endif
            write(6,'(4x,a)')'fit anchor:'
            select case(tail_form%def_anchor(itail))
            case('q')
              write(6,'(4x,a,es11.4)')'* q = ',tail_form%anchor(itail)
            case('A')
              write(6,'(4x,a,es11.4)')'* A = ',tail_form%anchor(itail)
            case default
              write(6,'(4x,a,es11.4)')'* (unset)'
            end select
          enddo ! itail
          if(tail_form%nconstraint<1)then
            write(6,'(2x,a)')'No constraints applied.'
          else
            write(6,'(2x,a)')trim(i2s(tail_form%nconstraint))//' constraints &
               &applied:'
            do iconstraint=1,tail_form%nconstraint
              write(6,'(2x,a)',advance='no')'* '
              itail=tail_form%lhs_itail(iconstraint)
              iparam=tail_form%lhs_iparam(iconstraint)
              nrhs=tail_form%nrhs(iconstraint)
              select case(itail)
              case(1)
                write(6,'(a)',advance='no')'l'
              case(2)
                write(6,'(a)',advance='no')'r'
              end select
              write(6,'(a)',advance='no')trim(i2s(iparam))//' ='
              do irhs=1,nrhs
                jtail=tail_form%rhs_itail(irhs,iconstraint)
                jparam=tail_form%rhs_iparam(irhs,iconstraint)
                t1=tail_form%rhs_coeff(irhs,iconstraint)
                write(6,'(1x)',advance='no')
                if(.not.are_equal(t1,1.d0))&
                   &write(6,'(es11.4,a)',advance='no')t1,'*'
                select case(jtail)
                case(1)
                  write(6,'(a)',advance='no')'l'
                case(2)
                  write(6,'(a)',advance='no')'r'
                end select
                write(6,'(a)',advance='no')trim(i2s(jparam))
              enddo ! irhs
              write(6,'()')
            enddo ! iconstraint
          endif ! nconstraint<1 or not
          write(6,'()')

        endif ! MPI_MASTER

      case('help')
        if(MPI_MASTER)then
          select case(trim(field(2,command)))
          case('')
            call pprint('TREAT is a toolbox for analyzing statistical samples &
               &drawn from distributions affected by Pareto-type heavy tails, &
               &P(A) ~ A^-mu at A->infty.  TREAT is particularly useful in &
               &that it can compute a reliable estimate of the mean for 2 < &
               &mu <= 3 (or even 1 < mu <= 2 for distributions satisfying &
               &certain symmetry conditions) and of the variance for 3 < mu &
               &<= 4.')
            call pprint('')
            call pprint('TREAT uses a command-line interface.  The list &
               &of available commands is:')
            call pprint('')
            call pprint('* inspect <file>',0,10)
            call pprint('* load <file> [type <type> using <columns>] &
               &[trim <first-line> <last-line>]',0,7)
            call pprint('* unload',0,9)
            call pprint('* generate <nsample> &
               &[<c1>*]{gauss(<lambda1>)|h(<mu1>,<lambda1>)} &
               &[[<c2>*]{gauss(<lambda2>)|h(<mu2>,<lambda2>)} [...]]',0,11)
            call pprint('* save <file>',0,7)
            call pprint('* report {stats|corr}',0,9)
            call pprint('* set [left|right] <variable> <value>',0,6)
            call pprint('* unset [left|right] <variable>',0,8)
            call pprint('* assess <variable>',0,9)
            call pprint('* plot {left|right} <subject> [at <grid-definition>] &
               &[to <file>]',0,7)
            call pprint('* evaluate TRE [plot {P|fit} [grid <ngrid>] &
               &[to <file>]]',0,11)
            call pprint('* status',0,9)
            call pprint('* help [<command> | set <variable>]',0,7)
            call pprint('')
            call pprint('Type help <command> for detailed information about &
               &<command>.  Also, refer to "help notation" for a summary of &
               &the notation used in the script, and to "help workflow" for &
               &suggestions of how to use this utility.',2,2)
            call pprint('')
          case('notation')
            call pprint('Help topic: notation',0,0)
            call pprint('')
            call pprint('The following symbols are used throughout TREAT:',2,2)
            call pprint('* M : number of samples in the dataset',2,4)
            call pprint('* A : variable sampled in the dataset',2,4)
            call pprint('* P(A) : probability distribution function followed &
               &by variable A',2,4)
            call pprint('* w : weight of each sample in the dataset',2,4)
            call pprint('* Ac : parameter which approximates the "centre" &
               &of the distribution.',2,4)
            call pprint('* mu : principal exponent of the tail of the &
               &probability distribution, P(A) ~ |A-Ac|^-mu for |A|->infty',2,4)
            call pprint('* dmu : exponent increment for beyond-leading &
               &order terms of the tail of the probability distribution, &
               &P(A) = c0 |A-Ac|^-mu + c1 + |A-Ac|^{-mu-dmu} + ... for &
               &|A|->infty',2,4)
            call pprint('* q : sample quantile of a statistic, e.g., the &
               &mth largest value of A in the sample, A_m, corresponds to &
               &q_m = (m-1/2)/M',2,4)
            call pprint('* mlogq : shorthand for -log(q)',2,4)
            call pprint('* logA : shorthand for log|A-Ac|',2,4)
            call pprint('* "qq" scale : (mlogq,logA) scale in which tail &
               &index estimation is performed, in which the asymptotic &
               &slope of logA as a function of mlogq is 1/(mu-1).',2,4)
            call pprint('* 1/A : shorthand for 1/|A-Ac|',2,4)
            call pprint('* x : shorthand for |A-Ac|^-dmu',2,4)
            call pprint('* y : shorthand for |A-Ac|^(mu-1)',2,4)
            call pprint('* "yx" scale : (x,y) scale in which the tail &
               &regression takes place, in which y is a polynomial in x &
               &near x=0.',2,4)
            call pprint('* "anchor" : point from which the fit to the &
               &distribution replaces the sampled data in the tail regression &
               &estimator.  The value of A at which this happens is referred &
               &to as the "anchor A", or AL (for the left tail) or AR (for &
               &the right tail).  The corresponding value of q is referred to &
               &as the "anchor q", and that of -log(q) as the "anchor &
               &mlogq".',2,4)
            call pprint('')
          case('workflow')
            call pprint('Help topic: workflow',0,0)
            call pprint('')
            call pprint('The following sequence of actions would be a &
               &reasonable way of using TREAT:',2,2)
            call pprint('* Load data from single-column file:',2,4)
            call pprint('TREAT> load file',4,6)
            call pprint('* Verify standard estimates of expectation value and &
               &and absence of serial correlation:',2,4)
            call pprint('TREAT> report stats',4,6)
            call pprint('TREAT> report corr',4,6)
            call pprint('* Verify expected principal exponents:',2,4)
            call pprint('TREAT> assess mu at mlogq=1:8:17',4,6)
            call pprint('* Set principal exponent for both tails:',2,4)
            call pprint('TREAT> set mu 4',4,6)
            call pprint('* Auto-set nparam and anchor (symmetrically &
               &sampled):',2,4)
            call pprint('TREAT> assess anchor,nparam at mlogq=1:8:17 apply',4,6)
            call pprint('* Plot e.g. left tail and tail fit and compare &
               &them:',2,4)
            call pprint('TREAT> plot P at left mlogq=1:12:100',4,6)
            call pprint('TREAT> plot fit at left mlogq=4:12:100',4,6)
            call pprint('* Evaluate the final estimate of the mean:',2,4)
            call pprint('TREAT> evaluate TRE',4,6)
            call pprint('')
          case('inspect')
            call pprint('Command: inspect <file>',0,9)
            call pprint('')
            call pprint('Report the number of data lines and columns &
              &detected in <file>.',2,2)
            call pprint('')
          case('load')
            call pprint('Command: load <file> [type <type> using <columns>] &
               &[trim <first-line> <last-line>]',0,9)
            call pprint('')
            call pprint('Load data from <file> into the working dataset, &
               &replacing any previously loaded data.  By default, data are &
               &assumed to be unweighted ("A" type) and all data on column 1 &
               &are loaded into a dataset.',2,2)
            call pprint('')
            call pprint('Unweighted data can be loaded from a different &
               &column with "type A using <A-column>".  Weighted data can be &
               &loaded with "type Aw using <A-column> <w-column>".  A subset &
               &of the data in <file> can be loaded using the "trim" &
               &subcommand.',2,2)
            call pprint('')
          case('unload')
            call pprint('Command: unload',0,9)
            call pprint('')
            call pprint('Unload the dataset from memory.',2,2)
            call pprint('')
          case('generate')
            call pprint('Command: generate <nsample> &
               &[<c1>*]{gauss([<lambda1>[,<centre1>]) | &
               &        h(<mu1>[,<lambda1>[,<centre1>]])} &
               &[<c2>*]{gauss([<lambda2>[,<centre2>]) | &
               &        h(<mu2>[,<lambda2>[,<centre2>]])} [...]]',0,9)
            call pprint('')
            call pprint('Generate dataset of <nsample> unweighted data &
               &distributed according to the specified linear combination of &
               &model probability distribution functions.',&
               &2,2)
            call pprint('Function "gauss" is a normalized Gaussian centred at &
               &<centre> (zero by default) of width <lambda> (one by &
               &default),',2,2)
            call pprint('')
            write(6,'(a)')'                  /   / A-centre \ 2 \'
            write(6,'(a)')'  gauss(A) = N exp| - | -------- |   | ,'
            write(6,'(a)')'                  \   \  lambda  /   /'
            call pprint('')
            call pprint('where N is the normalization factor, and function &
               &"h" is a normalized model heavy tail distribution,',2,2)
            call pprint('')
            write(6,'(a)')'                     1'
            write(6,'(a)')'  h(A) = N --------------------- ,'
            write(6,'(a)')'                | A-centre | mu'
            write(6,'(a)')'            1 + | -------- |'
            write(6,'(a)')'                |  lambda  |'
            call pprint('')
            call pprint('where N is the normalization factor.',2,2)
            call pprint('')
            call pprint('The analytic asymptote of the model distribution is &
               &reported during the generation process.',2,2)
            call pprint('')
            call pprint('The model distribution is generated by direct &
               &inversion of the CDF, and the resulting data is serially &
               &uncorrelated.',2,2)
            call pprint('')
          case('save')
            call pprint('Command: save <file>',0,9)
            call pprint('')
            call pprint('Saves the current dataset to <file>.  Unweighted &
               &datasets are written as a single-column file, and weighted &
               &datasets are written as two-column files, with A in the &
               &first column and weights in the second.',2,2)
            call pprint('')
            call pprint('The "save" command is useful to store model datasets &
               &created with the "generate" command (which are always &
               &unweighted).',2,2)
            call pprint('')
          case('report')
            call pprint('Command: report <report>',0,9)
            call pprint('')
            call pprint('Analyze dataset and report properties.  Available &
               &values of <report> are:',2,2)
            call pprint('* stats : report mean, standard error, variance, &
               &variance standard error, median, and max/min of the data.',2,4)
            call pprint('* corr  : report correlation length of whole dataset &
               &and each third, each fifth, and each tenth of the dataset. &
               &Also report these for datasets formed by picking one out of &
               &every 2, 4, 8, and 16 successive points of the original &
               &sample.  Note that the correlation time should be one for the &
               &TRE to be applicable.',2,4)
            call pprint('')
          case('unset')
            call pprint('Command: unset <variable>',0,9)
            call pprint('')
            call pprint('Sets <variable> to its default value.',2,2)
            call pprint('')
            call pprint('Type "help set <variable>" for detailed information &
               &on variables and their default values.',2,2)
            call pprint('')
          case('assess')
            call pprint('Command: assess <variable> <options> [to &
               &<filename>] [apply]',0,9)
            call pprint('')
            call pprint('Analyze dataset and suggest a value for <variable>, &
               &setting the suggested value if "apply" is provided, and &
               &optionally writing a plot of the data analyzed to <filename>. &
               &The following values of <variable> <options> can be used:',2,2)
            call pprint('* mu at <point-range> : perform a linear fit in qq &
               &scale at each point in <point-range> and report the value of &
               &mu that minimizes the least-squares function.',2,4)
            call pprint('* nparam,anchor at <anchor-range> for <nparam-range> &
               &: for each anchor in <anchor-range>, report the "optimal" &
               &value of nparam in <nparam-range> along with the resulting &
               &estimators of 0th, 1st and 2nd moment of the distribution. &
               &The smallest uncertainty determines the optimal anchor and &
               &nparam.',2,4)
            call pprint('')
          case('plot')
            call pprint('Command: plot <subject> at <point-range> [to &
               &<file>]',0,9)
            call pprint('')
            call pprint('Plot <subject> at the provided points to <file>.  &
               &The default value of <file> is "TREAT_<subject>_plot.dat".',&
               &2,2)
            call pprint('')
            call pprint('The plotting grid is defined using the "at &
               &<point-range>" option, where <point-range> can be one of:',2,2)
            call pprint('')
            call pprint('* <variable>=<comma-separated-list>',2,4)
            call pprint('* <variable>=<first>:<last>:<number>',2,4)
            call pprint('')
            call pprint('and <variable> is one of:',2,2)
            call pprint('')
            call pprint('* q',2,4)
            call pprint('* A',2,4)
            call pprint('* mlogq',2,4)
            call pprint('* 1/A',2,4)
            call pprint('* logA',2,4)
            call pprint('* x',2,4)
            call pprint('')
            call pprint('The following values of <subject> can be used:',2,2)
            call pprint('')
            call pprint('* model : plot the CDF of the current model &
               &probability.  This requires having run a "generate" command &
               &and having set mu.  The plot file contains two columns:',2,4)
            call pprint('1. -log(CDF)',4,7)
            call pprint('2. log|A-Ac|',4,7)
            call pprint('')
          case('evaluate')
            call pprint('Command: evaluate TRE [plot {P|fit} [grid <ngrid>] &
               &[to <file>]]',0,9)
            call pprint('')
            call pprint('Evaluate the tail-regression estimator of the mean &
               &of the dataset.',2,2)
            call pprint('')
            call pprint('Plots of the probability distribution and the fit &
               &can be generated using the same bootstrap sample if the &
               &respective "plot" subcommands are specified.  The number of &
               &points used to sample the target function in each tail can be &
               &set by specifying "grid" (100 by default), and the name of &
               &the plot file can be set with the "to" subcommand &
               &(TREAT_<function>_plot.dat by default).',2,2)
            call pprint('')
            call pprint('The plot files contain the following columns:')
            call pprint('* P :',2,4)
            call pprint('1. -log(q)',4,7)
            call pprint('2-5. -2,-1,+1,+2-sigma confidence interval limits &
               &for -log(q)',4,7)
            call pprint('6. log|A-Ac|',4,7)
            call pprint('7. 1/|A-Ac|',4,7)
            call pprint('8. q*|A-Ac|^{mu-1}',4,7)
            call pprint('9-12. -2,-1,+1,+2-sigma confidence interval limits &
               &for q*|A-Ac|^{mu-1}',4,7)
            call pprint('13. A',4,7)
            call pprint('14. P(A)',4,7)
            call pprint('15-18. -2,-1,+1,+2-sigma confidence interval limits &
               &for P(A)',4,7)
            call pprint('')
            call pprint('* fit :',2,4)
            call pprint('1. -log(q)',4,7)
            call pprint('2. log|A-Ac|',4,7)
            call pprint('3-6. -2,-1,+1,+2-sigma confidence interval limits &
               &for -log|A-Ac|',4,7)
            call pprint('7. 1/|A-Ac|',4,7)
            call pprint('8. q*|A-Ac|^{mu-1}',4,7)
            call pprint('9-12. -2,-1,+1,+2-sigma confidence interval limits &
               &for q*|A-Ac|^{mu-1}',4,7)
            call pprint('13. A',4,7)
            call pprint('14. P(A)',4,7)
            call pprint('15-18. -2,-1,+1,+2-sigma confidence interval limits &
               &for P(A)',4,7)
            call pprint('')
          case('status')
            call pprint('Command: status',0,9)
            call pprint('')
            call pprint('',2,2)
            call pprint('Report dataset details and values of internal &
               &variables.',2,2)
          case('set')
            if(nfield(command)==2)then
              call pprint('Command: set [left|right] <variable> <value>',0,9)
              call pprint('')
              call pprint('Sets <variable> to <value> (possibly only for the &
                 &left or right tail for tail-specific variables).  The list &
                 &of available variables is:',2,2)
              call pprint('')
              call pprint('- Ac',2,4)
              call pprint('- self-consistent',2,4)
              call pprint('- weight',2,4)
              call pprint('* mu',2,4)
              call pprint('* dmu',2,4)
              call pprint('* nparam',2,4)
              call pprint('* anchor',2,4)
              call pprint('- constraint',2,4)
              call pprint('- nbs',2,4)
              call pprint('- random',2,4)
              call pprint('- echo',2,4)
              call pprint('')
              call pprint('In the list above, asterisks * flag tail-specific &
                 &variables.  Omitting [left|right] for a tail-specific &
                 &variable will set the value of the variable for *both* &
                 &tails to <value>.',2,2)
              call pprint('')
              call pprint('Type "help set <variable>" for detailed &
                 &information on <variable>.',2,2)
              call pprint('')
            else
              select case(trim(field(3,command)))
              case('Ac')
                call pprint('Variable: Ac',0,10)
                call pprint('')
                call pprint('"Ac" defines the centre Ac of the distribution, &
                   &and can be set to "median" (default), "mean", or to a &
                   &real-valued number.  If "self-consistent" is set, setting &
                   &"Ac" determines the initial value of Ac, which is then &
                   &set to the mean when "evaluate mean" is run.  Ac is used &
                   &to divide the distribution in two, and in the fit formula &
                   &throughout the fitting process.',2,2)
                call pprint('')
              case('self-consistent','sc')
                call pprint('Variable: self-consistent',0,10)
                call pprint('')
                call pprint('"self-consistent" (or "sc") defines whether &
                   &variable Ac is to be self-consistently set to the &
                   &estimated mean of the data (false by default).  In this &
                   &case, variable "Ac" defines the initial value of Ac.',2,2)
                call pprint('')
              case('weight','weights')
                call pprint('Variable: weight',0,10)
                call pprint('')
                call pprint('"weight" sets the type of x-dependent fit weight &
                   &prefactor to use in the tail-fitting.  Allowed values &
                   &are:',2,2)
                call pprint('* "hill": use the Hill estimator weight, &
                   &{log[(M+1/2)/(k-1/2)]}^-1, as prefactor',2,4)
                call pprint('* "none": apply no prefactor.',2,4)
                call pprint('')
              case('mu')
                call pprint('Variable: mu',0,10)
                call pprint('')
                call pprint('"mu" is the (expected) real-valued principal &
                   &exponent of the left/right tail of the probability &
                   &distribution, P(A) ~ |A-Ac|^-mu as A->infty.  mu sets &
                   &the leading-order power of |A-Ac| in the fit form.  "mu" &
                   &is unset by default.',2,2)
                call pprint('')
              case('dmu')
                call pprint('Variable: dmu',0,10)
                call pprint('')
                call pprint('"dmu" is the increment to use between &
                   &consecutive exponents in the left/right tail fit.  "dmu" &
                   &is 1 by default.',2,2)
                call pprint('')
              case('nparam')
                call pprint('Variable: nparam',0,10)
                call pprint('')
                call pprint('"nparam" is the integer expansion order &
                   &of the fit form for each tail.  The expansion goes from &
                   &|A-Ac|^-mu to |A-Ac|^(-mu+(nparam-1)*dmu).  "nparam" is &
                   &unset by default.',2,2)
                call pprint('')
              case('anchor')
                call pprint('Variable: anchor',0,10)
                call pprint('')
                call pprint('"anchor" defines the point at which the fit &
                   &replaces the distribution for each tail.  "anchor" can &
                   &be specified in terms of the number of points in the tail &
                   &(M), the minus logarithm of the quantile function &
                   &(mlogq), the value of A (A), the value of 1/|A-Ac| &
                   &(1/A), or the value of log|A-Ac| (logA) at the desired &
                   &point, e.g., "set anchor mlogq=3.5".  "anchor" is unset &
                   &by default.',2,2)
                call pprint('')
              case('constraint')
                call pprint('Variable: constraint',0,10)
                call pprint('')
                call pprint('Setting "constraint" adds a linear, homogeneous &
                   &constraint on the parameters.  The constraint should be &
                   &specified as {l|r}<i> = [<coeff1>*]{l|r}<j1> &
                   &[[<coeff2>*{l|r}<j2> [...]] , where, e.g., l3 refers to &
                   &the third parameter in the left tail.  Unsetting &
                   &"constraint" clears *all* constraints.')
                call pprint('')
              case('nbs')
                call pprint('Variable: nbs',0,10)
                call pprint('')
                call pprint('"nbs" is the number of samples to use for &
                   &bootstrapping in "plot", "evaluate", and "assess" &
                   &commands.  By default "nbs" is '//&
                   &trim(i2s(bs_setup_default%nbs))//'.  Values &
                   &of "nbs" are rounded up to the nearest multiple of &
                   &the available number of MPI processes.  The uncertainty &
                   &in a bootstrapped uncertainty dx is approximately &
                   &dx / sqrt(2*nbs); e.g., nbs=5000 would yield &
                   &uncertainties which are accurate to 1% of their values.')
                call pprint('')
              case('random')
                call pprint('Variable: random',0,10)
                call pprint('')
                call pprint('"random" sets the random number seed.  The value &
                   &of "random" can be a positive integer or "timer".  The &
                   &default value of "random" is 1.  Operations performed &
                   &after setting "timer" to a value will use the same random &
                   &number sequence.  This is useful for reproducibility.')
                call pprint('')
              case('echo')
                call pprint('Variable: echo',0,10)
                call pprint('')
                call pprint('Setting "echo" causes input commands to be &
                   &echoed back to stdout, which is useful for scripts.')
                call pprint('')
              case default
                call master_msg('No help for variable "'//&
                   &trim(field(3,command))//'".')
              end select
            endif
          case default
            call master_msg('No help for command "'//trim(field(2,command))//&
               &'".')
          end select
        endif ! MPI_MASTER

      case('quit','exit')
        call quit()

      case('')
        continue

      case default
        call master_msg('Command "'//trim(field(1,command))//&
           &'" not recognized.')
      end select

    enddo user_loop

  END SUBROUTINE treat_main


  ! DATA READ/WRITE ROUTINES.


  SUBROUTINE parse_file(filename,get_nline,nline,ncolumn)
    !----------------------------------------------------------!
    ! Get number of data lines nline and number of columns per !
    ! line ncolumn in file filename.  If .not.GET_NLINE, just  !
    ! get number of columns in fist line of file and leave     !
    ! NLINE set to zero.                                       !
    !----------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: filename
    LOGICAL, INTENT(in) :: get_nline
    INTEGER, INTENT(inout) :: nline,ncolumn
    INTEGER, PARAMETER :: io=10
    CHARACTER(2048) line
    INTEGER j,ierr
    DOUBLE PRECISION t1

    ! Initialize.
    nline=-1 ! flag non-existing file
    ncolumn=0

    ! Open file.
    open(unit=io,file=trim(filename),status='old',iostat=ierr)
    if(ierr/=0)return
    nline=0

    ! Get first data line to count the columns, then count non-empty,
    ! non-comment lines.
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr<0)exit
      if(ierr>0)then
        nline=-2 ! flag reading error
        exit
      endif
      line=adjustl(line)
      ! Skip comments.
      if(line(1:1)=='#')cycle
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      if(ncolumn==0)then
        do
          read(line,*,iostat=ierr)(t1,j=1,ncolumn+1)
          if(ierr/=0)exit
          ncolumn=ncolumn+1
        enddo
        if(ncolumn==0)then
          nline=-3 ! flag column count problem
          exit
        endif
      endif
      if(.not.get_nline)exit
      nline=nline+1
    enddo

    ! Close file.
    close(io)

  END SUBROUTINE parse_file


  SUBROUTINE read_file(filename,icol_A,icol_w,Mstart,Mend,dataset,ierr)
    !------------------------------------------------------------!
    ! Read dataset from lines Mstart:Mend of file filename using !
    ! column icol_A for A and icol_w for w.                      !
    !------------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: Mstart,Mend,icol_A(:),icol_w
    TYPE(dataset_type),INTENT(inout) :: dataset
    INTEGER, INTENT(inout) :: ierr
    INTEGER, PARAMETER :: io=10
    CHARACTER(2048) line
    INTEGER i,j,k
    DOUBLE PRECISION t1(size(icol_A)),t2

    ! Initialize dataset.
    dataset%have_w=icol_w/=0
    dataset%M=Mend-Mstart+1
    dataset%description=trim(i2s(dataset%M))//' '//&
       &trim(type_string(dataset%have_w))//'-type data loaded from "'//&
       &trim(filename)//'", lines '//trim(i2s(Mstart))//':'//trim(i2s(Mend))
    if(dataset%have_w)then
      dataset%description=trim(dataset%description)//', columns '//&
         &trim(i2s(icol_A(1)))
    else
      dataset%description=trim(dataset%description)//', column '//&
         &trim(i2s(icol_A(1)))
    endif
    do i=2,size(icol_A)
      dataset%description=trim(dataset%description)//'+'//&
         &trim(i2s(icol_A(i)))
    enddo ! i
    if(dataset%have_w)dataset%description=trim(dataset%description)//' and '//&
       &trim(i2s(icol_w))
    if(allocated(dataset%A))deallocate(dataset%A)
    if(allocated(dataset%w))deallocate(dataset%w)
    if(allocated(dataset%indx))deallocate(dataset%indx)
    if(allocated(dataset%q))deallocate(dataset%q)
    allocate(dataset%A(dataset%M))
    if(dataset%have_w)allocate(dataset%w(dataset%M))

    ! Open file.
    open(unit=io,file=trim(filename),status='old')

    ! Read file.
    i=0
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr/=0)exit
      line=adjustl(line)
      if(line(1:1)=='#')cycle
      do k=1,size(icol_A)
        read(line,*,iostat=ierr)(t1(k),j=1,icol_A(k))
        if(ierr/=0)exit
      enddo ! k
      if(ierr/=0)exit
      if(icol_w>0)read(line,*,iostat=ierr)(t2,j=1,icol_w)
      if(ierr/=0)exit
      i=i+1
      if(i>=Mstart.and.i<=Mend)then
        dataset%A(i-Mstart+1)=sum(t1)
        if(icol_w>0)dataset%w(i-Mstart+1)=t2
      endif
      if(i==Mend)exit
    enddo

    ! Flag error.
    ierr=0
    if(i/=Mend)ierr=1
    if(ierr/=0)then
      ! Reset dataset.
      dataset%have_w=.false.
      dataset%M=0
      dataset%description=''
      if(allocated(dataset%A))deallocate(dataset%A)
      if(allocated(dataset%w))deallocate(dataset%w)
      if(allocated(dataset%indx))deallocate(dataset%indx)
      if(allocated(dataset%q))deallocate(dataset%q)
    endif

    ! Close file.
    close(io)

  END SUBROUTINE read_file


  SUBROUTINE dump_dataset(dataset,filename)
    !---------------------------!
    ! Dump dataset to filename. !
    !---------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    CHARACTER(*),INTENT(in) :: filename
    INTEGER i
    INTEGER, PARAMETER :: io=10
    open(unit=io,file=filename,status='replace')
    write(io,'(a)')'# '//trim(adjustl(dataset%description))
    do i=1,dataset%M
      if(dataset%have_w)then
        write(io,'(2(1x,es20.12))')dataset%A(i),dataset%w(i)
      else
        write(io,'(1(1x,es20.12))')dataset%A(i)
      endif
    enddo ! i
    close(io)
  END SUBROUTINE dump_dataset


  ! REPORT ROUTINES


  SUBROUTINE report_serial_correlation(dataset)
    !----------------------------------------------------!
    ! Perform serial-correlation analysis on portions of !
    ! {A,w}(1:M) with ascending-sort vector indx(1:M).   !
    !----------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    ! Parameters.
    INTEGER, PARAMETER :: NNDECORR=5,NNCHUNK=4
    INTEGER, PARAMETER :: NCHUNK_VECTOR(NNCHUNK)=(/1,3,5,10/)
    ! Local variables.
    INTEGER M_work,inchunk,nchunk,ichunk,k,ndecorr,indecorr,kk
    DOUBLE PRECISION t1
    DOUBLE PRECISION err_factor_array(NNDECORR,maxval(NCHUNK_VECTOR),NNCHUNK)
    LOGICAL mask_work(dataset%M)
    DOUBLE PRECISION, ALLOCATABLE :: A_work(:),w_work(:)

    ! Output header.
    write(6,'(a)')'Serial correlation analysis'
    write(6,'(a)')'==========================='

    ! Loop over decorrelation skip sizes.
    err_factor_array=-1.d0
    do indecorr=1,NNDECORR
      ndecorr=2**(indecorr-1)
      ! Loop over chunk sizes.
      do inchunk=1,NNCHUNK
        nchunk=NCHUNK_VECTOR(inchunk)
        M_work=dataset%M/(ndecorr*nchunk)
        if(M_work<10)cycle
        allocate(A_work(M_work))
        if(dataset%have_w)allocate(w_work(M_work))
        ! Loop over chunks.
        do ichunk=1,nchunk
          ! Construct contiguous version of data vector for this chunk.
          mask_work=.false.
          do k=1,M_work
            kk=(ichunk-1)*M_work*ndecorr+(k-1)*ndecorr+1
            mask_work(dataset%indx(kk))=.true.
          enddo ! k
          ! Analyze serial correlation for this chunk.
          A_work(1:M_work)=pack(dataset%A,mask_work)
          if(dataset%have_w)then
            w_work(1:M_work)=pack(dataset%w,mask_work)
            call reblock_w(M_work,A_work,w_work,err_factor=t1)
          else
            call reblock(M_work,A_work,err_factor=t1)
          endif
          err_factor_array(indecorr,ichunk,inchunk)=t1
        enddo ! ichunk
        deallocate(A_work)
        if(dataset%have_w)deallocate(w_work)
      enddo ! inchunk
    enddo ! iskip

    ! Report.
    do inchunk=1,NNCHUNK
      nchunk=NCHUNK_VECTOR(inchunk)
      write(6,'(a)')'* Error factor for '//trim(i2s(nchunk))//' chunks:'
      write(6,'(a9,'//trim(i2s(NNDECORR))//'(1x,i12))')'ndecorr',&
         &(/ (2**(indecorr-1), indecorr=1,NNDECORR) /)
      do ichunk=1,nchunk
        write(6,'(2x,a2,"-",a3,"%",'//trim(i2s(NNDECORR))//&
           &'(1x,f12.4))')&
           &trim(i2s(nint(dble(ichunk-1)*100.d0/dble(nchunk)))),&
           &trim(i2s(nint(dble(ichunk)*100.d0/dble(nchunk)))),&
           &err_factor_array(:,ichunk,inchunk)
      enddo ! idecorr
      write(6,'()')
    enddo ! inchunk

  END SUBROUTINE report_serial_correlation


  ! ASSESS ROUTINES


  SUBROUTINE assess_mu(dataset,tail_form,bs_setup,plot_grid,filename,&
     &apply_assessment)
    !--------------------------------------------------------!
    ! Perform tail-index estimation on both tails of A(1:M). !
    !--------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    TYPE(tail_form_type),INTENT(inout) :: tail_form
    TYPE(bs_setup_type),INTENT(in) :: bs_setup
    TYPE(plot_grid_type),INTENT(in) :: plot_grid
    CHARACTER(*),INTENT(in) :: filename
    LOGICAL, INTENT(in) :: apply_assessment
    ! Arguments to bootstrap routine.
    TYPE(assess_grid_type) assess_grid
    LOGICAL have_mom1
    DOUBLE PRECISION bs_A(bs_setup%nbs,plot_grid%npoint),&
       &bs_M(bs_setup%nbs,plot_grid%npoint),&
       &bs_TIE_param(bs_setup%nbs,2,plot_grid%npoint),&
       &bs_TIE_chi2(bs_setup%nbs,plot_grid%npoint)
    ! Variables for analysis.
    DOUBLE PRECISION x_qq,dx_qq,x_qq_CI(4),y_qq,dy_qq,y_qq_CI(4),&
       &mu,dmu,mu_CI(4),c,dc,c_CI(4),chi2,chi2_CI(4),t1,t2
    DOUBLE PRECISION bs_x_qq(bs_setup%nbs),bs_y_qq(bs_setup%nbs)
    ! Local variables.
    INTEGER ipoint,itail,imin,ierr
    DOUBLE PRECISION min_chi2,rounded_mu(2)
    ! Parameters.
    INTEGER, PARAMETER :: io=10

    ! Initialize.
    call initialize_assess_grid(assess_grid)

    ! Write header.
    if(MPI_MASTER)then
      write(6,'(a)')'Tail-index estimation'
      write(6,'(a)')'====================='
    endif ! MPI_MASTER

    ! Perform evaluation.
    call bootstrap_fit(dataset,tail_form,bs_setup,plot_grid,assess_grid,&
       &have_mom1,bs_A=bs_A,bs_M=bs_M,bs_TIE_param=bs_TIE_param,&
       &bs_TIE_chi2=bs_TIE_chi2)

    ! Perform assessment.
    rounded_mu=0.d0
    do itail=1,2
      imin=0
      min_chi2=0.d0
      do ipoint=1,plot_grid%npoint
        if(get_point_itail(dataset,tail_form,plot_grid,ipoint)/=itail)cycle
        call characterize_dist(bs_setup%nbs,bs_TIE_chi2(:,ipoint),mean=t1,&
           &var=t2)
        chi2=t1+3.d0*sqrt(t2)
        if(imin==0.or.chi2<min_chi2)then
          imin=ipoint
          min_chi2=chi2
        endif
      enddo ! ipoint
      ! Report.
      ipoint=imin
      if(ipoint==0)cycle
      if(itail==1)then
        bs_x_qq=-log((bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
      else
        bs_x_qq=-log(1.d0-(bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
      endif
      bs_y_qq=log(abs(bs_A(:,ipoint)-tail_form%A_centre))
      call characterize_dist(bs_setup%nbs,bs_x_qq,mean=t1,var=t2)
      x_qq=t1
      dx_qq=sqrt(t2)
      call characterize_dist(bs_setup%nbs,bs_y_qq,mean=t1,var=t2)
      y_qq=t1
      dy_qq=sqrt(t2)
      call characterize_dist(bs_setup%nbs,bs_TIE_param(:,1,ipoint),mean=t1,&
         &var=t2)
      mu=t1
      dmu=sqrt(t2)
      call characterize_dist(bs_setup%nbs,bs_TIE_param(:,2,ipoint),mean=t1,&
         &var=t2)
      c=t1
      dc=sqrt(t2)
      ! Round to nearest half-integer, and take large mu to mean the tail is
      ! not heavy.
      if(mu<8.d0)rounded_mu(itail)=0.5d0*nint(mu*2.d0)
      if(MPI_MASTER)then
        if(itail==1)then
          write(6,'(a)')'* Left tail:'
        else
          write(6,'(a)')'* Right tail:'
        endif
        write(6,'(2x,a,2(1x,es20.12)," (=~ ",f4.2,")")')'* mu    = ',mu,dmu,&
           &rounded_mu(itail)
        write(6,'(2x,a,2(1x,es20.12))')'* c     = ',c,dc
        write(6,'(2x,a,2(1x,es20.12))')'* mlogq = ',x_qq,dx_qq
        write(6,'(2x,a,2(1x,es20.12))')'* logA  = ',y_qq,dy_qq
      endif
    enddo ! itail
    if(MPI_MASTER)write(6,'()')

    ! Write plot if requested.
    if(MPI_MASTER.and.len_trim(filename)>0)then

      ! Open file.
      open(unit=io,file=trim(filename),status='replace',iostat=ierr)
      if(ierr/=0)then
        call master_msg('Could not open plot file.')
        return
      endif
      write(io,'(a)')'# x_qq|CIs  y_qq|CIs  mu|CIs  c|CIs  MSE|CIs'

      ! Loop over plot grid points.
      do ipoint=1,plot_grid%npoint
        ! Transform scales.
        itail=get_point_itail(dataset,tail_form,plot_grid,ipoint)
        if(itail==1)then
          bs_x_qq=-log((bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
        else
          bs_x_qq=-log(1.d0-(bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
        endif
        bs_y_qq=log(abs(bs_A(:,ipoint)-tail_form%A_centre))
        ! Evaluate statistics.
        call get_conf_intervals(bs_setup%nbs,bs_x_qq,x_qq,x_qq_CI)
        call get_conf_intervals(bs_setup%nbs,bs_y_qq,y_qq,y_qq_CI)
        call get_conf_intervals(bs_setup%nbs,bs_TIE_param(:,1,ipoint),mu,mu_CI)
        call get_conf_intervals(bs_setup%nbs,bs_TIE_param(:,2,ipoint),c,c_CI)
        call get_conf_intervals(bs_setup%nbs,bs_TIE_chi2(:,ipoint),chi2,chi2_CI)
        ! Write fit results.
        write(io,'(25(1x,es20.12))')x_qq,x_qq_CI,y_qq,y_qq_CI,mu,mu_CI,c,c_CI,&
           &chi2,chi2_CI
      enddo ! ipoint

      ! Close file.
      close(io)

      ! Report.
      call master_msg('Plot written to "'//trim(filename)//'".')

    endif ! MPI_MASTER

    ! Apply suggested values of mu.
    if(apply_assessment)then
      call master_msg('Setting mu to suggested values.')
      tail_form%mu=rounded_mu
    endif

  END SUBROUTINE assess_mu


  SUBROUTINE assess_fit(dataset,tail_form,bs_setup,assess_grid,filename,&
     &apply_assessment)
    !------------------------------------------------------!
    ! Make plot of integrals against anchor and/or nparam. !
    !------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    TYPE(tail_form_type),INTENT(inout) :: tail_form
    TYPE(bs_setup_type),INTENT(in) :: bs_setup
    TYPE(assess_grid_type),INTENT(in) :: assess_grid
    CHARACTER(*), INTENT(in) :: filename
    LOGICAL, INTENT(in) :: apply_assessment
    ! Arguments to bootstrap_fit.
    TYPE(plot_grid_type) plot_grid
    LOGICAL have_mom1
    DOUBLE PRECISION &
       &bs_chi2(bs_setup%nbs,2,assess_grid%nnparam,assess_grid%nanchor),&
       &bs_mom(bs_setup%nbs,0:2,assess_grid%nnparam,assess_grid%nanchor),&
       &bs_mom_std(bs_setup%nbs,2),&
       &bs_anchor_A(bs_setup%nbs,2,assess_grid%nanchor),&
       &bs_anchor_M(bs_setup%nbs,2,assess_grid%nanchor),&
       &bs_param(bs_setup%nbs,assess_grid%max_nparam,2,assess_grid%nnparam,&
       &   assess_grid%nanchor)
    ! Variables for analysis.
    DOUBLE PRECISION t1,t2,t3,t4,t5,t6,t7,t8,t9
    DOUBLE PRECISION chi2(assess_grid%nnparam,assess_grid%nanchor),&
       &dchi2(assess_grid%nnparam,assess_grid%nanchor),&
       &mom(0:2,assess_grid%nnparam,assess_grid%nanchor),&
       &dmom(0:2,assess_grid%nnparam,assess_grid%nanchor),&
       &mom_std(2),dmom_std(2),&
       &anchor_mlogq(2,assess_grid%nanchor),&
       &danchor_mlogq(2,assess_grid%nanchor),&
       &anchor_invA(2,assess_grid%nanchor),&
       &danchor_invA(2,assess_grid%nanchor)
    ! Local variables.
    CHARACTER(9) tstring
    LOGICAL pass_test(6,assess_grid%nnparam,assess_grid%nanchor),&
       &fail_test(6,assess_grid%nnparam,assess_grid%nanchor),&
       &ignore_point(3,assess_grid%nnparam,assess_grid%nanchor)
    INTEGER ianchor,janchor,inparam,jnparam,imom,itest,itail,ierr,ierror,&
       &nparam(2),ibest(2),neg_count,ibs,ipoint
    DOUBLE PRECISION dbest,x_rec,y_rec
    DOUBLE PRECISION, ALLOCATABLE :: delta_mu(:,:)
    ! Test tolerances.
    ! * Ignore-test 1 imposes M>nparam and has no associated tolerance.
    ! * Ignore-test 2: ignore if uncertainty in norm exceeds this:
    DOUBLE PRECISION, PARAMETER :: IGNORE_DNORM_ABS_TOL=0.05d0
    ! * Ignore-test 3: ignore if uncertainty in norm exceeds minimum
    !   uncertainty in norm at lower orders by over this factor:
    DOUBLE PRECISION, PARAMETER :: IGNORE_DNORM_REL_TOL=1.5d0
    ! * Test 1: does not pass if P(A)<0 at any A for a fraction of BS
    !   samples exceeding this value.  This test does not allow a "fail".
    DOUBLE PRECISION, PARAMETER :: PASS_NEG_FRAC=0.51d0
    INTEGER, PARAMETER :: TEST_NEG_NPOINT=101
    ! * Tests 2-6: these pass (and [nparam,anchor] becomes a candidate for
    !   selection) if the relevant estimator deviates by less than this many
    !   stderrs from its target:
    DOUBLE PRECISION, PARAMETER :: PASS_NORM1_TOL=1.5d0
    DOUBLE PRECISION, PARAMETER :: PASS_NORM_REL_TOL=1.5d0
    DOUBLE PRECISION, PARAMETER :: PASS_MEAN_REL_TOL=1.5d0
    DOUBLE PRECISION, PARAMETER :: PASS_VAR_REL_TOL=1.5d0
    DOUBLE PRECISION, PARAMETER :: PASS_CHI2_REL_TOL=3.d0 ! different meaning
    DOUBLE PRECISION, PARAMETER :: PASS_DCHI2_REL_TOL=2.d0
    !   ... and fail (which can be used to "un-pass" tests at other
    !   [nparam,anchor]) if the relevant estimator deviates by more than this
    !   many stderrs from its target:
    DOUBLE PRECISION, PARAMETER :: FAIL_NORM1_TOL=2.5d0
    DOUBLE PRECISION, PARAMETER :: FAIL_NORM_REL_TOL=2.5d0
    DOUBLE PRECISION, PARAMETER :: FAIL_MEAN_REL_TOL=2.5d0
    DOUBLE PRECISION, PARAMETER :: FAIL_VAR_REL_TOL=2.5d0
    DOUBLE PRECISION, PARAMETER :: FAIL_CHI2_REL_TOL=3.d0
    ! Parameters.
    INTEGER, PARAMETER :: io=10

    ! Initialize.
    call initialize_plot_grid(plot_grid)
    anchor_mlogq=0.d0
    danchor_mlogq=0.d0
    anchor_invA=0.d0
    danchor_invA=0.d0

    if(MPI_MASTER)then
      ! Write header.
      write(6,'(a)')'Fit assessment'
      write(6,'(a)')'=============='
      write(6,'(4(1x,a23),1x,a11,4(1x,a27),1x,a9)')&
         &'-mlogq of left anchor--','-mlogq of right anchor-',&
         &'--1/A of left anchor---','--1/A of right anchor--',&
         &'--nparam---','------Norm estimator-------',&
         &'------Mean estimator-------','----Variance estimator-----',&
         &'-----------Chi^2-----------','RDC+0U12X'
    endif ! MPI_MASTER

    ! Open file and write header.
    if(len_trim(filename)>0)then
      if(MPI_MASTER)open(unit=io,file=filename,status='replace',iostat=ierr)
      call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
      if(ierr/=0)then
        call master_msg('Could not open plot file.')
        return
      endif
      if(MPI_MASTER)write(io,'(a)')'# x_qqL  dx_qqL  x_qqR  dx_qqR  &
         &x_recL  dx_recL  x_recR  dx_recR  nparamL  nparamR &
         &mom0  dmom0  mom1  dmom1  mom2  dmom2  chi2  dchi2  #RDC+0U12X'
    endif

    ! Perform all fits at this anchor.
    call bootstrap_fit(dataset,tail_form,bs_setup,plot_grid,assess_grid,&
       &have_mom1,bs_param=bs_param,bs_mom=bs_mom,bs_chi2=bs_chi2,&
       &bs_mom_std=bs_mom_std,bs_anchor_A=bs_anchor_A,bs_anchor_M=bs_anchor_M)

    ! Compute standard estimator for this bootstrap sample.
    do imom=1,2
      call characterize_dist(bs_setup%nbs,bs_mom_std(1,imom),mean=t1,var=t2)
      mom_std(imom)=t1
      dmom_std(imom)=sqrt(t2)
    enddo ! imom

    ! Initialize test passes/failures.
    ignore_point=.false.
    pass_test=.false.
    fail_test=.false.

    ! Loop over anchors to analyze results.
    do ianchor=1,assess_grid%nanchor

      ! Catch out-of-range anchors.
      do itail=1,2
        if(tail_form%mu(itail)<=0.d0)cycle
        ! Check that anchor A is non-zero.
        if(any(are_equal(bs_anchor_A(:,itail,ianchor),0.d0)))&
           &ignore_point(1,:,ianchor)=.true.
        ! Check that anchor M exceeds number of parameters.
        do inparam=1,assess_grid%nnparam
          if(any(bs_anchor_M(:,itail,ianchor)<=&
             &dble(assess_grid%nparam(itail,inparam))))&
             &ignore_point(1,inparam,ianchor)=.true.
        enddo ! inparam
      enddo ! itail

      ! Statistical analysis.
      do inparam=1,assess_grid%nnparam
        chi2(inparam,ianchor)=0.d0
        dchi2(inparam,ianchor)=0.d0
        do itail=1,2
          if(tail_form%mu(itail)<=0.d0)cycle
          call characterize_dist(bs_setup%nbs,bs_chi2(1,itail,inparam,ianchor),&
             &mean=t1,var=t2)
          chi2(inparam,ianchor)=chi2(inparam,ianchor)+t1
          dchi2(inparam,ianchor)=sqrt(dchi2(inparam,ianchor)**2+t2)
        enddo
        do imom=0,2
          call characterize_dist(bs_setup%nbs,bs_mom(1,imom,inparam,ianchor),&
             &mean=t1,var=t2)
          mom(imom,inparam,ianchor)=t1
          dmom(imom,inparam,ianchor)=sqrt(t2)
        enddo ! imom
      enddo ! inparam
      do itail=1,2
        if(tail_form%mu(itail)<=0.d0)cycle
        call characterize_dist(bs_setup%nbs,&
           &-log((bs_anchor_M(:,itail,ianchor)-0.5d0)/dataset%M),&
           &mean=t1,var=t2)
        anchor_mlogq(itail,ianchor)=t1
        danchor_mlogq(itail,ianchor)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,1.d0/bs_anchor_A(:,itail,ianchor),&
           &mean=t1,var=t2)
        anchor_invA(itail,ianchor)=t1
        danchor_invA(itail,ianchor)=sqrt(t2)
      enddo ! itail

      ! Ignore point if norm uncertainty is > absolute tolerance.
      ignore_point(2,:,ianchor)=dmom(0,:,ianchor)>IGNORE_DNORM_ABS_TOL

      ! Ignore point if any lower order has a norm uncertainty which is lower
      ! by more than relative tolerance.
      do inparam=1,assess_grid%nnparam
        ! Skip if already ignoring.
        do jnparam=1,assess_grid%nnparam
          ! Do not compare against self.
          if(jnparam==inparam)cycle
          ! Skip if ignoring.
          ! Only compare against lower-order jnparam.
          if(.not.all(assess_grid%nparam(:,jnparam)<=&
             &assess_grid%nparam(:,inparam)))cycle
          if(dmom(0,inparam,ianchor)>IGNORE_DNORM_REL_TOL*&
             &dmom(0,jnparam,ianchor))ignore_point(3,inparam,ianchor)=.true.
        enddo ! jnparam
      enddo ! inparam

      ! Convergence tests.
      do inparam=1,assess_grid%nnparam

        ! Initialize.
        pass_test(:,inparam,ianchor)=.true.
        fail_test(:,inparam,ianchor)=.false.

        ! Test 1: fit is non-negative on all points in grid.
        ! NB, this test cannot be "failed" - i.e., negative fits do not
        ! invalidate non-negative fits at other points.
        call instantiate_delta_mu(tail_form,&
           &assess_grid%nparam(1:2,inparam),delta_mu)
        do itail=1,2
          if(tail_form%mu(itail)<=0.d0)cycle
          ! Loop over test grid points.
          do ipoint=1,TEST_NEG_NPOINT
            neg_count=0
            do ibs=1,bs_setup%nbs
              x_rec=(dble(ipoint-1)/dble(TEST_NEG_NPOINT-1))/&
                 &bs_anchor_A(ibs,itail,ianchor)
              y_rec=eval_fit(x_rec,assess_grid%nparam(itail,inparam),&
                 &delta_mu(1,itail),bs_param(ibs,&
                 &   1:assess_grid%nparam(itail,inparam),itail,inparam,&
                 &   ianchor))
              if(y_rec<0.d0)neg_count=neg_count+1
            enddo ! ibs
            if(dble(neg_count)>dble(bs_setup%nbs)*PASS_NEG_FRAC)then
              pass_test(1,inparam,ianchor)=.false.
              exit
            endif
          enddo ! ipoint
          if(.not.pass_test(1,inparam,ianchor))exit
        enddo ! itail
        deallocate(delta_mu)

        ! Test 2: norm deviates from unity by < relative tolerance.
        pass_test(2,inparam,ianchor)=abs(mom(0,inparam,ianchor)-1.d0)<&
           &PASS_NORM1_TOL*dmom(0,inparam,ianchor)
        fail_test(2,inparam,ianchor)=abs(mom(0,inparam,ianchor)-1.d0)>&
           &FAIL_NORM1_TOL*dmom(0,inparam,ianchor)

        ! Test 6 part 1: uncertainty in chi2 is less than tolerance times chi2.
        if(dchi2(inparam,ianchor)>PASS_DCHI2_REL_TOL*&
           &max(0.d0,chi2(inparam,ianchor)))&
           &pass_test(6,inparam,ianchor)=.false.

        ! Tests 3, 4, 5 & 6 involve comparing different inparam.
        do jnparam=1,assess_grid%nnparam
          ! Do not compare against self.
          if(jnparam==inparam)cycle
          ! Skip if ignoring.
          if(any(ignore_point(:,jnparam,ianchor)))cycle
          ! Only compare against higher-order jnparam.
          if(.not.all(assess_grid%nparam(:,jnparam)>=&
             &assess_grid%nparam(:,inparam)))cycle

          ! Test 3: norm deviates from higher-order by < relative tolerance.
          if(tail_form%momdef>=MOMDEF_MOM0)then
            if(abs(mom(0,inparam,ianchor)-mom(0,jnparam,ianchor))>&
               &PASS_NORM_REL_TOL*sqrt(dmom(0,inparam,ianchor)**2+&
               &dmom(0,jnparam,ianchor)**2))pass_test(3,inparam,ianchor)=.false.
            if(abs(mom(0,inparam,ianchor)-mom(0,jnparam,ianchor))>&
               &FAIL_NORM_REL_TOL*sqrt(dmom(0,inparam,ianchor)**2+&
               &dmom(0,jnparam,ianchor)**2))fail_test(3,inparam,ianchor)=.true.
          endif

          ! Test 4: mean deviates from higher-order by < relative tolerance.
          if(have_mom1)then
            if(abs(mom(1,inparam,ianchor)-mom(1,jnparam,ianchor))>&
               &PASS_MEAN_REL_TOL*sqrt(dmom(1,inparam,ianchor)**2+&
               &dmom(1,jnparam,ianchor)**2))pass_test(4,inparam,ianchor)=.false.
            if(abs(mom(1,inparam,ianchor)-mom(1,jnparam,ianchor))>&
               &FAIL_MEAN_REL_TOL*sqrt(dmom(1,inparam,ianchor)**2+&
               &dmom(1,jnparam,ianchor)**2))fail_test(4,inparam,ianchor)=.true.
          endif

          ! Test 5: variance deviates from higher-order by < rel. tolerance.
          if(tail_form%momdef>=MOMDEF_MOM2)then
            if(abs(mom(2,inparam,ianchor)-mom(2,jnparam,ianchor))>&
               &PASS_VAR_REL_TOL*sqrt(dmom(2,inparam,ianchor)**2+&
               &dmom(2,jnparam,ianchor)**2))pass_test(5,inparam,ianchor)=.false.
            if(abs(mom(2,inparam,ianchor)-mom(2,jnparam,ianchor))>&
               &FAIL_VAR_REL_TOL*sqrt(dmom(2,inparam,ianchor)**2+&
               &dmom(2,jnparam,ianchor)**2))fail_test(5,inparam,ianchor)=.true.
          endif

          ! Test 6: chi2 deviates from higher-order by < relative tolerance.
          if(chi2(inparam,ianchor)-chi2(jnparam,ianchor)>&
             &PASS_CHI2_REL_TOL*dchi2(jnparam,ianchor))&
             &pass_test(6,inparam,ianchor)=.false. ! this test is signed
          if(abs(chi2(inparam,ianchor)-chi2(jnparam,ianchor))>&
             &FAIL_CHI2_REL_TOL*sqrt(dchi2(inparam,ianchor)**2+&
             &dchi2(jnparam,ianchor)**2))fail_test(6,inparam,ianchor)=.true.
             ! this test is unsigned

        enddo ! jnparam
      enddo ! inparam

    enddo ! ianchor

    ! Dump to plot file prior to cross-nparam/anchor refinement.
    if(len_trim(filename)>0.and.MPI_MASTER)then
      write(io,'(8(1x,es20.13),2(1x,i5),8(1x,es20.13),1x,"#",a9)')&
         &0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0,0,1.d0,0.d0,&
         &mom_std(1),dmom_std(1),mom_std(2),dmom_std(2),&
         &0.d0,0.d0,'---------'
      do ianchor=1,assess_grid%nanchor
        do inparam=1,assess_grid%nnparam
          write(io,'(8(1x,es20.13),2(1x,i5),8(1x,es20.13),1x,"#",a9)')&
             &anchor_mlogq(1,ianchor),danchor_mlogq(1,ianchor),&
             &anchor_mlogq(2,ianchor),danchor_mlogq(2,ianchor),&
             &anchor_invA(1,ianchor),danchor_invA(1,ianchor),&
             &anchor_invA(2,ianchor),danchor_invA(2,ianchor),&
             &assess_grid%nparam(1:2,inparam),&
             &mom(0,inparam,ianchor),dmom(0,inparam,ianchor),&
             &mom(1,inparam,ianchor),dmom(1,inparam,ianchor),&
             &mom(2,inparam,ianchor),dmom(2,inparam,ianchor),&
             &chi2(inparam,ianchor),dchi2(inparam,ianchor),&
             &test_string(ignore_point(1,inparam,ianchor),&
             &            tail_form%momdef,have_mom1,&
             &            pass_test(1,inparam,ianchor),&
             &            fail_test(1,inparam,ianchor))
        enddo ! inparam
      enddo ! ianchor
    endif ! len_trim(filename)>0.and.MPI_MASTER

    ! Refine test results by comparing different inparam: tests failed at
    ! a certain inparam and ianchor should not pass at the same ianchor
    ! for lower expansion orders.
    do inparam=1,assess_grid%nnparam
      do jnparam=1,assess_grid%nnparam
        if(inparam==jnparam)cycle
        ! Check that inparam is lower than or equal to jnparam.
        if(any(assess_grid%nparam(:,inparam)>assess_grid%nparam(:,jnparam)))&
           &cycle
        do ianchor=1,assess_grid%nanchor
          if(any(ignore_point(:,inparam,ianchor).or.&
             &ignore_point(:,jnparam,ianchor)))cycle
          do itest=1,6
            if(fail_test(itest,jnparam,ianchor))&
               &pass_test(itest,inparam,ianchor)=.false.
          enddo ! itest
        enddo ! ianchor
      enddo ! jnparam
    enddo ! inparam

    ! Refine test results by comparing different ianchor: tests failed at
    ! a certain inparam and ianchor should not pass at the same inparam
    ! for more internal anchors.
    do ianchor=1,assess_grid%nanchor
      do janchor=1,assess_grid%nanchor
        if(janchor==ianchor)cycle
        ! Check that ianchor is more internal than janchor on both tails.
        do itail=1,2
          if(tail_form%mu(itail)<=0.d0)cycle
          ! Ignore test if anchors defined by different variables.
          if(assess_grid%def_anchor(itail,ianchor)/=&
             &assess_grid%def_anchor(itail,janchor))exit
          select case(assess_grid%def_anchor(itail,ianchor))
          case('q')
            if(assess_grid%anchor(itail,ianchor)<&
               &assess_grid%anchor(itail,janchor))exit
          case('A')
            if(assess_grid%anchor(itail,ianchor)>&
               &assess_grid%anchor(itail,janchor))exit
          end select
        enddo ! itail
        if(itail<=2)cycle
        do inparam=1,assess_grid%nnparam
          if(any(ignore_point(:,inparam,ianchor).or.&
             &ignore_point(:,inparam,janchor)))cycle
          do itest=1,6
            if(fail_test(itest,inparam,janchor))&
               &pass_test(itest,inparam,ianchor)=.false.
          enddo ! itest
        enddo ! inparam
      enddo ! janchor
    enddo ! ianchor

    ! Print standard estimator for this bootstrap sample.
    if(MPI_MASTER)write(6,'(8(1x,es11.4),2(1x,i5),8(1x,es13.6),1x,a9)')&
       &0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0,0,1.d0,0.d0,&
       &mom_std(1),dmom_std(1),mom_std(2),dmom_std(2),&
       &0.d0,0.d0,'---------'

    ! Loop over anchors and print best inparam.
    ibest=0
    dbest=0.d0
    do ianchor=1,assess_grid%nanchor
      inparam=0
      ! Pick best inparam.
      do jnparam=1,assess_grid%nnparam
        if(any(ignore_point(:,jnparam,ianchor)))cycle
        if(any(.not.pass_test(:,jnparam,ianchor)))cycle
        if(inparam==0)then
          inparam=jnparam
        elseif(all(assess_grid%nparam(:,jnparam)<=&
           &assess_grid%nparam(:,inparam)))then
          inparam=jnparam
        endif ! inparam==0 or jnparam "<" inparam
      enddo ! jnparam
      ! Prepare report.
      nparam=0
      t1=0.d0
      t2=0.d0
      t3=0.d0
      t4=0.d0
      t5=0.d0
      t6=0.d0
      t7=0.d0
      t8=0.d0
      tstring='---------'
      if(inparam>0)then
        nparam=assess_grid%nparam(1:2,inparam)
        t1=mom(0,inparam,ianchor)
        t2=dmom(0,inparam,ianchor)
        t3=mom(1,inparam,ianchor)
        t4=dmom(1,inparam,ianchor)
        t5=mom(2,inparam,ianchor)
        t6=dmom(2,inparam,ianchor)
        t7=chi2(inparam,ianchor)
        t8=dchi2(inparam,ianchor)
        ! Select convergence critetion depending on available quantities.
        if(tail_form%momdef>=MOMDEF_MOM2)then
          t9=t6
        elseif(have_mom1)then
          t9=t4
        elseif(tail_form%momdef>=MOMDEF_MOM0)then
          t9=t2
        else
          t9=t8
        endif
        ! Store "best" choice of ianchor,inparam.
        if(t9<dbest.or.all(ibest==0))then
          ibest=(/ianchor,inparam/)
          dbest=t9
        endif
        tstring=test_string(ignore_point(1,inparam,ianchor),&
           &tail_form%momdef,have_mom1,pass_test(1,inparam,ianchor),&
           &fail_test(1,inparam,ianchor))
      endif
      ! Report.
      if(MPI_MASTER)write(6,'(8(1x,es11.4),2(1x,i5),8(1x,es13.6),1x,a9)')&
         &anchor_mlogq(1,ianchor),danchor_mlogq(1,ianchor),&
         &anchor_mlogq(2,ianchor),danchor_mlogq(2,ianchor),&
         &anchor_invA(1,ianchor),danchor_invA(1,ianchor),&
         &anchor_invA(2,ianchor),danchor_invA(2,ianchor),&
         &nparam,t1,t2,t3,t4,t5,t6,t7,t8,tstring
    enddo ! ianchor
    if(MPI_MASTER)write(6,'()')

    ! Close plot file.
    if(len_trim(filename)>0.and.MPI_MASTER)then
      close(io)
      call master_msg('Plot written to "'//trim(filename)//'".')
    endif

    ! Suggest values for anchor,nparam.
    ianchor=ibest(1)
    inparam=ibest(2)
    if(all(ibest>0))then
      if(MPI_MASTER)then
        write(6,'(a)')'Suggested anchor and nparam:'
        write(6,'(a)',advance='no')'ASSESS  anchor'
        do itail=1,2
          write(6,'(a)',advance='no')'  '//&
             &trim(assess_grid%def_anchor(itail,ianchor))
          write(6,'(2x,es20.13)',advance='no')assess_grid%anchor(itail,ianchor)
        enddo ! itail
        write(6,'()')
        write(6,'(a)',advance='no')'ASSESS  anchor'
        do itail=1,2
          select case(assess_grid%def_anchor(itail,ianchor))
          case('q')
            t1=-log(assess_grid%anchor(itail,ianchor))
            write(6,'(2x,a,2x,es20.13)',advance='no')'mlogq',t1
          case('A')
            t1=1.d0/assess_grid%anchor(itail,ianchor)**&
               &tail_form%global_delta_mu(itail)
            write(6,'(2x,a,2x,es20.13)',advance='no')'x',t1
          end select
        enddo ! itail
        write(6,'()')
        write(6,'(a,2(1x,i3))')'ASSESS  nparam  ',assess_grid%nparam(:,inparam)
      endif
      if(apply_assessment)then
        call master_msg('Setting anchor and nparam to suggested values.')
        tail_form%def_anchor(:)=assess_grid%def_anchor(:,ianchor)
        tail_form%anchor(:)=assess_grid%anchor(:,ianchor)
        tail_form%nparam(:)=assess_grid%nparam(:,inparam)
      endif
    else
      if(MPI_MASTER)then
        write(6,'(a)')'No anchor or nparam to suggest in tested ranges.'
        write(6,'(a,2(1x,es20.13))')'ASSESS  anchor  ',0.d0,0.d0
        write(6,'(a)',advance='no')'ASSESS  anchor'
        do itail=1,2
          select case(assess_grid%def_anchor(itail,ianchor))
          case('q')
            write(6,'(2x,a,2x,es20.13)',advance='no')'mlogq',0.d0
          case('A')
            write(6,'(2x,a,2x,es20.13)',advance='no')'x',0.d0
          end select
        enddo ! itail
        write(6,'()')
        write(6,'(a,2(1x,i3))')'ASSESS  nparam  ',0,0
      endif
      if(apply_assessment)then
        call master_msg('Unsetting anchor and nparam.')
        tail_form%def_anchor(:)=''
        tail_form%anchor(:)=0.d0
        tail_form%nparam(:)=0
      endif
    endif

  END SUBROUTINE assess_fit


  CHARACTER(9) FUNCTION test_string(ignore_point,momdef,have_mom1,pass_test,&
     &fail_test)
    !------------------------------------------------!
    ! Generate a string detailing test pass/failure. !
    !------------------------------------------------!
    IMPLICIT NONE
    LOGICAL, INTENT(in) :: ignore_point(3),have_mom1,pass_test(6),fail_test(6)
    INTEGER, INTENT(in) :: momdef
    if(ignore_point(1))then
      test_string(1:1)='I'
    else
      test_string(1:1)='>'
    endif
    if(ignore_point(2))then
      test_string(2:2)='I'
    else
      test_string(2:2)='>'
    endif
    if(ignore_point(3))then
      test_string(3:3)='I'
    else
      test_string(3:3)='>'
    endif
    test_string(4:4)='~'
    if(pass_test(1))test_string(4:4)='P'
    if(fail_test(1))test_string(4:4)='X'
    test_string(5:6)='--'
    if(momdef>=MOMDEF_MOM0)then
      test_string(5:6)='~~'
      if(pass_test(2))test_string(5:5)='P'
      if(fail_test(2))test_string(5:5)='X'
      if(pass_test(3))test_string(6:6)='P'
      if(fail_test(3))test_string(6:6)='X'
    endif
    test_string(7:7)='-'
    if(have_mom1)then
      test_string(7:7)='~'
      if(pass_test(4))test_string(7:7)='P'
      if(fail_test(4))test_string(7:7)='X'
    endif
    test_string(8:8)='-'
    if(momdef>=MOMDEF_MOM2)then
      test_string(8:8)='~'
      if(pass_test(5))test_string(8:8)='P'
      if(fail_test(5))test_string(8:8)='X'
    endif
    test_string(9:9)='~'
    if(pass_test(6))test_string(9:9)='P'
    if(fail_test(6))test_string(9:9)='X'
  END FUNCTION test_string


  ! PLOT ROUTINES


  SUBROUTINE plot_model(model_P,plot_grid,filename)
    !-------------------------------------------!
    ! Generate qq-plot of model P, which should !
    ! roughly match qq-plot of data.            !
    !-------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    TYPE(plot_grid_type),INTENT(in) :: plot_grid
    CHARACTER(*),INTENT(in) :: filename
    ! Local variables.
    INTEGER i,ierr
    DOUBLE PRECISION q,A
    ! Parameters.
    INTEGER,PARAMETER :: io=10

    if(MPI_MASTER)then

      ! Open file.
      open(unit=io,file=trim(filename),status='replace',iostat=ierr)
      if(ierr/=0)return
      write(io,'(a)')'# q  A'

      ! Plot model CDF.
      do i=1,plot_grid%npoint
        select case(plot_grid%def_point(i))
        case('q')
          q=plot_grid%point(i)
          A=invert_model_P_CDF(model_P,q)
        case('A')
          A=plot_grid%point(i)
          q=model_P_CDF(model_P,A)
        end select
        write(io,'(2(1x,es20.12))')q,A
      enddo ! i

      ! Close file.
      close(io)

    endif ! MPI_MASTER

  END SUBROUTINE plot_model


  ! MAIN TRE ROUTINE


  SUBROUTINE evaluate_TRE(dataset,tail_form,bs_setup,plot_P_fname,&
     &plot_P_ngrid,plot_fit_fname,plot_fit_ngrid)
    !----------------------------------------------------------!
    ! Estimate the mean of the data in DATASET as per the tail !
    ! fit parameters in TAIL_FORM, and report the result.      !
    !----------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    TYPE(tail_form_type),INTENT(inout) :: tail_form
    TYPE(bs_setup_type),INTENT(in) :: bs_setup
    CHARACTER(*),INTENT(in) :: plot_P_fname,plot_fit_fname
    INTEGER,INTENT(in) :: plot_P_ngrid,plot_fit_ngrid
    ! Bootstrapping arrays,
    DOUBLE PRECISION &
       &tail_mom(0:2,2),dtail_mom(0:2,2),&
       &mom(0:2),dmom(0:2),&
       &centre_mom(0:2),dcentre_mom(0:2),&
       &param(maxval(tail_form%nparam),2),&
       &dparam(maxval(tail_form%nparam),2),&
       &anchor_A(2),danchor_A(2),&
       &anchor_logA(2),danchor_logA(2),&
       &anchor_invA(2),danchor_invA(2),&
       &anchor_mlogq(2),danchor_mlogq(2),&
       &anchor_M(2),danchor_M(2),&
       &mom_std(2),dmom_std(2)
    DOUBLE PRECISION &
       &bs_tail_mom(bs_setup%nbs,0:2,2),&
       &bs_mom(bs_setup%nbs,0:2),&
       &bs_centre_mom(bs_setup%nbs,0:2),&
       &bs_param(bs_setup%nbs,maxval(tail_form%nparam),2),&
       &bs_anchor_A(bs_setup%nbs,2),&
       &bs_anchor_M(bs_setup%nbs,2),&
       &bs_mom_std(bs_setup%nbs,2)
    DOUBLE PRECISION &
       &bs_A(bs_setup%nbs,2*plot_P_ngrid),&
       &bs_M(bs_setup%nbs,2*plot_P_ngrid),&
       &bs_P(bs_setup%nbs,2*plot_P_ngrid)
    ! Assessment grid for bootstrap routine.
    TYPE(assess_grid_type) :: assess_grid
    ! Variables for plots.
    TYPE(plot_grid_type) :: plot_grid
    INTEGER ibs,ipoint,k,k1,k2,ierr
    INTEGER, PARAMETER :: io=10
    DOUBLE PRECISION inv_A,x_rec,x_P,x_qq,y_rec,y_P,y_qq,x_qq_CI(4),&
       &y_rec_CI(4),y_P_CI(4),y_qq_CI(4),x1,x2,logx,logx1,logx2
    DOUBLE PRECISION, ALLOCATABLE :: bs_x_qq(:),bs_y_qq(:),bs_y_rec(:),&
       &bs_y_P(:),delta_mu(:,:)
    ! Local variables.
    LOGICAL have_mom1,is_htail(2),is_done
    INTEGER itail,iparam,imom,iiter
    DOUBLE PRECISION t1,t2
    ! Parameters.
    INTEGER :: MAX_NITER=10

    ! Print header.
    if(MPI_MASTER)then
      write(6,'(a)')'Tail-regression estimation'
      write(6,'(a)')'=========================='
    endif

    ! Decide which tails are to be regarded heavy in this routine.
    is_htail=tail_form%mu>0.d0.and.len_trim(tail_form%def_anchor)>0.and.&
       &tail_form%nparam>0

    ! Initialize things itail-dependent.
    tail_mom=0.d0
    dtail_mom=0.d0
    param=0.d0
    dparam=0.d0
    anchor_A=0.d0
    danchor_A=0.d0
    anchor_logA=0.d0
    danchor_logA=0.d0
    anchor_invA=0.d0
    danchor_invA=0.d0
    anchor_mlogq=0.d0
    danchor_mlogq=0.d0
    anchor_M=0.d0
    danchor_M=0.d0
    call set_assess_grid_from_tail_form(tail_form,assess_grid)

    ! Construct grid for plotting P if requested.
    plot_grid%npoint=0
    if(len_trim(plot_P_fname)>0)then
      if(allocated(plot_grid%point))deallocate(plot_grid%point)
      if(allocated(plot_grid%def_point))deallocate(plot_grid%def_point)
      plot_grid%npoint=2*plot_P_ngrid
      allocate(plot_grid%point(plot_grid%npoint),&
         &plot_grid%def_point(plot_grid%npoint))
      plot_grid%def_point='A'
      ! Construct range for left tail.
      k1=1
      k2=tail_form%Mhalf(1)
      x1=0.5d0*(dataset%A(dataset%indx(k1))+dataset%A(dataset%indx(k1+1)))
      x2=0.5d0*(dataset%A(dataset%indx(k2))+dataset%A(dataset%indx(k2+1)))
      logx1=log(abs(x1-tail_form%A_centre))
      logx2=log(abs(x2-tail_form%A_centre))
      do ipoint=1,plot_P_ngrid
        logx=logx1+(dble(ipoint-1)/dble(plot_P_ngrid-1))*(logx2-logx1)
        plot_grid%point(ipoint)=tail_form%A_centre-exp(logx)
      enddo ! ipoint
      ! Construct range for right tail.
      k1=dataset%M-tail_form%Mhalf(1)+2
      k2=dataset%M
      x1=0.5d0*(dataset%A(dataset%indx(k1))+dataset%A(dataset%indx(k1-1)))
      x2=0.5d0*(dataset%A(dataset%indx(k2))+dataset%A(dataset%indx(k2-1)))
      logx1=log(abs(x1-tail_form%A_centre))
      logx2=log(abs(x2-tail_form%A_centre))
      do ipoint=1,plot_P_ngrid
        logx=logx1+(dble(ipoint-1)/dble(plot_P_ngrid-1))*(logx2-logx1)
        plot_grid%point(plot_P_ngrid+ipoint)=tail_form%A_centre+exp(logx)
      enddo ! ipoint
    endif

    ! Loop over self-consistence loop iterations.
    do iiter=1,MAX_NITER

      ! Print iteration heading
      if(MPI_MASTER.and.tail_form%self_consistent_centre)write(6,'(a)')&
         &'Iteration '//trim(i2s(iiter))//':'

      ! Perform combined bootstrap.
      if(len_trim(plot_P_fname)>0)then
        call bootstrap_fit(dataset,tail_form,bs_setup,plot_grid,assess_grid,&
           &have_mom1,bs_param=bs_param,bs_tail_mom=bs_tail_mom,&
           &bs_centre_mom=bs_centre_mom,bs_mom=bs_mom,bs_anchor_A=bs_anchor_A,&
           &bs_anchor_M=bs_anchor_M,bs_mom_std=bs_mom_std,&
           &bs_A=bs_A,bs_M=bs_M,bs_P=bs_P)
      else
        call bootstrap_fit(dataset,tail_form,bs_setup,plot_grid,assess_grid,&
           &have_mom1,bs_param=bs_param,bs_tail_mom=bs_tail_mom,&
           &bs_centre_mom=bs_centre_mom,bs_mom=bs_mom,bs_anchor_A=bs_anchor_A,&
           &bs_anchor_M=bs_anchor_M,bs_mom_std=bs_mom_std)
      endif

      ! Compute statistics on parameters.
      do itail=1,2
        if(.not.is_htail(itail))cycle
        do iparam=1,tail_form%nparam(itail)
          call characterize_dist(bs_setup%nbs,bs_param(1,iparam,itail),&
             &mean=t1,var=t2)
          param(iparam,itail)=t1
          dparam(iparam,itail)=sqrt(t2)
        enddo ! iparam
      enddo ! itail

      ! Obtain moments from standard estimators.
      do imom=1,2
        call characterize_dist(bs_setup%nbs,bs_mom_std(1,imom),mean=t1,var=t2)
        mom_std(imom)=t1
        dmom_std(imom)=sqrt(t2)
      enddo ! imom

      ! Obtain partial and total contributions to moments from TRE.
      do imom=0,2
        do itail=1,2
          if(.not.is_htail(itail))cycle
          call characterize_dist(bs_setup%nbs,bs_tail_mom(1,imom,itail),&
             &mean=t1,var=t2)
          tail_mom(imom,itail)=t1
          dtail_mom(imom,itail)=sqrt(t2)
        enddo ! itail
        call characterize_dist(bs_setup%nbs,bs_centre_mom(1,imom),&
           &mean=t1,var=t2)
        centre_mom(imom)=t1
        dcentre_mom(imom)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,bs_mom(1,imom),mean=t1,var=t2)
        mom(imom)=t1
        dmom(imom)=sqrt(t2)
      enddo ! imom

      ! Obtain bootstrapped anchor values.
      do itail=1,2
        if(.not.is_htail(itail))cycle
        call characterize_dist(bs_setup%nbs,&
           &bs_anchor_A(:,itail)+tail_form%A_centre,mean=t1,var=t2)
        anchor_A(itail)=t1
        danchor_A(itail)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,log(bs_anchor_A(:,itail)),&
           &mean=t1,var=t2)
        anchor_logA(itail)=t1
        danchor_logA(itail)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,1.d0/bs_anchor_A(:,itail),&
           &mean=t1,var=t2)
        anchor_invA(itail)=t1
        danchor_invA(itail)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,&
           &-log((bs_anchor_M(:,itail)-0.5d0)/dataset%M),mean=t1,var=t2)
        anchor_mlogq(itail)=t1
        danchor_mlogq(itail)=sqrt(t2)
        call characterize_dist(bs_setup%nbs,bs_anchor_M(:,itail),&
           &mean=t1,var=t2)
        anchor_M(itail)=t1
        danchor_M(itail)=sqrt(t2)
      enddo ! itail

      ! Report.
      if(MPI_MASTER)then
        ! Report 1st and 2nd moments with standard estimator.
        write(6,'("* ",a,":")')'Standard estimator'
        write(6,'(2x,"* ",a,t25,":",2(1x,es20.12))')'Mean',&
           &mom_std(1),dmom_std(1)
        write(6,'(2x,"* ",a,t25,":",2(1x,es20.12))')'Variance',&
           &mom_std(2),dmom_std(2)
        ! Report Ac.
        write(6,'("* ",a,":")')'Tail-regression estimator'
        write(6,'(2x,"* ",a,t25,":",1x,es20.12)')'Ac',tail_form%A_centre
        ! Per tail report.
        do itail=1,2
          write(6,'(2x,"* ",a,":")')'Tail '//trim(i2s(itail))
          if(.not.is_htail(itail))then
            write(6,'(4x,"* ",a)')'Not a heavy tail'
            cycle
          endif
          ! Report anchor.
          write(6,'(4x,"* ",a,":")')'Fit anchor'
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'A',&
             &anchor_A(itail),danchor_A(itail)
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'log|A-Ac|',&
             &anchor_logA(itail),danchor_logA(itail)
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'|A-Ac|^-1',&
             &anchor_invA(itail),danchor_invA(itail)
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'M',&
             &anchor_M(itail),danchor_M(itail)
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'q',&
             &anchor_M(itail)/dble(dataset%M),danchor_M(itail)/dble(dataset%M)
          write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')'-logq',&
             &anchor_mlogq(itail),danchor_mlogq(itail)
          ! Report fit coefficients.
          write(6,'(4x,"* ",a,":")')'Coefficients'
          do iparam=1,tail_form%nparam(itail)
            write(6,'(6x,"* ",a,t25,":",2(1x,es20.12))')&
               &'c_'//trim(i2s(iparam)),&
               &param(iparam,itail)*(tail_form%mu(itail)+dble(iparam-2)),&
               &dparam(iparam,itail)*(tail_form%mu(itail)+dble(iparam-2))
          enddo ! iparam
          ! Report tail integrals.
          if(tail_form%momdef>=MOMDEF_MOM0)then
            write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to norm',&
               &tail_mom(0,itail),dtail_mom(0,itail)
            if(have_mom1)then
              write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to mean',&
                 &tail_mom(1,itail),dtail_mom(1,itail)
              if(tail_form%momdef>=MOMDEF_MOM2)then
                write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to var.',&
                   &tail_mom(2,itail),dtail_mom(2,itail)
              endif
            endif
          endif
        enddo ! itail
        ! Report central moments.
        write(6,'(2x,"* ",a,":")')'Centre'
        write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to norm',&
           &centre_mom(0),dcentre_mom(0)
        write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to mean',&
           &centre_mom(1),dcentre_mom(1)
        write(6,'(4x,"* ",a,t25,":",2(1x,es20.12))')'Contrib. to var.',&
           &centre_mom(2),dcentre_mom(2)
        ! Report total moments.
        if(tail_form%momdef>=MOMDEF_MOM0)then
          write(6,'(2x,"* ",a,t25,":",2(1x,es20.12))')'Norm',mom(0),dmom(0)
          if(have_mom1)then
            write(6,'(2x,"* ",a,t25,":",2(1x,es20.12))')'Mean',mom(1),dmom(1)
            if(tail_form%momdef>=MOMDEF_MOM2)then
              write(6,'(2x,"* ",a,t25,":",2(1x,es20.12))')'Variance',&
                 &mom(2),dmom(2)
            endif
          endif
        endif
        write(6,'()')
      endif ! MPI_MASTER

      ! Set result.
      tail_form%A_result_set=.true.
      tail_form%A_result=mom(1)
      tail_form%dA_result=dmom(1)

      ! Decide if we need to loop.
      if(.not.tail_form%self_consistent_centre)exit
      is_done=abs(tail_form%A_result-tail_form%A_centre)<&
         &2.d0*tail_form%dA_result
      call refresh_centre(dataset,tail_form)
      if(is_done)exit

    enddo ! iiter

    ! Dump P plot.
    if(len_trim(plot_P_fname)>0.and.MPI_MASTER)then

      ! Prepare for plot.
      open(unit=io,file=trim(plot_P_fname),status='replace',iostat=ierr)
      if(ierr/=0)then
        call master_msg('Could not open P plot file "'//trim(plot_P_fname)//&
           &'".')
        return
      endif
      write(io,'(a)')'# x_qq|CIs  y_qq  x_rec  y_rec|CIs  x_P  y_P|CIs'
      allocate(bs_x_qq(bs_setup%nbs),bs_y_rec(bs_setup%nbs))

      ! Loop over plot grid points.
      ipoint=0
      do itail=1,2
        do k=1,plot_P_ngrid
          ipoint=ipoint+1
          ! Transform scales.
          if(itail==1)then
            bs_x_qq=-log((bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
          else
            bs_x_qq=-log(1.d0-(bs_M(:,ipoint)-0.5d0)/dble(dataset%M))
          endif
          y_qq=log(abs(bs_A(1,ipoint)-tail_form%A_centre))
          x_rec=1.d0/abs(bs_A(1,ipoint)-tail_form%A_centre)**&
             &tail_form%global_delta_mu(itail)
          if(tail_form%mu(itail)>0.d0)then
            bs_y_rec=exp(-bs_x_qq)*abs(bs_A(:,ipoint)-tail_form%A_centre)**&
               &(tail_form%mu(itail)-1.d0)
          else
            bs_y_rec=0.d0
          endif
          x_P=bs_A(1,ipoint)
          ! Evaluate statistics.
          call get_conf_intervals(bs_setup%nbs,bs_x_qq,x_qq,x_qq_CI)
          call get_conf_intervals(bs_setup%nbs,bs_y_rec,y_rec,y_rec_CI)
          call get_conf_intervals(bs_setup%nbs,bs_P(:,ipoint),y_P,y_P_CI)
          ! Write P(A) in qq, rec, and linear scale.
          write(io,'(18(1x,es20.12))')x_qq,x_qq_CI,y_qq,x_rec,y_rec,y_rec_CI,&
             &x_P,y_P,y_P_CI
        enddo ! k
      enddo ! itail

      ! Tidy up.
      deallocate(bs_x_qq,bs_y_rec)
      close(io)
      write(6,'(a)')'Plotted P to "'//trim(plot_P_fname)//'".'
      write(6,'()')

    endif

    ! Dump fit plot.
    if(len_trim(plot_fit_fname)>0.and.MPI_MASTER)then

      ! Prepare for plot.
      open(unit=io,file=trim(plot_fit_fname),status='replace',iostat=ierr)
      if(ierr/=0)then
        call master_msg('Could not open fit plot file "'//&
           &trim(plot_fit_fname)//'".')
        return
      endif
      write(io,'(a)')'# x_qq  y_qq|CIs  x_rec  y_rec|CIs  x_P  y_P|CIs'
      allocate(bs_y_qq(bs_setup%nbs),bs_y_rec(bs_setup%nbs),&
         &bs_y_P(bs_setup%nbs))
      call instantiate_delta_mu(tail_form,tail_form%nparam,delta_mu)

      ! Loop over tails.
      do itail=1,2
        if(.not.is_htail(itail))cycle
        ! Loop over plot grid points.
        do ipoint=1,plot_fit_ngrid
          inv_A=(dble(ipoint-1)/dble(plot_fit_ngrid-1))**&
             &(1.d0/tail_form%global_delta_mu(itail))*anchor_invA(itail)
          if(ipoint==1)then
            ! Deal with x=0 separately.
            x_rec=0.d0
            x_P=0.d0
            x_qq=0.d0
          else
            x_rec=inv_A**tail_form%global_delta_mu(itail)
            x_P=tail_form%A_centre
            if(itail==1)then
              x_P=x_P-1.d0/inv_A
            else
              x_P=x_P+1.d0/inv_A
            endif
            k=index_of_xindx(x_P,dataset%M,dataset%A,dataset%indx)
            k=max(1,min(k,dataset%M))
            if(itail==1)then
              x_qq=-log((k-0.5d0)/dble(dataset%M))
            else
              x_qq=-log(1.d0-(k-0.5d0)/dble(dataset%M))
            endif
          endif
          ! Compute fit value at each BS sample.
          do ibs=1,bs_setup%nbs
            bs_y_rec(ibs)=eval_fit(inv_A,tail_form%nparam(itail),&
               &delta_mu(1,itail),&
               &bs_param(ibs,1:tail_form%nparam(itail),itail))
            bs_y_P(ibs)=eval_fit_P(tail_form%mu(itail),inv_A,&
               &tail_form%nparam(itail),delta_mu(1,itail),&
               &bs_param(ibs,1:tail_form%nparam(itail),itail))
          enddo ! ibs
          bs_y_qq=0.d0
          where(bs_y_rec>0.d0)bs_y_qq=(log(bs_y_rec)+x_qq)/&
             &(tail_form%mu(itail)-1.d0)
          ! Evaluate statistics.
          call get_conf_intervals(bs_setup%nbs,bs_y_qq,y_qq,y_qq_CI)
          call get_conf_intervals(bs_setup%nbs,bs_y_rec,y_rec,y_rec_CI)
          call get_conf_intervals(bs_setup%nbs,bs_y_P,y_P,y_P_CI)
          ! Write fit in qq, rec, and linear scale.
          write(io,'(30(1x,es20.12))')x_qq,y_qq,y_qq_CI,x_rec,&
             &y_rec,y_rec_CI,x_P,y_P,y_P_CI
        enddo ! ipoint
      enddo ! itail

      ! Tidy up.
      deallocate(bs_y_qq,bs_y_rec,bs_y_P,delta_mu)
      close(io)
      write(6,'(a)')'Plotted fit to "'//trim(plot_fit_fname)//'".'
      write(6,'()')

    endif

    ! Greppable line.
    if(MPI_MASTER)then
      write(6,'("EVAL  Std.mean ",2(1x,es20.12))')mom_std(1),dmom_std(1)
      write(6,'("EVAL  Std.var. ",2(1x,es20.12))')mom_std(2),dmom_std(2)
      if(tail_form%momdef>=MOMDEF_MOM0)then
        write(6,'("EVAL  Norm     ",2(1x,es20.12))')mom(0),dmom(0)
        if(have_mom1)then
          write(6,'("EVAL  Mean     ",2(1x,es20.12))')mom(1),dmom(1)
          if(tail_form%momdef>=MOMDEF_MOM2)then
            write(6,'("EVAL  Variance ",2(1x,es20.12))')mom(2),dmom(2)
          endif
        endif
      endif
      write(6,'()')
    endif ! MPI_MASTER

  END SUBROUTINE evaluate_TRE


  ! MAIN BOOSTSTRAPPED TRE ROUTINE.


  SUBROUTINE bootstrap_fit(dataset,tail_form,bs_setup,plot_grid,assess_grid,&
     &have_mom1,bs_param,bs_tail_mom,bs_centre_mom,bs_mom,bs_chi2,&
     &bs_anchor_A,bs_anchor_M,bs_mom_std,bs_A,bs_M,bs_P,bs_TIE_param,&
     &bs_TIE_chi2)
    !----------------------------------------------------------!
    ! Perform bootstrap on whole dataset, optionally returning !
    ! bootstrap samples of:                                    !
    ! * parameters bs_param at all requested expansion orders  !
    !   and anchors,                                           !
    ! * 0th, 1st, 2nd tail moments bs_tail_mom at all          !
    !   requested expansion orders and anchors,                !
    ! * 0th, 1st, 2nd central moments bs_centre_mom at all     !
    !   requested anchors,                                     !
    ! * 0th, 1st, 2nd total moments bs_mom at all requested    !
    !   exoansion orders and anchors,                          !
    ! * fit chi^2 bs_chi2 at all requested expansion orders    !
    !   and anchors,                                           !
    ! * value of A at the tail anchor bs_anchor_A at all       !
    !   requested anchors,                                     !
    ! * number of tail data bs_anchor_M at all requested       !
    !   anchors,                                               !
    ! * standard estimator of 1st and 2nd total moments        !
    !   bs_mom_std,                                            !
    ! * value of A at each requested point bs_A,               !
    ! * value of M at each requested point bs_M,               !
    ! * value of P at each requested point bs_P,               !
    ! * linear-fit parameters bs_TIE_param at each requested   !
    !   point,                                                 !
    ! * linear-fit chi^2 bs_TIE_chi2 at each requested point.  !
    !----------------------------------------------------------!
    IMPLICIT NONE
    ! Input arguments.
    TYPE(dataset_type), INTENT(in) :: dataset
    TYPE(tail_form_type), INTENT(in) :: tail_form
    TYPE(bs_setup_type), INTENT(in) :: bs_setup
    TYPE(plot_grid_type), INTENT(in) :: plot_grid
    TYPE(assess_grid_type), INTENT(in) :: assess_grid
    LOGICAL, INTENT(inout) :: have_mom1
    ! Output arguments.
    DOUBLE PRECISION, INTENT(inout), OPTIONAL :: &
       &bs_param(bs_setup%nbs,assess_grid%max_nparam,2,&
       &         assess_grid%nnparam,assess_grid%nanchor),&
       &bs_tail_mom(bs_setup%nbs,0:2,2,assess_grid%nnparam,&
       &            assess_grid%nanchor),&
       &bs_centre_mom(bs_setup%nbs,0:2,assess_grid%nanchor),&
       &bs_mom(bs_setup%nbs,0:2,assess_grid%nnparam,assess_grid%nanchor),&
       &bs_chi2(bs_setup%nbs,2,assess_grid%nnparam,assess_grid%nanchor),&
       &bs_anchor_A(bs_setup%nbs,2,assess_grid%nanchor),&
       &bs_anchor_M(bs_setup%nbs,2,assess_grid%nanchor),&
       &bs_mom_std(bs_setup%nbs,2),&
       &bs_A(bs_setup%nbs,plot_grid%npoint),&
       &bs_M(bs_setup%nbs,plot_grid%npoint),&
       &bs_P(bs_setup%nbs,plot_grid%npoint),&
       &bs_TIE_param(bs_setup%nbs,2,plot_grid%npoint),&
       &bs_TIE_chi2(bs_setup%nbs,plot_grid%npoint)
    ! Local variables.
    INTEGER ibs,itail,j,iparam,imom,ierror
    LOGICAL mom1def
    DOUBLE PRECISION t1,t2,t3
    DOUBLE PRECISION, ALLOCATABLE :: delta_mu(:,:)
    ! Anchor/nparam/plot-point loops.
    INTEGER ipoint,inparam,nparam(2),nparam_dim(2),ianchor,anchor_M(2)
    DOUBLE PRECISION anchor_A(2)
    ! Resample variables.
    INTEGER indx_work(dataset%M)
    ! Regularization.
    DOUBLE PRECISION :: regx(2),regy(2)
    ! Linear system.
    DOUBLE PRECISION param(assess_grid%max_nparam,2),chi2(2)
    DOUBLE PRECISION &
       &xx_expval((assess_grid%max_nparam*(assess_grid%max_nparam+1))/2,2),&
       &xy_expval(assess_grid%max_nparam,2),yy_expval(2)
    ! Constraints.
    TYPE(constraint_set_type) cset
    ! Central region.
    INTEGER M_work
    DOUBLE PRECISION cmom0,cmom1,cmom2
    ! Total moments.
    DOUBLE PRECISION tmom0,tmom1,tmom2
    DOUBLE PRECISION mom0(2),mom1(2),mom2(2)
    ! Internal flags.
    LOGICAL need_integrals,need_fit,need_TIE,is_htail(2)

    ! Initialize.
    if(present(bs_param))bs_param=0.d0
    if(present(bs_tail_mom))bs_tail_mom=0.d0
    if(present(bs_centre_mom))bs_centre_mom=0.d0
    if(present(bs_mom))bs_mom=0.d0
    if(present(bs_chi2))bs_chi2=0.d0
    if(present(bs_anchor_A))bs_anchor_A=0.d0
    if(present(bs_anchor_M))bs_anchor_M=0.d0
    if(present(bs_mom_std))bs_mom_std=0.d0
    if(present(bs_A))bs_A=0.d0
    if(present(bs_M))bs_M=0.d0
    if(present(bs_P))bs_P=0.d0
    if(present(bs_TIE_param))bs_TIE_param=0.d0
    if(present(bs_TIE_chi2))bs_TIE_chi2=0.d0
    have_mom1=.true.

    ! Set flags.
    need_integrals=present(bs_tail_mom).or.present(bs_mom).or.&
       &present(bs_centre_mom)
    need_fit=need_integrals.or.present(bs_param).or.present(bs_chi2)
    need_TIE=present(bs_TIE_param).or.present(bs_TIE_chi2)

    ! Loop over bootstrap samples.
    do ibs=1,bs_setup%nbs_proc

      ! Generate bootstrap index vector.  NB, the bootstrapped dataset is
      ! unweighted regardless of original dataset being weighted or not.
      call gen_bootstrap_resample(dataset,indx_work)

      ! Compute standard estimators.
      if(present(bs_mom_std))call characterize_dist_indx(dataset%M,dataset%A,&
         &indx_work,mean=bs_mom_std(ibs,1),var=bs_mom_std(ibs,2))

      ! Compute things depending solely on the resample.
      if(present(bs_A).or.present(bs_M).or.present(bs_P))then
        do ipoint=1,plot_grid%npoint
          call get_point_M_A(dataset,indx_work,plot_grid,ipoint,t1,t2)
          if(present(bs_A))bs_A(ibs,ipoint)=t2
          if(present(bs_M))bs_M(ibs,ipoint)=t1
          if(present(bs_P))bs_P(ibs,ipoint)=eval_KDE_P(dataset,indx_work,t2)
        enddo ! imlogq
      endif ! present(bs_A).or.present(bs_M).or.present(bs_P)

      if(need_TIE)then
        ! Perform tail-index estimation.
        do ipoint=1,plot_grid%npoint
          call get_point_iabs_M_itail(dataset,indx_work,tail_form,plot_grid,&
             &ipoint,j,itail)
          call perform_qq_fit(dataset,indx_work,itail==2,tail_form%A_centre,&
             &j,t1,t2,t3)
          if(present(bs_TIE_param))bs_TIE_param(ibs,1:2,ipoint)=&
             &(/1.d0+1.d0/t1,exp(t2)/)
          if(present(bs_TIE_chi2))bs_TIE_chi2(ibs,ipoint)=t3
        enddo ! ipoint
      endif ! need_TIE

      ! Loop over anchors.
      do ianchor=1,assess_grid%nanchor

        ! Decide which tails are heavy at this ianchor.
        is_htail=tail_form%mu>0.d0.and.&
           &len_trim(assess_grid%def_anchor(:,ianchor))>0

        ! Determine anchors for this resample.
        call get_anchor(dataset%M,dataset%A,indx_work,is_htail,&
           &tail_form%A_centre,assess_grid%def_anchor(1,ianchor),&
           &assess_grid%anchor(1,ianchor),anchor_M,anchor_A)
        if(present(bs_anchor_A))bs_anchor_A(ibs,1:2,ianchor)=anchor_A(1:2)
        if(present(bs_anchor_M))bs_anchor_M(ibs,1:2,ianchor)=anchor_M(1:2)
        if(any((anchor_M<=0.5d0.or.are_equal(anchor_A,0.d0)).and.is_htail))&
           &cycle

        ! Rest is for things depending on the tail fits - cycle if not needed.
        if(.not.need_fit)cycle

        ! Compute regularization factors.
        regx=1.d0
        regy=1.d0
        do itail=1,2
          if(.not.is_htail(itail))cycle
          regx(itail)=anchor_A(itail)**tail_form%global_delta_mu(itail)
          regy(itail)=(dataset%M/(anchor_M(itail)-0.5d0))*&
             &anchor_A(itail)**(1.d0-tail_form%mu(itail))
        enddo ! itail

        ! Generate arrays required for fit at all expansion orders.
        nparam_dim(1)=maxval(assess_grid%nparam(1,:))
        nparam_dim(2)=maxval(assess_grid%nparam(2,:))
        call instantiate_delta_mu(tail_form,nparam_dim,delta_mu)
        do itail=1,2
          if(.not.is_htail(itail))cycle
          call prepare_fit_data(dataset%M,dataset%A,indx_work,&
             &anchor_M(itail),itail==2,tail_form%mu(itail),nparam_dim(itail),&
             &delta_mu(1,itail),tail_form%A_centre,tail_form%fit_weight_type,&
             &regx(itail),regy(itail),xx_expval(1,itail),xy_expval(1,itail),&
             &yy_expval(itail))
        enddo ! itail

        ! Loop over relevant expansion orders.
        do inparam=1,assess_grid%nnparam

          ! Get nparam and adjust for is_htail.  nparam>0 is used within
          ! this loop to decide if each tail is heavy or not.
          nparam=assess_grid%nparam(:,inparam)
          where(.not.is_htail)nparam=0

          ! Perform fit.
          call instantiate_constraint_set(tail_form,nparam,nparam_dim,&
             &delta_mu,regx,regy,cset)
          call perform_fit(nparam,nparam_dim,xx_expval,xy_expval,yy_expval,&
             &cset,param,chi2)
          call unregularize_param(nparam,nparam_dim,delta_mu,regx,regy,param,&
             &chi2)
          if(any(param>sqrt(huge(1.d0))))cycle
          if(present(bs_param))bs_param(ibs,:,:,inparam,ianchor)=param
          if(present(bs_chi2))bs_chi2(ibs,1:2,inparam,ianchor)=chi2

          if(need_integrals)then
            ! Compute central integrals.
            M_work=dataset%M-sum(anchor_M,nparam>0)
            cmom0=dble(M_work)/dble(dataset%M)
            call characterize_dist_indx(M_work,dataset%A,&
               &indx_work(anchor_M(1)+1),mean=cmom1,var=cmom2)
            ! Adjust for norm.
            cmom1=cmom1*cmom0
            cmom2=cmom2*(dble(M_work-1)/dble(dataset%M-1))
            ! Recentre central variance around zero, which makes all
            ! contributions additive.
            cmom2=cmom2+(dble(dataset%M)/dble(dataset%M-1))*(cmom1*cmom1)/cmom0
            ! Obtain tail contribution to norm, mean, and variance.
            call integrate_fit(nparam,nparam_dim,tail_form%mu,delta_mu,&
               &tail_form%momdef,param,anchor_A,tail_form%A_centre,mom0,mom1,&
               &mom2,mom1def)
            have_mom1=have_mom1.and.mom1def
            ! Compute total moments, recentring total variance around total
            ! mean.
            tmom0=cmom0+sum(mom0)
            tmom1=cmom1+sum(mom1)
            tmom2=cmom2+sum(mom2)
            ! Recentre variances around total mean.
            cmom2=cmom2-2.d0*cmom1*tmom1+cmom0*tmom1*tmom1
            if(tail_form%momdef>=MOMDEF_MOM2)then
              mom2=mom2-2.d0*mom1*tmom1+mom0*tmom1*tmom1
              tmom2=tmom2-tmom1*tmom1*(2.d0-tmom0)
            endif
            ! Store tail contributions to moments.
            if(present(bs_tail_mom))then
              bs_tail_mom(ibs,0,1:2,inparam,ianchor)=mom0
              bs_tail_mom(ibs,1,1:2,inparam,ianchor)=mom1
              bs_tail_mom(ibs,2,1:2,inparam,ianchor)=mom2
            endif
            ! Store central contribution to moments.
            if(present(bs_centre_mom))then
              bs_centre_mom(ibs,0,ianchor)=cmom0
              bs_centre_mom(ibs,1,ianchor)=cmom1
              bs_centre_mom(ibs,2,ianchor)=cmom2
            endif
            ! Store total moments.
            if(present(bs_mom))then
              bs_mom(ibs,0,inparam,ianchor)=tmom0
              bs_mom(ibs,1,inparam,ianchor)=tmom1
              bs_mom(ibs,2,inparam,ianchor)=tmom2
            endif
          endif ! need_integrals

        enddo ! inparam

      enddo ! ianchor

    enddo ! ibs

    ! Gather+broadcast data over MPI - return if no MPI.
    if(MPI_NPROC<=1)return

    ! Global have_mom1 flag.
    mom1def=have_mom1
    call mpi_reduce(have_mom1,mom1def,1,mpi_logical,mpi_land,0,mpi_comm_world,&
       &ierror)

    ! bs_param array; indexed (iparam,itail,inparam,ianchor).
    if(present(bs_param))then
      do ianchor=1,assess_grid%nanchor
        do inparam=1,assess_grid%nnparam
          do itail=1,2
            do iparam=1,assess_grid%nparam(itail,inparam)
              if(MPI_MASTER)then
                call mpi_gather(mpi_in_place,&
                   &bs_setup%nbs_proc,mpi_double_precision,&
                   &bs_param(1,iparam,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
                   &ierror)
              else
                call mpi_gather(bs_param(1,iparam,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,&
                   &bs_param(1,iparam,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
                   &ierror)
              endif ! MPI_MASTER or not
            enddo ! iparam
          enddo ! itail
        enddo ! inparam
      enddo ! ianchor
    endif ! present(bs_param)

    ! bs_tail_mom array; indexed (imom,itail,inparam,ianchor)
    if(present(bs_tail_mom))then
      do ianchor=1,assess_grid%nanchor
        do inparam=1,assess_grid%nnparam
          do imom=0,2
            do itail=1,2
              if(MPI_MASTER)then
                call mpi_gather(mpi_in_place,&
                   &bs_setup%nbs_proc,mpi_double_precision,&
                   &bs_tail_mom(1,imom,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
                   &ierror)
              else
                call mpi_gather(bs_tail_mom(1,imom,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,&
                   &bs_tail_mom(1,imom,itail,inparam,ianchor),&
                   &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
                   &ierror)
              endif ! MPI_MASTER or not
            enddo ! itail
          enddo ! imom
        enddo ! inparam
      enddo ! ianchor
    endif ! present(bs_tail_mom)

    ! bs_centre_mom array; indexed (imom,ianchor).
    if(present(bs_centre_mom))then
      do ianchor=1,assess_grid%nanchor
        do imom=0,2
          if(MPI_MASTER)then
            call mpi_gather(mpi_in_place,bs_setup%nbs_proc,&
               &mpi_double_precision,bs_centre_mom(1,imom,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          else
            call mpi_gather(bs_centre_mom(1,imom,ianchor),bs_setup%nbs_proc,&
               &mpi_double_precision,bs_centre_mom(1,imom,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          endif ! MPI_MASTER or not
        enddo ! imom
      enddo ! ianchor
    endif ! present(bs_centre_mom)

    ! bs_mom array; indexed (imom,inparam,ianchor)
    if(present(bs_mom))then
      do ianchor=1,assess_grid%nanchor
        do inparam=1,assess_grid%nnparam
          do imom=0,2
            if(MPI_MASTER)then
              call mpi_gather(mpi_in_place,&
                 &bs_setup%nbs_proc,mpi_double_precision,&
                 &bs_mom(1,imom,inparam,ianchor),bs_setup%nbs_proc,&
                 &mpi_double_precision,0,mpi_comm_world,ierror)
            else
              call mpi_gather(bs_mom(1,imom,inparam,ianchor),&
                 &bs_setup%nbs_proc,mpi_double_precision,&
                 &bs_mom(1,imom,inparam,ianchor),bs_setup%nbs_proc,&
                 &mpi_double_precision,0,mpi_comm_world,ierror)
            endif ! MPI_MASTER or not
          enddo ! imom
        enddo ! inparam
      enddo ! ianchor
    endif ! present(bs_mom)

    ! bs_chi2 array; indexed (itail,inparam)
    if(present(bs_chi2))then
      do ianchor=1,assess_grid%nanchor
        do inparam=1,assess_grid%nnparam
          do itail=1,2
            if(MPI_MASTER)then
              call mpi_gather(mpi_in_place,&
                 &bs_setup%nbs_proc,mpi_double_precision,&
                 &bs_chi2(1,itail,inparam,ianchor),bs_setup%nbs_proc,&
                 &mpi_double_precision,0,mpi_comm_world,ierror)
            else
              call mpi_gather(bs_chi2(1,itail,inparam,ianchor),&
                 &bs_setup%nbs_proc,mpi_double_precision,&
                 &bs_chi2(1,itail,inparam,ianchor),bs_setup%nbs_proc,&
                 &mpi_double_precision,0,mpi_comm_world,ierror)
            endif ! MPI_MASTER or not
          enddo ! itail
        enddo ! inparam
      enddo ! ianchor
    endif ! present(bs_chi2)

    ! bs_mom_std array; indexed (imom[=1:2]).
    if(present(bs_mom_std))then
      do imom=1,2
        if(MPI_MASTER)then
          call mpi_gather(mpi_in_place,bs_setup%nbs_proc,&
             &mpi_double_precision,bs_mom_std(1,imom),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
        else
          call mpi_gather(bs_mom_std(1,imom),bs_setup%nbs_proc,&
             &mpi_double_precision,bs_mom_std(1,imom),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
        endif ! MPI_MASTER or not
      enddo ! imom
    endif ! present(bs_mom_std)

    ! bs_A, bs_M, and bs_P arrays; indexed (ipoint).
    if(present(bs_A).or.present(bs_P))then
      do ipoint=1,plot_grid%npoint
        if(MPI_MASTER)then
          if(present(bs_A))call mpi_gather(mpi_in_place,&
             &bs_setup%nbs_proc,mpi_double_precision,bs_A(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
          if(present(bs_M))call mpi_gather(mpi_in_place,&
             &bs_setup%nbs_proc,mpi_double_precision,bs_M(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
          if(present(bs_P))call mpi_gather(mpi_in_place,&
             &bs_setup%nbs_proc,mpi_double_precision,bs_P(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
        else
          if(present(bs_A))call mpi_gather(bs_A(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,bs_A(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
          if(present(bs_M))call mpi_gather(bs_M(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,bs_M(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
          if(present(bs_P))call mpi_gather(bs_P(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,bs_P(1,ipoint),&
             &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
             &ierror)
        endif
      enddo ! ipoint
    endif ! present(bs_A).or.present(bs_M).or.present(bs_P)

    ! bs_anchor_A and bs_anchor_M arrays; indexed (itail,ianchor).
    if(present(bs_anchor_A).or.present(bs_anchor_M))then
      do ianchor=1,assess_grid%nanchor
        do itail=1,2
          if(MPI_MASTER)then
            if(present(bs_anchor_A))call mpi_gather(&
               &mpi_in_place,bs_setup%nbs_proc,&
               &mpi_double_precision,bs_anchor_A(1,itail,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
            if(present(bs_anchor_M))call mpi_gather(&
               &mpi_in_place,bs_setup%nbs_proc,&
               &mpi_double_precision,bs_anchor_M(1,itail,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          else
            if(present(bs_anchor_A))call mpi_gather(&
               &bs_anchor_A(1,itail,ianchor),bs_setup%nbs_proc,&
               &mpi_double_precision,bs_anchor_A(1,itail,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
            if(present(bs_anchor_M))call mpi_gather(&
               &bs_anchor_M(1,itail,ianchor),bs_setup%nbs_proc,&
               &mpi_double_precision,bs_anchor_M(1,itail,ianchor),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          endif
        enddo ! itail
      enddo ! ianchor
    endif ! present(bs_anchor_A).or.present(bs_anchor_M)

    ! bs_TIE_param array; indexed (iparam[=1:2],ipoint).
    if(present(bs_TIE_param))then
      do ipoint=1,plot_grid%npoint
        do iparam=1,2
          if(MPI_MASTER)then
            call mpi_gather(mpi_in_place,bs_setup%nbs_proc,&
               &mpi_double_precision,bs_TIE_param(1,iparam,ipoint),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          else
            call mpi_gather(bs_TIE_param(1,iparam,ipoint),bs_setup%nbs_proc,&
               &mpi_double_precision,bs_TIE_param(1,iparam,ipoint),&
               &bs_setup%nbs_proc,mpi_double_precision,0,mpi_comm_world,&
               &ierror)
          endif
        enddo ! iparam
      enddo ! ipoint
    endif ! present(bs_TIE_param)

    ! bs_TIE_chi2 array; indexed (ipoint).
    if(present(bs_TIE_chi2))then
      do ipoint=1,plot_grid%npoint
        if(MPI_MASTER)then
          call mpi_gather(mpi_in_place,bs_setup%nbs_proc,&
             &mpi_double_precision,bs_TIE_chi2(1,ipoint),bs_setup%nbs_proc,&
             &mpi_double_precision,0,mpi_comm_world,ierror)
        else
          call mpi_gather(bs_TIE_chi2(1,ipoint),bs_setup%nbs_proc,&
             &mpi_double_precision,bs_TIE_chi2(1,ipoint),bs_setup%nbs_proc,&
             &mpi_double_precision,0,mpi_comm_world,ierror)
        endif
      enddo ! ipoint
    endif ! present(bs_TIE_chi2)

    ! Final broadcast from master process.
    call mpi_bcast(have_mom1,1,mpi_logical,0,mpi_comm_world,ierror)
    if(present(bs_param))call mpi_bcast(bs_param,size(bs_param),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_tail_mom))call mpi_bcast(bs_tail_mom,size(bs_tail_mom),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_centre_mom))call mpi_bcast(bs_centre_mom,&
       &size(bs_centre_mom),mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_mom))call mpi_bcast(bs_mom,size(bs_mom),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_chi2))call mpi_bcast(bs_chi2,size(bs_chi2),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_A))call mpi_bcast(bs_A,size(bs_A),mpi_double_precision,0,&
       &mpi_comm_world,ierror)
    if(present(bs_M))call mpi_bcast(bs_M,size(bs_M),mpi_double_precision,0,&
       &mpi_comm_world,ierror)
    if(present(bs_P))call mpi_bcast(bs_P,size(bs_P),mpi_double_precision,0,&
       &mpi_comm_world,ierror)
    if(present(bs_anchor_A))call mpi_bcast(bs_anchor_A,size(bs_anchor_A),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_anchor_M))call mpi_bcast(bs_anchor_M,size(bs_anchor_M),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_TIE_param))call mpi_bcast(bs_TIE_param,size(bs_TIE_param),&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    if(present(bs_TIE_chi2))call mpi_bcast(bs_TIE_chi2,size(bs_TIE_chi2),&
       &mpi_double_precision,0,mpi_comm_world,ierror)

  END SUBROUTINE bootstrap_fit


  ! QQ-FITTING.


  SUBROUTINE perform_qq_fit(dataset,indx,reverse,A_centre,k,afit,bfit,mse)
    !---------------------------------------------------------!
    ! Perform least-squares fit of qq-plot of dataset%A up to !
    ! the kth order statistic to a straight line, and return  !
    ! slope afit, offset bfit and mean-squared error mse.     !
    !---------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    INTEGER, INTENT(in) :: indx(dataset%M), k
    LOGICAL, INTENT(in) :: reverse
    DOUBLE PRECISION, INTENT(in) :: A_centre
    DOUBLE PRECISION, INTENT(inout) :: afit,bfit,mse
    ! Local variables.
    INTEGER i,ii
    DOUBLE PRECISION inv_M,x_k1,x,y,w
    DOUBLE PRECISION sum_w,sum_wx,sum_wy,sum_wx2,sum_wy2,sum_wxy

    ! Initialize.
    inv_M=1.d0/dble(dataset%M)
    x_k1=-log((dble(k)+0.5d0)*inv_M)

    ! Evaluate accumulators.
    sum_w=0.d0
    sum_wx=0.d0
    sum_wy=0.d0
    sum_wx2=0.d0
    sum_wy2=0.d0
    sum_wxy=0.d0
    do i=1,k
      ! Compute x,y,w at this point.
      if(.not.reverse)then
        ii=indx(i)
      else
        ii=indx(dataset%M+1-i)
      endif
      x=-log((dble(i)-0.5d0)*inv_M)
      y=log(abs(dataset%A(ii)-A_centre))
      w=1.d0/(x-x_k1)
      ! Accumulate.
      sum_w=sum_w+w
      sum_wx=sum_wx+w*x
      sum_wy=sum_wy+w*y
      sum_wx2=sum_wx2+w*x*x
      sum_wy2=sum_wy2+w*y*y
      sum_wxy=sum_wxy+w*x*y
    enddo ! i

    ! Perform least-squares fit.
    afit=(sum_wxy*sum_w-sum_wx*sum_wy)/(sum_wx2*sum_w-sum_wx*sum_wx)
    bfit=(sum_wy-afit*sum_wx)/sum_w
    mse=(afit*afit*sum_wx2+bfit*bfit*sum_w+sum_wy2+&
       &2.d0*afit*bfit*sum_wx-2.d0*afit*sum_wxy-&
       &2.d0*bfit*sum_wy)/sum_w

  END SUBROUTINE perform_qq_fit


  ! MODEL P PARSERS AND GENERATION ROUTINES


  SUBROUTINE generate_model_P(model_P,nsample,dataset)
    !--------------------------------------------------------!
    ! Generate M random numbers distributed according to the !
    ! model PDF defined by model_P by direct sampling.       !
    !--------------------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    INTEGER,INTENT(in) :: nsample
    TYPE(dataset_type),INTENT(inout) :: dataset
    ! Local variables.
    INTEGER i,ierr
    DOUBLE PRECISION xi,t1,t2
    ! Parallelization variables.
    INTEGER ierror,iproc,Mtot,Mthis
    INTEGER Mproc(0:MPI_NPROC-1),Moffset(0:MPI_NPROC-1)

    ! Set up dataset.
    dataset%M=nsample
    dataset%have_w=.false.
    dataset%description=trim(i2s(dataset%M))//' '//&
       &trim(type_string(dataset%have_w))//'-type data from &
       &model with asymptote P(A) ~ '//trim(model_asymptote(model_P))
    call get_model_P_mean_var(model_P,t1,t2,ierr)
    if(allocated(dataset%A))deallocate(dataset%A)
    if(allocated(dataset%w))deallocate(dataset%w)
    allocate(dataset%A(dataset%M))

    ! Report.
    if(MPI_MASTER)then
      write(6,'(a)')'Generating sample drawn from model distribution'
      write(6,'(a)')'==============================================='
      if(ierr>=2)then
        write(6,'(a)')'Mean              : (undefined)'
        write(6,'(a)')'Variance          : (undefined)'
      elseif(ierr==1)then
        write(6,'(a,es20.12)')'Mean              : ',t1
        write(6,'(a)')'Variance          : (undefined)'
      else
        write(6,'(a,es20.12)')'Mean              : ',t1
        write(6,'(a,es20.12)')'Variance          : ',t2
      endif
      write(6,'(a)')'Asymptote         : '//trim(model_asymptote(model_P))
      write(6,'(a)')'Number of samples : '//trim(i2s(nsample))
      write(6,'()')
    endif

    ! Compute number of steps per process.
    Mtot=0
    do iproc=0,MPI_NPROC-1
      Mthis=nsample/MPI_NPROC
      if(iproc<mod(nsample,MPI_NPROC))Mthis=Mthis+1
      Mproc(iproc)=Mthis
      Moffset(iproc)=Mtot
      Mtot=Mtot+Mthis
    enddo ! iproc

    ! Generate data.
    do i=1,Mproc(MPI_IPROC)
      call random_number(xi)
      dataset%A(i)=invert_model_P_CDF(model_P,xi)
    enddo ! i

    ! Gather data on master.
    if(MPI_MASTER)then
      call mpi_gatherv(mpi_in_place,Mproc(MPI_IPROC),mpi_double_precision,&
         &dataset%A,Mproc,Moffset,mpi_double_precision,0,MPI_COMM_WORLD,ierror)
    else ! .not.MPI_MASTER
      call mpi_gatherv(dataset%A,Mproc(MPI_IPROC),mpi_double_precision,&
         &dataset%A,Mproc,Moffset,mpi_double_precision,0,mpi_comm_world,ierror)
    endif ! MPI_MASTER or not
    ! Broadcast.
    call mpi_bcast(dataset%A,dataset%M,mpi_double_precision,0,mpi_comm_world,&
       &ierror)

  END SUBROUTINE generate_model_P


  SUBROUTINE get_model_P_mean_var(model_P,mean,var,ierr)
    !-----------------------------------------------!
    ! Compute the mean and variance of the provided !
    ! model distribution.                           !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    DOUBLE PRECISION,INTENT(inout) :: mean,var
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION t0,t1,t2
    INTEGER iterm
    ierr=0
    mean=0.d0
    var=0.d0
    do iterm=1,model_P%nterm
      t1=0.d0
      t2=0.d0
      if(model_P%muvec(iterm)<=0.d0)then
        t1=model_P%centrevec(iterm)
        t2=0.5d0
      elseif(model_P%muvec(iterm)>3.d0)then
        t0=pi/model_P%muvec(iterm)
        t1=model_P%centrevec(iterm)
        t2=sin(t0)/sin(3.d0*t0)
      elseif(model_P%muvec(iterm)>2.d0)then
        ierr=max(ierr,1)
        t1=model_P%centrevec(iterm)
      else
        ierr=max(ierr,2)
      endif
      mean=mean+model_P%avec(iterm)*t1
      var=var+model_P%avec(iterm)*(model_P%centrevec(iterm)**2+&
         &model_P%lambdavec(iterm)**2*t2)
    enddo ! iterm
    var=var-mean**2
  END SUBROUTINE get_model_P_mean_var


  DOUBLE PRECISION FUNCTION invert_model_P_CDF(model_P,cdf)
    !-----------------------------------------------------------!
    ! Find the value of A such that the cumulative distribution !
    ! function of the model PDF defined in model_P equals cdf.  !
    ! This is done by bisection since the inversion is not      !
    ! analytic.                                                 !
    !-----------------------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    DOUBLE PRECISION,INTENT(in) :: cdf
    ! Local variables.
    DOUBLE PRECISION AL,AR,AC,cL,cR,cC
    ! Tolerances.
    DOUBLE PRECISION, PARAMETER :: HUGE_TOL=sqrt(huge(1.d0))
    DOUBLE PRECISION, PARAMETER :: CONV_TOL=1.d-8

    ! Get initial bracket.
    AL=0.d0
    cL=model_P_CDF(model_P,AL)
    if(cdf<cL)then
      AR=AL
      cR=cL
      AL=-1.d0
      do
        cL=model_P_CDF(model_P,AL)
        if(cdf<cL.and.-HUGE_TOL<AL)then
          AL=2.d0*AL
        elseif(cL<cdf)then
          exit
        else
          invert_model_P_CDF=AL
          return
        endif
      enddo
    elseif(cL<cdf)then
      AR=1.d0
      do
        cR=model_P_CDF(model_P,AR)
        if(cR<cdf.and.AR<HUGE_TOL)then
          AR=2.d0*AR
        elseif(cdf<cR)then
          exit
        else
          invert_model_P_CDF=AR
          return
        endif
      enddo
    else
      invert_model_P_CDF=AL
      return
    endif

    ! Perform bisection.
    do
      AC=0.5d0*(AL+AR)
      if(AR-AL<CONV_TOL*min(abs(AR),abs(AL)))exit
      cC=model_P_CDF(model_P,AC)
      if(cC<cdf)then
        AL=AC
        cL=cC
      elseif(cdf<cC)then
        AR=AC
        cR=cC
      else
        exit
      endif
    enddo

    ! Set output value.
    invert_model_P_CDF=AC

  END FUNCTION invert_model_P_CDF


  DOUBLE PRECISION FUNCTION model_P_CDF(model_P,A)
    !--------------------------------------------------!
    ! Evaluate the cumulative distribution function of !
    ! the model distribution defined in model_P at A.  !
    !--------------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    DOUBLE PRECISION,INTENT(in) :: A
    ! Local variables.
    INTEGER iterm
    DOUBLE PRECISION x,absx,absx_mu,xarg,t1

    ! Initialize.
    model_P_CDF=0.d0

    ! Loop over terms.
    do iterm=1,model_P%nterm
      ! Get reduced variable and its absolute value.
      x=(A-model_P%centrevec(iterm))/model_P%lambdavec(iterm)
      absx=abs(x)
      if(model_P%muvec(iterm)<=0.d0)then ! gauss(lambda;A).
        ! CDF difference between zero and |x|.
        t1=half_sqrt_pi*erf(absx)
        ! Prefactors.
        t1=t1*model_P%lambdavec(iterm)*model_P%normvec(iterm)
        ! Opposite sign for negative x.
        if(x<0.d0)t1=-t1
        t1=0.5d0+t1
        ! Add contribution.
        model_P_CDF=model_P_CDF+model_P%avec(iterm)*t1
      else ! h(mu,lambda;A).
        absx_mu=absx**model_P%muvec(iterm)
        ! Decide if it is more efficient to compute the CDF of h(A) from the
        ! tail end or from the centre.  Tests yield turnover points of
        ! absx=1.060,1.100,1.135,1.160,1.175,1.210,1.225,1.250 at
        ! mu=2.5,3,3.5,4,4.5,5,5.5,6, which is approximately 1+(mu-1)*0.05.
        if(absx>1.d0+(model_P%muvec(iterm)-1.d0)*0.05d0)then
          xarg=1.d0/(1.d0+absx_mu)
          ! CDF at -|x|.
          t1=(absx/((model_P%muvec(iterm)-1.d0)*(1.d0+absx_mu)))*&
             &hyper2f1_11(2.d0-1.d0/model_P%muvec(iterm),xarg)
          t1=t1*model_P%lambdavec(iterm)*model_P%normvec(iterm)
          ! Reverse for positive x.
          if(x>0.d0)t1=1.d0-t1
        else
          xarg=absx_mu/(1.d0+absx_mu)
          ! CDF difference between zero and |x|.
          t1=(absx/(absx_mu+1.d0))*&
             &hyper2f1_11(1.d0+1.d0/model_P%muvec(iterm),xarg)
          t1=t1*model_P%lambdavec(iterm)*model_P%normvec(iterm)
          ! Change sign for negative x.
          if(x<0.d0)t1=-t1
          ! Add CDF at zero.
          t1=0.5d0+t1
        endif
        ! Add contribution.
        model_P_CDF=model_P_CDF+model_P%avec(iterm)*t1
      endif
    enddo ! iterm

  END FUNCTION model_P_CDF


  SUBROUTINE parse_model_P_token(token,iterm,model_P,ierr)
    !-------------------------------------------!
    ! Parse a token of one of the forms:        !
    ! - [a*]h(mu[,lambda[,A0]])                 !
    ! - [a*]gauss([lambda[,A0]])                !
    ! and return the parameters as the iterm-th !
    ! entry in model_P.                         !
    !-------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: token
    INTEGER,INTENT(in) :: iterm
    TYPE(model_P_type),INTENT(inout) :: model_P
    INTEGER,INTENT(inout) :: ierr
    CHARACTER(len(token)) remainder
    INTEGER ipos

    ! Initialize.
    ierr=0
    remainder=token

    ! Parse a.
    ipos=scan(remainder,'*')
    if(ipos<1)then
      model_P%avec(iterm)=1.d0
    else
      read(remainder(1:ipos-1),*,iostat=ierr)model_P%avec(iterm)
      if(ierr/=0)return
      if(model_P%avec(iterm)<=0.d0)then
        ierr=1
        return
      endif
      remainder=remainder(ipos+1:)
    endif

    ! Parse function.
    ipos=scan(remainder,'(')
    if(ipos<1)then
      ierr=1
      return
    endif
    select case(remainder(1:ipos-1))
    case('gauss')
      ! Initialize default Gaussian.
      model_P%muvec(iterm)=0.d0 ! flags use of Gaussian
      model_P%lambdavec(iterm)=1.d0
      model_P%centrevec(iterm)=0.d0
    case('h')
      ! Parse mu, which is a mandatory argument.
      remainder=remainder(ipos+1:)
      ipos=scan(remainder,',')
      if(ipos<1)then
        ipos=scan(remainder,')')
        if(ipos<1)then
          ierr=1
          return
        endif
      endif
      read(remainder(1:ipos-1),*,iostat=ierr)model_P%muvec(iterm)
      if(ierr/=0)return
      if(model_P%muvec(iterm)<=1.d0)then
        ierr=2
        return
      endif
      ! Initialize default model function.
      model_P%lambdavec(iterm)=1.d0
      model_P%centrevec(iterm)=0.d0
    case default
      ierr=1
      return
    end select

    ! Parse lambda if present.
    remainder=remainder(ipos+1:)
    if(len_trim(remainder)<1)return
    ipos=scan(remainder,',')
    if(ipos<1)then
      ipos=scan(remainder,')')
      if(ipos<1)then
        ierr=1
        return
      endif
    endif
    read(remainder(1:ipos-1),*,iostat=ierr)model_P%lambdavec(iterm)
    if(ierr/=0)return
    if(model_P%lambdavec(iterm)<=0.d0)then
      ierr=1
      return
    endif

    ! Parse A_centre if present.
    remainder=remainder(ipos+1:)
    if(len_trim(remainder)<1)return
    ipos=scan(remainder,')')
    if(ipos<1)then
      ierr=1
      return
    endif
    read(remainder(1:ipos-1),*,iostat=ierr)model_P%centrevec(iterm)
    if(ierr/=0)return

    ! Ensure there is no trailing junk.
    remainder=remainder(ipos+1:)
    if(len_trim(remainder)>0)then
      ierr=1
      return
    endif

  END SUBROUTINE parse_model_P_token


  CHARACTER(16384) FUNCTION model_asymptote(model_P)
    !-----------------------------------------------!
    ! Return a string specifying the asymptote of a !
    ! model distribution (assuming all centrevec(:) !
    ! are the same; this routine does not take into !
    ! account the corresponding change in the       !
    ! beyond-leading-order asymptotic coefficients. !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(in) :: model_P
    ! Local variables.
    CHARACTER(80) char1,char2
    LOGICAL isfirst
    INTEGER iterm
    DOUBLE PRECISION t1
    model_asymptote=''
    isfirst=.true.
    do iterm=1,model_P%nterm
      if(model_P%muvec(iterm)<=0.d0)cycle
      t1=model_P%avec(iterm)*model_P%normvec(iterm)*&
         &model_P%lambdavec(iterm)**model_P%muvec(iterm)
      if(.not.isfirst)model_asymptote=trim(model_asymptote)//' +'
      char1=dble2human(t1)
      char2=dble2human(model_P%muvec(iterm))
      model_asymptote=trim(model_asymptote)//' '//trim(char1)//' * A^-'//&
         &trim(char2)
      isfirst=.false.
    enddo ! iterm
    if(isfirst)model_asymptote='0'
    model_asymptote=adjustl(model_asymptote)
  END FUNCTION model_asymptote


  SUBROUTINE normalize_model_P(model_P)
    !---------------------------------------------!
    ! Rescale terms in model distribution so that !
    ! they become normalized prefactors.          !
    !---------------------------------------------!
    IMPLICIT NONE
    TYPE(model_P_type),INTENT(inout) :: model_P
    INTEGER iterm
    DOUBLE PRECISION t0,t1

    ! Compute individual term norms.
    do iterm=1,model_P%nterm
      if(model_P%muvec(iterm)>0.d0)then ! h
        t0=pi/model_P%muvec(iterm)
        t1=sin(t0)/(2.d0*t0*model_P%lambdavec(iterm))
      else ! gauss
        t1=inv_sqrt_pi/model_P%lambdavec(iterm)
      endif ! h or gauss
      model_P%normvec(iterm)=t1
    enddo ! iterm

    ! Normalize prefactors.
    model_P%avec=model_P%avec/sum(model_P%avec)

  END SUBROUTINE normalize_model_P


  ! QUANTILE GENERATORS.


  SUBROUTINE gen_quantiles_w(M,lda,w,inv_total_w,q,indx)
    !--------------------------------------------------!
    ! Evaluate the quantile function for weighted data !
    ! given the ascending sort indices INDX(:).        !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M,lda,indx(M)
    DOUBLE PRECISION, INTENT(in) :: w(lda),inv_total_w
    DOUBLE PRECISION, INTENT(inout) :: q(M)
    INTEGER k,kk
    DOUBLE PRECISION cum_w,w_prev
    kk=indx(1)
    cum_w=0.5d0*w(kk)
    q(1)=cum_w*inv_total_w
    w_prev=w(kk)
    do k=2,M
      kk=indx(k)
      cum_w=cum_w+0.5d0*(w(kk)+w_prev)
      q(k)=cum_w*inv_total_w
      w_prev=w(kk)
    enddo ! k
  END SUBROUTINE gen_quantiles_w


  ! TAIL FITTING.


  SUBROUTINE prepare_fit_data(M,A,indx,Mfit,reverse,mu,nparam,delta_mu,&
     &A_centre,fit_weight_type,regx,regy,xx_expval,xy_expval,yy_expval)
    !-----------------------------------------------------!
    ! Construct the x/y/w vectors containing the data and !
    ! least-squares weights to be fitted to a polynomial, !
    ! and arrays containing unique entries of <x*x> and   !
    ! <x*y> expectation values for later use in fits.     !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M,indx(M),Mfit,nparam
    LOGICAL, INTENT(in) :: reverse
    CHARACTER(*), INTENT(in) :: fit_weight_type
    DOUBLE PRECISION, INTENT(in) :: A(M),mu,delta_mu(nparam),A_centre,&
       &regx,regy
    DOUBLE PRECISION, INTENT(inout) :: xx_expval((nparam*(nparam+1))/2),&
       &xy_expval(nparam),yy_expval
    ! Local variables.
    INTEGER iparam,jparam,k,ijindx
    DOUBLE PRECISION log_qMfit1,sum_w,inv_M,Ak,q,x,y,w,xpow(nparam)

    ! Initialize.
    inv_M=1.d0/dble(M)
    log_qMfit1=log((dble(Mfit)+0.5d0)*inv_M)

    ! Construct fit data vectors.
    xy_expval=0.d0
    xx_expval=0.d0
    yy_expval=0.d0
    sum_w=0.d0
    do k=1,Mfit
      q=(dble(k)-0.5d0)*inv_M
      if(.not.reverse)then
        Ak=A(indx(k))
      else
        Ak=A(indx(M+1-k))
      endif
      x=regx/abs(Ak-A_centre)
      y=regy*q*abs(Ak-A_centre)**(mu-1.d0)
      w=x**(mu-1.d0)
      select case(trim(fit_weight_type))
      case('hill')
        w=w/(log_qMfit1-log(q))
      end select
      ! Accumulate.
      sum_w=sum_w+w
      xpow(:)=x**delta_mu(:)
      xy_expval(:)=xy_expval(:)+w*y*xpow(:)
      yy_expval=yy_expval+w*y*y
      ijindx=0
      do iparam=1,nparam
        do jparam=iparam,nparam
          ijindx=ijindx+1
          xx_expval(ijindx)=xx_expval(ijindx)+w*xpow(iparam)*xpow(jparam)
        enddo ! jparam
      enddo ! iparam
    enddo ! k

    ! Renormalize by total weight.
    xy_expval=xy_expval/sum_w
    xx_expval=xx_expval/sum_w
    yy_expval=yy_expval/sum_w

  END SUBROUTINE prepare_fit_data


  SUBROUTINE perform_fit(nparam,nparam_dim,xx_expval,xy_expval,yy_expval,&
     &cset,param,chi2)
    !--------------------------------------------------------------!
    ! Perform least-squares fit to a polynomial from pre-processed !
    ! {xx,xy}_expval, and return the coefficients param.  This     !
    ! version simultaneously fits both tails constraining the      !
    ! requested parameters to be linear combinations of others.    !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nparam(2),nparam_dim(2)
    TYPE(constraint_set_type), INTENT(in) :: cset
    DOUBLE PRECISION, INTENT(in) :: &
       &xx_expval((maxval(nparam_dim)*(maxval(nparam_dim)+1))/2,2),&
       &xy_expval(maxval(nparam_dim),2),yy_expval(2)
    DOUBLE PRECISION, INTENT(inout) :: param(maxval(nparam_dim),2),chi2(2)
    ! Local variables.
    INTEGER ieq,irhs,iparam,itail,jp,jparam,jtail,k,kparam,ktail,l,lparam,&
      &ltail,i,ijindx,klindx,ierr
    DOUBLE PRECISION t1,t2
    ! Linear system.
    INTEGER np
    DOUBLE PRECISION,ALLOCATABLE :: Xmat(:,:),Xinv(:,:),pvec(:)
    INTEGER,ALLOCATABLE :: piv(:)

    ! Initialize.
    param=0.d0
    chi2=0.d0
    if(all(nparam==0))return

    ! See which parameters are free.
    np=count(cset%ndependant>0)
    if(np==0)return
    allocate(Xmat(np,np),Xinv(np,np),pvec(np),piv(np))
    Xmat=0.d0
    pvec=0.d0
    piv=0

    ! Initialize equation counter.
    ieq=0
    ! Loop over tails.
    do itail=1,2
      ! Loop over parameters in tail.
      do iparam=1,nparam(itail)
        ! Skip if this parameter is determined.
        if(cset%ndependant(iparam,itail)<1)cycle
        ieq=ieq+1
        ! Loop over parameters depending on this parameter (which includes
        ! itself with coefficient 1).
        do k=1,cset%ndependant(iparam,itail)
          ktail=cset%dependant_itail(k,iparam,itail)
          kparam=cset%dependant_iparam(k,iparam,itail)
          if(kparam>nparam(ktail))cycle
          t1=cset%constraint_coeff(kparam,ktail,iparam,itail)
          pvec(ieq)=pvec(ieq)+t1*xy_expval(kparam,ktail)
          ! Loop over all free parameters.
          jp=0
          do jtail=1,2
            do jparam=1,nparam(jtail)
              if(cset%ndependant(jparam,jtail)<1)cycle
              jp=jp+1
              do l=1,cset%ndependant(jparam,jtail)
                ltail=cset%dependant_itail(l,jparam,jtail)
                if(ltail/=ktail)cycle
                lparam=cset%dependant_iparam(l,jparam,jtail)
                if(lparam>nparam(ltail))cycle
                t2=cset%constraint_coeff(lparam,ltail,jparam,jtail)
                klindx=ltriang_index(nparam_dim(jtail),kparam,lparam)
                Xmat(jp,ieq)=Xmat(jp,ieq)+t1*t2*xx_expval(klindx,ktail)
              enddo ! l
            enddo ! jparam
          enddo ! jtail
          if(jp/=np)return
        enddo ! k
      enddo ! iparam
    enddo ! itail
    if(ieq/=np)return

    ! Flush.
    where(abs(Xmat)<epsilon(1.d0)**2*maxval(abs(Xmat)))Xmat=0.d0

    ! Invert Xmat.
    call dgetrf(np,np,Xmat,np,piv,ierr)
    if(ierr/=0)return
    Xinv=0.d0
    do i=1,np
      Xinv(i,i)=1.d0
    enddo ! i
    call dgetrs('N',np,np,Xmat,np,piv,Xinv,np,ierr)
    if(ierr/=0)return

    ! Evaluate fit coefficients.
    pvec=matmul(Xinv,pvec)

    ! Guard against very large parameter values.
    if(any(abs(pvec)>sqrt(huge(1.d0))))return

    ! Put parameters in output arrays.
    jp=0
    do itail=1,2
      ! Loop over parameters in tail.
      do iparam=1,nparam(itail)
        ! Skip if this parameter is determined.
        if(cset%ndependant(iparam,itail)<1)cycle
        jp=jp+1
        param(iparam,itail)=pvec(jp)
      enddo ! iparam
    enddo ! itail
    ! Reconstruct determined parameters.
    do itail=1,2
      ! Loop over parameters in tail.
      do iparam=1,nparam(itail)
        ! Skip if this parameter is free.
        if(cset%ndependant(iparam,itail)>0)cycle
        ! Loop over terms in linear combination.
        do irhs=1,cset%ncombine(iparam,itail)
          ktail=cset%combine_itail(irhs,iparam,itail)
          kparam=cset%combine_iparam(irhs,iparam,itail)
          if(kparam>nparam(ktail))cycle
          t1=cset%constraint_coeff(iparam,itail,kparam,ktail)
          param(iparam,itail)=param(iparam,itail)+t1*param(kparam,ktail)
        enddo
      enddo ! iparam
    enddo ! itail

    ! Evaluate chi2 for each tail.
    do itail=1,2
      t1=yy_expval(itail)
      do iparam=1,nparam(itail)
        do jparam=1,nparam(itail)
          ijindx=ltriang_index(nparam_dim(itail),iparam,jparam)
          t1=t1+param(iparam,itail)*param(jparam,itail)*&
             &xx_expval(ijindx,itail)
        enddo ! jparam
        t1=t1-2.d0*param(iparam,itail)*xy_expval(iparam,itail)
      enddo ! iparam
      chi2(itail)=t1
    enddo ! itail

  END SUBROUTINE perform_fit


  INTEGER FUNCTION ltriang_index(n,i,j)
    !---------------------------------------------------------------!
    ! Return the index of the (i,j)th element of an N x N symmetric !
    ! matrix in a lower-triangular packed representation.           !
    !---------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n,i,j
    INTEGER ii,jj
    ii=max(i,j)
    jj=min(i,j)
    ltriang_index=ii+((jj-1)*(2*n-jj))/2
  END FUNCTION ltriang_index


  ! KERNEL DENSITY ESTIMATION.


  DOUBLE PRECISION FUNCTION eval_KDE_P(dataset,indx,A)
    !------------------------------------------------!
    ! Evaluate the kernel density estimator of P(A). !
    !------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    INTEGER,INTENT(in) :: indx(dataset%M)
    DOUBLE PRECISION,INTENT(in) :: A
    INTEGER ibackforth,k,k0,k1,k2,kk,dk
    DOUBLE PRECISION x,t1,Ak,h
    DOUBLE PRECISION,PARAMETER :: norm=1.d0/sqrt(8.d0*atan(1.d0))

    ! Initialize.
    eval_KDE_P=0.d0

    ! Locate closest datum.
    k0=index_of_xindx(A,dataset%M,dataset%A,indx)

    ! Loop over directions from this point.
    do ibackforth=0,2
      select case(ibackforth)
      case(0)
        k1=k0
        k2=k0
        dk=1
      case(1)
        k1=k0-1
        k2=1
        dk=-1
      case(2)
        k1=k0+1
        k2=dataset%M
        dk=1
      end select
      ! Loop over points in chosen direction.
      do k=k1,k2,dk
        kk=indx(k)
        Ak=dataset%A(kk)
        ! Get value of h at this A.
        if(Ak<dataset%A_KDE_onset(1))then
          h=dataset%h_KDE_centre+(dataset%h_KDE_end(1)-dataset%h_KDE_centre)*&
             &abs(Ak-dataset%A_KDE_onset(1))/&
             &abs(dataset%A_KDE_end(1)-dataset%A_KDE_onset(1))
        elseif(Ak>dataset%A_KDE_onset(2))then
          h=dataset%h_KDE_centre+(dataset%h_KDE_end(2)-dataset%h_KDE_centre)*&
             &abs(Ak-dataset%A_KDE_onset(2))/&
             &abs(dataset%A_KDE_end(2)-dataset%A_KDE_onset(2))
        else
          h=dataset%h_KDE_centre
        endif
        ! Evaluate Gaussian kernel; exit if negligible.
        x=(A-Ak)/h
        if(abs(x)>dataset%max_x)exit
        t1=exp(-0.5d0*x*x)
        ! Add contribution.
        if(dataset%have_w)t1=t1*dataset%w(kk)
        eval_KDE_P=eval_KDE_P+t1/h
      enddo ! k
    enddo ! ibackforth
    ! Apply norm.
    eval_KDE_P=norm*dataset%inv_total_w*eval_KDE_P
  END FUNCTION eval_KDE_P


  ! EVALUATION OF FIT VALUE AND INTEGRALS.


  DOUBLE PRECISION FUNCTION eval_fit(x,nparam,delta_mu,param)
    !--------------------------------------------!
    ! Evaluate polynomial of order nparam-1 with !
    ! coefficients param(1:nparam) at x.         !
    !--------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x
    INTEGER, INTENT(in) :: nparam
    DOUBLE PRECISION, INTENT(in) :: delta_mu(nparam),param(nparam)
    INTEGER iparam
    eval_fit=0.d0
    do iparam=1,nparam
      eval_fit=eval_fit+param(iparam)*x**delta_mu(iparam)
    enddo ! iparam
  END FUNCTION eval_fit


  DOUBLE PRECISION FUNCTION eval_fit_P(mu,inv_A,nparam,delta_mu,param)
    !-----------------------------------!
    ! Evaluate P(A) from parameters for !
    ! reciprocal-scale fit.             !
    !-----------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: mu,inv_A
    INTEGER, INTENT(in) :: nparam
    DOUBLE PRECISION, INTENT(in) :: delta_mu(nparam),param(nparam)
    INTEGER iparam
    eval_fit_P=0.d0
    do iparam=1,nparam
      eval_fit_P=eval_fit_P+param(iparam)*(mu+delta_mu(iparam)-1.d0)*&
         &inv_A**(mu+delta_mu(iparam))
    enddo ! iparam
  END FUNCTION eval_fit_P


  SUBROUTINE integrate_fit(nparam,nparam_dim,mu,delta_mu,momdef,param,&
     &anchor_A,A_centre,mom0,mom1,mom2,have_mom1)
    !-----------------------------------------------------!
    ! Obtain analytical contributions to 0th, 1st and 2nd !
    ! moment of fitted tails of the provided parameters.  !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nparam(2),nparam_dim(2),momdef
    DOUBLE PRECISION, INTENT(in) :: mu(2),delta_mu(maxval(nparam_dim),2),&
       &param(maxval(nparam_dim),2),anchor_A(2),A_centre
    DOUBLE PRECISION, INTENT(inout) :: mom0(2),mom1(2),mom2(2)
    LOGICAL, INTENT(inout) :: have_mom1
    DOUBLE PRECISION cvec(maxval(nparam),2),x,t0,t1,t2,rexp
    INTEGER itail,jtail,iparam,jparam

    ! Initialize.
    mom0=0.d0
    mom1=0.d0
    mom2=0.d0
    if(all(nparam==0))return

    ! Evaluate expansion parameters from fit parameters.
    do itail=1,2
      do iparam=1,nparam(itail)
        cvec(iparam,itail)=(mu(itail)+delta_mu(iparam,itail)-1.d0)*&
           &param(iparam,itail)
      enddo ! iparam
    enddo ! itail

    ! Check for divergent contributions with unmatched coefficients.
    if(momdef==MOMDEF_MOM1_MAYBE)then
      have_mom1=.true.
      itail=1
      jtail=2
      do iparam=1,nparam(itail)
        rexp=mu(itail)+delta_mu(iparam,itail)
        if(rexp>2.d0)cycle
        if(mu(jtail)<=0.d0)then
          have_mom1=.false.
          exit
        endif
        ! Find matching exponent in tail 2.  NB, we assume all exponents
        ! for a given tail are different, as enforced outside this routine.
        do jparam=1,nparam(jtail)
          if(are_equal(mu(jtail)+delta_mu(jparam,jtail),rexp))exit
        enddo ! jparam
        if(jparam>nparam(jtail))then
          have_mom1=.false.
          exit
        elseif(.not.are_equal(cvec(iparam,itail),cvec(jparam,jtail)))then
          have_mom1=.false.
          exit
        endif
      enddo ! iparam
    else
      have_mom1=momdef>=MOMDEF_MOM1
    endif

    ! Evaluate contribution to mom0, mom1, and mom2.
    ! Loop over tails.
    do itail=1,2
      ! Skip if non-heavy tail.
      if(mu(itail)<=0.d0.or.nparam(itail)<1)cycle
      x=1.d0/anchor_A(itail)
      ! Loop over contributions in this tail.
      do iparam=1,nparam(itail)
        rexp=mu(itail)+delta_mu(iparam,itail)
        ! Variance.
        t2=0.d0
        if(momdef>=MOMDEF_MOM2)then
          if(rexp>3.d0)t2=x**(rexp-3.d0)/(rexp-3.d0)
          t2=cvec(iparam,itail)*t2
        endif
        ! Mean.
        t1=0.d0
        if(have_mom1)then
          if(are_equal(rexp,2.d0))then
            t1=log(x)
          elseif(rexp>1.d0)then
            t1=x**(rexp-2.d0)/(rexp-2.d0)
          endif
          t1=cvec(iparam,itail)*t1
          if(itail==1)t1=-t1
        endif
        ! Norm.
        t0=0.d0
        if(momdef>=MOMDEF_MOM0)then
          if(rexp>1.d0)t0=x**(rexp-1.d0)/(rexp-1.d0)
          t0=cvec(iparam,itail)*t0
        endif
        ! Add contributions.
        mom2(itail)=mom2(itail)+t2
        mom1(itail)=mom1(itail)+t1
        mom0(itail)=mom0(itail)+t0
      enddo ! iparam
      ! Recentre mom2 and mom1 around A_centre.
      if(momdef>=MOMDEF_MOM2)mom2(itail)=mom2(itail)+&
         &(2.d0*mom1(itail)+mom0(itail)*A_centre)*A_centre
      if(have_mom1)mom1(itail)=mom1(itail)+mom0(itail)*A_centre
    enddo ! itail

  END SUBROUTINE integrate_fit


  ! BOOTSTRAP RESAMPLE GENERATOR


  SUBROUTINE gen_bootstrap_resample(dataset,indx_work)
    !-----------------------------------------------!
    ! Generate a bootstrap index vector for a whole !
    ! dataset resample.                             !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    INTEGER, INTENT(inout) :: indx_work(dataset%M)
    INTEGER i,it1
    DOUBLE PRECISION t1
    ! Loop over samples.
    do i=1,dataset%M
      call random_number(t1)
      if(dataset%have_w)then
        t1=t1*(dataset%q(dataset%M)+dataset%q(1))
        it1=index_of_x_floor(t1-dataset%q(1),dataset%M,dataset%q)
      else
        it1=floor(t1*dble(dataset%M))+1
      endif
      indx_work(i)=it1
    enddo ! i
    ! In-place sort.
    call quicksort3_int_inplace(dataset%M,indx_work)
    ! Map.
    do i=1,dataset%M
      indx_work(i)=dataset%indx(indx_work(i))
    enddo ! i
  END SUBROUTINE gen_bootstrap_resample


  ! REFRESH ROUTINES


  SUBROUTINE refresh_dataset(dataset)
    !--------------------------------------------------!
    ! Generate sort indices and quantile functions for !
    ! dataset, both left-to-right and right-to-left,   !
    ! and compute the inverse of the total weight.     !
    !--------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(inout) :: dataset
    INTEGER k,k0,sgn,itail,nKDE_min
    DOUBLE PRECISION t1,t2,t3,h,dA

    ! Initialize.
    if(allocated(dataset%indx))deallocate(dataset%indx)
    allocate(dataset%indx(dataset%M))
    if(allocated(dataset%q))deallocate(dataset%q)
    if(dataset%have_w)allocate(dataset%q(dataset%M))

    ! Get sort indices.
    call quicksort(dataset%M,dataset%A,dataset%indx)

    ! Get inverse of sum of weights using robust sum routine.
    if(dataset%have_w)then
      dataset%inv_total_w=1.d0/sum(dataset%w)
    else
      dataset%inv_total_w=1.d0/dble(dataset%M)
    endif

    ! Evaluate (left-to-right) quantiles from sorted weights.
    if(dataset%have_w)call gen_quantiles_w(dataset%M,dataset%M,dataset%w,&
       &dataset%inv_total_w,dataset%q,dataset%indx)

    ! Compute mean, stderr.
    if(dataset%have_w)then
      call characterize_dist_w(dataset%M,dataset%A,dataset%w,&
         &mean=dataset%A_mean,stderr=dataset%A_stderr,&
         &var=dataset%A_var,err_var=dataset%A_err_var)
    else
      call characterize_dist(dataset%M,dataset%A,mean=dataset%A_mean,&
         &stderr=dataset%A_stderr,var=dataset%A_var,&
         &err_var=dataset%A_err_var)
    endif

    ! Compute h for Gaussian kernel density estimation using a quantile-based
    ! measure of the standard deviation to avoid being thrown way off by
    ! heavy tails.  NB, 0.5*(t2-t1) should equal sqrt(var) for a normal
    ! distribution.
    if(dataset%have_w)then
      t1=dataset%A(dataset%indx(index_of_x(0.158655254d0,dataset%M,dataset%q)))
      t2=dataset%A(dataset%indx(index_of_x(1.d0-0.158655254d0,dataset%M,&
         &dataset%q)))
    else
      t1=dataset%A(dataset%indx(nint(0.158655254d0*dble(dataset%M))))
      t2=dataset%A(dataset%indx(nint((1.d0-0.158655254d0)*dble(dataset%M))))
    endif
    t3=0.5d0*(t2-t1)
    h=t3/(0.75d0*dble(dataset%M))**0.2d0
    dataset%h_KDE_centre=h

    ! Get data density at end of tail.
    nKDE_min=nint(sqrt(dble(dataset%M)))
    do itail=1,2
      if(itail==1)then
        k0=0
        sgn=1
      else
        k0=dataset%M+1
        sgn=-1
      endif
      do k=1,dataset%M/2
        dA=abs(dataset%A(dataset%indx(k0+sgn*(k+nKDE_min)))-&
         &dataset%A(dataset%indx(k0+sgn*k)))
        if(dA<dataset%max_x/dble(nKDE_min))exit
      enddo ! k
      ! This k defines linear h onset.
      dataset%A_KDE_end(itail)=dataset%A(dataset%indx(k0+sgn))
      if(k==1)then
        ! Just use single value.
        dataset%A_KDE_onset(itail)=dataset%A(dataset%indx(k0+sgn))
        dataset%h_KDE_end(itail)=h
      else
        ! Linear increase.
        dataset%A_KDE_onset(itail)=dataset%A(dataset%indx(k0+sgn*k))
        dataset%h_KDE_end(itail)=&
           &abs(dataset%A(dataset%indx(k0+sgn*(1+nKDE_min)))-&
           &dataset%A(dataset%indx(k0+sgn)))/dataset%max_x
      endif
    enddo ! itail

    ! Compute median.
    if(dataset%have_w)then
      dataset%A_median=dataset%A(dataset%indx(index_of_x(0.5d0,dataset%M,&
         &                                               dataset%q)))
    else
      dataset%A_median=dataset%A(dataset%indx(nint(0.5d0*dble(dataset%M))))
    endif

  END SUBROUTINE refresh_dataset


  SUBROUTINE refresh_centre(dataset,tail_form)
    !------------------------------!
    ! Evaluate A_centre and Mhalf. !
    !------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(in) :: dataset
    TYPE(tail_form_type),INTENT(inout) :: tail_form

    ! Do not attempt to do anything unless we have data loaded.
    if(dataset%M<1)return

    ! Reset A_centre, Mhalf.
    if(tail_form%self_consistent_centre.and.tail_form%A_result_set)then
      tail_form%A_centre=tail_form%A_result
    else
      select case(trim(tail_form%A_centre_def))
      case('median')
        tail_form%A_centre=dataset%A_median
      case('mean')
        tail_form%A_centre=dataset%A_mean
      case default
        read(tail_form%A_centre_def,*)tail_form%A_centre
      end select
    endif
    tail_form%Mhalf(1)=count(dataset%A<tail_form%A_centre)
    tail_form%Mhalf(2)=count(dataset%A>tail_form%A_centre)

  END SUBROUTINE refresh_centre


  SUBROUTINE refresh_momdef(tail_form)
    !------------------!
    ! Evaluate momdef. !
    !------------------!
    IMPLICIT NONE
    TYPE(tail_form_type),INTENT(inout) :: tail_form
    if(all(tail_form%mu<=0.d0))then
      tail_form%momdef=MOMDEF_MOM2
    elseif(any(tail_form%mu<=1.d0.and..not.tail_form%mu<=0.d0))then
      tail_form%momdef=MOMDEF_NONE
    elseif(any(tail_form%mu<=2.d0.and..not.tail_form%mu<=0.d0))then
      tail_form%momdef=MOMDEF_MOM0
      if(are_equal(tail_form%mu(1),tail_form%mu(2)).and.&
         &tail_form%nconstraint>0)tail_form%momdef=MOMDEF_MOM1_MAYBE
    elseif(any(tail_form%mu<=3.d0.and..not.tail_form%mu<=0.d0))then
      tail_form%momdef=MOMDEF_MOM1
    else
      tail_form%momdef=MOMDEF_MOM2
    endif
  END SUBROUTINE refresh_momdef


  SUBROUTINE refresh_bs_setup(bs_setup)
    !-----------------------------------!
    ! Adjust nbs to be divisible by the !
    ! number of MPI processes.          !
    !-----------------------------------!
    IMPLICIT NONE
    TYPE(bs_setup_type), INTENT(inout) :: bs_setup
    bs_setup%nbs=max(bs_setup%nbs,1)
    bs_setup%nbs_proc=bs_setup%nbs/MPI_NPROC
    if(mod(bs_setup%nbs,MPI_NPROC)/=0)bs_setup%nbs_proc=bs_setup%nbs_proc+1
    bs_setup%nbs=bs_setup%nbs_proc*MPI_NPROC
  END SUBROUTINE refresh_bs_setup


  ! INSTANTIATION ROUTINES


  SUBROUTINE instantiate_constraint_set(tail_form,nparam,nparam_dim,delta_mu,&
     &regx,regy,cset)
    !--------------------------------------------!
    ! Construct constraint arrays for a specific !
    ! value of nparam(:).                        !
    !--------------------------------------------!
    IMPLICIT NONE
    TYPE(tail_form_type), INTENT(in) :: tail_form
    INTEGER, INTENT(in) :: nparam(2), nparam_dim(2)
    DOUBLE PRECISION, INTENT(in) :: delta_mu(maxval(nparam_dim),2),regx(2),&
       &regy(2)
    TYPE(constraint_set_type), INTENT(inout) :: cset
    ! Local variables
    INTEGER itail, iparam, jtail, jparam, iconstraint, nrhs, irhs, it1
    DOUBLE PRECISION t1,t2,t3

    ! De-, re-allocate and initialize arrays.
    if(allocated(cset%ndependant))deallocate(cset%ndependant)
    if(allocated(cset%dependant_itail))deallocate(cset%dependant_itail)
    if(allocated(cset%dependant_iparam))deallocate(cset%dependant_iparam)
    if(allocated(cset%ncombine))deallocate(cset%ncombine)
    if(allocated(cset%combine_itail))deallocate(cset%combine_itail)
    if(allocated(cset%combine_iparam))deallocate(cset%combine_iparam)
    if(allocated(cset%constraint_coeff))deallocate(cset%constraint_coeff)
    allocate(cset%ndependant(maxval(nparam),2),&
       &cset%dependant_itail(sum(nparam),maxval(nparam),2),&
       &cset%dependant_iparam(sum(nparam),maxval(nparam),2),&
       &cset%ncombine(maxval(nparam),2),&
       &cset%combine_itail(sum(nparam),maxval(nparam),2),&
       &cset%combine_iparam(sum(nparam),maxval(nparam),2),&
       &cset%constraint_coeff(maxval(nparam),2,maxval(nparam),2))
    cset%ndependant=0
    cset%dependant_itail=0
    cset%dependant_iparam=0
    cset%ncombine=0
    cset%combine_itail=0
    cset%combine_iparam=0
    cset%constraint_coeff=0.d0
    do itail=1,2
      if(tail_form%mu(itail)<=0.d0)cycle
      do iparam=1,nparam(itail)
        cset%ndependant(iparam,itail)=1
        cset%dependant_itail(1,iparam,itail)=itail
        cset%dependant_iparam(1,iparam,itail)=iparam
        cset%constraint_coeff(iparam,itail,iparam,itail)=1.d0
      enddo ! iparam
    enddo ! itail

    ! Store applicable constraints.
    do iconstraint=1,tail_form%nconstraint
      ! Figure out if this constraint is applicable.
      itail=tail_form%lhs_itail(iconstraint)
      iparam=tail_form%lhs_iparam(iconstraint)
      if(iparam>nparam(itail))cycle
      nrhs=tail_form%nrhs(iconstraint)
      do irhs=1,nrhs
        jtail=tail_form%rhs_itail(irhs,iconstraint)
        jparam=tail_form%rhs_iparam(irhs,iconstraint)
        if(jparam>nparam(jtail))exit
      enddo ! irhs
      if(irhs<=nrhs)cycle
      ! Apply constraint.
      cset%ndependant(iparam,itail)=0
      cset%dependant_itail(:,iparam,itail)=0
      cset%dependant_iparam(:,iparam,itail)=0
      do irhs=1,nrhs
        jtail=tail_form%rhs_itail(irhs,iconstraint)
        jparam=tail_form%rhs_iparam(irhs,iconstraint)
        t1=tail_form%rhs_coeff(irhs,iconstraint)
        ! Update constrained->free parameter map.
        cset%ncombine(iparam,itail)=cset%ncombine(iparam,itail)+1
        cset%combine_itail(irhs,iparam,itail)=jtail
        cset%combine_iparam(irhs,iparam,itail)=jparam
        ! Update free->constrained parameter map.
        it1=cset%ndependant(jparam,jtail)+1
        cset%ndependant(jparam,jtail)=it1
        cset%dependant_itail(it1,jparam,jtail)=itail
        cset%dependant_iparam(it1,jparam,jtail)=iparam
        ! Construct regularization coefficient.
        if(itail==jtail)then
          ! Simplified for better numerics.
          t2=regx(itail)**(delta_mu(jparam,jtail)-delta_mu(iparam,itail))
        else
          ! Do by pieces for better numerics.
          t2=(regx(jtail)/regx(itail))**min(delta_mu(iparam,itail),&
             &delta_mu(jparam,jtail))
          t3=delta_mu(jparam,jtail)-delta_mu(iparam,itail)
          if(are_equal(t3,0.d0))then
            continue
          elseif(t3>0.d0)then
            t2=t2*regx(jtail)**t3
          else
            t2=t2*regx(itail)**t3
          endif
          t2=t2*(regy(itail)/regy(jtail))
        endif
        ! Set coefficient.
        cset%constraint_coeff(iparam,itail,jparam,jtail)=t1*t2
      enddo ! irhs
    enddo ! iconstraint

  END SUBROUTINE instantiate_constraint_set


  SUBROUTINE instantiate_delta_mu(tail_form,nparam,delta_mu)
    !-------------------------------------------------------!
    ! Construct delta_mu for a specific value of nparam(:). !
    ! NB, we are assuming that the exponent increment is    !
    ! constant, a constraint we might want to lift at some  !
    ! point.                                                !
    !-------------------------------------------------------!
    IMPLICIT NONE
    TYPE(tail_form_type), INTENT(in) :: tail_form
    INTEGER, INTENT(in) :: nparam(2)
    DOUBLE PRECISION, ALLOCATABLE :: delta_mu(:,:)
    INTEGER itail, iparam
    if(allocated(delta_mu))deallocate(delta_mu)
    allocate(delta_mu(maxval(nparam),2))
    delta_mu=0.d0
    do itail=1,2
      do iparam=1,nparam(itail)
        delta_mu(iparam,itail)=dble(iparam-1)*tail_form%global_delta_mu(itail)
      enddo ! iparam
    enddo ! itail
  END SUBROUTINE instantiate_delta_mu


  SUBROUTINE unregularize_param(nparam,nparam_dim,delta_mu,regx,regy,param,chi2)
    !-------------------------------------------------------!
    ! Remove regularization factors from param(:) and chi2. !
    !-------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nparam(2), nparam_dim(2)
    DOUBLE PRECISION, INTENT(in) :: delta_mu(maxval(nparam_dim),2),&
       &regx(2),regy(2)
    DOUBLE PRECISION, INTENT(inout) :: param(maxval(nparam_dim),2),chi2(2)
    INTEGER itail, iparam
    DOUBLE PRECISION t1, t2
    ! Unregularize.
    do itail=1,2
      do iparam=1,nparam(itail)
        ! Prefactor is regx**delta_mu/regy; improve numerics by case.
        t2=delta_mu(iparam,itail)
        if(are_equal(t2,0.d0))then
          t1=1.d0/regy(itail)
        elseif(are_equal(t2,1.d0))then
          t1=regx(itail)/regy(itail)
        elseif(t2>1.d0)then
          t1=regy(itail)**(1.d0/t2)
          t1=regx(itail)/t1
          if(t1>huge(1.d0)**(1.d0/t2))then
            param=0.d0
            chi2=0.d0
            return
          endif
          t1=t1**t2
        else
          t1=regx(itail)**t2/regy(itail)
        endif
        if(t1>1.d0)then
          if(abs(param(iparam,itail))>huge(1.d0)/t1)then
            param=0.d0
            chi2=0.d0
            return
          endif
        endif
        param(iparam,itail)=param(iparam,itail)*t1
      enddo ! iparam
      chi2(itail)=chi2(itail)/regy(itail)**2
    enddo ! itail
  END SUBROUTINE unregularize_param


  ! CONSTRAINT HANDLING


  SUBROUTINE parse_constraint_token(token,iparam,itail,coeff)
    !-------------------------------------------------------!
    ! Parse a string of the form "[<coeff>*]{l|r}<iparam>". !
    !-------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: token
    INTEGER,INTENT(inout) :: iparam,itail
    DOUBLE PRECISION,INTENT(inout) :: coeff
    INTEGER ipos,ierr
    iparam=0
    itail=0
    coeff=1.d0
    ipos=scan(token,'*')
    if(ipos>0)then
      read(token(:ipos-1),*,iostat=ierr)coeff
      if(ierr/=0)return
    endif
    select case(token(ipos+1:ipos+1))
    case('l')
      itail=1
    case('r')
      itail=2
    case default
      return
    end select
    read(token(ipos+2:),*,iostat=ierr)iparam
    if(ierr/=0)iparam=0
  END SUBROUTINE parse_constraint_token


  SUBROUTINE increment_constraint_storage(nrhs,tail_form)
    !-------------------------------------------------!
    ! Add space for one more constraint of up to NRHS !
    ! terms on the right-hand side.                   !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nrhs
    TYPE(tail_form_type), INTENT(inout) :: tail_form
    INTEGER n,m
    INTEGER, ALLOCATABLE :: lhs_iparam(:),lhs_itail(:),cnrhs(:),&
       &rhs_iparam(:,:),rhs_itail(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: rhs_coeff(:,:)
    n=tail_form%nconstraint
    m=0
    if(n>0)then
      m=maxval(tail_form%nrhs)
      allocate(lhs_iparam(n),lhs_itail(n),cnrhs(n),&
         &rhs_iparam(m,n),rhs_itail(m,n),rhs_coeff(m,n))
      lhs_iparam=tail_form%lhs_iparam
      lhs_itail=tail_form%lhs_itail
      cnrhs=tail_form%nrhs
      rhs_iparam=tail_form%rhs_iparam
      rhs_itail=tail_form%rhs_itail
      rhs_coeff=tail_form%rhs_coeff
      deallocate(tail_form%lhs_iparam,tail_form%lhs_itail,tail_form%nrhs,&
         &tail_form%rhs_iparam,tail_form%rhs_itail,tail_form%rhs_coeff)
    endif
    allocate(tail_form%lhs_iparam(n+1),tail_form%lhs_itail(n+1),&
       &tail_form%nrhs(n+1),tail_form%rhs_iparam(nrhs,n+1),&
       &tail_form%rhs_itail(nrhs,n+1),tail_form%rhs_coeff(nrhs,n+1))
    tail_form%lhs_iparam=0
    tail_form%lhs_itail=0
    tail_form%nrhs=0
    tail_form%rhs_iparam=0
    tail_form%rhs_itail=0
    tail_form%rhs_coeff=0.d0
    if(n>0)then
      tail_form%lhs_iparam(1:n)=lhs_iparam(1:n)
      tail_form%lhs_itail(1:n)=lhs_itail(1:n)
      tail_form%nrhs(1:n)=cnrhs(1:n)
      tail_form%rhs_iparam(1:m,1:n)=rhs_iparam(1:m,1:n)
      tail_form%rhs_itail(1:m,1:n)=rhs_itail(1:m,1:n)
      tail_form%rhs_coeff(1:m,1:n)=rhs_coeff(1:m,1:n)
      deallocate(lhs_iparam,lhs_itail,cnrhs,rhs_iparam,rhs_itail,rhs_coeff)
    endif
  END SUBROUTINE increment_constraint_storage


  ! ROUTINES FOR MANIPULATING PLOT GRIDS.


  SUBROUTINE initialize_plot_grid(plot_grid)
    !-------------------------------------------------!
    ! Empty and deallocate the contents of plot_grid. !
    !-------------------------------------------------!
    IMPLICIT NONE
    TYPE(plot_grid_type), INTENT(inout) :: plot_grid
    plot_grid%npoint=0
    if(allocated(plot_grid%def_point))deallocate(plot_grid%def_point)
    if(allocated(plot_grid%point))deallocate(plot_grid%point)
  END SUBROUTINE initialize_plot_grid


  SUBROUTINE parse_plot_grid_points(dataset,tail_form,itail,string,&
     &plot_grid,ierr)
    !----------------------------------------------------!
    ! Parse points from token and add them to plot_grid. !
    !----------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    TYPE(tail_form_type), INTENT(in) :: tail_form
    INTEGER, INTENT(in) :: itail
    CHARACTER(*), INTENT(in) :: string
    TYPE(plot_grid_type), INTENT(inout) :: plot_grid
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    CHARACTER var
    CHARACTER(5) invar
    CHARACTER, ALLOCATABLE :: def_point_save(:)
    INTEGER i,ii,n,npoint_save,ipos,nmult
    DOUBLE PRECISION, ALLOCATABLE :: xl(:),xr(:),point_save(:)

    ! Parse string into range.
    call parse_range_string(string,dataset%M,tail_form%A_centre,&
       &tail_form%global_delta_mu(1),var,n,xl,ierr)
    if(ierr/=0)return
    ! Repeat with right-hand side delta_mu, in case we are setting "x"
    ! which depends on it.
    call parse_range_string(string,dataset%M,tail_form%A_centre,&
       &tail_form%global_delta_mu(2),var,n,xr,ierr)

    ! Deal with two-sidedness.
    nmult=1
    if(itail==0)then
      ipos=scan(string,'=')
      invar=trim(adjustl(string(1:ipos-1)))
      select case(trim(invar))
      case('mlogq','1/A','logA','x')
        nmult=2
      end select
    endif

    ! Store current points.
    npoint_save=plot_grid%npoint
    if(npoint_save>0)then
      allocate(def_point_save(npoint_save),point_save(npoint_save))
      def_point_save=plot_grid%def_point
      point_save=plot_grid%point
      deallocate(plot_grid%def_point,plot_grid%point)
    endif

    ! Allocate point and restore saved contents.
    plot_grid%npoint=nmult*n+npoint_save
    allocate(plot_grid%point(plot_grid%npoint),&
       &plot_grid%def_point(plot_grid%npoint))
    if(npoint_save>0)then
      plot_grid%def_point(1:npoint_save)=def_point_save
      plot_grid%point(1:npoint_save)=point_save
      deallocate(def_point_save,point_save)
    endif

    ! Add new range to point list.
    ii=npoint_save
    do i=1,n
      ii=ii+1
      plot_grid%def_point(ii)=var
      if(nmult==1)then
        select case(var)
        case('q') ! absolute
          if(itail==2)then
            plot_grid%point(ii)=1.d0-xr(i)
          else
            plot_grid%point(ii)=xl(i)
          endif
        case('A') ! absolute
          if(itail==2)then
            plot_grid%point(ii)=xr(i)+tail_form%A_centre
          else
            plot_grid%point(ii)=tail_form%A_centre-xl(i)
          endif
        end select
      else ! nmult/=1 (implies itail==0)
        plot_grid%def_point(ii+n)=var
        select case(var)
        case('q') ! absolute
          plot_grid%point(ii)=xl(i)
          plot_grid%point(ii+n)=1.d0-xr(i)
        case('A') ! absolute
          plot_grid%point(ii)=tail_form%A_centre-xl(i)
          plot_grid%point(ii+n)=xr(i)+tail_form%A_centre
        end select
      endif ! nmult==1 or not
    enddo ! i

  END SUBROUTINE parse_plot_grid_points


  ! ROUTINES FOR MANIPULATING ASSESSMENT GRIDS.


  SUBROUTINE initialize_assess_grid(assess_grid)
    !---------------------------------------------------!
    ! Empty and deallocate the contents of assess_grid. !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(assess_grid_type), INTENT(inout) :: assess_grid
    assess_grid%nanchor=0
    if(allocated(assess_grid%def_anchor))deallocate(assess_grid%def_anchor)
    if(allocated(assess_grid%anchor))deallocate(assess_grid%anchor)
    assess_grid%nnparam=0
    if(allocated(assess_grid%nparam))deallocate(assess_grid%nparam)
    assess_grid%max_nparam=0
  END SUBROUTINE initialize_assess_grid


  SUBROUTINE parse_assess_grid_anchor(dataset,tail_form,string,assess_grid,&
     &ierr)
    !--------------------------------------------!
    ! Parse anchor pairs from token and add them !
    ! to assess_grid.                            !
    !--------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    TYPE(tail_form_type), INTENT(in) :: tail_form
    CHARACTER(*), INTENT(in) :: string
    TYPE(assess_grid_type), INTENT(inout) :: assess_grid
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    CHARACTER var
    CHARACTER, ALLOCATABLE :: def_anchor_save(:,:)
    INTEGER i,ii,n,nanchor_save
    DOUBLE PRECISION, ALLOCATABLE :: xl(:),xr(:),anchor_save(:,:)

    ! Parse string into range.
    call parse_range_string(string,dataset%M,tail_form%A_centre,&
       &tail_form%global_delta_mu(1),var,n,xl,ierr)
    if(ierr/=0)return

    ! Repeat with right-hand side delta_mu, in case we are setting "x"
    ! which depends on it.
    call parse_range_string(string,dataset%M,tail_form%A_centre,&
       &tail_form%global_delta_mu(2),var,n,xr,ierr)

    ! Store current anchors.
    nanchor_save=assess_grid%nanchor
    if(nanchor_save>0)then
      allocate(def_anchor_save(2,nanchor_save),anchor_save(2,nanchor_save))
      def_anchor_save=assess_grid%def_anchor
      anchor_save=assess_grid%anchor
      deallocate(assess_grid%def_anchor,assess_grid%anchor)
    endif

    ! Allocate anchor and restore saved contents.
    assess_grid%nanchor=n+nanchor_save
    allocate(assess_grid%anchor(2,assess_grid%nanchor),&
       &assess_grid%def_anchor(2,assess_grid%nanchor))
    if(nanchor_save>0)then
      assess_grid%def_anchor(1:2,1:nanchor_save)=def_anchor_save
      assess_grid%anchor(1:2,1:nanchor_save)=anchor_save
      deallocate(def_anchor_save,anchor_save)
    endif

    ! Add new range to anchor list.
    do i=1,n
      ii=nanchor_save+i
      assess_grid%def_anchor(1:2,ii)=var
      select case(var)
      case('q') ! each anchor from its tail end
        assess_grid%anchor(1:2,ii)=(/xl(i),xr(i)/)
      case('A') ! magnitude of distance to A_centre
        assess_grid%anchor(1:2,ii)=(/xl(i),xr(i)/)
      end select
    enddo ! i

  END SUBROUTINE parse_assess_grid_anchor


  SUBROUTINE parse_assess_grid_nparam(string,assess_grid,ierr)
    !-------------------------------------------------------!
    ! Parse nparams from token and add them to assess_grid. !
    !-------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: string
    TYPE(assess_grid_type), INTENT(inout) :: assess_grid
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    INTEGER i,ii,n,nnparam_save
    INTEGER, ALLOCATABLE :: x(:),nparam_save(:,:)

    ! Parse string into range.
    call parse_int_range_string(string,n,x,ierr)
    if(ierr/=0)return

    ! Store current nparam.
    nnparam_save=assess_grid%nnparam
    if(nnparam_save>0)then
      allocate(nparam_save(2,nnparam_save))
      nparam_save=assess_grid%nparam
      deallocate(assess_grid%nparam)
    endif

    ! Allocate nparam and restore saved contents.
    assess_grid%nnparam=n+nnparam_save
    allocate(assess_grid%nparam(2,assess_grid%nnparam))
    if(nnparam_save>0)then
      assess_grid%nparam(1:2,1:nnparam_save)=nparam_save
      deallocate(nparam_save)
    endif

    ! Add new range to nparam list.
    do i=1,n
      ii=nnparam_save+i
      assess_grid%nparam(1:2,ii)=x(i)
    enddo ! i
    assess_grid%max_nparam=maxval(assess_grid%nparam)

  END SUBROUTINE parse_assess_grid_nparam


  SUBROUTINE set_assess_grid_from_tail_form(tail_form,assess_grid)
    !-------------------------------------------------!
    ! Set a single nparam and a single anchor pair in !
    ! assess_grid following contents of tail_form.     !
    !-------------------------------------------------!
    IMPLICIT NONE
    TYPE(tail_form_type), INTENT(in) :: tail_form
    TYPE(assess_grid_type), INTENT(inout) :: assess_grid
    if(allocated(assess_grid%nparam))deallocate(assess_grid%nparam)
    if(allocated(assess_grid%anchor))deallocate(assess_grid%anchor)
    assess_grid%nnparam=1
    allocate(assess_grid%nparam(2,1))
    assess_grid%nparam(1:2,1)=tail_form%nparam(1:2)
    assess_grid%max_nparam=maxval(assess_grid%nparam)
    assess_grid%nanchor=1
    allocate(assess_grid%def_anchor(2,1),assess_grid%anchor(2,1))
    assess_grid%def_anchor(1:2,1)=tail_form%def_anchor(1:2)
    assess_grid%anchor(1:2,1)=tail_form%anchor(1:2)
  END SUBROUTINE set_assess_grid_from_tail_form


  SUBROUTINE parse_range_string(string,M,A_centre,delta_mu,var,n,x,ierr)
    !----------------------------------------------!
    ! Parse a range string into a vector X(1:N) of !
    ! values of variable VAR.  NB, for A-based     !
    ! variables we return A-A_centre, and for      !
    ! q-based variables we return the fraction of  !
    ! data from the tail end.                      !
    !----------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: string
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A_centre,delta_mu
    CHARACTER, INTENT(inout) :: var
    INTEGER, INTENT(inout) :: n
    DOUBLE PRECISION, ALLOCATABLE, INTENT(inout) :: x(:)
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    CHARACTER(len_trim(string)) remainder
    CHARACTER(5) invar
    INTEGER ipos,jpos,i,jerr
    DOUBLE PRECISION t1,t2

    ! Initialize error flag.
    ierr=1

    ! Extract variable.
    ipos=scan(string,'=')
    if(ipos<1)return
    invar=trim(adjustl(string(1:ipos-1)))
    remainder=adjustl(string(ipos+1:))

    ! Get whether this is a list or a range.
    ipos=scan(remainder,':')
    if(ipos>0)then
      ! Range (x0:x1:n).
      ! Get start point.
      t1=parse_dble(remainder(1:ipos-1),jerr)
      if(jerr/=0)return
      remainder=remainder(ipos+1:)
      ipos=scan(remainder,':')
      if(ipos<1)return
      ! Get end point.
      t2=parse_dble(remainder(1:ipos-1),jerr)
      if(jerr/=0)return
      remainder=remainder(ipos+1:)
      ! Get number of points.
      n=parse_int(remainder,jerr)
      if(jerr/=0)return
      ! Checks.
      if(n<1)return
      if(n==1.and..not.are_equal(t1,t2))return
      ! Generate values.
      allocate(x(n))
      x(1)=t1
      do i=2,n
        x(i)=t1+(t2-t1)*(dble(i-1)/dble(n-1))
      enddo ! i
    else
      ! Comma-separated list.
      ! Get number of values.
      n=1
      ipos=0
      do
        jpos=scan(remainder(ipos+1:),',')
        if(jpos<1)exit
        n=n+1
        ipos=ipos+jpos
      enddo
      ! Parse values.
      allocate(x(n))
      i=0
      do while(len_trim(remainder)>0)
        ipos=scan(remainder,',')
        if(ipos<1)ipos=len_trim(remainder)+1
        t1=parse_dble(remainder(1:ipos-1),jerr)
        if(jerr/=0)return
        i=i+1
        x(i)=t1
        remainder=remainder(ipos+1:)
      enddo
    endif

    ! Convert "invar" to "var".
    select case(invar)
    case('q')
      var='q'
    case('mlogq')
      var='q'
      x=exp(-x)
    case('M')
      var='q'
      x=(x-0.5d0)/dble(M)
    case('A')
      var='A'
      x=x-A_centre ! make relative to Ac
    case('1/A')
      var='A'
      if(any(x<=0.d0))return
      x=1.d0/x
    case('logA')
      var='A'
      x=exp(x)
    case('x')
      var='A'
      if(any(x<=0.d0))return
      x=1.d0/x**delta_mu
    case default
      return
    end select

    ! Unset error flag.
    ierr=0

  END SUBROUTINE parse_range_string


  SUBROUTINE parse_int_range_string(string,n,x,ierr)
    !-------------------------------------------!
    ! Parse a range string into a vector X(1:N) !
    ! of integer values.                        !
    !-------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: string
    INTEGER, INTENT(inout) :: n
    INTEGER, ALLOCATABLE, INTENT(inout) :: x(:)
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    CHARACTER(len_trim(string)) remainder
    CHARACTER(6) invar
    INTEGER ipos,jpos,i,j,jerr
    INTEGER it1,it2,it3

    ! Initialize error flag.
    ierr=1

    ! Extract variable.
    ipos=scan(string,'=')
    if(ipos<1)return
    invar=trim(adjustl(string(1:ipos-1)))
    if(invar/='nparam')return
    remainder=adjustl(string(ipos+1:))

    ! Get whether this is a list or a range.
    ipos=scan(remainder,':')
    if(ipos>0)then
      ! Range (x0:x1[:stride]).
      ! Get start point.
      it1=parse_int(remainder(1:ipos-1),jerr)
      if(jerr/=0)return
      remainder=remainder(ipos+1:)
      ipos=scan(remainder,':')
      if(ipos<1)ipos=len_trim(remainder)+1
      ! Get end point.
      it2=parse_int(remainder(1:ipos-1),jerr)
      if(jerr/=0)return
      remainder=remainder(ipos+1:)
      ! Get stride.
      if(len_trim(remainder)==0)then
        if(it1==it2)then
          it3=1
        else
          it3=sign(1,it2-it1)
        endif
      else
        it3=parse_int(remainder,jerr)
        if(jerr/=0)return
      endif
      ! Count number of steps.
      n=0
      do i=it1,it2,it3
        n=n+1
      enddo ! i
      ! Checks.
      if(n<1)return
      ! Generate values.
      allocate(x(n))
      j=0
      do i=it1,it2,it3
        j=j+1
        x(j)=i
      enddo ! i
    else
      ! Comma-separated list.
      ! Get number of values.
      n=1
      ipos=0
      do
        jpos=scan(remainder(ipos+1:),',')
        if(jpos<1)exit
        n=n+1
        ipos=ipos+jpos
      enddo
      ! Parse values.
      allocate(x(n))
      i=0
      do while(len_trim(remainder)>0)
        ipos=scan(remainder,',')
        if(ipos<1)ipos=len_trim(remainder)+1
        it1=parse_int(remainder(1:ipos-1),jerr)
        if(jerr/=0)return
        i=i+1
        x(i)=it1
        remainder=remainder(ipos+1:)
      enddo
    endif

    ! Unset error flag.
    ierr=0

  END SUBROUTINE parse_int_range_string


  SUBROUTINE parse_sum(string,ilist)
    !-------------------------------------------------!
    ! Given string "n1+n2+n3+...+nN" return a list of !
    ! integers ILIST(1:N)=(/n1,n2,...,nN/).           !
    !-------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: string
    INTEGER, ALLOCATABLE, INTENT(inout) :: ilist(:)
    ! Local variables.
    CHARACTER(len(string)) remainder
    INTEGER n, ipos, i, ierr

    ! Ensure array is deallocated.
    if(allocated(ilist))deallocate(ilist)

    ! Count number of terms and check string is valid.
    remainder=string
    n=0
    do
      ipos=scan(trim(remainder),'+')
      if(ipos<1)ipos=len_trim(remainder)+1
      i=parse_int(remainder(1:ipos-1),ierr)
      if(ierr/=0)return
      n=n+1
      remainder=trim(remainder(ipos+1:))
      if(len_trim(remainder)==0)exit
    enddo
    if(n<1)return

    ! Allocate array and read entries.
    allocate(ilist(n))
    remainder=string
    n=0
    do
      ipos=scan(trim(remainder),'+')
      if(ipos<1)ipos=len_trim(remainder)+1
      n=n+1
      ilist(n)=parse_int(remainder(1:ipos-1),ierr)
      remainder=trim(remainder(ipos+1:))
      if(len_trim(remainder)==0)exit
    enddo

  END SUBROUTINE parse_sum


  SUBROUTINE get_anchor(M,A,indx,is_htail,A_centre,def_anchor,anchor,anchor_M,&
     &anchor_A)
    !-----------------------------------------!
    ! Get the location of the anchors for the !
    ! sorted dataset A(indx(1:M)).            !
    !-----------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M, indx(M)
    DOUBLE PRECISION, INTENT(in) :: A(M), A_centre, anchor(2)
    LOGICAL, INTENT(in) :: is_htail(2)
    CHARACTER, INTENT(in) :: def_anchor(2)
    INTEGER, INTENT(inout) :: anchor_M(2)
    DOUBLE PRECISION, INTENT(inout) :: anchor_A(2)
    ! Local variables.
    INTEGER itail, k, k1
    DOUBLE PRECISION sgn, t1

    ! Intialize.
    anchor_M=0
    anchor_A=0.d0

    ! Loop over tails.
    do itail=1,2
      if(.not.is_htail(itail))cycle
      ! Set M as per the chosen selection method.
      select case(def_anchor(itail))
      case('q')
        ! Exploiting the fact that q is unweighted.
        anchor_M(itail)=nint(dble(M)*anchor(itail)+0.5d0)
      case('A')
        sgn=1.d0
        if(itail==1)sgn=-1.d0
        anchor_M(itail)=count(sgn*(A(indx(1:M))-(A_centre+sgn*anchor(itail)))&
           &>=0.d0)
      end select
      ! Set anchor_A from value of anchor_M.
      if(itail==1)then
        k=anchor_M(itail)
        k1=k+1
      else
        k=M+1-anchor_M(itail)
        k1=k-1
      endif
      ! Guard against anchor crossing over to wrong tail.
      t1=0.5d0*(A(indx(k))+A(indx(k1)))-A_centre
      if(itail==1)t1=-t1
      if(t1<0.d0)then
        anchor_M(itail)=0
        if(itail==1)then
          t1=A_centre-A(indx(1))
        else
          t1=A(indx(M))-A_centre
        endif
      endif
      anchor_A(itail)=t1
    enddo ! itail

  END SUBROUTINE get_anchor


  INTEGER FUNCTION get_point_itail(dataset,tail_form,plot_grid,ipoint)
    !---------------------------------------------!
    ! Decide whether point IPOINT in plot_grid is !
    ! in the left or the right tail.              !
    !---------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    TYPE(tail_form_type), INTENT(in) :: tail_form
    TYPE(plot_grid_type), INTENT(in) :: plot_grid
    INTEGER, INTENT(in) :: ipoint
    if(plot_grid%def_point(ipoint)=='q')then
      if(nint(plot_grid%point(ipoint)*dble(dataset%M))<=&
         &tail_form%Mhalf(1))then
        get_point_itail=1
      else
        get_point_itail=2
      endif
    else
      if(nint(plot_grid%point(ipoint))<tail_form%A_centre)then
        get_point_itail=1
      else
        get_point_itail=2
      endif
    endif
  END FUNCTION get_point_itail


  SUBROUTINE get_point_M_A(dataset,indx,plot_grid,ipoint,M,A)
    !-----------------------------------------------!
    ! Get the value of M and A corresponding to the !
    ! IPOINTth point in plot_grid for DATASET under !
    ! bootstrapping index INDX.                     !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    INTEGER, INTENT(in) :: indx(dataset%M), ipoint
    TYPE(plot_grid_type), INTENT(in) :: plot_grid
    DOUBLE PRECISION, INTENT(inout) :: M, A
    INTEGER j, j1, j2
    DOUBLE PRECISION t1, t2
    select case(plot_grid%def_point(ipoint))
    case('q')
      M=plot_grid%point(ipoint)*dble(dataset%M)
      j=nint(M)
      if(are_equal(dble(j),M))then
        A=dataset%A(indx(j))
      else
        j=floor(M)
        A=dataset%A(indx(j))*(dble(j+1)-M)+dataset%A(indx(j+1))*(M-dble(j))
      endif
    case('A')
      A=plot_grid%point(ipoint)
      j1=index_of_xindx(A,dataset%M,dataset%A,indx)
      t1=dataset%A(indx(j1))
      if(t1<=plot_grid%point(ipoint))then
        j2=min(j1+1,dataset%M)
      else
        j2=max(j1-1,1)
      endif
      t2=dataset%A(indx(j2))
      if(are_equal(t1,t2))then
        M=dble(j1)
      else
        M=dble(j1)+((A-t1)/(t2-t1))*dble(j2-j1)
      endif
    end select
  END SUBROUTINE get_point_M_A


  SUBROUTINE get_point_iabs_M_itail(dataset,indx,tail_form,plot_grid,ipoint,&
     &M,itail)
    !---------------------------------------------------!
    ! Get the nearest integer value of M from the tail  !
    ! and tail identifier ITAIL corresponding to the    !
    ! IPOINTth point in PLOT_GRID for DATASET under     !
    ! bootstrapping index INDX.                         !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type), INTENT(in) :: dataset
    INTEGER, INTENT(in) :: indx(dataset%M), ipoint
    TYPE(tail_form_type), INTENT(in) :: tail_form
    TYPE(plot_grid_type), INTENT(in) :: plot_grid
    INTEGER, INTENT(inout) :: M, itail
    DOUBLE PRECISION t1
    select case(plot_grid%def_point(ipoint))
    case('q')
      M=nint(plot_grid%point(ipoint)*dble(dataset%M))
      t1=dataset%A(indx(M))
      itail=1
      if(t1>tail_form%A_centre)then
        itail=2
        M=dataset%M+1-M
      endif
    case('A')
      itail=1
      M=index_of_xindx(plot_grid%point(ipoint),dataset%M,dataset%A,indx)
      if(plot_grid%point(ipoint)>tail_form%A_centre)then
        itail=2
        M=dataset%M+1-M
      endif
    end select
  END SUBROUTINE get_point_iabs_M_itail


  ! NUMERICAL UTILITIES


  INTEGER FUNCTION index_of_x(xk,M,x)
    !---------------------------------------------------!
    ! Given a sorted vector x(1:M), return index k such !
    ! that x(k) is the closest value to xk.  Works for  !
    ! both ascending and descending sorts.              !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: xk,x(M)
    INTEGER k1,k2

    ! Handle out-of-range values.
    if(xk<=x(1).and.xk<=x(M))then
      if(x(1)<=x(M))then
        index_of_x=1
      else
        index_of_x=M
      endif
      return
    elseif(xk>=x(1).and.xk>=x(M))then
      if(x(1)<=x(M))then
        index_of_x=M
      else
        index_of_x=1
      endif
      return
    endif

    ! Locate by bisection.  NB, s=1 (-1) for ascending (descending) sorts.
    k1=1 ! k1 : s*x(k1) < s*xk (smallest candidate)
    k2=M ! k2 : s*x(k2) > s*xk (largest candidate)
    do
      index_of_x=(k1+k2)/2
      if(index_of_x==k1)then
        ! Converged - select closest of the two points.
        if(abs(xk-x(k1))>abs(xk-x(k2)))index_of_x=k2
        exit
      endif
      if(.not.(x(index_of_x)<xk.or.x(index_of_x)>xk))exit
      if(x(k1)<xk.eqv.xk<x(index_of_x))then
        k2=index_of_x
      else
        k1=index_of_x
      endif
    enddo ! bisection iterations

  END FUNCTION index_of_x


  INTEGER FUNCTION index_of_xindx(xk,M,x,indx)
    !-------------------------------------------------!
    ! Given a vector x(1:M) such that x(indx(1:M)) is !
    ! sorted, return index k such that x(indx(k)) is  !
    ! the closest value to xk.  Works for both        !
    ! ascending and descending sorts.                 !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M,indx(M)
    DOUBLE PRECISION, INTENT(in) :: xk,x(*)
    INTEGER k1,k2

    ! Handle out-of-range values.
    if(xk<=x(indx(1)).and.xk<=x(indx(M)))then
      if(x(indx(1))<=x(indx(M)))then
        index_of_xindx=1
      else
        index_of_xindx=M
      endif
      return
    elseif(xk>=x(indx(1)).and.xk>=x(indx(M)))then
      if(x(indx(1))<=x(indx(M)))then
        index_of_xindx=M
      else
        index_of_xindx=1
      endif
      return
    endif

    ! Locate by bisection.  NB, s=1 (-1) for ascending (descending) sorts.
    k1=1 ! k1 : s*x(k1) < s*xk (smallest candidate)
    k2=M ! k2 : s*x(k2) > s*xk (largest candidate)
    do
      index_of_xindx=(k1+k2)/2
      if(index_of_xindx==k1)then
        ! Converged - select closest of the two points.
        if(abs(xk-x(indx(k1)))>abs(xk-x(indx(k2))))index_of_xindx=k2
        exit
      endif
      if(.not.(x(indx(index_of_xindx))<xk.or.x(indx(index_of_xindx))>xk))exit
      if(x(indx(k1))<xk.eqv.xk<x(indx(index_of_xindx)))then
        k2=index_of_xindx
      else
        k1=index_of_xindx
      endif
    enddo ! bisection iterations

  END FUNCTION index_of_xindx


  INTEGER FUNCTION index_of_x_floor(xk,M,x)
    !--------------------------------------------------!
    ! Given a vector x(1:M) sorted in ascending order, !
    ! return smallest index k such that x(k) > xk.     !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: xk,x(M)
    INTEGER k1,k2

    ! Handle first quantile and out-of-range values.
    if(xk<=x(1))then
      index_of_x_floor=1
      return
    elseif(xk>x(M))then
      index_of_x_floor=M
      return
    endif

    ! Locate by bisection.
    k1=1 ! k1 : x(k1) < xk  (not a candidate)
    k2=M ! k2 : x(k2) >= xk (largest candidate)
    do
      index_of_x_floor=(k1+k2)/2
      ! Detect k2=k1+1.
      if(index_of_x_floor==k1)then
        index_of_x_floor=k2
        exit
      endif
      ! See if we have found the correct point.
      if(x(index_of_x_floor-1)<xk.and.xk<=x(index_of_x_floor))exit
      if(x(index_of_x_floor)<xk)then
        k1=index_of_x_floor
      else
        k2=index_of_x_floor
      endif
    enddo ! bisection iterations

  END FUNCTION index_of_x_floor


  DOUBLE PRECISION FUNCTION hyper2f1_11(c,x)
    !--------------------------------------------------------!
    ! Evaluate the hypergeometric 2F1 function at (1,1,c,x). !
    !--------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: c,x
    INTEGER k
    DOUBLE PRECISION r
    DOUBLE PRECISION,PARAMETER :: CONV_TOL=1.d-10
    ! Compute hypergeometric series.
    hyper2f1_11=1.d0
    r=1.d0
    k=1
    do
      r=r*(c+dble(k-2))**2/(dble(k)*(c+dble(k-1)))*x
      hyper2f1_11=hyper2f1_11+r
      if(abs(r)<=abs(hyper2f1_11)*CONV_TOL)exit
      k=k+1
    enddo
    hyper2f1_11=(1.d0-x)**(c-2.d0)*hyper2f1_11
  END FUNCTION hyper2f1_11


  SUBROUTINE reblock_w(M,A,w,mean,stderr,var,err_factor)
    !------------------------------------------------------------!
    ! Given a series of M data A(1:M) and weights W(1:M), return !
    ! the mean and reblocked errorbar.                           !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M),w(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,var,err_factor
    INTEGER nrtn,rtn,block_m,nblock,k,i
    DOUBLE PRECISION Aw(M),sum_w,ave_A,var_A,var0,sum_block_w2,&
       &block_sum_w,block_ave_A,eff_nblock,err_A,delta_err_A,ncorr,cc_err_A,&
       &err_A_vector(floor(log(dble(M))/log(2.d0)))

    ! Quick return.
    if(M==1)then
      if(present(mean))mean=A(1)
      if(present(stderr))stderr=0.d0
      if(present(var))var=0.d0
      return
    endif

    ! Basic quantities.
    Aw=A*w
    sum_w=sum(w)
    ave_A=sum(Aw)/sum_w

    ! Loop over reblocking transformation numbers (RTNs).
    nrtn=floor(log(dble(M))/log(2.d0))
    block_m=1
    do rtn=1,nrtn

      ! Number of blocks at this RTN.
      nblock=int(M/block_m)

      ! Evaluate the sum of the squares of the deviations from the average.
      ! Last, incomplete block has fewer data points and hence a smaller
      ! weight.
      var_A=0.d0
      sum_block_w2=0.d0
      k=0
      do i=1,nblock
        block_sum_w=sum(w(k+1:k+block_m))
        block_ave_A=sum(Aw(k+1:k+block_m))/block_sum_w
        var_A=var_A+block_sum_w*(block_ave_A-ave_A)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        k=k+block_m
      enddo ! i
      eff_nblock=dble(nblock)
      if(M>k)then
        block_sum_w=sum(w(k+1:M))
        block_ave_A=sum(Aw(k+1:M))/block_sum_w
        var_A=var_A+block_sum_w*(block_ave_A-ave_A)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        eff_nblock=eff_nblock+dble(M-k)/dble(block_m)
      endif

      ! Evaluate variance, standard error in mean and error in standard error.
      var_A=var_A/(sum_w-sum_block_w2/sum_w)
      if(rtn==1)var0=var_A
      err_A=sqrt(var_A/eff_nblock)
      delta_err_A=err_A/sqrt(2.d0*(eff_nblock-1.d0))
      err_A_vector(rtn)=err_A

      ! Double block length for next reblock.
      block_m=block_m*2

    enddo ! rtn

    ! Analyze reblock plot to obtain correlation-corrected errorbar.
    cc_err_A=maxval(err_A_vector)
    if(err_A_vector(1)>0.d0)then
      block_m=1
      do rtn=1,nrtn
        ncorr=(err_A_vector(rtn)/err_A_vector(1))**2
        if(dble(block_m**3)>=2.d0*dble(M)*ncorr**2)then
          cc_err_A=err_A_vector(rtn)
          exit
        endif
        block_m=block_m*2
      enddo ! rtn
    endif

    ! Return result.
    if(present(mean))mean=ave_A
    if(present(stderr))stderr=cc_err_A
    if(present(var))var=var0
    if(present(err_factor))err_factor=cc_err_A/err_A_vector(1)

  END SUBROUTINE reblock_w


  SUBROUTINE reblock(M,A,mean,stderr,var,err_factor)
    !----------------------------------------------------!
    ! Given a series of M unweighted data A(1:M), return !
    ! the mean and reblocked errorbar.                   !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,var,err_factor
    INTEGER nrtn,rtn,block_m,nblock,k,i
    DOUBLE PRECISION ave_A,var_A,var0,sum_block_w2,block_sum_w,block_ave_A,&
       &eff_nblock,err_A,delta_err_A,ncorr,cc_err_A,&
       &err_A_vector(floor(log(dble(M))/log(2.d0)))

    ! Quick return.
    if(M==1)then
      if(present(mean))mean=A(1)
      if(present(stderr))stderr=0.d0
      if(present(var))var=0.d0
      return
    endif

    ! Basic quantities.
    ave_A=sum(A)/dble(M)

    ! Loop over reblocking transformation numbers (RTNs).
    nrtn=floor(log(dble(M))/log(2.d0))
    block_m=1
    do rtn=1,nrtn

      ! Number of blocks at this RTN.
      nblock=int(M/block_m)

      ! Evaluate the sum of the squares of the deviations from the average.
      ! Last, incomplete block has fewer data points and hence a smaller
      ! weight.
      var_A=0.d0
      sum_block_w2=0.d0
      k=0
      do i=1,nblock
        block_sum_w=dble(block_m)
        block_ave_A=sum(A(k+1:k+block_m))/block_sum_w
        var_A=var_A+block_sum_w*(block_ave_A-ave_A)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        k=k+block_m
      enddo ! i
      eff_nblock=dble(nblock)
      if(M>k)then
        block_sum_w=dble(M-k)
        block_ave_A=sum(A(k+1:M))/block_sum_w
        var_A=var_A+block_sum_w*(block_ave_A-ave_A)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        eff_nblock=eff_nblock+dble(M-k)/dble(block_m)
      endif

      ! Evaluate variance, standard error in mean and error in standard error.
      var_A=var_A/(dble(M)-sum_block_w2/dble(M))
      if(rtn==1)var0=var_A
      err_A=sqrt(var_A/eff_nblock)
      delta_err_A=err_A/sqrt(2.d0*(eff_nblock-1.d0))
      err_A_vector(rtn)=err_A

      ! Double block length for next reblock.
      block_m=block_m*2

    enddo ! rtn

    ! Analyze reblock plot to obtain correlation-corrected errorbar.
    cc_err_A=maxval(err_A_vector)
    if(err_A_vector(1)>0.d0)then
      block_m=1
      do rtn=1,nrtn
        ncorr=(err_A_vector(rtn)/err_A_vector(1))**2
        if(dble(block_m**3)>=2.d0*dble(M)*ncorr**2)then
          cc_err_A=err_A_vector(rtn)
          exit
        endif
        block_m=block_m*2
      enddo ! rtn
    endif

    ! Return result.
    if(present(mean))mean=ave_A
    if(present(stderr))stderr=cc_err_A
    if(present(var))var=var0
    if(present(err_factor))err_factor=cc_err_A/err_A_vector(1)

  END SUBROUTINE reblock


  SUBROUTINE characterize_dist_w(M,A,w,mean,stderr,err_stderr,var,err_var,&
     &skew,kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! a weighted data set.                            !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M),w(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,err_stderr,var,&
       &err_var,skew,kurt
    DOUBLE PRECISION V1,V2,V3,V4,m1,m2,m3,m4,K2,K3,K4,num1,num2,denom
    ! Initialize.
    if(present(mean))mean=0.d0
    if(present(stderr))stderr=0.d0
    if(present(var))var=0.d0
    if(present(err_var))err_var=0.d0
    if(present(skew))skew=0.d0
    if(present(kurt))kurt=0.d0
    ! Compute mean.
    if(present(mean).or.present(stderr).or.present(err_stderr).or.&
       &present(var).or.present(err_var).or.present(skew).or.present(kurt))then
      V1=sum(w)
      m1=sum(w*A)/V1
      if(present(mean))mean=m1
    endif
    if(M<2)return
    ! Compute variance.
    if(present(stderr).or.present(err_stderr).or.present(var).or.&
       &present(err_var).or.present(skew).or.present(kurt))then
      V2=sum(w**2)
      m2=sum(w*(A-m1)**2)/V1
      K2=m2*(V1*V1/(V1*V1-V2))
      if(present(stderr))stderr=sqrt(max(0.d0,K2/dble(M)))
      ! NB, err_stderr not properly weighted, and only valid for Gaussian
      ! distributions (use err_var below if non-Gauss).
      if(present(err_stderr))err_stderr=sqrt(max(0.d0,&
         &0.5d0*K2/(dble(M)*dble(M-1))))
      if(present(var))var=K2
    endif
    if(K2<=0.d0)return
    if(M<3)return
    ! Compute skewness.
    if(present(err_var).or.present(skew).or.present(kurt))then
      V3=sum(w**3)
      m3=sum(w*(A-m1)**3)/V1
      K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
      if(present(skew))skew=K3/K2**1.5d0
    endif
    if(M<4)return
    ! Compute kurtosis.
    if(present(err_var).or.present(kurt))then
      V4=sum(w**4)
      m4=sum(w*(A-m1)**4)/V1
      num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
      num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
      denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
      K4=(m4*(num1/denom)-m2*m2*(num2/denom))
      if(present(kurt))kurt=K4/(K2*K2)
      ! NB, err_var not properly weighted.
      if(present(err_var))err_var=&
         &sqrt((K4/(K2*K2)+dble(2*M)/dble(M-1))/dble(M))*K2
    endif
  END SUBROUTINE characterize_dist_w


  SUBROUTINE characterize_dist(M,A,mean,stderr,err_stderr,var,err_var,skew,&
     &kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! an unweighted data set.  This routine uses the  !
    ! same formulae as above with w(:)=1.             !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,err_stderr,var,&
       &err_var,skew,kurt
    DOUBLE PRECISION V1,V2,V3,V4,m1,m2,m3,m4,K2,K3,K4,num1,num2,denom
    ! Initialize.
    if(present(mean))mean=0.d0
    if(present(stderr))stderr=0.d0
    if(present(var))var=0.d0
    if(present(err_var))err_var=0.d0
    if(present(skew))skew=0.d0
    if(present(kurt))kurt=0.d0
    ! Compute mean.
    if(present(mean).or.present(stderr).or.present(err_stderr).or.&
       &present(var).or.present(err_var).or.present(skew).or.present(kurt))then
      V1=dble(M)
      m1=sum(A)/V1
      if(present(mean))mean=m1
    endif
    if(M<2)return
    ! Compute variance.
    if(present(stderr).or.present(err_stderr).or.present(var).or.&
       &present(err_var).or.present(skew).or.present(kurt))then
      V2=dble(M)
      m2=sum((A-m1)**2)/V1
      K2=m2*(V1*V1/(V1*V1-V2))
      if(present(stderr))stderr=sqrt(max(0.d0,K2/dble(M)))
      ! NB, err_stderr not properly weighted, and only valid for Gaussian
      ! distributions (use err_var below if non-Gauss).
      if(present(err_stderr))err_stderr=sqrt(max(0.d0,&
         &0.5d0*K2/(dble(M)*dble(M-1))))
      if(present(var))var=K2
    endif
    if(K2<=0.d0)return
    if(M<3)return
    ! Compute skewness.
    if(present(err_var).or.present(skew).or.present(kurt))then
      V3=dble(M)
      m3=sum((A-m1)**3)/V1
      K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
      if(present(skew))skew=K3/K2**1.5d0
    endif
    if(M<4)return
    ! Compute kurtosis.
    if(present(err_var).or.present(kurt))then
      V4=dble(M)
      m4=sum((A-m1)**4)/V1
      num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
      num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
      denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
      K4=(m4*(num1/denom)-m2*m2*(num2/denom))
      if(present(kurt))kurt=K4/(K2*K2)
      if(present(err_var))err_var=&
         &sqrt((K4/(K2*K2)+dble(2*M)/dble(M-1))/dble(M))*K2
    endif
  END SUBROUTINE characterize_dist


  SUBROUTINE characterize_dist_indx(M,A,indx,mean,stderr,err_stderr,var,&
     &err_var,skew,kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! an indexed unweighted data set.  This routine   !
    ! uses the same formulae as above.                !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M,indx(M)
    DOUBLE PRECISION, INTENT(in) :: A(*)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,err_stderr,var,&
       &err_var,skew,kurt
    DOUBLE PRECISION V1,V2,V3,V4,m1,m2,m3,m4,K2,K3,K4,num1,num2,denom
    ! Initialize.
    if(present(mean))mean=0.d0
    if(present(stderr))stderr=0.d0
    if(present(var))var=0.d0
    if(present(err_var))err_var=0.d0
    if(present(skew))skew=0.d0
    if(present(kurt))kurt=0.d0
    ! Compute mean.
    if(present(mean).or.present(stderr).or.present(err_stderr).or.&
       &present(var).or.present(err_var).or.present(skew).or.present(kurt))then
      V1=dble(M)
      m1=sum(A(indx(1:M)))/V1
      if(present(mean))mean=m1
    endif
    if(M<2)return
    ! Compute variance.
    if(present(stderr).or.present(err_stderr).or.present(var).or.&
       &present(err_var).or.present(skew).or.present(kurt))then
      V2=dble(M)
      m2=sum((A(indx(1:M))-m1)**2)/V1
      K2=m2*(V1*V1/(V1*V1-V2))
      if(present(stderr))stderr=sqrt(max(0.d0,K2/dble(M)))
      ! NB, err_stderr not properly weighted, and only valid for Gaussian
      ! distributions (use err_var below if non-Gauss).
      if(present(err_stderr))err_stderr=sqrt(max(0.d0,&
         &0.5d0*K2/(dble(M)*dble(M-1))))
      if(present(var))var=K2
    endif
    if(K2<=0.d0)return
    if(M<3)return
    ! Compute skewness.
    if(present(err_var).or.present(skew).or.present(kurt))then
      V3=dble(M)
      m3=sum((A(indx(1:M))-m1)**3)/V1
      K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
      if(present(skew))skew=K3/K2**1.5d0
    endif
    if(M<4)return
    ! Compute kurtosis.
    if(present(err_var).or.present(kurt))then
      V4=dble(M)
      m4=sum((A(indx(1:M))-m1)**4)/V1
      num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
      num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
      denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
      K4=(m4*(num1/denom)-m2*m2*(num2/denom))
      if(present(kurt))kurt=K4/(K2*K2)
      ! NB, err_var not properly weighted.
      if(present(err_var))err_var=&
         &sqrt((K4/(K2*K2)+dble(2*M)/dble(M-1))/dble(M))*K2
    endif
  END SUBROUTINE characterize_dist_indx


  SUBROUTINE get_conf_intervals(M,A,mean,CI)
    !------------------------------------------------------------------!
    ! Given A(1:M), return the values of A from the dataset closest to !
    ! the 2.3% | 15.9% | 84.1% | 97.7% quantiles.  These correspond to !
    ! the one-sigma and two-sigma confidence intervals for a normal    !
    ! distribution (hence the significance of these particular         !
    ! quantiles).                                                      !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M)
    DOUBLE PRECISION, INTENT(inout) :: mean, CI(4)
    INTEGER indx(M)
    mean=sum(A)/dble(M)
    call quicksort(M,A,indx)
    CI(1)=A(indx(floor(1.d0+0.022750131948179d0*dble(M))))
    CI(2)=A(indx(floor(1.d0+0.158655253931457d0*dble(M))))
    CI(3)=A(indx(floor(1.d0+0.841344746068543d0*dble(M))))
    CI(4)=A(indx(floor(1.d0+0.977249868051821d0*dble(M))))
  END SUBROUTINE get_conf_intervals


  SUBROUTINE quicksort(n,x,indx,preinitialized)
    !------------------------------------------------------------------!
    ! Apply the quicksort algorithm to create an integer index vector  !
    ! INDX such that the values X(INDX(I)) increase with increasing I. !
    !                                                                  !
    ! This is Leonard J. Moss' 1986 SORTX routine, see                 !
    !  http://www.fortran.com/quick_sort2.f                            !
    ! de-goto-ed by PLR, 06.2009.                                      !
    !                                                                  !
    ! In the comments the square-bracket notation 'X[I]' is used as    !
    ! shorthand for 'X(INDX(I))'.                                      !
    !                                                                  !
    ! For information about quicksort, see:                            !
    !  http://en.wikipedia.org/wiki/Quicksort                          !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(inout) :: indx(n)
    DOUBLE PRECISION, INTENT(in) :: x(n)
    LOGICAL, INTENT(in), OPTIONAL :: preinitialized
    ! Size of subsequence below which we stop QuickSort and move to the final
    ! straight-insertion sorting stage.  Optimal value by 1986 standards is
    ! m =~ 9.  A simple test for array sizes ranging from 500 to 100000 on an
    ! Intel Core 2 processor with the ifort compiler seem to point to m =~ 30
    ! as the optimal value in 2009 -- the difference in timing is small,
    ! though.
    INTEGER, PARAMETER :: m=30
    ! Maximum stack size
    INTEGER, PARAMETER :: stack_size=128
    INTEGER l,r,i,j,p,lstk(stack_size),rstk(stack_size),istk
    LOGICAL already_init
    DOUBLE PRECISION xp

    already_init=.false.
    if(present(preinitialized))already_init=preinitialized

    if(.not.already_init)then
      ! Initialize INDX.
      do i=1,n
        indx(i)=i
      enddo ! i
    endif

    if(n>m)then ! Skip QuickSort for small vectors.

      ! Initialize left/right region boundaries and boundary stack, and start
      ! main loop.
      istk=0
      l=1
      r=n
      do

        ! Sort the subsequence X[L]..X[R].
        ! At this point, X[KL] <= X[KM] <= X[KR] for all KL < L, KR > R, and
        ! L <= KM <= R.  This is not applicable to the first iteration where
        ! there is no data for KL < L or KR > R.
        i=l
        j=r

        ! Let the pivot, P, be the midpoint of this subsequence, P=(L+R)/2.
        ! Then rearrange INDX(L), INDX(P), and INDX(R) so the corresponding X
        ! values are in increasing order.  The pivot key is then X[P].
        p=(l+r)/2
        if(x(indx(l))>x(indx(p)))call swap1(indx(p),indx(l))
        if(x(indx(p))>x(indx(r)))call swap1(indx(p),indx(r))
        if(x(indx(l))>x(indx(p)))call swap1(indx(p),indx(l))
        xp=x(indx(p))

        ! Inner loop over elements to sort in L..R.  We want to swap values
        ! between the right and left sides and/or move X[P] until all smaller
        ! values are left of P and all larger values are right of P.  At the
        ! end of this process neither the left or right side will be
        ! internally ordered yet, but X[P] will be in its final position.
        do
          ! Search for datum on left >= X[P] and datum on right <= X[P].
          ! At this point X[L] <= X[P] <= X[R], therefore we can start
          ! scanning up from L and down from R to find the required elements.
          i=i+1
          do while(x(indx(i))<xp)
            i=i+1
          enddo
          j=j-1
          do while(x(indx(j))>xp)
            j=j-1
          enddo
          if(i>=j)exit ! exit when the two scans collide
          call swap1(indx(i),indx(j))
        enddo ! loop over elements to sort in L..R.

        ! Select next subsequence to sort.  At this point, I >= J.  The
        ! elements in the left subsequence {X[KL], L <= KL < I} and right
        ! subsequence {X[KR], J < KR <= R} verify that
        ! X[KL] <= X[I] == X[P] <= X[KR].
        if(r-j>=i-l.and.i-l>m)then
          ! Both subsequences are more than M elements long. Push longer (left)
          ! on stack and QuickSort the shorter (right).
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort: stack size too small <1>.'
            stop
          endif
          lstk(istk)=j+1
          rstk(istk)=r
          r=i-1
        elseif(i-l>r-j.and.r-j>m)then
          ! Both subsequences are more than M elements long. Push longer
          ! (right) on stack and QuickSort the shorter (left).
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort: stack size too small <2>.'
            stop
          endif
          lstk(istk)=l
          rstk(istk)=i-1
          l=j+1
        elseif(r-j>m)then
          ! Only right subsequence is more than M elements long. QuickSort it.
          l=j+1
        elseif(i-l>m)then
          ! Only left subsequence is more than M elements long. QuickSort it.
          r=i-1
        else
          ! Both subsequences are less than M elements long. Pop the stack, or
          ! terminate QuickSort if empty
          if(istk<1)exit
          l=lstk(istk)
          r=rstk(istk)
          istk=istk-1
        endif

      enddo ! main loop

    endif ! n>m

    ! Final straight-insertion sorting stage.
    do i=2,n
      if(x(indx(i-1))<=x(indx(i)))cycle
      call swap1(indx(i-1),indx(i))
      do j=i-1,2,-1
        if(x(indx(j-1))<=x(indx(j)))exit
        call swap1(indx(j-1),indx(j))
      enddo ! j
    enddo ! i

  END SUBROUTINE quicksort


  SUBROUTINE quicksort3_int_inplace(n,x)
    !-------------------------------------------------------!
    ! 3-way QuickSort (=QuickSort3) for an integer vector   !
    ! X(1:N), resulting in X(1:N) sorted in ascending       !
    ! order (sorting performed in place in this version).   !
    !                                                       !
    ! QuickSort3 is the best method for integers because it !
    ! takes repetitions into account, and is therefore      !
    ! a stable sorting algorithm.  Repetitions are rare     !
    ! or non-exisiting in many real cases, but are to be    !
    ! expected in integer problems.  There is additional    !
    ! overhead due to checking for equalities, so 2-way     !
    ! QuickSort (=QuickSort) should be used when no         !
    ! repetitions are expected.                             !
    !                                                       !
    ! PLR 02.2012                                           !
    !-------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(inout) :: x(n)
    ! Size of subsequence below which we stop QuickSort and move to the final
    ! straight-insertion sorting stage.  Optimal value by 1986 standards is
    ! m =~ 9.  A simple test for array sizes ranging from 500 to 100000 on an
    ! Intel Core 2 processor with the ifort compiler seem to point to m =~ 30
    ! as the optimal value in 2009 -- the difference in timing is small,
    ! though.
    INTEGER, PARAMETER :: m=30
    ! Maximum stack size
    INTEGER, PARAMETER :: stack_size=128
    INTEGER l,r,i,j,k,p,q,xr,lstk(stack_size),rstk(stack_size),istk

    if(n>m)then ! Skip QuickSort for small vectors.

      ! Initialize left/right region boundaries and boundary stack, and start
      ! main loop.
      istk=0
      l=1
      r=n
      do
        i=l-1
        j=r
        xr=x(r)
        p=l-1
        q=r
        do
          i=i+1
          do while(x(i)<xr)
            i=i+1
          enddo
          j=j-1
          do while(x(j)>xr)
            if(j==l)exit
            j=j-1
          enddo
          if(i>=j)exit
          call swap1(x(i),x(j))
          if(x(i)==xr)then
            p=p+1
            call swap1(x(p),x(i))
          endif
          if(xr==x(j))then
            q=q-1
            call swap1(x(j),x(q))
          endif
        enddo
        call swap1(x(i),x(r))
        j=i+1
        i=i-1
        do k=l,p-1
          call swap1(x(k),x(i))
          i=i-1
        enddo
        do k=r-1,q+1,-1
          call swap1(x(j),x(k))
          j=j+1
        enddo
        if(r-j>=i-l.and.i-l>=m)then
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort3_int: stack size too small <1>.'
            stop
          endif
          lstk(istk)=j
          rstk(istk)=r
          r=i
        elseif(i-l>r-j.and.r-j>=m)then
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort3_int: stack size too small <2>.'
            stop
          endif
          lstk(istk)=l
          rstk(istk)=i
          l=j
        elseif(i-l>=m)then
          r=i
        elseif(r-j>=m)then
          l=j
        else
          if(istk<1)exit
          l=lstk(istk)
          r=rstk(istk)
          istk=istk-1
        endif
      enddo

    endif

    ! Final straight-insertion sorting stage.
    do i=2,n
      if(x(i-1)<=x(i))cycle
      call swap1(x(i-1),x(i))
      do j=i-1,2,-1
        if(x(j-1)<=x(j))exit
        call swap1(x(j-1),x(j))
      enddo ! j
    enddo ! i

  END SUBROUTINE quicksort3_int_inplace


  SUBROUTINE swap1(x,y)
    !------------------------!
    ! Swap integers X and Y. !
    !------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: x,y
    INTEGER z
    z=x
    x=y
    y=z
  END SUBROUTINE swap1


  LOGICAL ELEMENTAL FUNCTION are_equal(x,y,tol)
    !------------------------------------------------------!
    ! Check if two floating-point numbers are equal within !
    ! a reasonable tolerance.                              !
    !------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION abs_x,abs_y,big,small
    ! Parameters.
    DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-12
    DOUBLE PRECISION,PARAMETER :: tol_rel=1.d-9
    if(present(tol))then
      are_equal=abs(x-y)<=tol
    else
      abs_x=abs(x)
      abs_y=abs(y)
      if(abs_x<=tol_zero.and.abs_y<=tol_zero)then
        are_equal=.true.
      elseif(x>0.d0.eqv.y>0.d0)then
        big=max(abs_x,abs_y)
        small=min(abs_x,abs_y)
        are_equal=big-small<=big*tol_rel
      else
        are_equal=.false.
      endif
    endif
  END FUNCTION are_equal


  ! FIELD PARSERS.


  FUNCTION field(n,line)
    !--------------------------------------------------------!
    ! Return the N-th field in string LINE, where the fields !
    ! are separated by one or more spaces.  An empty string  !
    ! is returned if N<1.                                    !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    CHARACTER(*),INTENT(in) :: line
    CHARACTER(len(line)) :: field
    CHARACTER(len(line)) remainder
    INTEGER i,k
    ! Initialize.
    field=''
    if(n<1)return
    remainder=adjustl(line)
    ! Loop over fields.
    i=0
    do
      i=i+1
      if(remainder(1:1)=='"')then
        ! Potentially start of a multi-word field.
        ! Locate end of field (quote followed by <space> or <EOL>).
        k=index(trim(remainder(2:)),'" ')+1
        if(k==1)then
          ! quote-space not found, so see if there is a quote at EOL.
          k=len_trim(remainder)
          if(remainder(k:k)/='"')k=1
        endif
        if(k>1)then
          ! Found end of field.
          if(i==n)then
            field=trim(remainder(2:k-1))
          else
            remainder=adjustl(remainder(k+1:))
          endif
          cycle
        endif
      endif
      ! Single-word field.
      ! Locate end of field.
      k=scan(trim(remainder),' ')
      if(k==0)then
        ! End is EOL.
        if(i==n)field=trim(remainder)
        return
      elseif(i==n)then
        field=trim(remainder(1:k-1))
        return
      else
        remainder=adjustl(remainder(k:))
      endif
    enddo
  END FUNCTION field


  INTEGER FUNCTION nfield(line)
    !--------------------------------------!
    ! Return the number of fields in LINE. !
    !--------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: line
    nfield=0
    do
      if(len_trim(field(nfield+1,line))==0)exit
      nfield=nfield+1
    enddo
  END FUNCTION nfield


  INTEGER FUNCTION int_field(ifield,command,ierr)
    !----------------------------------------------------!
    ! Like field, but returning the value as an integer. !
    ! ierr is set to a non-zero value if the requested   !
    ! field could not be parsed as an integer.           !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ifield
    CHARACTER(*),INTENT(in) :: command
    INTEGER,INTENT(inout) :: ierr
    int_field=parse_int(field(ifield,command),ierr)
  END FUNCTION int_field


  INTEGER FUNCTION parse_int(string,ierr)
    !--------------------------------------!
    ! Parse a string to obtain an integer. !
    !--------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ierr
    read(string,*,iostat=ierr)parse_int
  END FUNCTION parse_int


  DOUBLE PRECISION FUNCTION dble_field(ifield,command,ierr)
    !--------------------------------------------------------!
    ! Like field, but returning the value as an real number. !
    ! ierr is set to a non-zero value if the requested field !
    ! could not be parsed as a real number.                  !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ifield
    CHARACTER(*),INTENT(in) :: command
    INTEGER,INTENT(inout) :: ierr
    dble_field=parse_dble(field(ifield,command),ierr)
  END FUNCTION dble_field


  DOUBLE PRECISION FUNCTION parse_dble(string,ierr)
    !-----------------------------------------------------!
    ! Parse a string to obtain a double-precision number. !
    ! Supports fractions, e.g., string="1/3".             !
    !-----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ierr
    INTEGER ipos
    DOUBLE PRECISION t1,t2
    ipos=scan(string,'/')
    if(ipos==0)then
      read(string,*,iostat=ierr)parse_dble
    else
      ierr=-1
      if(ipos<2)return
      if(ipos>=len_trim(string))return
      read(string(1:ipos-1),*,iostat=ierr)t1
      if(ierr/=0)return
      read(string(ipos+1:),*,iostat=ierr)t2
      if(ierr/=0)return
      if(are_equal(t2,0.d0))then
        ierr=-1
        return
      endif
      parse_dble=t1/t2
    endif
  END FUNCTION parse_dble


  LOGICAL FUNCTION are_equal_string(cx,cy,tol)
    !-------------------------------------------------------!
    ! Check if two strings are equal, or if their numerical !
    ! values are equal within a reasonable tolerance.       !
    !-------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: cx,cy
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION x,y
    INTEGER ierr
    are_equal_string=trim(cx)==trim(cy)
    if(are_equal_string)return
    x=parse_dble(cx,ierr)
    if(ierr/=0)return
    y=parse_dble(cy,ierr)
    if(ierr/=0)return
    are_equal_string=are_equal(x,y,tol)
  END FUNCTION are_equal_string


  CHARACTER(2) FUNCTION type_string(have_w)
    !---------------------------------------!
    ! Produce a dataset type string "A[w]". !
    !---------------------------------------!
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: have_w
    type_string='A'
    if(have_w)type_string=trim(type_string)//'w'
  END FUNCTION type_string


  SUBROUTINE parse_type_string(string,ipos_A,ipos_w,ierr)
    !------------------------------------------------!
    ! Parse a dataset type string "A[w]" into flags. !
    !------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ipos_A,ipos_w,ierr
    ! Local variables.
    CHARACTER(len_trim(string)) remainder
    INTEGER ipos
    ierr=0
    ipos_A=0
    ipos_w=0
    remainder=adjustl(trim(string))
    ipos=0
    do while(len_trim(remainder)>0)
      ipos=ipos+1
      if(remainder(1:1)=='a'.or.remainder(1:1)=='A')then
        if(ipos_A>0)ierr=1
        ipos_A=ipos
        remainder=adjustl(remainder(2:))
      elseif(remainder(1:1)=='w')then
        if(ipos_w>0)ierr=1
        ipos_w=ipos
        remainder=adjustl(remainder(2:))
      else
        ierr=1
      endif
      if(ierr/=0)return
    enddo
  END SUBROUTINE parse_type_string


  ! STRING UTILITIES.


  CHARACTER(11) FUNCTION i2s(n)
    !---------------------------------------------!
    ! Convert integers to left-justified strings. !
    !---------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    INTEGER i,j
    INTEGER,PARAMETER :: ichar0=ichar('0')
    i2s=''
    i=abs(n)
    do j=len(i2s),1,-1
      i2s(j:j)=achar(ichar0+mod(i,10))
      i=i/10
      if(i==0)exit
    enddo ! j
    if(n<0)then
      i2s='-'//adjustl(i2s(2:11))
    else
      i2s=adjustl(i2s)
    endif ! n<0
  END FUNCTION i2s


  SUBROUTINE pprint(text,indent1,indent)
    !-------------------------------------------------------------!
    ! Print TEXT to stdout, folding lines at column 79 and using  !
    ! an indentation of INDENT1 spaces on the first line and      !
    ! INDENT on the rest.  INDENT1 defaults to INDENT, and INDENT !
    ! defaults to zero.                                           !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: text
    INTEGER,INTENT(in),OPTIONAL :: indent1,indent
    INTEGER,PARAMETER :: line_width=79
    CHARACTER(line_width) line
    CHARACTER(len(text)) remainder,word
    INTEGER ind1,indn,ind,ipos
    if(len_trim(text)==0)then
      write(6,'()')
      return
    endif
    indn=0
    if(present(indent))indn=max(indent,0)
    ind1=indn
    if(present(indent1))ind1=max(indent1,0)
    ind=ind1
    line=''
    remainder=adjustl(text)
    do while(len_trim(remainder)>0)
      ipos=scan(trim(remainder),' ')
      if(ipos==0)then
        ipos=len_trim(remainder)+1
        word=trim(remainder)
        remainder=''
      else
        word=remainder(1:ipos-1)
        remainder=adjustl(remainder(ipos:))
      endif
      if(len_trim(line)==0)then
        ! Ensure overlong words get flushed without passing through buffer.
        do while(ind+len_trim(word)>line_width)
          write(6,'(a)')repeat(' ',ind)//word(1:line_width-ind)
          word=word(line_width-ind+1:)
          ind=indn
        enddo
        ! Start new line buffer.
        line=repeat(' ',ind)//trim(word)
        ind=indn
      else
        if(len_trim(line)+1+len_trim(word)>line_width)then
          ! Flush current line buffer.
          write(6,'(a)')trim(line)
          ! Ensure overlong words get flushed without passing through buffer.
          do while(ind+len_trim(word)>line_width)
            write(6,'(a)')repeat(' ',ind)//word(1:line_width-ind)
            word=word(line_width-ind+1:)
            ind=indn
          enddo
          ! Start new line buffer.
          line=repeat(' ',ind)//trim(word)
          ind=indn
        else
          ! Add line to buffer.
          line=trim(line)//' '//trim(word)
        endif
      endif
    enddo
    if(len_trim(line)>0)write(6,'(a)')trim(line)
  END SUBROUTINE pprint


  CHARACTER(80) FUNCTION dble2human(x)
    !--------------------------------------------------------------!
    ! Print x to 4 significant figures and remove trailing zeroes. !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x
    INTEGER m,i,j
    if(abs(x)<=0.d0)then
      dble2human='0'
      return
    endif
    m=floor(log10(abs(x)))
    write(dble2human,'(f'//trim(i2s(6-m))//'.'//trim(i2s(3-m))//')')x
    dble2human=adjustl(dble2human)
    if(m-3<0)then
      j=len_trim(dble2human)
      do i=m-3,min(m,-1)
        if(dble2human(j:j)/='0')exit
        dble2human(j:j)=''
        j=j-1
      enddo ! i
      if(i==0)dble2human(j:j)=''
    endif
  END FUNCTION dble2human


  ! RUN CONTROL UTILITIES.


  SUBROUTINE init_mpi(nproc,iproc,master)
    !---------------------------------------------!
    ! Initialize MPI, get number of processes and !
    ! current process index.                      !
    !---------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: nproc,iproc
    LOGICAL, INTENT(inout) :: master
    INTEGER ierror
    call mpi_init(ierror)
    call mpi_comm_size(mpi_comm_world,nproc,ierror)
    call mpi_comm_rank(mpi_comm_world,iproc,ierror)
    master=iproc==0
  END SUBROUTINE init_mpi


  SUBROUTINE init_random(irandom)
    !-----------------------------------------------!
    ! Initialize built-in random seed so that       !
    ! different processess use different sequences. !
    !-----------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: irandom
    INTEGER i,j,k,l,n,iskip
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER, PARAMETER :: NSKIP=100
    DOUBLE PRECISION t1
    call random_seed(size=n)
    allocate(seed(n))
    do i=1,n
      if(irandom<1)then
        call system_clock(j,k,l)
        seed(i)=MPI_IPROC*n+j
      else
        seed(i)=((irandom-1)*MPI_NPROC+MPI_IPROC)*n+i
      endif
    enddo ! i
    call random_seed(put=seed)
    ! Skip a few random numbers (ifort has trouble with first).
    do iskip=1,NSKIP
      call random_number(t1)
    enddo ! iskip
    deallocate(seed)
  END SUBROUTINE init_random


  SUBROUTINE finish_mpi()
    !---------------!
    ! Finalize MPI. !
    !---------------!
    IMPLICIT NONE
    INTEGER ierror
    call mpi_finalize(ierror)
  END SUBROUTINE finish_mpi


  SUBROUTINE master_msg(msg)
    !----------------------------------------------!
    ! Write a message from the master MPI process. !
    !----------------------------------------------!
    IMPlICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if(MPI_MASTER)then
      write(6,'(a)')msg
      write(6,'()')
    endif
  END SUBROUTINE master_msg


  SUBROUTINE quit(msg)
    !---------------------!
    ! Quit with an error. !
    !---------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if(MPI_MASTER)then
      if(present(msg))then
        write(6,'(a)')'ERROR : '//msg
      else
        write(6,'(a)')'Quitting.'
      endif
    endif ! MPI_MASTER
    call finish_mpi()
    stop
  END SUBROUTINE quit


END PROGRAM treat
