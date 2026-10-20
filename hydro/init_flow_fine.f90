!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_flow_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  use amr_parameters, only: flen
  use input_hydro_condinit_module, only: r_input_hydro_condinit
  use input_hydro_grafic_module, only: r_input_hydro_grafic
  use input_hydro_gadget_module, only: r_input_hydro_gadget
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to input a given initial
  ! condition into the exisiting AMR structure, either from a file
  ! of from an analytical model given as an external routine.
  !--------------------------------------------------------------------
  character(len=flen)::filename
  logical::ok_file1,ok_file2,ok_file

  associate(s=>pst%s)

  if(s%m%noct_tot(ilevel)==0)return
  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering init_flow_fine for level ',I2)

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  if(s%r%hydro)then
     filename=TRIM(s%r%initfile(ilevel))//'/ic_d'
     INQUIRE(file=filename,exist=ok_file1)
     filename=TRIM(s%r%initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
     ok_file=ok_file1.or.ok_file2
     if(ok_file)then
        ! Read external grafic files
        if(s%r%verbose)write(*,*)'Reading initial conditions from grafic file'
        call r_input_hydro_grafic(pst,ilevel,1)
     else if (s%r%filetype=='gadget')then
        ! Read external gadget files
        if(s%r%verbose)write(*,*)'Reading initial conditions from gadget file'
        call r_input_hydro_gadget(pst,ilevel,1)
     else
        ! Use internal-defined or user-defined functions
        if(s%r%verbose)write(*,*)'Computing initial conditions from analytical model'
        call r_input_hydro_condinit(pst,ilevel,1)
     endif
  endif

  end associate

end subroutine m_init_flow_fine
!###############################################
!###############################################
!###############################################
!###############################################
