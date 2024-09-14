!===========================================================================
!===========================================================================
!This file is part of ElVibRot-TnumTana.
! MIT License
!
! Copyright (c) 2022 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
!================================================================
!     Module for "one-the-fly" (OTF) calculation of PES (pot, gradrient, hessian)
!================================================================
  MODULE mod_OTF_def
  USE TnumTana_system_m
   IMPLICIT NONE

   PRIVATE

        TYPE param_OTF
!         for the ab initio calculation
           integer :: charge       = 0
           integer :: multiplicity = -1

          character (len=Name_len)      :: ab_initio_prog     = ' '  ! gaussian, gamess
          character (len=Name_longlen)  :: ab_initio_meth     = ' '
          character (len=Name_longlen)  :: ab_initio_basis    = ' '

          character (len=Line_len)      :: commande_unix      = ' '
          logical                       :: header             = .FALSE.
          logical                       :: footer             = .FALSE.
          character (len=Name_len)      :: file_name          = ' '
          TYPE (File_t)   :: file_header
          TYPE (File_t)   :: file_footer
          TYPE (File_t)   :: file_data
          TYPE (File_t)   :: file_log
          TYPE (File_t)   :: file_FChk
          TYPE (File_t)   :: file_pun
        CONTAINS
          PROCEDURE, PRIVATE, PASS(para_OTF1) :: OTF2_TO_OTF1
          GENERIC,   PUBLIC  :: assignment(=) => OTF2_TO_OTF1
        END TYPE param_OTF


      PUBLIC :: param_OTF, dealloc_OTF, write_OTF, init_OTF

  CONTAINS

      SUBROUTINE init_OTF(para_OTF)

      TYPE (param_OTF)    :: para_OTF
!     - for the files -----------------------------------------------


!     write(out_unit,*) 'init_OTF'

      para_OTF%ab_initio_meth     = 'hf '
      para_OTF%ab_initio_basis    = ' sto-3g'
      para_OTF%ab_initio_prog     = 'g03'
      para_OTF%header             = .FALSE.
      para_OTF%footer               = .FALSE.
      para_OTF%file_name          = 'xx'
      para_OTF%commande_unix      = './g03.run ' //                     &
                          trim(adjustl(para_OTF%file_name)) // ' 2>err'

      para_OTF%file_header%name      ='xx.header'
      para_OTF%file_header%unit      = 0
      para_OTF%file_header%formatted = .TRUE.
      para_OTF%file_header%append    = .FALSE.
      para_OTF%file_header%old       = .TRUE.

      para_OTF%file_footer%name      ='xx.footer'
      para_OTF%file_footer%unit      = 0
      para_OTF%file_footer%formatted = .TRUE.
      para_OTF%file_footer%append    = .FALSE.
      para_OTF%file_footer%old       = .TRUE.

      para_OTF%file_data%name      ='xx.com'
      para_OTF%file_data%unit      = 0
      para_OTF%file_data%formatted = .TRUE.
      para_OTF%file_data%append    = .FALSE.
      para_OTF%file_data%old       = .TRUE.

      para_OTF%file_log%name      ='xx.log'
      para_OTF%file_log%unit      = 0
      para_OTF%file_log%formatted = .TRUE.
      para_OTF%file_log%append    = .FALSE.
      para_OTF%file_log%old       = .TRUE.

      para_OTF%file_FChk%name      ='Test.FChk'
      para_OTF%file_FChk%unit      = 0
      para_OTF%file_FChk%formatted = .TRUE.
      para_OTF%file_FChk%append    = .FALSE.
      para_OTF%file_FChk%old       = .TRUE.


      para_OTF%file_pun%name      ='xx.dat'
      para_OTF%file_pun%unit      = 0
      para_OTF%file_pun%formatted = .TRUE.
      para_OTF%file_pun%append    = .FALSE.
      para_OTF%file_pun%old       = .TRUE.

!     write(out_unit,*) 'END init_OTF'

      END SUBROUTINE init_OTF

      SUBROUTINE write_OTF(para_OTF)

      TYPE (param_OTF)    :: para_OTF

      write(out_unit,*) 'write_OTF'

      write(out_unit,*) 'para_OTF%charge',para_OTF%charge
      write(out_unit,*) 'para_OTF%multiplicity',para_OTF%multiplicity

      write(out_unit,*) 'para_OTF%ab_initio_meth',trim(adjustl(para_OTF%ab_initio_meth))
      write(out_unit,*) 'para_OTF%ab_initio_basis',trim(adjustl(para_OTF%ab_initio_basis))
      write(out_unit,*) 'para_OTF%ab_initio_prog',trim(adjustl(para_OTF%ab_initio_prog))
      write(out_unit,*) 'para_OTF%header',para_OTF%header
      write(out_unit,*) 'para_OTF%footer',para_OTF%footer
      write(out_unit,*) 'para_OTF%file_name ',trim(adjustl(para_OTF%file_name))
      write(out_unit,*) 'para_OTF%commande_unix',trim(adjustl(para_OTF%commande_unix))

      write(out_unit,*) 'para_OTF%file_header'
      CALL file_Write(para_OTF%file_header)

      write(out_unit,*) 'para_OTF%file_footer'
      CALL file_Write(para_OTF%file_footer)

      write(out_unit,*) 'para_OTF%file_data'
      CALL file_Write(para_OTF%file_data)

      write(out_unit,*) 'para_OTF%file_log'
      CALL file_Write(para_OTF%file_log)

      write(out_unit,*) 'para_OTF%file_FChk'
      CALL file_Write(para_OTF%file_FChk)

      write(out_unit,*) 'para_OTF%file_pun'
      CALL file_Write(para_OTF%file_pun)

      write(out_unit,*) 'END write_OTF'
      flush(out_unit)

      END SUBROUTINE write_OTF

      SUBROUTINE OTF2_TO_OTF1(para_OTF1,para_OTF2)

      CLASS (param_OTF), intent(inout)    :: para_OTF1
      TYPE (param_OTF),  intent(in)       :: para_OTF2

       para_OTF1%charge          = para_OTF2%charge
       para_OTF1%multiplicity    = para_OTF2%multiplicity

       para_OTF1%ab_initio_prog  = para_OTF2%ab_initio_prog
       para_OTF1%ab_initio_meth  = para_OTF2%ab_initio_meth
       para_OTF1%ab_initio_basis = para_OTF2%ab_initio_basis
       para_OTF1%commande_unix   = para_OTF2%commande_unix
       para_OTF1%header          = para_OTF2%header
       para_OTF1%footer          = para_OTF2%footer

       para_OTF1%file_name       = para_OTF2%file_name
       para_OTF1%file_header     = para_OTF2%file_header
       para_OTF1%file_footer     = para_OTF2%file_footer
       para_OTF1%file_data       = para_OTF2%file_data
       para_OTF1%file_log        = para_OTF2%file_log
       para_OTF1%file_FChk       = para_OTF2%file_FChk
       para_OTF1%file_pun        = para_OTF2%file_pun

      END SUBROUTINE OTF2_TO_OTF1
      SUBROUTINE dealloc_OTF(para_OTF)

      CLASS (param_OTF), intent(inout)    :: para_OTF

       para_OTF%charge          = 0
       para_OTF%multiplicity    = -1

       para_OTF%ab_initio_prog  = ''
       para_OTF%ab_initio_meth  = ''
       para_OTF%ab_initio_basis = ''
       para_OTF%commande_unix   = ''
       para_OTF%header          = .FALSE.
       para_OTF%footer          = .FALSE.

       para_OTF%file_name       = ''

       CALL file_dealloc(para_OTF%file_header)
       CALL file_dealloc(para_OTF%file_footer)
       CALL file_dealloc(para_OTF%file_data)
       CALL file_dealloc(para_OTF%file_log)
       CALL file_dealloc(para_OTF%file_FChk)
       CALL file_dealloc(para_OTF%file_pun)

      END SUBROUTINE dealloc_OTF
  END MODULE mod_OTF_def

