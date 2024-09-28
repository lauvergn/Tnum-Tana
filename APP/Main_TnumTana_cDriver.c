/*==========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!===============================================================================
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
!===========================================================================*/
#include <stdio.h>
#include <math.h>


void Qact_TO_Qcart_TnumTanaDriver_FOR_c(double *, int *, double *, int *);
void Qcart_TO_Qact_TnumTanaDriver_FOR_c(double *, int *, double *, int *);
void get_Qact0_TnumTanaDriver_FOR_c(double *, int *);
void Init_TnumTana_FOR_Driver_FOR_c(int *, int *, int *);

int main(void)
{
  int        init,nb_act,nb_cart;
  double     Qact[6];
  double     Qcart[12];

  int        j;


  init = 0;
  Init_TnumTana_FOR_Driver_FOR_c(&nb_act,&nb_cart,&init);

  printf("nb_act = %i\n",nb_act);
  printf("nb_cart = %i\n",nb_cart);
  printf("---------------------------------------------\n");
  printf("---------------------------------------------\n");

  get_Qact0_TnumTanaDriver_FOR_c(Qact,&nb_act);

  for (j=0;j<nb_act;j++)
    printf(" Qact0[:] = %30.22le \n",Qact[j]);

  printf("---------------------------------------------\n");
  printf("---------------------------------------------\n");


  Qact[0]=2.2791110577179996;
  Qact[1]=2.0824950811500056;
  Qact[2]=2.1237280405061139;
  Qact[3]=2.0824950811500056;
  Qact[4]=2.1237280405061139;
  Qact[5]=3.1415926535897931;

  for (j=0;j<nb_act;j++)
    printf(" Qact[:] = %30.22le \n",Qact[j]);

  printf("---------------------------------------------\n");
  printf("---------------------------------------------\n");


   Qact_TO_Qcart_TnumTanaDriver_FOR_c(Qact,&nb_act,Qcart,&nb_cart);

  for (j=0;j<nb_cart;j=j+3)
    printf(" X %30.22le %30.22le %30.22le \n",Qcart[j+0],Qcart[j+2],Qcart[j+2]);

  printf("---------------------------------------------\n");
  printf("---------------------------------------------\n");

   Qcart_TO_Qact_TnumTanaDriver_FOR_c(Qact,&nb_act,Qcart,&nb_cart);

  for (j=0;j<nb_act;j++)
    printf(" Qact[:] = %30.22le \n",Qact[j]);

  printf("---------------------------------------------");
  printf("---------------------------------------------");

}
