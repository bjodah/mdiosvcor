!******************************************************************************
!*                    Code generated with sympy 0.7.1.rc1                     *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                      This file is part of 'autowrap'                       *
!******************************************************************************

REAL*8 function autofunc(DUMMY, P, T)
implicit none
REAL*8, intent(in) :: DUMMY
REAL*8, intent(in) :: P
REAL*8, intent(in) :: T

autofunc = -0.0031055900621118d0*DUMMY*(7.50574125177766d-19*DUMMY**17* &
      exp(-9.30200381340186d-11*DUMMY**4)/T**10 + 9.77191149429517d-45* &
      DUMMY**15*exp(-0.0031055900621118d0*DUMMY)/T - &
      4.71983325174457d-41*DUMMY**14*exp(-0.0031055900621118d0*DUMMY)/T &
      + 2.51880416167177d-15*DUMMY**13*exp(-9.64468963388758d-6*DUMMY** &
      2)/T**8 - 2.82413283290351d-8*DUMMY**13*exp(-9.30200381340186d-11 &
      *DUMMY**4)/T**10 + 1.15537115850702d-8*DUMMY**13*exp( &
      -0.0031055900621118d0*DUMMY)/T**13 - 4.8363836695104d-5*DUMMY**12 &
      *exp(-0.0031055900621118d0*DUMMY)/T**13 - 1.66237760020362d-23* &
      DUMMY**11*exp(-0.0031055900621118d0*DUMMY)/T**4 - &
      7.59437221367534d-14*DUMMY**11*exp(-9.64468963388758d-6*DUMMY**2) &
      /T**6 - 1.56695814419266d-9*DUMMY**11*exp(-9.64468963388758d-6* &
      DUMMY**2)/T**8 - 1.636512053919d-6*DUMMY**11*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**9 + 1.49052061805185d+96*DUMMY &
      **11*exp(-8.97149397534996d-16*DUMMY**6)/T**44 - &
      9.96362226323932d+101*DUMMY**11*exp(-8.97149397534996d-16*DUMMY** &
      6)/T**46 + 6.50979775821756d+112*DUMMY**11*exp( &
      -8.97149397534996d-16*DUMMY**6)/T**50 - 8.99959902004564d-25* &
      DUMMY**10*exp(-9.64468963388758d-6*DUMMY**2)/T + &
      2.03043816323944d-21*DUMMY**10*exp(-9.64468963388758d-6*DUMMY**2) &
      /T**2 - 1.74636261650203d-18*DUMMY**10*exp(-9.64468963388758d-6* &
      DUMMY**2)/T**3 + 5.96812100205062d-16*DUMMY**10*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**4 + 5.8881414599212d-20*DUMMY** &
      10*exp(-0.0031055900621118d0*DUMMY)/T**4 + 8.49425442369174d-6* &
      DUMMY**10*exp(-9.64468963388758d-6*DUMMY**2)/T**8 - &
      0.000802098216797691d0*DUMMY**10*exp(-0.0031055900621118d0*DUMMY) &
      /T**11 + 7.35423360355437d-17*DUMMY**9*exp(-0.0031055900621118d0* &
      DUMMY)/T**4 + 3.93707444301357d-8*DUMMY**9*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**6 + 0.848400578992687d0*DUMMY** &
      9*exp(-9.64468963388758d-6*DUMMY**2)/T**9 + 2.58275625808856d0* &
      DUMMY**9*exp(-0.0031055900621118d0*DUMMY)/T**11 + &
      4.19901491157486d-19*DUMMY**8*exp(-9.64468963388758d-6*DUMMY**2)/ &
      T - 9.47357777327933d-16*DUMMY**8*exp(-9.64468963388758d-6*DUMMY &
      **2)/T**2 + 8.14814376882287d-13*DUMMY**8*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**3 - 2.78459396089478d-10*DUMMY &
      **8*exp(-9.64468963388758d-6*DUMMY**2)/T**4 - &
      2.13125689831005d-13*DUMMY**8*exp(-0.0031055900621118d0*DUMMY)/T &
      **4 - 3.96323224049725d0*DUMMY**8*exp(-9.64468963388758d-6*DUMMY &
      **2)/T**8 - 3159117.90507434d0*DUMMY**8*exp(-9.64468963388758d-6* &
      DUMMY**2)/T**10 + 1.02254379546168d+112*DUMMY**8*exp( &
      -8.97149397534996d-16*DUMMY**6)/T**50 - 4.24663692802952d-13* &
      DUMMY**7*exp(-0.0031055900621118d0*DUMMY)/T**3 - &
      0.00839840035799165d0*DUMMY**7*exp(-9.64468963388758d-6*DUMMY**2) &
      /T**6 + 9781982838.47122d0*DUMMY**7*exp(-9.64468963388758d-6* &
      DUMMY**2)/T**10 - 8.42470937778289d+45*DUMMY**7*exp( &
      -2.9952452279154d-8*DUMMY**3)/T**23 + 9.57191963577854d-10*DUMMY &
      **6*exp(-0.0031055900621118d0*DUMMY)/T**3 - 9896303873621.27d0* &
      DUMMY**6*exp(-9.64468963388758d-6*DUMMY**2)/T**10 - &
      6.53568956117047d+45*DUMMY**6*exp(-2.9952452279154d-8*DUMMY**3)/T &
      **22 + 9.28027724040571d+48*DUMMY**6*exp(-2.9952452279154d-8* &
      DUMMY**3)/T**23 - 7.82273680553224d-7*DUMMY**5*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**3 + 2612.33922815402d0*DUMMY**5 &
      *exp(-9.64468963388758d-6*DUMMY**2)/T**6 + 860775.117270037d0* &
      DUMMY**5*exp(-9.64468963388758d-6*DUMMY**2)/T**7 + &
      43643296.5373839d0*DUMMY**5*exp(-0.0031055900621118d0*DUMMY)/T**9 &
      + 2.48467953918655d+15*DUMMY**5*exp(-9.64468963388758d-6*DUMMY**2 &
      )/T**10 - 3.57242757573539d+31*DUMMY**5*exp(-2.9952452279154d-8* &
      DUMMY**3)/T**16 - 1.66139622023622d+111*DUMMY**5*exp( &
      -8.97149397534996d-16*DUMMY**6)/T**44 + 1.11058674180859d+117* &
      DUMMY**5*exp(-8.97149397534996d-16*DUMMY**6)/T**46 - &
      7.25609109932399d+127*DUMMY**5*exp(-8.97149397534996d-16*DUMMY**6 &
      )/T**50 + 1.56163899349298d-6*DUMMY**4*exp(-0.0031055900621118d0* &
      DUMMY)/T**2 - 70265707425.1881d0*DUMMY**4*exp( &
      -0.0031055900621118d0*DUMMY)/T**9 + 1.47132628100921d+18*DUMMY**4 &
      *exp(-9.64468963388758d-6*DUMMY**2)/T**10 - 1.27656210848006d+20* &
      DUMMY**4*exp(-0.0031055900621118d0*DUMMY)/T**13 + &
      4.68782394357642d+53*DUMMY**4*exp(-2.9952452279154d-8*DUMMY**3)/T &
      **23 - 9.37699268376458d-7*DUMMY**3*(-0.124223602484472d0*DUMMY + &
      40)*exp(-150*(-1.21d0 + 647.096d0/T)**2 - 20*( &
      0.0031055900621118d0*DUMMY - 1)**2) + 0.000611430825150138d0* &
      DUMMY**3*(-0.124223602484472d0*DUMMY + 40)*exp(-150*(-1.21d0 + &
      647.096d0/T)**2 - 20*(0.0031055900621118d0*DUMMY - 1)**2)/T + &
      6.82944170646771d-7*DUMMY**3/T - 6.82695784551734d-7*DUMMY**3*exp &
      (-9.64468963388758d-6*DUMMY**2)/T - 0.00201139102361896d0*DUMMY** &
      3*exp(-0.0031055900621118d0*DUMMY)/T**2 + 0.162218528588961d0* &
      DUMMY**3*exp(-9.64468963388758d-6*DUMMY**2)/T**3 - &
      13241415.9529575d0*DUMMY**3*(-0.124223602484472d0*DUMMY + 40)*exp &
      (-250*(-1.25d0 + 647.096d0/T)**2 - 20*(0.0031055900621118d0*DUMMY &
      - 1)**2)/T**4 - 844.218129675036d0*DUMMY**3*exp( &
      -0.0031055900621118d0*DUMMY)/T**4 - 178497214518.053d0*DUMMY**3* &
      exp(-9.64468963388758d-6*DUMMY**2)/T**7 - 2.63758990189279d+17* &
      DUMMY**3*exp(-9.64468963388758d-6*DUMMY**2)/T**9 - &
      8.35825213726398d+20*DUMMY**3*exp(-9.64468963388758d-6*DUMMY**2)/ &
      T**10 + 1.64421199572232d+23*DUMMY**3*exp(-0.0031055900621118d0* &
      DUMMY)/T**13 + 2.90936203386998d+53*DUMMY**3*exp( &
      -2.9952452279154d-8*DUMMY**3)/T**22 - 4.13111516609254d+56*DUMMY &
      **3*exp(-2.9952452279154d-8*DUMMY**3)/T**23 - 2.26263700343296d-7 &
      *DUMMY**2*(647.096d0/T)**0.375d0 - 0.000905817493251659d0*DUMMY** &
      2*exp(-150*(-1.21d0 + 647.096d0/T)**2 - 20*(0.0031055900621118d0* &
      DUMMY - 1)**2) + 0.590642177095033d0*DUMMY**2*exp(-150*(-1.21d0 + &
      647.096d0/T)**2 - 20*(0.0031055900621118d0*DUMMY - 1)**2)/T + &
      0.0012003218364612d0*DUMMY**2*exp(-0.0031055900621118d0*DUMMY)/T &
      - 12791207810.5569d0*DUMMY**2*exp(-250*(-1.25d0 + 647.096d0/T)**2 &
      - 20*(0.0031055900621118d0*DUMMY - 1)**2)/T**4 + &
      815514.713266085d0*DUMMY**2*exp(-0.0031055900621118d0*DUMMY)/T**4 &
      + 281331099.548702d0*DUMMY**2*exp(-0.0031055900621118d0*DUMMY)/T &
      **5 + 98915788009978.5d0*DUMMY**2*exp(-9.64468963388758d-6*DUMMY &
      **2)/T**7 + 2.16345912482069d+23*DUMMY**2*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**10 + 1.1926995300554d+39*DUMMY &
      **2*exp(-2.9952452279154d-8*DUMMY**3)/T**16 - &
      5.69884903379088d+126*DUMMY**2*exp(-8.97149397534996d-16*DUMMY**6 &
      )/T**50 + 0.00197531113946696d0*DUMMY*sqrt(647.096d0/T) - &
      0.0016239462024446d0*DUMMY*(647.096d0/T)**0.75d0 + &
      0.000987767418585217d0*DUMMY*(-0.198757763975155d0*DUMMY + 64.0d0 &
      )*((0.32d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      1.66666666666667d0 + 1 - 647.096d0/T)**2 + 0.2d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**3.5d0)**0.95d0*exp(-800.0d0* &
      (-1 + 647.096d0/T)**2 - 32.0d0*(0.0031055900621118d0*DUMMY - 1)** &
      2) - 0.000461945368221242d0*DUMMY*(-0.173913043478261d0*DUMMY + &
      56.0d0)*((0.32d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      1.66666666666667d0 + 1 - 647.096d0/T)**2 + 0.2d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**3.5d0)**0.85d0*exp(-700.0d0* &
      (-1 + 647.096d0/T)**2 - 28.0d0*(0.0031055900621118d0*DUMMY - 1)** &
      2) - 0.000461945368221242d0*DUMMY*(0.544d0*(0.010351966873706d0* &
      DUMMY - 3.33333333333333d0)*(0.32d0*((0.0031055900621118d0*DUMMY &
      - 1)**2)**1.66666666666667d0 + 1 - 647.096d0/T)*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0/( &
      0.0031055900621118d0*DUMMY - 1)**2 + 0.17d0*(0.0217391304347826d0 &
      *DUMMY - 7.0d0)*((0.0031055900621118d0*DUMMY - 1)**2)**3.5d0/( &
      0.0031055900621118d0*DUMMY - 1)**2)*((0.32d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      3.5d0)**0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 28.0d0*( &
      0.0031055900621118d0*DUMMY - 1)**2)/((0.32d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      3.5d0) + 0.000987767418585217d0*DUMMY*(0.608d0*( &
      0.010351966873706d0*DUMMY - 3.33333333333333d0)*(0.32d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)*((0.0031055900621118d0*DUMMY - 1)**2)** &
      1.66666666666667d0/(0.0031055900621118d0*DUMMY - 1)**2 + 0.19d0*( &
      0.0217391304347826d0*DUMMY - 7.0d0)*((0.0031055900621118d0*DUMMY &
      - 1)**2)**3.5d0/(0.0031055900621118d0*DUMMY - 1)**2)*((0.32d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      3.5d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - 32.0d0*( &
      0.0031055900621118d0*DUMMY - 1)**2)/((0.32d0*(( &
      0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      3.5d0) + 0.070784629725462d0*DUMMY*exp(-9.64468963388758d-6*DUMMY &
      **2)/T - 0.77300726268101d0*DUMMY*exp(-0.0031055900621118d0*DUMMY &
      )/T + 364051727.461536d0*DUMMY*exp(-0.0031055900621118d0*DUMMY)/T &
      **4 - 181177228109.364d0*DUMMY*exp(-0.0031055900621118d0*DUMMY)/T &
      **5 - 46591403258930.8d0*DUMMY*exp(-0.0031055900621118d0*DUMMY)/T &
      **6 + 2.73475871387852d+22*DUMMY*exp(-9.64468963388758d-6*DUMMY** &
      2)/T**9 - 3.2181183921207d+25*DUMMY*exp(-9.64468963388758d-6* &
      DUMMY**2)/T**10 + 1.10843284498229d+27*DUMMY*exp( &
      -0.0031055900621118d0*DUMMY)/T**12 + 0.012533547935523d0*( &
      647.096d0/T)**-0.5d0 + 7.8957634722828d0*(647.096d0/T)**0.875d0 - &
      8.7803203303561d0*(647.096d0/T)**1.0d0 - 0.14874640856724d0*(( &
      0.32d0*((0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 &
      + 1 - 647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)** &
      2)**3.5d0)**0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 28.0d0*( &
      0.0031055900621118d0*DUMMY - 1)**2) + 0.31806110878444d0*((0.32d0 &
      *((0.0031055900621118d0*DUMMY - 1)**2)**1.66666666666667d0 + 1 - &
      647.096d0/T)**2 + 0.2d0*((0.0031055900621118d0*DUMMY - 1)**2)** &
      3.5d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - 32.0d0*( &
      0.0031055900621118d0*DUMMY - 1)**2) - 117224656242.615d0*exp( &
      -0.0031055900621118d0*DUMMY)/T**4 + 1.50024318493757d+16*exp( &
      -0.0031055900621118d0*DUMMY)/T**6 - 5.1279922820133d+18*exp( &
      -9.64468963388758d-6*DUMMY**2)/T**7 - 3.56915376084297d+29*exp( &
      -0.0031055900621118d0*DUMMY)/T**12) - 1 + 0.00216676249173786d0*P &
      /(DUMMY*T)

end function
