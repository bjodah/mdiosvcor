!******************************************************************************
!*                      Code generated with sympy 0.7.1                       *
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

autofunc = -2.17961356419406d-9*DUMMY*P/T + 0.00139707852303047d0*DUMMY* &
      (647.096d0/T)**0.375d0 - 0.00611639528321138d0*DUMMY*sqrt( &
      647.096d0/T) + 0.00502842144428018d0*DUMMY*(647.096d0/T)**0.75d0 &
      - 0.518538380298064d0*DUMMY*(149.711627154954d0 - &
      6600.28910432124d0/T)*((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - &
      123.059491404039d0)/((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0) + 1.26761291979494d0*DUMMY*( &
      167.324759761419d0 - 7376.79370482962d0/T)*((4.77350664551323d0 - &
      647.096d0/T)**2 + 35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + &
      647.096d0/T)**2 - 140.639418747473d0)/((4.77350664551323d0 - &
      647.096d0/T)**2 + 35.5942454927211d0) + 60.4602391362888d0*DUMMY* &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 123.059491404039d0) - &
      169.060355714069d0*DUMMY*((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - &
      140.639418747473d0) - 0.0031055900621118d0*DUMMY*( &
      0.012533547935523d0*(647.096d0/T)**-0.5d0 - 0.224929642207906d0*( &
      647.096d0/T)**0.375d0 + 1.96947928119407d0*sqrt(647.096d0/T) - &
      1.61915170505822d0*(647.096d0/T)**0.75d0 + 7.8957634722828d0*( &
      647.096d0/T)**0.875d0 - 8.7803203303561d0*(647.096d0/T)**1.0d0 - &
      0.46058153248749d0*(149.711627154954d0 - 6600.28910432124d0/T)*(( &
      4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 123.059491404039d0)/ &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0) + &
      0.984851159229072d0*(167.324759761419d0 - 7376.79370482962d0/T)* &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - 140.639418747473d0)/ &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0) + &
      53.9233353038905d0*((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - &
      123.059491404039d0) - 131.820361050913d0*((4.77350664551323d0 - &
      647.096d0/T)**2 + 35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + &
      647.096d0/T)**2 - 140.639418747473d0) + 77037.6245166147d0*exp( &
      -150*(-1.21d0 + 647.096d0/T)**2 - 87.8996367171704d0) - &
      50232713.1035893d0*exp(-150*(-1.21d0 + 647.096d0/T)**2 - &
      87.8996367171704d0)/T + 662.561833699109d0/T + 51356.2794065008d0 &
      /T**2 - 38016624.6826219d0/T**3 + 1.08786181738036d+18*exp(-250*( &
      -1.25d0 + 647.096d0/T)**2 - 87.8996367171704d0)/T**4 + &
      26696814223.2155d0/T**4 + 4477218632261.18d0/T**5 - &
      4.22110838645506d+15/T**6 + 5.24079092675672d+16/T**7 + &
      3.62208019609034d+20/T**8 - 6.92177444909024d+22/T**9 + &
      2.11381160675712d+24/T**10 + 7.84950375199669d+25/T**11 + &
      3.3828688472416d+28/T**12 + 5.66937219839112d+28/T**13 - &
      4.34863326759208d+33/T**16 - 7.84008793023439d+50/T**22 + &
      1.17461541620322d+53/T**23 + 2.41487073879849d-254/T**44 - &
      1.61425877405078d-248/T**46 + 1.05468652589692d-237/T**50) + &
      61043.9831933747d0*DUMMY*exp(-150*(-1.21d0 + 647.096d0/T)**2 - &
      87.8996367171704d0) - 39803990.7602783d0*DUMMY*exp(-150*(-1.21d0 &
      + 647.096d0/T)**2 - 87.8996367171704d0)/T - 7.02934873970537d0* &
      DUMMY/T + 1426.86725389069d0*DUMMY/T**2 - 1576003.96153247d0* &
      DUMMY/T**3 + 8.62012801063964d+17*DUMMY*exp(-250*(-1.25d0 + &
      647.096d0/T)**2 - 87.8996367171704d0)/T**4 + 700632183.499809d0* &
      DUMMY/T**4 - 10118170751.9452d0*DUMMY/T**5 - 65615594990210.5d0* &
      DUMMY/T**6 + 2.28920498567896d+15*DUMMY/T**7 + &
      7.4479517812196d+18*DUMMY/T**8 - 1.76524587558945d+21*DUMMY/T**9 &
      + 1.20500931343597d+23*DUMMY/T**10 - 1.32979728779013d+24*DUMMY/T &
      **11 + 1.7013285960784d+26*DUMMY/T**12 + 6.11042716176489d+28* &
      DUMMY/T**13 - 1.13387577850419d+33*DUMMY/T**16 - &
      2.01901148593555d+50*DUMMY/T**22 + 3.34572211489632d+52*DUMMY/T** &
      23 + 3.95771564239841d-253*DUMMY/T**44 - 2.64559800170437d-247* &
      DUMMY/T**46 + 1.72851875436044d-236*DUMMY/T**50 + &
      1.42615442020995d0*(149.711627154954d0 - 6600.28910432124d0/T)*( &
      -0.0372665720468512d0*DUMMY*(4.77350664551323d0 - 647.096d0/T) - &
      0.369099466531815d0*DUMMY)*((4.77350664551323d0 - 647.096d0/T)**2 &
      + 35.5942454927211d0)**0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 &
      - 123.059491404039d0)/((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**2 + 1.42615442020995d0*(149.711627154954d0 - &
      6600.28910432124d0/T)*(0.0316765862398235d0*DUMMY*( &
      4.77350664551323d0 - 647.096d0/T) + 0.313734546552043d0*DUMMY)*(( &
      4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 123.059491404039d0)/ &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)**2 - &
      3.04951400547434d0*(167.324759761419d0 - 7376.79370482962d0/T)*( &
      -0.0372665720468512d0*DUMMY*(4.77350664551323d0 - 647.096d0/T) - &
      0.369099466531815d0*DUMMY)*((4.77350664551323d0 - 647.096d0/T)**2 &
      + 35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 &
      - 140.639418747473d0)/((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**2 - 3.04951400547434d0*(167.324759761419d0 - &
      7376.79370482962d0/T)*(0.0354032434445086d0*DUMMY*( &
      4.77350664551323d0 - 647.096d0/T) + 0.350644493205224d0*DUMMY)*(( &
      4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - 140.639418747473d0)/ &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)**2 - &
      166.969358455976d0*(0.0316765862398235d0*DUMMY*( &
      4.77350664551323d0 - 647.096d0/T) + 0.313734546552043d0*DUMMY)*(( &
      4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0)** &
      0.85d0*exp(-700.0d0*(-1 + 647.096d0/T)**2 - 123.059491404039d0)/ &
      ((4.77350664551323d0 - 647.096d0/T)**2 + 35.5942454927211d0) + &
      1.42615442020995d0*(0.035256290415342d0*DUMMY*(4.77350664551323d0 &
      - 647.096d0/T) + 1.08797156062459d0*DUMMY)*((4.77350664551323d0 - &
      647.096d0/T)**2 + 35.5942454927211d0)**0.85d0*exp(-700.0d0*(-1 + &
      647.096d0/T)**2 - 123.059491404039d0)/((4.77350664551323d0 - &
      647.096d0/T)**2 + 35.5942454927211d0) + 408.171360173971d0*( &
      0.0354032434445086d0*DUMMY*(4.77350664551323d0 - 647.096d0/T) + &
      0.350644493205224d0*DUMMY)*((4.77350664551323d0 - 647.096d0/T)**2 &
      + 35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 &
      - 140.639418747473d0)/((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0) - 3.04951400547434d0*(0.0394040892877351d0* &
      DUMMY*(4.77350664551323d0 - 647.096d0/T) + 1.21596821481572d0* &
      DUMMY)*((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0)**0.95d0*exp(-800.0d0*(-1 + 647.096d0/T)**2 - &
      140.639418747473d0)/((4.77350664551323d0 - 647.096d0/T)**2 + &
      35.5942454927211d0) + 2.17317852864848d-6/T

end function
