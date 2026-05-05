c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(13,9) = 0
         pd(18,21) = -km*kmax4*Ya-1/delt
         pd(19,11) = 0
         pd(4,22) = 0
         pd(18,4) = 1/fcn*faa/foa*ktrans*(1-swAquaDiva)+1/fcn*faa/foa/de
     +lt
         pd(21,3) = 0
         pd(11,18) = 1/fcn*faa*foo/fdoco/foa*kpd-1/delt
         pd(16,23) = 0
         pd(17,6) = 0
         pd(5,5) = -fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)
     ++1)**2*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kind
     +ox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3
     +,j))+khyd*sp(3,j))/st/amming*exp((amming-sp(5,j))/st/amming)-fcn*k
     +hyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j
     +))/st/amming)+1)**2*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((d
     +oxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khy
     +d*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))/st/amming*exp((amming-sp(5
     +,j))/st/amming)-fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**
     +2/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(2,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swA
     +quaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))/st/am
     +ming*exp((amming-sp(5,j))/st/amming)-fcn*kmax5/(ksamm+sp(5,j))/(ex
     +p(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)+fcn*kmax5*sp(5,j)/(ks
     +amm+sp(5,j))**2/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)-0.
     +1D2*fcn*kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(
     +5,j))/st/amming)+1)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/
     +st/amming)-1/delt
         pd(7,19) = 0
         pd(13,3) = 2*khyd*sp(3,j)/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +s*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4
     ++sp(15,j))/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)
     ++1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp
     +(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))*(-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*km
     +ax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(
     +13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**
     +2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*km
     +ax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))+khyd)+kdeac*sp(14,j)/(exp((so4min-kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(1,12) = 0
         pd(4,16) = 0
         pd(5,23) = -0.1D-1*fcn*km*kmax4*Ya-fcn/delt
         pd(17,18) = 0
         pd(11,13) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(12,19) = 0
         pd(16,17) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+kreac/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+0.1
     +D-1*km*kmax3*Ys+1/delt
         pd(20,13) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(23,5) = kdeac*sp(20,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,
     +j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,
     +j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kdeac*sp(21,j)/(exp((ammin-
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1
     +)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)
     +/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kr
     +eac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +oxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(kso
     +xa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j
     +)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +xa+sp(4,j)))/st/ammin)+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))
     +**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(4,13) = 0
         pd(13,21) = 0
         pd(4,4) = -foo/fdoco*khyd*sp(3,j)**3/(exp((no3min-kmax2*kindox/
     +(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)
     +))/st/no3min)+1)**2*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**3
     +*sp(10,j)**2/(ksno3+sp(10,j))**2/(ksndoc+sp(3,j))**2/(1/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,
     +j)/fdocn)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-foo/fdo
     +co*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax
     +2*sp(8,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))/(ks
     +ndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*
     +sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no
     +3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,
     +j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(
     +3,j)/fdocn)**2*(-1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2)
     +*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**3*sp(10,j)**2/(ksno3
     ++sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/st/no3min*exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)-1/(exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j)))-foo/fdoco*khyd*sp(3,j)**3/(exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2**2*sp(9,j)*kindox**2
     +/(kindox+sp(4,j))**3*sp(10,j)**2/(ksno3+sp(10,j))**2/(ksndoc+sp(3,
     +j))**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j))+khyd*sp(3,j)/fdocn)/st/no3min*exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/
     +(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco
     +*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*
     +sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+
     +sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,
     +j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j))+khyd*sp(3,j)/fdocn)**2*(-1/((exp((no3min-kmax2*kindox/
     +(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)
     +))/st/no3min)+1)**2)*kmax2**2*sp(9,j)*kindox**2/(kindox+sp(4,j))**
     +3*sp(10,j)**2/(ksno3+sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/s
     +t/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-1/(exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))**
     +2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
         pd(9,14) = -1/delt
         pd(19,23) = 0
         pd(18,14) = -1/delt
         pd(21,15) = 0
         pd(8,18) = 0
         pd(12,14) = -km*kmax3*Ys
         pd(16,11) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(13,15) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(13,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(exp((ammi
     +ng-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,
     +j))**2/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khy
     +d*sp(3,j))-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*
     +kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))*(-kmax3/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4m
     +in*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
     +-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(1
     +3,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(13,j)*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+k
     +max3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp
     +(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))+kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+kreac*s
     +p(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso
     +4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min)
         pd(6,15) = 0
         pd(9,19) = 0
         pd(8,13) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(11,23) = -0.1D-1*km*kmax4*Ya-1/delt
         pd(12,8) = -km*kmax2*Yn
         pd(16,5) = 0
         pd(17,12) = 0
         pd(19,4) = 0
         pd(5,18) = (fcn*fdoco*foa+faa*foo)/fdoco/foa*kpd
         pd(6,2) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(22,21) = -kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1))
         pd(23,10) = 0
         pd(2,11) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(7,14) = 0
         pd(1,7) = -katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(7,5) = khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp(
     +(amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))/st/amming*ex
     +p((amming-sp(5,j))/st/amming)
         pd(12,2) = -khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Yo*kmax1/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAq
     +uaDiva)*faz)**3*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)+1)*Yo*
     +kmax1**2*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)**2/(ksox+sp(4,j))**2/
     +(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j))**2+1/delt
         pd(21,8) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)*sp(21,j)
         pd(18,16) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(4,7) = 0
         pd(19,16) = 0
         pd(3,15) = 1/fsulph*fdocs/delt
         pd(8,8) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)*sp(9,j)
         pd(3,7) = 0
         pd(5,10) = -fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming
     +)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*s
     +p(3,j))+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yn
     +*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**
     +2/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(
     +exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+k
     +hyd*sp(3,j))+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+
     +1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2
     +/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j))+khyd*sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,
     +j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*s
     +p(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j))+khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
     +/st/no3min)+1)**2)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+s
     +p(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)
     +))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j
     +))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-1/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))-fcn*kdeac*sp(9,
     +j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(k
     +indox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*kre
     +ac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*
     +kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+faa*foo/fdoco/foa/fnitra*fdocn/delt
         pd(14,16) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(23,16) = 0
         pd(6,10) = 0
         pd(11,7) = -1/delt
         pd(3,2) = -fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**
     +2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**3*sp(3,j)**3/(
     +exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1**2*sp(2,j)/(k
     +sodoc+sp(3,j))**2*sp(4,j)**2/(ksox+sp(4,j))**2/(1/(exp((doxmin-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j))**2
         pd(3,10) = 1/fnitra*fdocn/delt
         pd(6,6) = -0.1D-1*km*kmax1*Yo-kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-kl*vel**0.58D0-kdet/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8
     +,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax
     +)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j
     +)-1/delt
         pd(8,2) = 0
         pd(21,20) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)*sp(21,j)
         pd(14,22) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(1,21) = 0
         pd(11,1) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(20,1) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(22,10) = 0
         pd(2,9) = 0
         pd(3,19) = 0
         pd(23,23) = -kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp
     +(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-0.1D-1*km*kmax4*Ya-1/delt
         pd(4,8) = foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,
     +j)**3/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2**2*sp(8,
     +j)*kindox**2/(kindox+sp(4,j))**2*sp(10,j)**2/(ksno3+sp(10,j))**2/(
     +ksndoc+sp(3,j))**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*k
     +max2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)**2
         pd(10,21) = 0
         pd(13,22) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(9,2) = -1/delt
         pd(14,15) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(14,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(exp((ammi
     +ng-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,
     +j))**2/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khy
     +d*sp(3,j))-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*
     +kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))*(-kmax3/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4m
     +in*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
     +-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(1
     +4,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(14,j)*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+k
     +max3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp
     +(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))-kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-kreac*s
     +p(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso
     +4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min)
         pd(22,16) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(6,19) = 0
         pd(12,1) = -khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Yo*kmax1/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAq
     +uaDiva)*faz)**3*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)+1)*Yo*
     +kmax1**2*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)**2/(ksox+sp(4,j))**2/
     +(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j))**2+1/delt
         pd(22,3) = 0
         pd(2,1) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(6,1) = -km*kmax1*Yo+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(
     +3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0
     ++(1.D0-swAquaDiva)*faz)**3*sp(3,j)**3/(exp((amming-sp(5,j))/st/amm
     +ing)+1)*Yo*kmax1**2*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)**2/(ksox+s
     +p(4,j))**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+
     +khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2-kl*vel**0.58D0-kdet/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8
     +,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax
     +)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j
     +)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)*sp(7,j)-1/delt
         pd(22,22) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-k
     +det*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)*sp(23,j)
         pd(5,14) = -fcn*km*kmax3*Ys-fcn/delt
         pd(14,21) = 0
         pd(23,11) = 0
         pd(1,22) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(15,14) = 0
         pd(20,2) = 0
         pd(22,9) = 0
         pd(7,2) = -km*kmax1*Yo+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(
     +3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0
     ++(1.D0-swAquaDiva)*faz)**3*sp(3,j)**3/(exp((amming-sp(5,j))/st/amm
     +ing)+1)*Yo*kmax1**2*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)**2/(ksox+s
     +p(4,j))**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+
     +khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2-katt*(1-1/(exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)+1))-1/delt
         pd(14,3) = 2*khyd*sp(3,j)/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +s*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4
     ++sp(15,j))/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j))+khyd*sp(3,j))-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)
     ++1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp
     +(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))*(-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)-khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*km
     +ax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**
     +2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*km
     +ax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))+khyd)-kdeac*sp(14,j)/(exp((so4min-kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)-kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(21,21) = Ya*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-kdeac*
     +(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4
     +,j)))/st/ammin)+1))-km*kmax4*Ya-1/delt
         pd(5,9) = -fcn*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-fcn*kdea
     +c*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))-fcn*km*kmax2*Yn
     +-fcn/delt
         pd(8,23) = 0
         pd(23,17) = 0
         pd(1,16) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(11,6) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(22,15) = 0
         pd(8,7) = 0
         pd(2,21) = 0
         pd(10,6) = 0
         pd(4,17) = 0
         pd(17,19) = 0
         pd(15,2) = 0
         pd(20,14) = 0
         pd(5,6) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1
     +)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*kdet*sp(20,j)/(
     +exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(1
     +6,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
         pd(2,4) = -kdeac*sp(2,j)/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)-kreac*sp(7,j)/(exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(18,5) = -1/fcn/delt
         pd(12,23) = 0
         pd(14,10) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min
     +-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(1
     +4,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)*
     +*3/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(14,j)*sp(15,
     +j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp
     +(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-khyd*sp(3,
     +j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,
     +j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(14,j)*sp(15,j)**
     +2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(ki
     +ndox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*
     +sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2)-kdeac*sp(14,j)/(exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2
     +*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-kreac*sp(17,j)
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(16,22) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(17,5) = 0
         pd(20,8) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(4,2) = 0
         pd(13,8) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(10,16) = 0
         pd(15,21) = -fsulph/fdocs/foo*fdoco*foa*kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)
         pd(18,9) = -1/delt
         pd(8,14) = 0
         pd(7,10) = 0
         pd(8,19) = 0
         pd(18,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(9,9) = -kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +))-1/delt
         pd(12,13) = -km*kmax3*Ys
         pd(16,10) = -kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)-kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)
         pd(17,17) = -kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)-0.1D-1*km*kmax3*Ys-1/delt
         pd(20,20) = Ya*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)*sp(21,j)-km*kmax4*Ya-1/delt
         pd(13,20) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(4,3) = 2*foo/fdoco*khyd*sp(3,j)/(exp((no3min-kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
     +/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khy
     +d*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*s
     +p(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+s
     +p(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j
     +)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j))+khyd*sp(3,j)/fdocn)*(-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3m
     +in*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-foo/fdoco*khyd*sp(3,j)**
     +2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/
     +(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))**2/(1/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+
     +khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)**2*(-
     +1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2)*kmax2*sp(8,j)*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(
     +ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kma
     +x2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksnd
     +oc+sp(3,j))-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(
     +8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))**2+khyd/fdocn)+2*foo/fdoco*khyd*sp(3,j)/(exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax
     +2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdoc
     +n)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
     ++1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)*(-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kind
     +ox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3
     +,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-foo/fdoc
     +o*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2
     +*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc
     ++sp(3,j))**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*s
     +p(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3
     +,j)/fdocn)**2*(-1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2)*
     +kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))/(ksndoc+sp(3,j))-1/(exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j))**2+khyd/fdocn)
         pd(15,9) = 0
         pd(5,11) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*kdet*sp(20,j)/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
         pd(14,17) = kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+1)
         pd(1,18) = 0
         pd(2,15) = 0
         pd(6,9) = 0
         pd(15,10) = -fsulph/fnitra*fdocn/fdocs/delt
         pd(22,13) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(8,9) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*
     +Yn*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc
     ++sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,
     +j))+khyd*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2**
     +2*sp(9,j)*kindox**2/(kindox+sp(4,j))**2*sp(10,j)**2/(ksno3+sp(10,j
     +))**2/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j))+khyd*sp(3,j))**2+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+
     +kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))+km*kmax2*Yn
     ++1/delt
         pd(3,6) = 0
         pd(10,4) = -fnitra/fdocn*khyd*sp(3,j)**3/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(13,j)
     +*sp(15,j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(ki
     +ndox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/(1/(exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j
     +)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)/st/so4min*
     +exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fn
     +itra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,
     +j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(
     +10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*
     +sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15
     +,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+khyd*sp(3,j)/fdocs)**2*(-1/((exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(13,j)*sp(15
     +,j)**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**
     +2/(kindox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/st/so4min*exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(ex
     +p((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*k
     +max3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j)))-fnitra/fdocn*k
     +hyd*sp(3,j)**3/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)**2*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+sp(15,j))
     +**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**3*kinno3**2/(k
     +inno3+sp(10,j))**2/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))+khyd*sp(3,j)/fdocs)/st/so4min*exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)-fnitra/fdocn*khyd*sp(3,j)**2/(e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*
     +kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*k
     +hyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kma
     +x3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)*
     +*2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)**2)*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+sp(15,j))**2*sp
     +(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**3*kinno3*
     +*2/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kin
     +no3/(kinno3+sp(10,j)))
         pd(20,19) = 0
         pd(8,20) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(14,23) = 0
         pd(23,13) = 0
         pd(15,16) = 0
         pd(6,7) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(8,3) = -2*khyd*sp(3,j)/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +n*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/
     +(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khy
     +d*sp(3,j))+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yn*
     +kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(k
     +sndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kh
     +yd*sp(3,j))+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yn
     +*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(
     +ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2/(1/(
     +exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+k
     +hyd*sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/a
     +mming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))+khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)**2)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/n
     +o3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))-1/(exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2+khyd)+kdeac*sp(9,j)/(e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/
     +(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)
     +)**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+kreac*sp(12
     +,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
         pd(14,5) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*
     +*2*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/amming*exp((
     +amming-sp(5,j))/st/amming)
         pd(21,23) = kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)
         pd(2,17) = 0
         pd(3,23) = 0
         pd(9,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(1,4) = -kreac*sp(6,j)/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)-kdeac*sp(1,j)/(exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(1,6) = kreac/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+0.
     +1D-1*km*kmax1*Yo+kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+
     +kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)*sp(7,j)+1/delt
         pd(6,21) = 0
         pd(10,17) = 0
         pd(5,17) = -0.1D-1*fcn*km*kmax3*Ys-fcn/delt
         pd(6,3) = 2*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)/(exp((
     +amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/
     +(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-
     +swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0
     +-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3
     +,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*
     +(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(
     +exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((ex
     +p((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*f
     +az)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
     +*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j
     +))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))
         pd(3,17) = 0
         pd(9,20) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kd
     +et*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)-km*kmax4*Ya-1/delt
         pd(10,23) = 0
         pd(19,17) = 0
         pd(21,9) = 0
         pd(3,11) = 0
         pd(9,4) = -kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)-kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*
     +*2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3
     +,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)-kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*k
     +max2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-
     +kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2
     +*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/fc
     +n*faa/foa*ktrans*(1-swAquaDiva)+1/fcn*faa/foa/delt
         pd(15,1) = 0
         pd(22,18) = 0
         pd(18,13) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(21,14) = 0
         pd(8,17) = 0
         pd(9,7) = -1/delt
         pd(12,15) = 0
         pd(16,12) = 0
         pd(4,14) = 0
         pd(19,22) = 0
         pd(8,12) = -kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
         pd(9,13) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(12,9) = -km*kmax2*Yn
         pd(16,6) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(1,1) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax
     +*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)-kdeac*(1-1/(exp((doxmin-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1))
         pd(5,15) = faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(4,21) = 0
         pd(14,8) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(17,7) = 0
         pd(20,10) = 0
         pd(13,10) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min
     +-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(1
     +3,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)*
     +*3/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(13,j)*sp(15,
     +j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp
     +(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-khyd*sp(3,
     +j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,
     +j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(13,j)*sp(15,j)**
     +2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(ki
     +ndox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*
     +sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2)+kdeac*sp(14,j)/(exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2
     +*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+kreac*sp(17,j)
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(18,22) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-0.1D-1*km*kmax4*Ya-1/delt
         pd(2,2) = -kdeac*(1-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1))
         pd(9,23) = -0.1D-1*km*kmax4*Ya-1/delt
         pd(15,19) = 0
         pd(22,4) = -kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(
     +5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-kreac*sp(23,j)/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+
     +1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(4,15) = 0
         pd(13,23) = 0
         pd(18,1) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(1,5) = 0
         pd(10,1) = 0
         pd(11,14) = -1/delt
         pd(14,14) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**3/(exp((ammi
     +ng-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+
     +sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**2*kin
     +no3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))+khyd*sp(3,j))**2-katt*(1-1/(exp((Bfmax-sp(1,j
     +)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st
     +/Bfmax)+1))-kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j)))/st/so4min)+1))-km*kmax3*Ys-1/delt
         pd(17,1) = 0
         pd(23,4) = kdeac*sp(20,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kdeac*sp(21,j)/(exp((ammin-
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1
     +)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)
     +/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kr
     +eac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +oxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(kso
     +xa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))*
     +*2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +xa+sp(4,j)))/st/ammin)+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))
     +*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(5,4) = fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+
     +1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(
     +10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     +/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j))+khyd*sp(3,j))+fcn*khyd*sp(3,j)**3/(exp((amming-sp(5,j))/st/amm
     +ing)+1)*Yn*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**3*sp(10,j)
     +**2/(ksno3+sp(10,j))**2/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1)**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
     ++1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j))/st/no3min*exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/
     +st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min)+1)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,
     +j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j))+khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)**2)*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**3*sp(1
     +0,j)**2/(ksno3+sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/st/no3m
     +in*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-1/(exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))**2*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))-fcn*kdeac*sp(9,j)/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kindox
     ++sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st
     +/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*kreac*sp(12,j)/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kindox
     ++sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st
     +/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*khyd*(1.D0+(1.D
     +0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+
     +1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-k
     +max1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2
     +/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*kh
     +yd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((do
     +xmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+fcn*khyd*(
     +1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st
     +/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxmin-kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*
     +faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))-1/(
     +exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +))**2)-fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((a
     +mming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))/(ks
     +ox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((
     +doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kh
     +yd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAq
     +uaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*
     +kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((dox
     +min-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,
     +j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st
     +/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(
     +-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox
     ++sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**
     +2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/dox
     +min)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D
     +0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/(
     +(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,
     +j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+s
     +p(3,j))/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))**2)+faa/foa*ktrans*(1-swAquaDiva)+faa
     +/foa/delt
         pd(7,20) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(13,4) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13
     +,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**
     +3/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(13,j)*sp(15,j
     +)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(
     +4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-khyd*sp(3,j
     +)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j
     +)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(13,j)*sp(15,j)**2
     +/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kin
     +dox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*s
     +p(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j)))+kdeac*sp(14,j)/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+kreac*sp(17,j)/
     +(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1
     +)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(23,22) = -kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp
     +(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-0.1D-1*km*kmax4*Ya-1/delt
         pd(7,9) = 0
         pd(1,3) = -kreac*sp(6,j)/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)-kdeac*sp(1,j)/(exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(21,19) = 0
         pd(12,10) = 0
         pd(16,7) = 0
         pd(2,16) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(7,4) = khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp(
     +(amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))/(
     +ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+
     +khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAqua
     +Diva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*km
     +ax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmi
     +n-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/
     +(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(-kmax1
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,
     +j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((am
     +ming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxm
     +in-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-
     +kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+
     +sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))/(k
     +sox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(
     +3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))**2)
         pd(13,19) = 0
         pd(19,21) = 0
         pd(3,16) = 0
         pd(21,13) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)*sp(21,j)
         pd(9,8) = -kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +))-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)-1/delt
         pd(17,16) = -kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)-0.1D-1*km*kmax3*Ys-1/delt
         pd(7,21) = 0
         pd(13,5) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*
     +*2*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/amming*exp((
     +amming-sp(5,j))/st/amming)
         pd(18,17) = -1/delt
         pd(2,3) = -kdeac*sp(2,j)/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)-kreac*sp(7,j)/(exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)**2*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(7,3) = 2*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)/(exp((
     +amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/
     +(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-
     +swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0
     +-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3
     +,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*
     +(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(
     +exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((ex
     +p((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*f
     +az)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
     +*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j
     +))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))
         pd(17,20) = 0
         pd(11,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(17,2) = 0
         pd(20,15) = 0
         pd(23,3) = 0
         pd(7,15) = 0
         pd(15,20) = -fsulph/fdocs/foo*fdoco*foa*kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)
         pd(4,20) = 0
         pd(12,22) = 0
         pd(14,9) = 0
         pd(16,21) = 0
         pd(20,9) = 0
         pd(23,9) = 0
         pd(8,4) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +n*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j
     +))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+
     +khyd*sp(3,j))+khyd*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)+1)*
     +Yn*kmax2**2*sp(9,j)*kindox**2/(kindox+sp(4,j))**3*sp(10,j)**2/(ksn
     +o3+sp(10,j))**2/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax
     +2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j))+khyd*sp(3,j))/st/no3min*exp((no3min-kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j)))/st/no3min)+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+
     +khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**
     +2)*kmax2**2*sp(9,j)*kindox**2/(kindox+sp(4,j))**3*sp(10,j)**2/(ksn
     +o3+sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/st/no3min*exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)-1/(exp((no3min-kmax2*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))+kdeac*sp(9,j)/(exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+kreac*sp(12,j)/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j)))/st/no3min)
         pd(2,5) = 0
         pd(11,22) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-0.1D-1*k
     +m*kmax4*Ya-1/delt
         pd(14,4) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14
     +,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**
     +3/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(14,j)*sp(15,j
     +)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(
     +4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-khyd*sp(3,j
     +)**2/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j
     +)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+khyd*sp(3,j))**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(14,j)*sp(15,j)**2
     +/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kin
     +dox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*s
     +p(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j)))-kdeac*sp(14,j)/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-kreac*sp(17,j)/
     +(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1
     +)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(17,11) = 0
         pd(19,3) = 0
         pd(21,22) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)*sp(21,j)
         pd(13,14) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+kdeac*(1-1/
     +(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1
     +))
         pd(23,12) = 0
         pd(1,23) = 0
         pd(6,14) = 0
         pd(15,15) = -1/delt
         pd(22,8) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(3,5) = 0
         pd(10,5) = 0
         pd(14,18) = 0
         pd(23,18) = 0
         pd(1,17) = 0
         pd(6,8) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(11,5) = -1/fcn/delt
         pd(22,14) = 0
         pd(5,22) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*kdet*sp(20,j)/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-
     +0.1D-1*fcn*km*kmax4*Ya-fcn/delt
         pd(21,10) = 0
         pd(9,3) = -kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/
     +(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)-kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)-kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)-kreac*sp(12,j)/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)-1/fcn*faa*foo/fdoco/foa/delt
         pd(12,18) = kpd+1/delt
         pd(16,16) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-k
     +det*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)*sp(17,j)
         pd(18,8) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(22,17) = 0
         pd(3,18) = kpd
         pd(1,11) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(4,9) = foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,
     +j)**3/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2**2*sp(9,
     +j)*kindox**2/(kindox+sp(4,j))**2*sp(10,j)**2/(ksno3+sp(10,j))**2/(
     +ksndoc+sp(3,j))**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*k
     +max2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)**2
         pd(10,22) = 0
         pd(19,18) = 0
         pd(5,16) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-0.1D-1*fcn*km*kmax
     +3*Ys-fcn*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)-fcn/delt
         pd(21,4) = Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-Ya*sp(21,j)*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/a
     +mmin)-kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+s
     +p(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(22,23) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+kreac/(exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/a
     +mmin)+1)+0.1D-1*km*kmax4*Ya+1/delt
         pd(9,16) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(6,20) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(19,12) = 0
         pd(16,20) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(17,3) = kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +)+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+kdeac*sp(14,
     +j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)+kreac*sp(16,j)/(exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4
     +min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(20,6) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(13,6) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(6,22) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(4,19) = 0
         pd(17,21) = 0
         pd(11,10) = -kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)-kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min)+1/fcn*faa*foo/fdoco/foa/fnitra*fdocn/delt
         pd(16,14) = -kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1))
         pd(22,19) = 0
         pd(2,13) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(12,11) = -0.1D-1*km*kmax2*Yn
         pd(16,8) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(17,15) = kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +)+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+kdeac*sp(14,
     +j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(
     +15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+kreac*sp(16,j)/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4mi
     +n*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j
     +)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(20,18) = 0
         pd(13,18) = 0
         pd(1,19) = 0
         pd(22,12) = 0
         pd(12,5) = -khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))/st/amming*
     +exp((amming-sp(5,j))/st/amming)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)*
     +*2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(2
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j))/st/amming*exp((amming-sp(5,j))/st/amming)-kmax5/(ksamm+sp(5
     +,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)+kmax5*sp(5,j)
     +/(ksamm+sp(5,j))**2/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1
     +)-0.1D2*kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(
     +5,j))/st/amming)+1)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/
     +st/amming)
         pd(16,2) = 0
         pd(17,9) = 0
         pd(7,12) = 0
         pd(13,12) = 0
         pd(23,14) = 0
         pd(7,17) = 0
         pd(1,10) = 0
         pd(7,8) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax
     +*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1
     +)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(10,13) = fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(15,j)/(ksso
     +4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+
     +khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**3/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(
     +13,j)*sp(15,j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**
     +2/(kindox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**2/(1/(exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp
     +(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)**2
         pd(15,18) = -fsulph/fdocs*kpd
         pd(5,21) = -fcn*km*kmax4*Ya-fcn/delt
         pd(8,16) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(2,23) = 0
         pd(14,7) = 0
         pd(23,7) = 0
         pd(3,21) = 0
         pd(18,19) = 0
         pd(10,19) = 0
         pd(19,13) = 0
         pd(8,11) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(21,5) = Ya*sp(21,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*s
     +p(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-Ya*sp(21,j)*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +*2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4
     +*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/a
     +mmin)-kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox
     +a+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(2,6) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)*sp(7,j)
         pd(3,8) = 0
         pd(10,2) = 0
         pd(14,13) = kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kd
     +et*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(14,j)
         pd(15,5) = -fsulph/fdocs/foo*fdoco*foa*sp(20,j)*kmax4/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoc
     +o*foa*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*sp(20,j)*kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kma
     +x4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp
     +(5,j))**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,
     +j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-fsulph/fdocs
     +/foo*fdoco*foa*sp(21,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,
     +j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,
     +j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j
     +)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)+fsulph/
     +fdocs/foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/
     +(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp
     +(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksoxa+sp(4,j)))/st/ammin)
         pd(5,7) = 0
         pd(23,19) = 0
         pd(6,13) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(15,6) = 0
         pd(18,23) = -0.1D-1*km*kmax4*Ya-1/delt
         pd(8,5) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*
     +*2*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1
     +/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     ++khyd*sp(3,j))/st/amming*exp((amming-sp(5,j))/st/amming)
         pd(9,6) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(10,8) = 0
         pd(19,2) = 0
         pd(15,12) = 0
         pd(20,4) = Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-Ya*sp(20,j)*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/a
     +mmin)+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+s
     +p(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(22,7) = 0
         pd(12,17) = -0.1D-1*km*kmax3*Ys
         pd(21,17) = 0
         pd(9,12) = kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
         pd(10,12) = 0
         pd(14,1) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(21,6) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)*sp(21,j)
         pd(14,12) = 0
         pd(23,2) = 0
         pd(3,22) = 0
         pd(7,22) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(9,21) = -km*kmax4*Ya-1/delt
         pd(18,18) = 1/fcn*faa*foo/fdoco/foa*kpd-1/delt
         pd(1,15) = 0
         pd(8,21) = 0
         pd(10,18) = 0
         pd(15,23) = 0
         pd(19,14) = 0
         pd(3,1) = -fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**
     +2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**3*sp(3,j)**3/(
     +exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1**2*sp(1,j)/(k
     +sodoc+sp(3,j))**2*sp(4,j)**2/(ksox+sp(4,j))**2/(1/(exp((doxmin-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j))**2
         pd(5,20) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*(kl*vel**0.58D
     +0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(20,j)/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-fcn*km*kmax4*Ya
     +-fcn/delt
         pd(3,12) = 0
         pd(23,8) = 0
         pd(7,16) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(1,9) = 0
         pd(7,7) = -0.1D-1*km*kmax1*Yo-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1))-1/delt
         pd(12,4) = -khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAq
     +uaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*
     +kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((dox
     +min-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**
     +2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/dox
     +min)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(
     +4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((
     +amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/
     +(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-
     +swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((do
     +xmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*
     +(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(kso
     +x+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/
     +(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))**2)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3
     +,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodo
     +c+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+
     +(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/ammi
     +ng)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(
     +2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2
     +*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))
     +/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)
     +**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1
     +/((exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.
     +D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/dox
     +min)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc
     ++sp(3,j))/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))**2)
         pd(16,1) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(19,8) = 0
         pd(6,4) = khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp(
     +(amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))/(
     +ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+
     +khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAqua
     +Diva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*km
     +ax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmi
     +n-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j))-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/
     +(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(-kmax1
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,
     +j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)-khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((am
     +ming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxm
     +in-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-
     +kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+
     +sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(k
     +sox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(
     +3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))**2)
         pd(21,18) = 0
         pd(9,11) = kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+0.1
     +D-1*km*kmax2*Yn-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)
         pd(19,7) = 0
         pd(4,1) = 0
         pd(6,18) = 0
         pd(15,11) = 0
         pd(2,18) = 0
         pd(12,21) = 0
         pd(21,12) = 0
         pd(10,9) = 0
         pd(19,1) = 0
         pd(5,1) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1
     +)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*khyd*(1.D0+(1.D
     +0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+
     +1)*Yo*kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-k
     +max1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**3*sp(3,j)**3
     +/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1**2*sp(1,j)/(ksodoc+s
     +p(3,j))**2*sp(4,j)**2/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(
     +3,j))**2-fcn*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)
         pd(6,12) = 0
         pd(11,9) = -kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1))-1/delt
         pd(19,20) = 0
         pd(5,3) = -2*fcn*khyd*sp(3,j)/(exp((amming-sp(5,j))/st/amming)+
     +1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1
     +/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     ++khyd*sp(3,j))+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming
     +)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j))+khyd*sp(3,j))+fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/a
     +mming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j
     +)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j))+khyd*sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)+fcn*khyd*sp(3,j)**2/(exp((ammi
     +ng-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1)/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     +*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)**2)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))-1/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2+khyd)-
     +fcn*kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp
     +(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)-fcn*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/
     +(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)-2*fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3
     +,j)/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*kh
     +yd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox
     ++sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((do
     +xmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAqua
     +Diva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*km
     +ax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-k
     +max1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((d
     +oxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+fcn*khyd*(1.D0+(1.D0-swAquaDiva
     +)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j))**2*(-1/((exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2
     +)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*
     +sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D
     +0-swAquaDiva)*faz))-2*fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(
     +3,j)/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swA
     +quaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*k
     +hyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j
     +))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(kso
     +x+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((d
     +oxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khy
     +d*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.D0-swAqu
     +aDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*k
     +max1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-
     +kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.
     +D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((
     +doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+fcn*khyd*(1.D0+(1.D0-swAquaDiv
     +a)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1
     +*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j))**2*(-1/((exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**
     +2)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/dox
     +min*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)
     +*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.
     +D0-swAquaDiva)*faz))-faa*foo/fdoco/foa/delt
         pd(13,1) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(5,12) = fcn*kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
         pd(11,11) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-k
     +det*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(16,15) = -kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-kreac*sp(17
     +,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4mi
     +n)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp
     +(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)
         pd(20,11) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(4,10) = -foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1)**2*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)*(-kmax2
     +*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     ++kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))/(ksnd
     +oc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp
     +(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2/(ksndoc+sp(3,j))/(1/(exp((n
     +o3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3
     +,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp
     +(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)**2*(-1/((exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2)*kmax2*sp(8,j)*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*
     +(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,
     +j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +))-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kind
     +ox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+s
     +p(3,j)))-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)*(-kmax2*kindo
     +x/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+foo
     +/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*
     +kmax2*sp(9,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))/(ksndoc+sp(
     +3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j))+khyd*sp(3,j)/fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3mi
     +n-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))**2/(ksndoc+sp(3,j))/(1/(exp((no3min-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/
     +fdocn)-foo/fdoco*khyd*sp(3,j)**2/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))/(ksndoc+sp(3,j))/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,j)/fdocn)**2*(-1/((exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)**2)*kmax2*sp(9,j)*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax
     +2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)
     +)+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,
     +j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kin
     +dox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-1/(
     +exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)
     +))+foo/fdoco/fnitra*fdocn/delt
         pd(18,2) = -1/delt
         pd(2,14) = 0
         pd(20,5) = Ya*sp(20,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)-Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*s
     +p(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-Ya*sp(20,j)*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +*2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4
     +*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/a
     +mmin)+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox
     +a+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(7,13) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(13,13) = khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Ys*kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)/(1/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j))-khyd*sp(3,j)**3/(exp((ammi
     +ng-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(13,j)*sp(15,j)**2/(ksso4+
     +sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**2*kin
     +no3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)**2/(1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))+khyd*sp(3,j))**2-kl*vel**0.58D0-kdet/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)-kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-km*kmax3*Ys-1/del
     +t
         pd(11,4) = -kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)-kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min)+1/fcn*faa/foa*ktrans*(1-swAquaDiva)+1/fcn*faa/foa/delt
         pd(18,7) = -1/delt
         pd(11,21) = -km*kmax4*Ya-1/delt
         pd(12,6) = 1/delt
         pd(16,3) = -kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-kreac*sp(17
     +,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4mi
     +n)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)
         pd(17,10) = kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+kreac*sp(17,j)/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)
         pd(20,23) = -kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp
     +(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)-0.1D-1*km*kmax4*Ya-1/delt
         pd(14,19) = 0
         pd(2,22) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(18,12) = -0.1D-1*km*kmax2*Yn-1/delt
         pd(20,17) = 0
         pd(3,3) = -2*fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,
     +j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-
     +swAquaDiva)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+fdoco*khyd*(1.D0+(1.D0-
     +swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)*kmax1*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))
     +/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1.D
     +0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1)*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-s
     +wAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodo
     +c+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*
     +faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,
     +j))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))-2
     +*fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)/(exp((doxmin-k
     +max1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*f
     +az)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
     +*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j
     +))+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((dox
     +min-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(2,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-k
     +max1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz
     +)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1
     +*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmi
     +n-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1
     +.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-swAquaD
     +iva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((
     +doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kh
     +yd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +)+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-swA
     +quaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0
     ++(1.D0-swAquaDiva)*faz)*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/
     +(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))-1/delt
         pd(10,7) = 0
         pd(3,4) = fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/
     +(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*khyd*(1.D0+(1.D0-swAqu
     +aDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)*kmax1*sp(1,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(1/(exp((doxmin
     +-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.
     +D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-swAquaDi
     +va)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*
     +faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +)*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(1/(exp
     +((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+
     +khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1.D0+(1.D0-
     +swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)*kmax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxm
     +in-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-
     +kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+
     +sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(k
     +sox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(
     +3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))**2)+fdoco*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*
     +sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*s
     +p(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(1/(exp((doxmin-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp
     +(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(
     +1.D0-swAquaDiva)*faz)*sp(3,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAqua
     +Diva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/
     +doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*khyd*(1.D0+(
     +1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(1/(e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAqu
     +aDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)
     +)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1.D0+(1.D
     +0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st
     +/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*
     +*2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.
     +D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fdoco*khyd*(1
     +.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((doxmin-kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j))/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D
     +0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/(
     +(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,
     +j))/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+s
     +p(3,j))/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))**2)
         pd(14,20) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(23,20) = kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1))
         pd(15,7) = 0
         pd(21,16) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)*sp(21,j)
         pd(10,11) = 0
         pd(14,2) = 0
         pd(19,5) = 0
         pd(6,16) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(15,13) = 0
         pd(22,6) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(8,15) = 0
         pd(18,3) = -1/fcn*faa*foo/fdoco/foa/delt
         pd(21,2) = 0
         pd(4,5) = 0
         pd(23,6) = 0
         pd(7,18) = 0
         pd(10,14) = fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(15,j)/(ksso
     +4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+
     +khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**3/(exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(
     +14,j)*sp(15,j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**
     +2/(kindox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**2/(1/(exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp
     +(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)**2
         pd(19,10) = 0
         pd(8,10) = -khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)
     +*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))/(ksndoc
     ++sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3,
     +j))+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*s
     +p(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2/(ksndo
     +c+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(3
     +,j))+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*
     +sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+
     +sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2/(1/(exp((no
     +3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,
     +j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*sp(
     +3,j))*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j)))/st/no3min)+khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1
     +)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+
     +khyd*sp(3,j))**2*(-1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**
     +2)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3min-kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
     +/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j))-1/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))**2*sp(3,j)/(ksndoc+sp(3,j)))+kdeac*sp(9,j)/(exp((no3min-kmax
     +2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*
     +exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+kreac*sp(12,j)/(exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
         pd(11,19) = 0
         pd(3,9) = 0
         pd(9,1) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(15,4) = -fsulph/fdocs/foo*fdoco*foa*sp(20,j)*kmax4*sp(5,j)/(
     +ksamm+sp(5,j))/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoc
     +o*foa*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j
     +))**2/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(
     +4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*sp(20,j)*kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,
     +j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-fsulph/fdocs
     +/foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksoxa+sp(4,
     +j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,
     +j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))**2/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)+fsulph/
     +fdocs/foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksoxa+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/
     +(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,
     +j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksoxa+sp(4,j)))/st/ammin)-fsulph/fdocs/foo*fdoco*ktrans*(1-swAqua
     +Diva)-fsulph/fdocs/foo*fdoco/delt
         pd(20,12) = 0
         pd(3,20) = 0
         pd(18,20) = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-k
     +det*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)
     +-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)-km*kmax4*Ya-1/delt
         pd(1,13) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(5,13) = fcn*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*km*kmax3*Ys-fc
     +n*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)-fcn/delt
         pd(10,20) = 0
         pd(2,19) = 0
         pd(13,17) = -kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)-0.1D-1*km*kmax3*Ys-1/delt
         pd(9,17) = -1/delt
         pd(11,8) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(19,19) = 1
         pd(3,14) = 0
         pd(11,17) = -1/delt
         pd(18,10) = 1/fcn*faa*foo/fdoco/foa/fnitra*fdocn/delt
         pd(21,11) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)*sp(21,j)
         pd(17,14) = kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+1))
         pd(2,8) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(7,11) = kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(13,11) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(9,22) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-0.1D-1*km*kmax4*Ya-1/delt
         pd(11,2) = -1/delt
         pd(17,8) = 0
         pd(20,21) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+kdeac*(1-1/
     +(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))
     +/st/ammin)+1))
         pd(4,12) = 0
         pd(5,2) = -fcn*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/
     +(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1
     +/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+fcn*khyd*(1.D0+(1.
     +D0-swAquaDiva)*faz)**3*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)
     ++1)*Yo*kmax1**2*sp(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)**2/(ksox+sp(4,
     +j))**2/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2/(1/(exp((do
     +xmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2
         pd(12,3) = -2*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)/(exp
     +((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*
     +faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D
     +0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.
     +D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)
     ++1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(e
     +xp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j
     +)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp
     +(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0
     +-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+1)**2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*
     +kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +)*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp
     +(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D
     +0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j)))/st/doxmin)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2
     +/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0
     ++(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((
     +exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-
     +swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))*(-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksod
     +oc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1
     +*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiv
     +a)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)
     +*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3
     +,j))**2*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))-
     +2*khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)/(exp((amming-sp(5,j
     +))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+s
     +p(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxm
     +in-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*f
     +az)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(
     +1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))+khyd*(1.D0+(1.D0-swAquaDiva)*
     +faz)**2*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(2,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax
     +1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAq
     +uaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((
     +amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*fa
     +z)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*
     +*2/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.
     +D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))*(-kmax1*(1.D0
     ++(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+km
     +ax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-swAquaD
     +iva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp((amming-s
     +p(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/(1/(exp(
     +(doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1.D0-swAquaDi
     +va)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+k
     +hyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))**2*(-1/((exp((doxmin-kma
     +x1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*(1.D0+(1.D0-swAquaDiva)*faz
     +)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*
     +(1.D0+(1.D0-swAquaDiva)*faz)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +))+kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))**2*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*(1.D0+(1.D0-sw
     +AquaDiva)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st
     +/doxmin)+1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D
     +0+(1.D0-swAquaDiva)*faz)*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))-1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+
     +(1.D0-swAquaDiva)*faz)*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)
     +/(ksox+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz))
         pd(22,1) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(16,19) = 0
         pd(20,7) = 0
         pd(7,23) = 0
         pd(4,6) = 0
         pd(8,22) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(17,22) = 0
         pd(16,13) = -kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(23,1) = 0
         pd(12,16) = -0.1D-1*km*kmax3*Ys
         pd(7,1) = kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet
     +*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(12,7) = 1/delt
         pd(16,4) = -kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno
     +3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)-kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(
     +kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)
         pd(20,22) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)-0.1D-1*km*kmax4*Ya-1/d
     +elt
         pd(11,3) = -kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)-kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-k
     +max2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j)))/st/no3min)-1/fcn*faa*foo/fdoco/foa/delt
         pd(20,3) = 0
         pd(3,13) = 0
         pd(8,6) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(11,16) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kdet*sp(20,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(18,11) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(9,5) = -1/fcn/delt
         pd(17,13) = kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+1))
         pd(20,16) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(5,8) = -fcn*khyd*sp(3,j)**2/(exp((amming-sp(5,j))/st/amming)
     ++1)*Yn*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/(1/(exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+khyd*s
     +p(3,j))+fcn*khyd*sp(3,j)**3/(exp((amming-sp(5,j))/st/amming)+1)*Yn
     +*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**2*sp(10,j)**2/(ksno3
     ++sp(10,j))**2/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)**2/(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*
     +sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j))+khyd*sp(3,j))**2+fcn*(kl*vel**0.58D0+kdet/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax))+fcn*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*kdet*sp(2
     +0,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)
         pd(13,16) = -kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-0.1D-1*km*kmax3*Ys-1/d
     +elt
         pd(9,18) = 1/fcn*faa*foo/fdoco/foa*kpd-1/delt
         pd(2,7) = kreac/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+0.
     +1D-1*km*kmax1*Yo+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+1/delt
         pd(17,23) = 0
         pd(10,3) = 2*fnitra/fdocn*khyd*sp(3,j)/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j
     +)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3
     +*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)*(-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*
     +exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fn
     +itra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,
     +j))/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*
     +sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15
     +,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+khyd*sp(3,j)/fdocs)**2*(-1/((exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(13,j)*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*s
     +p(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))+khyd/fdocs)+2*fnitra/fdocn*khyd*sp(3,j)/(exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)*
     +*2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)**2*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)*(-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/
     +(ksso4+sp(15,j))/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*s
     +p(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)**2*(-1/((exp((so4min-kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(1
     +4,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kma
     +x3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))+khyd/fdocs)
         pd(11,12) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+kreac/(exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
         pd(13,2) = 0
         pd(2,12) = 0
         pd(4,23) = 0
         pd(14,6) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(16,18) = 0
         pd(4,11) = 0
         pd(15,17) = 0
         pd(22,2) = 0
         pd(23,15) = 0
         pd(6,17) = 0
         pd(22,5) = -kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(
     +5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5
     +,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4
     +,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)-kreac*sp(23,j)/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+
     +1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j))+kmax4*sp(5,j
     +)/(ksamm+sp(5,j))**2*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin*exp((ammin-
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)
         pd(6,5) = khyd*(1.D0+(1.D0-swAquaDiva)*faz)**2*sp(3,j)**2/(exp(
     +(amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(1,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva
     +)*faz)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
     ++1)/(1/(exp((doxmin-kmax1*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*(1.D0+(1
     +.D0-swAquaDiva)*faz)*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))+khyd*(1.D0+(1.D0-swAquaDiva)*faz)*sp(3,j))/st/amming*ex
     +p((amming-sp(5,j))/st/amming)
         pd(8,1) = -kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(9,10) = -kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)-kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3
     +,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax
     +2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j)))/st/no3min)-kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
     +/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)-kreac*sp(12,j)/(exp((no3min-kmax2*kind
     +ox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3
     +,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((n
     +o3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3
     +,j)/(ksndoc+sp(3,j)))/st/no3min)+1/fcn*faa*foo/fdoco/foa/fnitra*fd
     +ocn/delt
         pd(19,6) = 0
         pd(23,21) = kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksoxa+sp(4,j)))/st/ammin)+1))
         pd(1,20) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(6,11) = -kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(15,8) = 0
         pd(22,11) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(12,20) = 0
         pd(12,12) = 0
         pd(16,9) = 0
         pd(1,14) = 0
         pd(15,22) = 0
         pd(19,15) = 0
         pd(4,18) = 0
         pd(11,20) = -kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-kl*vel**0.58D0-kdet/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-km*kmax4*Ya-1/delt
         pd(21,7) = 0
         pd(14,11) = kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(15,3) = fsulph/fdocs/delt
         pd(18,6) = -kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-1/delt
         pd(22,20) = -kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(1,2) = 0
         pd(2,10) = 0
         pd(13,7) = 0
         pd(1,8) = kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax
     +*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(2,20) = -kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfm
     +ax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(6,23) = 0
         pd(7,6) = kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet
     +*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(10,15) = -fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(13,j)*s
     +p(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j))+khyd*sp(3,j)/fdocs)*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+
     +kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+fnitra/fdocn*k
     +hyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(13,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fd
     +ocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))**2/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/
     +(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)
     +/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ks
     +so4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +)+khyd*sp(3,j)/fdocs)**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(13,j)*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)
     +/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(14,j)*sp(15,j)/(
     +ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j))+khyd*sp(3,j)/fdocs)*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(
     +15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+fnitra/fdocn*khyd*sp(3,
     +j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4
     +min)+1)*kmax3*sp(14,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(1
     +5,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*khyd*
     +sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))**2/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax
     +3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-f
     +nitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15
     +,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp
     +(3,j)/fdocs)**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+1)**2)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+s
     +p(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/
     +so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))+fn
     +itra/fdocn/fsulph*fdocs/delt
         pd(19,9) = 0
         pd(5,19) = 0
         pd(21,1) = katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)*sp(21,j)
         pd(10,10) = -fnitra/fdocn*khyd*sp(3,j)**3/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(13,j
     +)*sp(15,j)**2/(ksso4+sp(15,j))**2/(kssdoc+sp(3,j))**2*kindox**2/(k
     +indox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/(1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)/st/so4min
     +*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-f
     +nitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15
     +,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j))**2/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd
     +*sp(3,j)/fdocs)-fnitra/fdocn*khyd*sp(3,j)**2/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(1
     +5,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j))/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))+khyd*sp(3,j)/fdocs)**2*(-1/((exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2)*kmax3**2*sp(13,j)*sp(1
     +5,j)**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox*
     +*2/(kindox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/st/so4min*ex
     +p((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-1/(e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*
     +kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2)-fnitra/fdocn*
     +khyd*sp(3,j)**3/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)**2*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+sp(15,j)
     +)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**2*kinno3**2/(
     +kinno3+sp(10,j))**3/(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+khyd*sp(3,j)/fdocs)/st/so4min*exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)-fnitra/fdocn*khyd*sp(3,j)**2/(
     +exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)
     +*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/(1/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)-fnitra/fdocn*
     +khyd*sp(3,j)**2/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(1/(exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*km
     +ax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+khyd*sp(3,j)/fdocs)
     +**2*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/
     +so4min)+1)**2)*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+sp(15,j))**2*s
     +p(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j))**2*kinno3
     +**2/(kinno3+sp(10,j))**3/st/so4min*exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)-1/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))**2)
         pd(17,4) = kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3
     ++sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min)+kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(k
     +inno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j)))/st/so4min)+kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinn
     +o3/(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j)))/st/so4min)+kreac*sp(17,j)/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2
     +*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)
      end
