        implicit real*8(a-h,o-z) 
        parameter (nsolid=8,ndiss=15,ncomp=23,nreac=68,nx=1361)
        common/kinetics/kdet,kmax1,kl,Yo,st,ksodoc,ksox,ksamm,doxmin,Bfm
     +ax,Yn,ksndoc,kmax2,ksno3,kindox,no3min,katt,kdeac,kmax3,ksso4,kssd
     +oc,kinno3,so4min,Ys,kmax4,kreac,km,kpd,ammin,swAquaDiva,sw2,ksnit,
     +fdoco,fdocn,fdocs,foo,foa,fnitra,fsulph,fcn,faa,Ya,amming,kmax5,Ym
     +,ws,vel,khyd,ksoxa,ktrans,faz,kmax6
        real*8 kdet,kmax1,kl,Yo,st,ksodoc,ksox,ksamm,doxmin,Bfmax,Yn,ksn
     +doc,kmax2,ksno3,kindox,no3min,katt,kdeac,kmax3,ksso4,kssdoc,kinno3
     +,so4min,Ys,kmax4,kreac,km,kpd,ammin,swAquaDiva,sw2,ksnit,fdoco,fdo
     +cn,fdocs,foo,foa,fnitra,fsulph,fcn,faa,Ya,amming,kmax5,Ym,ws,vel,k
     +hyd,ksoxa,ktrans,faz,kmax6
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,endt
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,end
     +t
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic
        integer iq,iw,iDb,ipor,igrid,iarea,ic
