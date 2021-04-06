*deck,user01       USERDISTRIB                                          ANSYS,INC
      function  user09()
      
c     TUTORIAL // DEMONSTRATOR
c#################################################     
c     SIMP structural + displacement constraint
c#################################################
c     Reserved APDL parameter names: LF1, NALL, NCON, v2nd, v2el, vels
c     user09 Input: [radel,volobj,ITMAX,esel,nsel,penal,ndcon,dispobj,ndis]
c         radel   (dp)        = filter radius in units of elements
c         volobj  (dp)        = percentaje of initial volume to retain
c         ITMAX   (int)       = maximum number of iterations for the algorithm
c         esel    (CHAR*4)    = ELEMENT NAMED SELECTION of the elements involved in the optimization
c         nsel    (CHAR*4)    =   NODAL NAMED SELECTION of the elements involved in the optimization
c         penal   (dp)        = penalization factor [E=Exmin+xx**penal*(Ex-Exmin)]
c         ndcon   (int)       = NODE NUMBER of the node to which the displacement constraint is applied
c         dispobj (dp)        = maximum displacement alowed in ndcon
c         ddir    (int)       = direction of the constrained displacement (1=-x,2=-y,3=-z)
c     Known limitations: 
c         * mesh required to be linear with hexahedral elements
              ! Mesh -> Element Order -> Linear
c         * No extra nodes allowed including the use of 'remote loads'
c         * 99 iterations maximum due to file naming: max(ITMAX)=99
c         * Required Sparse Direct Solver for .full file writing and assembled matrix extraction
              ! Analysis Settings -> Solver Controls -> Solver Type -> Direct
#include "impcom.inc"
#include "ansysdef.inc"
c  /*******************************************************************\
c  | Variable declaration            |
c  \*******************************************************************/
c     !!! EXTERNAL
      external  TrackBegin, TrackEnd, elmget,getnod ,
     + wrinqr, dpinfun, intinfun , erhandler,  RunCommand  ,     
     + pardef, pardim, parevl , dspget, ch4infun,mpget, elmput , mpput,
     + KS9, KCond8nod, CentrVol, CPout,VTKmesh_hexa9,VTK_CELLDATA9
      character string,int2str1*(1),int2str2*(2),cName(4),
     + LAB*(PARMSIZE), chValue*(STRING_MAX_LENG), ch4infun*(4),
     + ch8infun*(8), cNameIn*(PARMSIZE),cNameIn1*(PARMSIZE),esel*(4),
     + nsel*(4), nobj*(4),desellistchar*(4), eli*(10),kxminchar*(10),
     + file*(15), valw*(10), ndis*(4), Exminchar*(10)
      integer i, j, l, f, ii, it, tt, flag, flag1, m0, user09,
     + mi, iprop, ntab, nnod,outtxt, ntot, ntob, eltot, nelo, dt(8),
     + ndcon, KERR, ITMAX , nnel , nall, nodes(8),filtsizeestim,
     + ndof, ddir,ntemp,elmat(EL_MAT), tempmat, elmdat(EL_DIM),intinfun,
     + iott, wrinqr, isub3(3)
      double precision dummy, tempt, kcel, one, v(6), ab(8,3), lx, ly,
     + lz, length(3), vpos(8,3), Ce(3), kc, et, K185(8,8), nuxy, Ex,
     + K185m(24,24), adj(8), dispi(8), force(8), adjm(24), dispim(24),
     + forcem(24), kxmin, xmin, Exmin,vol, vt, rmin, rmin2, radel, sf,
     + distt, dist, etm, lamv, lamm, qv, qd, dres, dobj,penal, volobj,
     + vchange, voltot, dispobj, ll1, ll2, av, Tol, vs, v0, emin, emax,
     + vr, atarg, evol, pi, k1, k2,dpinfun, value2, dpValue3(3),
     + dpValue1, dpValue,SUBC3(3) , subc(2), subc1 ,et4(4),et5    
      double precision, allocatable::Lf(:), connel(:),conndt(:),ei(:),
     + eiav(:),ef(:),eiold(:), e0(:), xold(:), xx(:),centrx(:), 
     + centry(:) ,centrz(:),conn(:),conn1(:),connd(:),connd1(:),
     + volel(:), xel(:), xi(:), aux1(:), eic(:), eim(:), efm(:), 
     + dL(:), dadj(:), ein(:), eimn(:), eivn(:), ecm(:), disptot(:),
     + disptotadj(:),dst(:)
      integer, allocatable::ellist(:),aux2(:),timestore(:)

c     !!! MMA declarations
      external MMASUB
      integer m, n
      integer, allocatable:: iyfree(:)
      double precision, allocatable:: xval(:), dfdx(:), p(:), q(:),
     + fval(:), fmax(:), a(:), c(:), hessf(:), xminMMA(:), xmaxMMA(:),
     + xold1(:), xold2(:), xmma(:), alfa(:), beta(:), p0(:), q0(:),
     + dsrch(:), uu(:), df0dx(:), ULAM(:), xlow(:), xupp(:), Y(:),
     + gradf(:), b(:), volel1(:),cenx1(:),ceny1(:),cenz1(:)
      double precision f0val, geps, ghinit(1), ghmove(1), Z
      iott = wrinqr(WR_OUTPUT) 

      write (iott,'(A)')('*** user01 START')
c  /*******************************************************************\
c  | 0. Read Data           |
c  \*******************************************************************/
      pi = 3.1415927
c     Input: radel,volobj,ITMAX,esel,nsel,nobj,ndis,penal,ndcon,dispobj
      radel=dpinfun(2);volobj=dpinfun(3);ITMAX=intinfun(4)      
      esel=ch4infun(5);nsel=ch4infun(6)
      
      penal=dpinfun(7);ndcon=intinfun(8)
      dispobj=dpinfun(9);ddir=intinfun(10)
      ndof=3
      write (iott,'(A)')('*** Input data to user01 read')
      write (iott,'(F8.4,I10,F8.4,2A)')radel,ITMAX,volobj,esel,nsel
c  /*******************************************************************\
c  | 0. Element selection           |
c  \*******************************************************************/
      !!!!!!!!!!!!!!!!!!!!!!!!
      ! 0.1 Creating mesh related parameters in APDL
      !!!!!!!!!!!!!!!!!!!!!!!!
      call RunCommand(29,"allsel                  ")
      call RunCommand(29,"*GET, nall, NODE, 0, COUNT   ") ! nall = total number of nodes in the model
      call RunCommand(29,"allsel                       ")
      
      call RunCommand(29,"cmsel,s,"//nsel//",NODE")
      call RunCommand(29,"*GET, v2nd, NODE, 0, COUNT")    !v2nd = number of nodes to use for the optimization
      call RunCommand(29,"allsel")

      call RunCommand(29,"cmsel,s,"//esel//",ELEM")        
      call RunCommand(29,"*GET, v2el, ELEM, 0, COUNT   ") !v2el = number of elements to work on
      call RunCommand(29,"*DIM, vels,ARRAY -- ,v2el,1,1")
      call RunCommand(29,"*VGET, vels, ELEM,n, ELIST   ") !vels = element list to work on
      call RunCommand(29,"*VGET, vole, ELEM,n, GEOM   ") !vole = element list to work on
      call RunCommand(29,"*VGET, cenx, ELEM,n,CENT, X   ") !vole = element list to work on
      call RunCommand(29,"*VGET, ceny, ELEM,n,CENT, Y   ") !vole = element list to work on
      call RunCommand(29,"*VGET, cenz, ELEM,n,CENT, Z   ") !vole = element list to work on
      call RunCommand(29,"allsel                       ")     
      
      !!!!!!!!!!!!!!!!!!!!!!!!
      ! 0.2 Reading mesh related values from APDL
      !!!!!!!!!!!!!!!!!!!!!!!!
      LAB="V2EL"
      call PAREVL(LAB,0,SUBC(1),2,dpValue,chValue,kerr)
      eltot=int(dpValue)

      LAB="V2ND"
      call PAREVL(LAB,0,SUBC(1),2,dpValue,chValue,kerr)
      ntot=int(dpValue)
      SUBC3=(/1.,1.,1./)
      allocate(ellist(eltot)); LAB="VELS"
          do i=1,eltot
              SUBC3(1)=DBLE(i)
              call PAREVL(LAB,3,SUBC3,2,dpValue,chValue,kerr)
              ellist(i)=int(dpValue)
          enddo       
c  /*******************************************************************\
c  | 5. Reading mesh/material characteristics           |
c  \*******************************************************************/
          allocate(volel1(eltot));  allocate(connel(eltot+1))
      allocate(conndt,cenx1,ceny1,cenz1,xx,xel,ei, ef, ein, 
     + eim, efm, eivn, eimn, eiold, dL,eiav, xi, e0, eic,dadj,dst,ecm,
     + mold=volel1)       
      
            SUBC3=(/1.,1.,1./);kerr=0
      do i=1,eltot
          LAB="VOLE";SUBC3(1)=DBLE(ellist(i))
          call parevl(LAB,3,SUBC3,2,dpValue,chValue,kerr)
          volel1(i)=dpValue; LAB="CENX"
          call parevl(LAB,3,SUBC3,2,dpValue,chValue,kerr)
          cenx1(i)=dpValue;LAB="CENY"
          call parevl(LAB,3,SUBC3,2,dpValue,chValue,kerr)
          ceny1(i)=dpValue;LAB="CENZ"
          call parevl(LAB,3,SUBC3,2,dpValue,chValue,kerr)
          cenz1(i)=dpValue
      enddo
          kerr=0; nnel=8
      write(iott,'(A)')('Characteristic length read')
      voltot=sum(volel1)
           
c  /*******************************************************************\
c  | 2. Filtering matrix            |
c  \*******************************************************************/

      rmin=radel*volel1(1)**(1./3.); sf=3.;kerr=0
      filtsizeestim= nint(eltot*(4./3.*rmin**3.*3.1415)/volel1(1)*sf)
      allocate(conn1(filtsizeestim)); allocate(connd1(filtsizeestim))
      rmin2=(rmin)**2; f=0;  connel(1)=0;conn1=1
      do i=1,eltot
          distt=0.
          do j=1,eltot
                  dist=       ((cenx1(i)-cenx1(j))**2+(ceny1(i)-
     + ceny1(j))**2+(cenz1(i)-cenz1(j))**2)**0.5
                  if(dist<rmin)THEN
                  f=f+1
                  conn1(f)=j
                  connd1(f)=dist
                  distt=distt+dist
                  endif
           enddo
          conndt(i)=distt; connel(i+1)=f
      enddo
      allocate(conn(f)); allocate(connd(f))
      do i=1,f
          conn(i)= conn1(i); connd(i)=connd1(i)
      enddo
      deallocate(conn1); deallocate(connd1)
                  
c  /*******************************************************************\
c  | 2. Material creation            |
c  \*******************************************************************/
      m0=100;      nnod=8;    l=1
      call mpget(l,16,tempmat,kc)
      call mpget(l,1,tempmat,Ex);call mpget(l,4,tempmat,nuxy)
           write (iott,'(A,3E16.8)')'MatProp',kc,Ex,nuxy
           kxmin=1e-9; Exmin=Ex*1e-6
      do i=1,eltot
          mi=i+m0;ntab=1;tempt=25.
          iprop=16;kcel=(kc-kxmin)*1.**(penal)
          call mpput(mi,iprop,ntab,tempt,kcel)
          iprop=1;kcel=Exmin+1.**(penal)*(Ex-Exmin)
          call mpput(mi,iprop,ntab,tempt,kcel)
          iprop=4;kcel=0.3
          call mpput(mi,iprop,ntab,tempt,kcel)
      enddo
          do i=1,eltot
              call elmget(ellist(i),elmdat,nodes)
              elmdat(EL_MAT)=m0+i
              call elmput(ellist(i),elmdat,nnod,nodes)
          enddo
c  /*******************************************************************\
c  | 2. Lf ADJ initialization           |
c  \*******************************************************************/
      write(kxminchar,'(I10)')(ndof)
      call RunCommand(55,'*DIM,Lf1,,(NALL)*'//kxminchar)   
      cNameIn='LF1'; dpValue=-1.
      subc1=DBLE((ndcon-1)*ndof+ddir)
      call pardef(cNameIn,0,1,subc1,dpValue,chValue,kerr)
      call RunCommand(55,'*DMAT,Lfm1,D,IMPORT,APDL,LF1')
c  /*******************************************************************\
c  | 2.INITIALIZE MMA OPTIMIZER           |
c  \*******************************************************************/
          xx=volobj
        m     = 2                !% The number of general constraints.
        n     = eltot             !% The number of design variables x_j.
       xmin=0.00001
      allocate(xmaxMMA,xminMMA, xmma, alfa, beta, p0 , q0, df0dx,
     + mold=volel1)
      allocate(xupp, xlow,xval, xold1, xold2, source=xx)
          xmaxMMA=1.
          xminMMA=xmin
      allocate(a(m));a=0.    
      allocate(fmax, b, gradf, fval, c, Y, dsrch, uu, ULAM, source=a)
      c=1.e4 
      allocate(iyfree(m)); allocate(dfdx(m*n))
      allocate(p,q, mold=dfdx); allocate(hessf(m*(m+1)/2))
      ghinit=0.5;ghmove=0.7; geps=1e-8
c  /*******************************************************************\
c  | 2. ITERATION start & file openings           |
c  \*******************************************************************/
      et4=0.
      do ii=1,ITMAX !!!! Iteration START
           if(ii>5)THEN
               if(et4(4)/(sum(et4)/4.)-1.<0.01)THEN   !!! Convergence criteria!!!
                   EXIT
               endif
           endif
           
              write(iott,'(A,2I10)')'Iteration',ITMAX,ii
              call RunCommand(22,"FINISH")
      if (ii<10) THEN
          write(int2str1,'(I1)')(ii)
          call RunCommand(22,"/FILNAME,'TestUPF"//int2str1//"',0")
          file="TestUPF"//int2str1//".full"
      else
          write(int2str2,'(I2)')(ii)
          call RunCommand(23,"/FILNAME,'TestUPF"//int2str2//"',0")
          file="TestUPF"//int2str2//".full"
      endif
c  /*******************************************************************\
c  | Solve FEM THERMAL+Struct          |
c  \*******************************************************************/
      call RunCommand(5,"/SOLU")
      call RunCommand(50,"OUTRES, BASIC")
      call RunCommand(5,"SOLVE")  
      call RunCommand(5,'/AUX2')
      call RunCommand(17,"COMBINE, FULL, ALL")
c  /*******************************************************************\
c  | Solve  Adjoint           |
c  \*******************************************************************/
            
       if (ii .EQ. 1)THEN
      call RunCommand(55,'*SMAT,U2S,D,IMPORT,FULL,'//file//',USR2SOLV')
      call RunCommand(55,'*MULT,U2S,,Lfm1,,Lfs1')
      call RunCommand(55,'*DMAT,VecX1,D,COPY,Lfs1')
       endif
       
      call RunCommand(55,'*SMAT,MatK,D,IMPORT,FULL,'//file//',STIFF')
      call RunCommand(55,'*LSENGINE,BCS,MyBcsSolver1,MatK')
      call RunCommand(55,'*LSFACTOR,MyBcsSolver1')
      call RunCommand(55,'*LSBAC,MyBcsSolver1,Lfs1,VecX1')
      call RunCommand(55,'*MULT,U2S,T,VecX1,,XU1')
      cNameIn1='XU1'
c  /*******************************************************************\
c  | Calculating sensitivities            |
c  \*******************************************************************/
       
      et=0.;  isub3=(/1,2,3/);  adjm=0.
      call dspget(ndcon,1,ddir,dpValue)
      dres=dpValue

      do i=1,eltot
          call elmget(ellist(i),elmdat,nodes)
          do j=1,nnel
              l=nodes(j)
          call dspget(l,3,isub3,dpValue3)
              dispim((j-1)*3+1:(j-1)*3+3)=dpValue3
              do f=1,3
                  subc1=DBLE((nodes(j)-1)*ndof+f)
              call parevl(cNameIn1,1,subc1,2,dpValue,chValue,kerr)
                  adjm((j-1)*3+f)=dpValue
              enddo
          enddo

          do j=1,nnel
               call getnod(nodes(j),v,kerr,0)
               ab(j,1:3)=v(1:3)
          enddo
          call KS9(K185m,nuxy,ab)
          eim(i) = DOT_PRODUCT (MATMUL(K185m,dispim),adjm)
          eic(i)=DOT_PRODUCT (MATMUL(K185m,dispim),dispim)
          etm=etm+eic(i)*(Exmin+xx(i)**penal*(Ex-Exmin))
          eic(i)=-eic(i)*(penal*xx(i)**(penal-1.)*(Ex-Exmin))
          eim(i) =-eim(i)*(penal*xx(i)**(penal-1.)*(Ex-Exmin))
      enddo
      
      do i=1,3    !!! Update convergence criteria!!!
      et4(i)=et4(i+1)    
      enddo
      et4(4)=etm
      
      call RunCommand(55,'*FREE,MatK')
      call RunCommand(55,'*FREE,MyBcsSolver1')
      call RunCommand(55,'*FREE,XU1')
      write (iott,'(A)')('Sensitivities calculated')      
c  /*******************************************************************\
c  | Filtering sensitivities            |
c  \*******************************************************************/
      write (iott,'(A)')('SUBSTEP: Filtered sensitivities initialized')
      efm=0.; ecm=0.
      do i=1,eltot
              do j=connel(i)+1,connel(i+1)
                  efm(i)=efm(i)+eim(conn(j))*(rmin-connd(j))
                  ecm(i)=ecm(i)+eic(conn(j))*(rmin-connd(j))
              enddo
          efm(i)=efm(i)/conndt(i)
          ecm(i)=ecm(i)/conndt(i)
      enddo
c  /*******************************************************************\
c  | Elements selection , MMA           |
c  \*******************************************************************/
             write (iott,'(A)')('SUBSTEP: MMA')
       
      call dspget(ndcon,1,3,dpValue)
      xval  = xx; f0val = etm;  df0dx = ecm
      do i=1,n
          xi(i)=volel1(i)*xx(i)
      enddo
      fval(1)  = sum(xi)/voltot - volobj
      fval(2)  = abs(dpValue) - dispobj
      do i=1,n
      dfdx((i-1)*m+1)  = volel1(i) / voltot
      dfdx((i-1)*m+2)  = efm(i)
      enddo
      kerr=ii
      call MMASUB(kerr,m,n,geps,iyfree,xval,xmma, 
     +            xminMMA,xmaxMMA,xold1,xold2,xlow,xupp, 
     +            alfa,beta,a,b,c,Y,Z,ULAM, 
     +            f0val,fval,fmax,df0dx,dfdx, 
     +            p,q,p0,q0,uu,gradf,dsrch,hessf, 
     +            ghinit,ghmove)
      xold2    = xold1;  xold1    = xx;    xx=xmma
c  /*******************************************************************\
c  | ANSYS Material Property Change            |
c  \*******************************************************************
       write (iott,'(A)')('SUBSTEP: Material Change')
      m0=100                                              ! Initial material number for SIMP variables
      do i=1,eltot
          mi=i+m0;ntab=1;tempt=25.                        ! material number, number of temperature values & temperature values
          iprop=16;kcel=kxmin+xmma(i)**(penal)*(kc-kxmin) ! number of material property & property value
          call mpput(mi,iprop,ntab,tempt,kcel)            ! change thermal conductivity
          iprop=1;kcel=Exmin+xmma(i)**(penal)*(Ex-Exmin)  ! number of material property & property value
          call mpput(mi,iprop,ntab,tempt,kcel)            ! change Young modulus
      enddo   
      write (iott,'(A)')('ANSYS element selection/deselection finished')
c  /*******************************************************************\
c  | ANSYS Plotting           |
c  \*******************************************************************/
      write (iott,'(A)')('SUBSTEP: Opt. Output')
      if(ii .EQ. 1)THEN
      open(101, file='TOPOutput.txt', status = 'replace') 
      write(101,'(4A)')'Iteration,','Compliance,','Vol.,','Displacement'
      write(101,'(I10,4E16.8)')ii,etm,sum(xi)/voltot,dpValue
      ELSE
      open(101, file='TOPOutput.txt', status = 'old', access='append')  
      write(101,'(I10,4E16.8)')ii,etm,sum(xi)/voltot,dpValue
      endif; close(101)
      
      write (iott,'(A)')('SUBSTEP: VTK Generation')
      LAB="NALL"
      call PAREVL(LAB,0,SUBC(1),2,dpValue,chValue,kerr)
      nall=int(dpValue)
      if(ii < 10)THEN     
      open(101, file='VTKxMMA_'//int2str1//'.vtk', status = 'replace')  
      call VTKmesh_hexa9(nall,eltot,ellist,101)   
      call VTK_CELLDATA9(xmma,eltot,101,'xMMA00','xxt00',1)  
      else
      open(101, file='VTKxMMA_'//int2str2//'.vtk', status = 'replace')  
      call VTKmesh_hexa9(nall,eltot,ellist,101)   
      call VTK_CELLDATA9(xmma,eltot,101,'xMMA00','xxt00',1)  
      endif;  close(101)
      write (iott,'(A)')('SUBSTEP: VTK written')
      enddo
c  /*******************************************************************\
c  | END of subroutine write            |
c  \*******************************************************************/
      call erhandler ('user09',3000,
     x             2,'NODE OFFSET COMPLETE 2',0.0d0,' ')
      user09 = 0
      call TrackEnd ('user09')
      return
      end
c  /*******************************************************************\
c  | Stiffness matrix subroutine            |
c  \*******************************************************************/
        subroutine KS9(K185,nuv,coordinates)
        integer i,j,k,l,m
        double precision coordinates(8,3), FindDet7,C0(6,6), 
     + GaussPoint(2), xi1,xi2,xi3,dShape(3,8) ,JM(3,3) ,Jacinv(3,3), 
     + Bt(24,6),auxiliar(3,8), B(6,24), Bii(6,3), K185(24,24),detJ,
     + length(3), ff(24,24), f(6,24),length_x, length_y,length_z,Ai,Bi,
     + Ci,Di,Ei,Fi,Gi,Hi,Ii,nuv
    
        K185=0.; C0(:,:)=0.
        C0(1,1:3)=(/1.-nuv,nuv,nuv/) ;C0(2,1:3) = (/nuv, 1.-nuv, nuv/)
        C0(3,1:3) = (/nuv, nuv, 1.-nuv/);  C0(4,4) = (1.-2.*nuv)/2.
        C0(5,5) =(1.-2.*nuv)/2.;           C0(6,6) = (1.-2.*nuv)/2.
        C0=1./((1.+nuv)*(1.-2.*nuv))*C0
        GaussPoint = (/-1./sqrt(3.), 1./sqrt(3.)/)
        do i=1,2
            xi1=GaussPoint(i)
            do j=1,2
                xi2=GaussPoint(j)
                do k=1,2
                    xi3=GaussPoint(k)
                ! Compute shape functions derivatives
          dShape(1,1:8) = (/-(1.-xi2)*(1.-xi3),(1.-xi2)*(1.-xi3),
     + (1.+xi2)*(1.-xi3),-(1.+xi2)*(1.-xi3),-(1.-xi2)*(1.+xi3),
     + (1.-xi2)*(1.+xi3),(1.+xi2)*(1.+xi3),-(1.+xi2)*(1.+xi3)/)
          dShape(2,1:8) =  (/-(1.-xi1)*(1.-xi3),-(1.+xi1)*(1.-xi3),
     + (1.+xi1)*(1.-xi3),(1.-xi1)*(1.-xi3),-(1.-xi1)*(1.+xi3),
     + -(1.+xi1)*(1.+xi3),(1.+xi1)*(1.+xi3),(1.-xi1)*(1.+xi3)/)
          dShape(3,1:8) = (/-(1.-xi1)*(1.-xi2),-(1.+xi1)*(1.-xi2),
     + -(1.+xi1)*(1.+xi2),-(1.-xi1)*(1.+xi2),(1.-xi1)*(1.-xi2),
     + (1.+xi1)*(1.-xi2),(1.+xi1)*(1.+xi2),(1.-xi1)*(1.+xi2)/)
          dShape   =      dShape/8.       
              ! Compute Jacobian matrix
              do l=1,3
              do m=1,3
              JM(l,m) = DOT_PRODUCT(dShape(l,1:8),coordinates(1:8,m))
              enddo
              enddo
          Ai= (JM(2,2)*JM(3,3)-JM(3,2)*JM(2,3));
          Bi=-(JM(2,1)*JM(3,3)-JM(3,1)*JM(2,3))
          Ci= (JM(2,1)*JM(3,2)-JM(3,1)*JM(2,2))
          Di=-(JM(1,2)*JM(3,3)-JM(3,2)*JM(1,3))
          Ei= (JM(1,1)*JM(3,3)-JM(3,1)*JM(1,3))
          Fi=-(JM(1,1)*JM(3,2)-JM(3,1)*JM(1,2))
          Gi= (JM(1,2)*JM(2,3)-JM(2,2)*JM(1,3))
          Hi=-(JM(1,1)*JM(2,3)-JM(2,1)*JM(1,3))
          Ii= (JM(1,1)*JM(2,2)-JM(2,1)*JM(1,2))
          detJ=(JM(1,1)*Ai+JM(1,2)*Bi+JM(1,3)*Ci)
          Jacinv(1,1:3)=(/Ai,Di,Gi/)/detJ
          Jacinv(2,1:3)=(/Bi,Ei,Hi/)/detJ
          Jacinv(3,1:3)=(/Ci,Fi,Ii/)/detJ
                auxiliar = MATMUL(Jacinv,dShape)
                ! Preallocate memory for B-Operator
                B(:,:) = 0.; Bii(:,:) = 0.
                do l=1,8
                    Bii(1,1) = auxiliar(1,l);
                    Bii(2,2) = auxiliar(2,l);
                    Bii(3,3) = auxiliar(3,l);
                    Bii(4,1:2) = (/auxiliar(2,l),auxiliar(1,l)/)
                    Bii(5,2:3) = (/auxiliar(3,l),auxiliar(2,l)/)
                    Bii(6,1) = auxiliar(3,l)
                    Bii(6,3) = auxiliar(1,l)
                
                    B(:,(l-1)*3+1:(l)*3)=Bii
                enddo
                ff=MATMUL(transpose(B),MATMUL(C0,B))
                 K185 = K185 + ff*detJ
            enddo
        enddo
        enddo
        return
      end      
c  /*******************************************************************\
c  | VTKmesh         |
c  \*******************************************************************/
            subroutine VTKmesh_hexa9(nnode,nel,ellist,file)  
      external elmget, ndgxyz, getnod
      integer nnode, nel, file, zz, elmat(10), nodes(8), 
     + ndgxyz, kerr, zz1, cells, ellist(nel)
      double precision xyz(3), v(6)
      
      write(file,'(A)')('# vtk DataFile Version 2.0')
      write(file,'(A)')('Unstructured Hexa Mesh')
      write(file,'(A)')('ASCII')
      write(file,'(A)')('DATASET UNSTRUCTURED_GRID')
      write(file,"(A,I10,A)") 'POINTS ', nnode, ' float'
      do zz=1,nnode
          kerr = ndgxyz(zz,xyz)
          write(file,'(3F16.8)') xyz(1:3)
      enddo
      cells=nel*(8+1)
      write(file,'(A,I15,A,I15)') 'CELLS ',nel,' ',cells
      do zz=1,nel
          call elmget(ellist(zz),elmat,nodes)
          nodes=nodes-1
          write(file,'(I1,8I10)')8,nodes(1:8)
      enddo
      write(file,'(A,I10)')'CELL_TYPES ',nel
          write(file,'(I2)')(12, zz=1,nel)
      endsubroutine
c  /*******************************************************************\
c  | VTK Cell data         |
c  \*******************************************************************/
            subroutine VTK_CELLDATA9(datavec,nel,file,dataname,tablename
     +, comp)  
      integer  nel, file, zz, comp
      character*5 dataname, tablename, datatype
      double precision datavec(nel)

      write(file,'(A,I10)')'CELL_DATA ',nel
      write(file,'(A,A,A,I1)')'SCALARS ',dataname,' float ', comp
      write(file,'(A,A)')'LOOKUP_TABLE ', tablename
      write(file,'(E16.8)') (datavec(zz),zz=1,nel)
      endsubroutine