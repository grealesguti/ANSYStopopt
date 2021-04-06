C********+*********+*********+*********+*********+*********+*********+
C
      subroutine KStiff8nod(K185,E,nuv,coordinates)
        double precision nuv,E
        integer i,j,k,l,m
        double precision coordinates(8,3), FindDet7
        double precision C0(6,6), GaussPoint(2), xi1,xi2,xi3
        double precision dShape(3,8) ,JM(3,3) ,Jacinv(3,3), Bt(24,6)
        double precision auxiliar(3,8), B(6,24), Bii(6,3), K185(24,24)
        double precision detJ,    length(3), ff(24,24), f(6,24)
        double precision length_x, length_y,length_z
        double precision Ai,Bi,Ci,Di,Ei,Fi,Gi,Hi,Ii
    
      K185=0.
    
        C0(:,:)=0.
        C0(1,1:3)=(/1.-nuv,nuv,nuv/) 
        C0(2,1:3) = (/nuv, 1.-nuv, nuv/)
        C0(3,1:3) = (/nuv, nuv, 1.-nuv/)
        C0(4,4) = (1.-2.*nuv)/2.
        C0(5,5) =(1.-2.*nuv)/2.
        C0(6,6) = (1.-2.*nuv)/2.
    
        C0=E/((1.+nuv)*(1.-2.*nuv))*C0

        GaussPoint = (/-1./sqrt(3.), 1./sqrt(3.)/)

        do i=1,2
            xi1=GaussPoint(i)
            do j=1,2
                xi2=GaussPoint(j)
                do k=1,2
                    xi3=GaussPoint(k)

                ! Compute shape functions derivatives
          dShape(1,1:2) = (/-(1.-xi2)*(1.-xi3),(1.-xi2)*(1.-xi3)/)
          dShape(1,3:4) = (/(1.+xi2)*(1.-xi3),-(1.+xi2)*(1.-xi3)/)
          dShape(1,5:6) = (/-(1.-xi2)*(1.+xi3),(1.-xi2)*(1.+xi3)/)
          dShape(1,7:8) = (/(1.+xi2)*(1.+xi3),-(1.+xi2)*(1.+xi3)/)
      
          dShape(2,1:2) =  (/-(1.-xi1)*(1.-xi3),-(1.+xi1)*(1.-xi3)/)
          dShape(2,3:4) = (/(1.+xi1)*(1.-xi3),(1.-xi1)*(1.-xi3)/)
          dShape(2,5:6) = (/-(1.-xi1)*(1.+xi3),-(1.+xi1)*(1.+xi3)/)
          dShape(2,7:8) = (/(1.+xi1)*(1.+xi3),(1.-xi1)*(1.+xi3)/)
      
          dShape(3,1:2) = (/-(1.-xi1)*(1.-xi2),-(1.+xi1)*(1.-xi2)/)
          dShape(3,3:4) = (/-(1.+xi1)*(1.+xi2),-(1.-xi1)*(1.+xi2)/)
          dShape(3,5:6) = (/(1.-xi1)*(1.-xi2),(1.+xi1)*(1.-xi2)/)
          dShape(3,7:8) = (/(1.+xi1)*(1.+xi2),(1.-xi1)*(1.+xi2)/)
      
          dShape   =      dShape/8.       

              ! Compute Jacobian matrix
              do l=1,3
              do m=1,3
              JM(l,m) = DOT_PRODUCT(dShape(l,1:8),coordinates(1:8,m))
              enddo
              enddo
          
          Ai= (JM(2,2)*JM(3,3)-JM(3,2)*JM(2,3))
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
                B(:,:) = 0.
                Bii(:,:) = 0.
                !% Construct first three rows
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
                
                f=MATMUL(C0,B)
                
                ff=MATMUL(transpose(B),f)
                
                K185 = K185 + ff*detJ

          
            enddo
        enddo
        enddo
        return
      end      
      
!  /*******************************************************************\
!  | Centroid element 185            |
!  \*******************************************************************/

        subroutine CentrVol(nodeposh,cmf,vt)
    
        integer  facei(3,3), file
        integer th(6,4), faces2(4,3), f, k ,i ,j,l,m
        double precision vt, cmf(3), vpos(8,3), fvertex(3,3)
        double precision nodeposh(8,3), xi(3), yi(3), zi(3)
        double precision beta(2), gam(2), d0, V
        double precision cmv(6,3), vv(6), nodepos(4,3)

        th(1,:)=(/4, 1, 2, 5/)
        th(2,:)=(/4, 2, 6, 5/)
        th(3,:)=(/4, 6, 8, 5/)
        th(4,:)=(/4, 2, 3, 7/)
        th(5,:)=(/4, 6, 2, 7/)
        th(6,:)=(/4, 8, 6, 7/)
    
        faces2(1,:)=(/1, 3, 2/)
        faces2(2,:)=(/2, 4, 1/)
        faces2(3,:)=(/3, 4, 2/)
        faces2(4,:)=(/3, 1, 4/)
        
        
        !write(file,'(A)')('nodeposh')
       ! do i=1,8
       !     write(file,'(F10.5)')(nodeposh(i,j),j=1,3)
       ! enddo
        
        !write(file,'(A)')('faces')
       ! do i=1,4
        !    write(file,'(I10)')(faces2(i,j),j=1,3)
       ! enddo        
    
        f = 4

        vt=0.
        do k=1,6
            do i=1,4
                nodepos(i,:)=nodeposh(th(k,i),:)  
            enddo

                        !     write(file,'(A)')('nodepos tetra')
                  !do i =1,4
               ! write(file,'(E20.10)')(nodepos(i,j),j=1,3)
                 ! enddo
            
              V=0.
            do i =1,f
             ! write(file,'(A)')('facei')                
      
                    
                    fvertex(1,:)=nodepos(faces2(i,1),:)
                      !write(file,'(I10)')(faces2(i,1))
                      !write(file,'(E20.10)')(fvertex(1,m),m=1,3)
               !write(file,'(E20.10)')(nodepos(faces2(i,1),m),m=1,3)

                    fvertex(2,:)=nodepos(faces2(i,2),:)
                      !write(file,'(I10)')(faces2(i,2))

                      !write(file,'(E20.10)')(fvertex(2,m),m=1,3)
              !write(file,'(E20.10)')(nodepos(faces2(i,2),m),m=1,3)

                    fvertex(3,:)=nodepos(faces2(i,3),:)
                      !write(file,'(I10)')(faces2(i,3))
                !write(file,'(E20.10)')(fvertex(3,m),m=1,3)
             ! write(file,'(E20.10)')(nodepos(faces2(i,3),m),m=1,3)

                xi=fvertex(:,1)
                yi=fvertex(:,2)
                zi=fvertex(:,3)
                beta=(/yi(2)-yi(1),yi(3)-yi(1)/)
                gam=(/zi(2)-zi(1),zi(3)-zi(1)/)
                d0=beta(1)*gam(2)-beta(2)*gam(1)
                
               ! write(file,'(A)')('xyz')
               ! write(file,'(E20.10)')(xi(j),j=1,3)
                !write(file,'(E20.10)')(yi(j),j=1,3)
                !write(file,'(E20.10)')(zi(j),j=1,3)

                !write(file,'(A)')('g0')
                !write(file,'(E20.10)')(gam)
                !write(file,'(A)')('b0')
               ! write(file,'(E20.10)')(beta)
               ! write(file,'(A)')('d0')
                !write(file,'(E20.10)')(d0)

                V=V+d0*sum(xi)/6.
               ! write(file,'(A)')('V')
                !write(file,'(E20.10)')(V)
                
            enddo

            cmv(k,1)=sum(nodeposh(:,1))/4.;
            cmv(k,2)=sum(nodeposh(:,2))/4.;
            cmv(k,3)=sum(nodeposh(:,3))/4.;
            vv(k)=V;
            vt=V+vt;
        enddo
        
       ! do i=1,6
       ! write(file,'(A,I5)')('tetra',i)
        !write(file,'(E20.10)')(cmv(i,j),j=1,3)
        !write(file,'(A)')('volume')
        !write(file,'(E20.10)')(vv(i))
       ! enddo
        

        
        !write(file,'(A)')('hexahedron volume')
        !write(file,'(E20.10)')(vt)
        
        cmf(1)=DOT_PRODUCT(cmv(:,1),vv)/vt/2.
        cmf(2)=DOT_PRODUCT(cmv(:,2),vv)/vt/2.
        cmf(3)=DOT_PRODUCT(cmv(:,3),vv)/vt/2.
        
        !write(file,'(A)')('hexahedron cm')
        !write(file,'(E20.10)')(cmf(j),j=1,3)
        end

c  /*******************************************************************\
c  | Thermal conductivity matrix            |
c  \*******************************************************************/
      subroutine KCond8nod(K185,kc,coordinates)
        double precision nuv,E
        integer i,j,k,l,m
        double precision coordinates(8,3), FindDet06,kc
        double precision C0(6,6), GaussPoint(2), xi1,xi2,xi3
        double precision dShape(3,8) ,JacobianMatrix(3,3) ,Jacinv(3,3)
        double precision auxiliar(3,8), B(6,8), K185(8,8)
        double precision detJ,    length(3), ff(8,8), f(3,8)
        double precision length_x, length_y,length_z, JM(3,3)
        double precision Ai,Bi,Ci,Di,Ei,Fi,Gi,Hi,Ii    
      K185=0.
    
        
       ! open(9, file = 'F9_C0.txt', status = 'replace')  
        !write(9,'(E20.10)')(E)  
       ! write(9,'(E20.10)')(nuv)  
       ! do i=1,6
        !              write(9,'(A,I2)')('row',i)  
       !   write(9,'(E20.10)')(C0(l,m),m=i,6)  
       !     enddo
       ! close(9)
    
        GaussPoint = (/-1./sqrt(3.), 1./sqrt(3.)/)
               !open(1, file = 'F1_DHAPEv2.txt', status = 'replace')  
               !open(2, file = 'F2_Jinvv2.txt', status = 'replace')  
               !open(3, file = 'F3_Bv2.txt', status = 'replace')  
               !open(4, file = 'F4_detJv2.txt', status = 'replace')  
               !open(5, file = 'F5_fv2.txt', status = 'replace')  
               !open(6, file = 'F6_ffv2.txt', status = 'replace')  
               !open(7, file = 'F7_Jv2.txt', status = 'replace')  
             !open(8, file = 'F8_kv2.txt', status = 'replace')  
                              
              ! open(9, file = 'F9_coords.txt', status = 'replace')  


        do i=1,2
            xi1=GaussPoint(i)
            do j=1,2
                xi2=GaussPoint(j)
                do k=1,2
                    xi3=GaussPoint(k)

                ! Compute shape functions derivatives
          dShape(1,1:2) = (/-(1.-xi2)*(1.-xi3),(1.-xi2)*(1.-xi3)/)
          dShape(1,3:4) = (/(1.+xi2)*(1.-xi3),-(1.+xi2)*(1.-xi3)/)
          dShape(1,5:6) = (/-(1.-xi2)*(1.+xi3),(1.-xi2)*(1.+xi3)/)
          dShape(1,7:8) = (/(1.+xi2)*(1.+xi3),-(1.+xi2)*(1.+xi3)/)
      
          dShape(2,1:2) =  (/-(1.-xi1)*(1.-xi3),-(1.+xi1)*(1.-xi3)/)
          dShape(2,3:4) = (/(1.+xi1)*(1.-xi3),(1.-xi1)*(1.-xi3)/)
          dShape(2,5:6) = (/-(1.-xi1)*(1.+xi3),-(1.+xi1)*(1.+xi3)/)
          dShape(2,7:8) = (/(1.+xi1)*(1.+xi3),(1.-xi1)*(1.+xi3)/)
      
          dShape(3,1:2) = (/-(1.-xi1)*(1.-xi2),-(1.+xi1)*(1.-xi2)/)
          dShape(3,3:4) = (/-(1.+xi1)*(1.+xi2),-(1.-xi1)*(1.+xi2)/)
          dShape(3,5:6) = (/(1.-xi1)*(1.-xi2),(1.+xi1)*(1.-xi2)/)
          dShape(3,7:8) = (/(1.+xi1)*(1.+xi2),(1.-xi1)*(1.+xi2)/)
      
          dShape   =      dShape/8.       
          !do l=1,3
          !write(1,'(A,I2)')('row',l)  
          !write(1,'(E20.10)')(dShape(l,m),m=1,8)  
          !enddo

          !            do l=1,8
          !write(9,'(A,I2)')('row',l)  
          !write(9,'(E20.10)')(coordinates(l,m),m=1,3)  
          !enddo

            
              ! Compute Jacobian matrix
              do l=1,3
              do m=1,3
              JM(l,m) = DOT_PRODUCT(dShape(l,1:8),coordinates(1:8,m))
              enddo
              enddo
              !JacobianMatrix = MATMUL(dShape,coordinates)
               ! Compute auxiliar matrix for construction of B-Operator
          !do l=1,3
          !write(7,'(A,I2)')('row',l)  
          !write(7,'(E20.10)')(JM(l,m),m=1,3)  
          !enddo
          
          
          Ai= (JM(2,2)*JM(3,3)-JM(3,2)*JM(2,3))
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
          
          !do l=1,3
          !write(2,'(A,I2)')('row',l)  
          !write(2,'(E20.10)')(Jacinv(l,m),m=1,3)  
          !enddo
                
                auxiliar = MATMUL(Jacinv,dShape)
            
                
          !write(3,'(A,I2)')('it',i)
          !!do l=1,3
          !write(3,'(A,I2)')('row',l)  
          !write(3,'(E20.10)')(auxiliar(l,m),m=1,8)  
          !enddo
          
         ! do l=1,6
         !     do m=1,24
          !f(l,m) = DOT_PRODUCT(C0(l,1:6),B(1:6,m))
          !    enddo
          !enddo
                 f=kc*auxiliar
                
          !write(5,'(A,I2)')('it',i)
      !do l=1,6
          !write(5,'(A,I2)')('row',l)  
          !write(5,'(E20.10)')(f(l,m),m=1,8)  
          !enddo      
                ff=MATMUL(transpose(auxiliar),f)
                
          !write(6,'(A,I2)')('it',i)
          !do l=1,8
          !rite(6,'(A,I2)')('row',l)  
          !write(6,'(E20.10)')(ff(l,m),m=1,8)  
          !enddo 
                ! Add to stiffness matrix
                
                K185 = K185 + ff*detJ
          !do l=1,8
          !write(8,'(A,I2)')('row',l)  
          !write(8,'(E20.10)')(K185(l,m),m=1,8)  
          !enddo
          
                    write(4,'(E16.8)')(detJ)  

          
            enddo
        enddo
        enddo
       !close(1)
        
          !close(2)
          
          !close(3)
          
          !close(4)
          !close(5)

           !close(6)

          !close(7)
          !close(8)
          !close(9)
      end      
      
c  /*******************************************************************\
c  | ANSYS Plotting density based         |
c  \*******************************************************************/
      subroutine PLTel(ellist,xx,elmtot,dof,PARMSIZE)
      
      external RunCommand,elmget,elmput, parevl,pardef
      integer i, kerr, elmtot
      double precision SUBC3(3), nodes(8), dpValue, ellist(elmtot),
     + xx(elmtot)
      character dof*(4), cNameIn*(PARMSIZE),LAB*(PARMSIZE)
      
          call RunCommand(6,'/post1')  
          call RunCommand(15,"/SHOW, PNG")
          call RunCommand(25,'/RGB,INDEX,100,100,100, 0')
		call RunCommand(25,'/RGB,INDEX, 80, 80, 80,13')
		call RunCommand(25,'/RGB,INDEX, 60, 60, 60,14')
		call RunCommand(25,'/RGB,INDEX, 0, 0, 0,15   ')
          call RunCommand(9,'/win,1,on')
          call RunCommand(13,'/VIEW,1,1,1,1')          
      call RunCommand(55,'*DIM,XRHO,,NALL')
      SUBC3=(/1.,1.,1./)
      cNameIn='XRHO'
      LAB='XRHO'
      do i=1,eltot
          call elmget(ellist(i),elmdat,nodes)
          do j=1,8
          SUBC3(1)=DBLE(nodes(j))
          call PAREVL(LAB,3,SUBC3(1),2,dpValue,chValue,kerr)
          subc1=DBLE(nodes(j))
          dpValue=xx(i)/8.+dpValue
          call pardef(cNameIn,0,1,subc1,dpValue,chValue,kerr)
          enddo
      enddo
      call RunCommand(55,'*VPUT,XRHO(1),NODE,1,TEMP')
      
          call RunCommand(25,'allsel            ')
          call RunCommand(25,'PLNSOL , '//dof) ! Plt MAG
                    
          write (iott,'(A)')('ANSYS Density Plotting finished')
      endsubroutine
      
c  /*******************************************************************\
c  | ANSYS Plotting DOF         |
c  \*******************************************************************/
      subroutine PLTDOF(dof)     
      character*4 dof
            external RunCommand

          call RunCommand(6,'/post1')  
          call RunCommand(15,"/SHOW, PNG")
          call RunCommand(25,'/RGB,INDEX,100,100,100, 0')
		call RunCommand(25,'/RGB,INDEX, 80, 80, 80,13')
		call RunCommand(25,'/RGB,INDEX, 60, 60, 60,14')
		call RunCommand(25,'/RGB,INDEX, 0, 0, 0,15   ')
          call RunCommand(9,'/win,1,on')
          call RunCommand(13,'/VIEW,1,1,1,1')          
          call RunCommand(25,'allsel            ')
          call RunCommand(25,'PLNSOL ,'//dof) ! Plt dof
                    
          write (iott,'(A)')('ANSYS '//dof//' Plotting finished')
      endsubroutine
      
      
c  /*******************************************************************\
c  | Time write         |
c  \*******************************************************************/
      subroutine CPout()     
      external RunCommand
           call RunCommand(50,'*get,TCPU,active,,time,CPU ')
           call RunCommand(50,'/com TCPU ')
      endsubroutine
      
c  /*******************************************************************\
c  | VTKmesh         |
c  \*******************************************************************/
      subroutine VTKmesh_hexa(nnode,nel,ellist,file)  
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
          write(file,'(F16.8,F16.8,F16.8)') xyz(1),xyz(2),xyz(3)
      enddo
      cells=nel*(8+1)
      write(file,'(A,I15,A,I15)') 'CELLS ',nel,' ',cells
      !write(file,'(A)') 'CELLS done'

      do zz=1,nel
          call elmget(ellist(zz),elmat,nodes)
          nodes=nodes-1
      !    write(file,'(A)') 'elem called'
          write(file,'(I1,I10,I10,I10,I10,I10,I10,I10,I10)')
     +  8,
     +  nodes(1),nodes(2),nodes(3),nodes(4),nodes(5),nodes(6),
     +  nodes(7),nodes(8)
      enddo
      write(file,'(A,I10)')'CELL_TYPES ',nel
      do zz=1,nel
          write(file,'(I2)')(12)
      enddo
      endsubroutine

      
c  /*******************************************************************\
c  | VTKmesh CELL data        |
c  \*******************************************************************/
      subroutine VTK_CELLDATA(datavec,nel,file,dataname,tablename
     +, comp)  
      integer  nel, file, zz, comp
      character*5 dataname, tablename, datatype
      double precision datavec(nel)

      
      write(file,'(A,I10)')'CELL_DATA ',nel
      write(file,'(A,A,A,I1)')'SCALARS ',dataname,' float ', comp
      write(file,'(A,A)')'LOOKUP_TABLE ', tablename
      do zz=1,nel
          write(file,'(E16.8)') datavec(zz)
      enddo
      endsubroutine
c  /*******************************************************************\
c  | VTKmesh POINT data        |
c  \*******************************************************************/
      subroutine VTK_POINTDATA()  

      endsubroutine

C********+*********+*********+*********+*********+*********+*********+
