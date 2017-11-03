!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
      
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ie
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt reference coords
      double precision  ::  dNdy(20,3)                        ! Derivative of shape functions wrt deformed coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  stress(6)                        ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  F(3,3),sigma(3,3), q(3,3) ,CauchyStress(3,3)                           ! Deformation gradient
      double precision  ::  Finv(3,3)                         ! Inverse of deformation gradient
      double precision  ::  B(3,3) , id(3,3), kab(NNODE, NNODE)  ! C-G deformation tensor
      double precision  ::  JJ                                ! det(F)
      double precision  ::  G(6,9),H(6,9)
      double precision  ::  D(6,6)                            ! Material tangent
      double precision  ::  Bbar(6,60)                        ! strain = Bbar*(dof_total)
      double precision  ::  Bstar(9,60)                       ! F = Bstar*(dof_total)
      double precision  ::  Pvec(3*NNODE),qvec(9),cauchy(6)
      double precision  ::  Pmat(3*NNODE,3*NNODE), Y(3*NNODE,3*NNODE)
      double precision  ::  P(3*NNODE,3*NNODE)     !
      double precision  ::  S(3,NNODE)
      double precision  ::  Svec(3*NNODE)
      double precision  ::  Smat(3*NNODE,3*NNODE)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant

    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      kab=0.d0
      H=0.d0
      ENERGY(1:8) = 0.d0
      q=0.d0
      qvec=0.d0
      id=0.d0
      Y=0.d0
      Bstar=0.d0
      C=0.d0
      F=0.d0
      D=0.d0
      JJ=0.d0
    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)

        ! Caculate the deformation gradient
        do i = 1,3
            ie = 3*(NNODE-1)+i
            F(i,1:3) = matmul(U(i:ie:3),dNdx(1:NNODE,1:3))        
            F(i,i) = F(i,i) + 1.d0                           ! F=Del[u] + I
        end do
       ! B = matmul(F,transpose(F))                           ! Left cauchy green tensor.
        
        call abq_UEL_invert3d(F,Finv,JJ)                     ! Givens Finv and Jacobian
       
        !dNdy(1:NNODE,1:3) = matmul(dNdx(1:NNODE,1:3),Finv)   NOT needed.   
        
        call fung(PROPS(1:NPROPS),NPROPS,F,JJ,stress,sigma,D) ! Subroutine gives stress vector, sigma matrix and tangent stiffness
        
        q=matmul(sigma,transpose(F))                             !MAtrix
        
        qvec(1)=q(1,1)
        qvec(2)=q(2,2)
        qvec(3)=q(3,3)
        qvec(4)=q(2,1)
        qvec(5)=q(1,2)
        qvec(6)=q(3,1)
        qvec(7)=q(1,3)
        qvec(8)=q(3,2)
        qvec(9)=q(2,3)
        
        
        H(1,1)=F(1,1)
        H(2,2)=F(2,2)
        H(3,3)=F(3,3)
        H(4,1)=F(1,2)
        H(4,2)=F(2,1)
        H(5,1)=F(1,3)
        H(5,3)=F(3,1)
        H(6,2)=F(2,3)
        H(6,3)=F(3,2)
        
        H(1,5)=F(2,1)
        H(2,4)=F(1,2)
        H(3,6)=F(1,3)
        H(4,4)=F(1,1)
        H(4,5)=F(2,2)
        H(5,5)=F(2,3)
        H(5,6)=F(1,1)
        H(6,4)=F(1,3)
        H(6,6)=F(1,2)
        
        H(1,7)=F(3,1)
        H(2,9)=F(3,2)
        H(3,8)=F(2,3)
        H(4,7)=F(3,2)
        H(4,9)=F(3,1)
        H(5,7)=F(3,3)
        H(5,8)=F(2,1)
        H(6,8)=F(2,2)
        H(6,9)=F(3,3)
        
        
        !
        !Bbar = 0.d0
        !Bbar(1,1:3*NNODE-2:3) = dNdy(1:NNODE,1)
        !Bbar(2,2:3*NNODE-1:3) = dNdy(1:NNODE,2)
        !Bbar(3,3:3*NNODE:3)   = dNdy(1:NNODE,3)       !used for geometrical stiffness
        !Bbar(4,1:3*NNODE-2:3) = dNdy(1:NNODE,2)
        !Bbar(4,2:3*NNODE-1:3) = dNdy(1:NNODE,1)
        !Bbar(5,1:3*NNODE-2:3) = dNdy(1:NNODE,3)
        !Bbar(5,3:3*NNODE:3)   = dNdy(1:NNODE,1)       !Not needed.
        !Bbar(6,2:3*NNODE-1:3) = dNdy(1:NNODE,3)
        !Bbar(6,3:3*NNODE:3)   = dNdy(1:NNODE,2)        
        !
        Bstar = 0.d0
        Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)       !CHECK?
        Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)
   
     !!
     !!   G(1,1:9) = [B(1,1),0.d0,0.d0,B(1,2),0.d0,B(1,3),0.d0,0.d0,0.d0]
     !!   G(2,1:9) = [0.d0,B(2,2),0.d0,0.d0,B(1,2),0.d0,0.d0,B(2,3),0.d0]
     !!   G(3,1:9) = [0.d0,0.d0,B(3,3),0.d0,0.d0,0.d0,B(1,3),0.d0,B(2,3)]
     !!   G(4,1:9) = [B(1,2),B(1,2),0.d0,
     !!1                          B(2,2),B(1,1),B(2,3),0.d0,B(1,3), 0.d0]
     !!   G(5,1:9) = [B(1,3),0.d0,B(1,3),
     !!1                          B(2,3),0.d0,B(3,3),B(1,1),0.d0,B(1,2)]
     !!   G(6,1:9) = [0.d0,B(2,3),B(2,3),
     !!1                           0.d0,B(1,3),0.d0,B(1,2),B(3,3),B(2,2)]
     !!
     !!   G = 2.d0*G
     !!
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     2   - matmul(transpose(Bstar(1:9,1:3*NNODE)),qvec)*    !CHECK POSITIVE SIGN & dimension of Bstar?
     3                                          w(kint)*determinant
     !!
     !!
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + matmul(transpose(Bstar(1:9,1:3*NNODE)),
     2     matmul(transpose(H),matmul(D,matmul(H,Bstar(1:9,            ! I think, done.
     3   1:3*NNODE))))) *w(kint)*determinant

!       Geometric stiffness
        
     !!   S = reshape(matmul(transpose(Bbar),stress),(/3,3*NNODE/3/))
     !!   do i = 1,NNODE
     !!     Pvec(1:3*NNODE) = reshape(spread(transpose(dNdy(i:i,1:3)),       !Not needed to 
     !!1                              dim=2,ncopies=NNODE),(/3*NNODE/))
     !!     Pmat(3*i-2:3*i,1:3*NNODE) = spread(Pvec,dim=1,ncopies=3)
     !!     Svec(1:3*NNODE) = reshape(spread(S(1:3,i:i),
     !!1                              dim=2,ncopies=NNODE),(/3*NNODE/))
     !!     Smat(3*i-2:3*i,1:3*NNODE) = spread(Svec,dim=1,ncopies=3)
     !!   end do
     !!   
        
        kab=matmul(dNdx(1:NNODE,1:3), matmul(sigma, transpose(dNdx
     1   (1:NNODE,1:3))))
        
        do i=1,3
          id(i,i)=1.d0           ! Identity matrix   !initialize zero.
        end do
        
        do i=1,NNODE
            do j=1,NNODE
          Y(3*i-2:3*i,3*j-2:3*j)=kab(i,j)*id                ! Geometric stiffness assembled
            end do
        end do
        
        
        
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1              + Y(1:3*NNODE,1:3*NNODE)**w(kint)*determinant 
        
        CauchyStress=1.d0/JJ *matmul(F,matmul(sigma,transpose(F)))
        cauchy(1)=CauchyStress(1,1)
        cauchy(2)=CauchyStress(2,2)
        cauchy(3)=CauchyStress(3,3)
        cauchy(4)=CauchyStress(1,2)
        cauchy(5)=CauchyStress(1,3)
        cauchy(6)=CauchyStress(2,3)
        
        if (NSVARS>=n_points*6) then   ! Store Cauchy stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) =  Cauchy(1:6)
        endif
      end do


      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention.
    !
    !   This is coded to apply nominal tractions to the element face (the residual force does not change as the element deforms)
    !
    !
     ! ! do j = 1,NDLOAD
     !!
     !!   call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     !!1                                     face_node_list,nfacenodes)
     !!
     !!   do i = 1,nfacenodes
     !!       face_coords(1:3,i) = coords(1:3,face_node_list(i))
     !!   end do
     !!
     !!   if (nfacenodes == 3) n_points = 3
     !!   if (nfacenodes == 6) n_points = 4
     !!   if (nfacenodes == 4) n_points = 4
     !!   if (nfacenodes == 8) n_points = 9
     !!
     !!   call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)
     !!
     !!   do kint = 1,n_points
     !!       call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     !!1                        nfacenodes,N2,dNdxi2)
     !!       dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     !!1                           dNdxi2(1:nfacenodes,1:2))
     !!       norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
     !!       norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
     !!       norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))
     !!
     !!       do i = 1,nfacenodes
     !!           ipoin = 3*face_node_list(i)-2
     !!           RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     !!1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
     !!       end do
     !!   end do
     !! end do

      return

      END SUBROUTINE UEL


      subroutine fung(element_properties,n_properties,F,J,stress,
     1    sigma,D)

       implicit none

       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: element_properties(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: J
       double precision, intent(out) :: stress(6)
       double precision, intent(out) :: D(6,6)

       double precision :: C(3,3), sigma(3,3)
       double precision :: Cinv(3,3),G(6,6),Omega(6,6),Omegaa(6,6)
       double precision :: Cvec(6),Cbarvec(6),Cstarvec(6),Cbarstarvec(6)
       double precision :: Cinvvec(6),Pvec(6)
       double precision :: eyevec(6)
       double precision :: mu
       double precision :: K,Q
       double precision :: ss
       double precision :: trC

       integer :: i

       !  This subroutine calculates the Kirchhoff stress tau = J*cauchy_stress (stored as a vector stress(i) = [tau_11, tau_22, tau_33, etc]
       !  and the tangent matrix D[I,J] = [dtau_11/dB_11, dtau_11/dB_22,... 
       !                                   dtau_22/dB_11, dtau_22/dB_22,
       !                                    etc
       G=0.d0
       Stress=0.d0
       Omegaa=0.d0
       Omega=0.d0
       sigma=0.d0
       mu = element_properties(1)          ! unchanged
       K  = element_properties(2)
       D=0.d0
       C = matmul(transpose(F),F)            !Right Cauchy green tensor
       !call abq_UEL_invert3d(C,Cinv,ss)      ! ss is just a dummy variable here
       call abq_UEL_invert3d(C,Cinv,ss)
       ss = J**(-2.d0/3.d0)
       
       !do i= 1,6
       !    G(i,i)=1.d0                    ! Need correction
       !end do
       G(1,1)=element_properties(3)
       G(2,2)=element_properties(4)
       G(3,3)=element_properties(5)
       G(4,4)=element_properties(6)
       G(5,5)=G(4,4)
       G(6,6)=G(4,4)
       
       !do i = 1,3
       !  stress(i) = mu*B(i,i)*ss
       !end do
       !
       !trB = sum(stress(1:3))/3.d0
       !stress(1:3) = stress(1:3) - trB + K*J*(J-1.d0)
       !stress(4) = mu*B(1,2)*ss
       !stress(5) = mu*B(1,3)*ss
       !stress(6) = mu*B(2,3)*ss
       !D = 0.d0
       !D(1,1) = 1.d0
       !D(2,2) = 1.d0
       !D(3,3) = 1.d0
       !D(4,4) = 0.5d0
       !D(5,5) = 0.5d0
       !D(6,6) = 0.5d0
       !D = D*mu*ss

       eyevec(1:3) = 1.d0                     !CHECK?
       eyevec(4:6) = 0.d0
       Cvec(1) = C(1,1)
       Cvec(2) = C(2,2)
       Cvec(3) = C(3,3)
       Cvec(4) = C(1,2)
       Cvec(5) = C(1,3)
       Cvec(6) = C(2,3)
       
       Cbarvec(1:6)=Cvec(1:6)*ss
       
       Cstarvec(1) = C(1,1)
       Cstarvec(2) = C(2,2)
       Cstarvec(3) = C(3,3)
       Cstarvec(4) = 2*C(1,2)
       Cstarvec(5) = 2*C(1,3)
       Cstarvec(6) = 2*C(2,3)
       
       Cbarstarvec(1:6)=Cstarvec(1:6)*ss
       
       Cinvvec(1) = Cinv(1,1)
       Cinvvec(2) = Cinv(2,2)
       Cinvvec(3) = Cinv(3,3)
       Cinvvec(4) = Cinv(1,2)
       Cinvvec(5) = Cinv(1,3)
       Cinvvec(6) = Cinv(2,3)
       
       Pvec=ss/2 * (matmul(G,(Cbarstarvec-eyevec)) 
     1 - 1.d0/3*(dot_product(Cstarvec,matmul(G,(Cbarstarvec-eyevec)))
     2  *Cinvvec))                                                   ! P vector defined 

       Q=1.d0/4 *dot_product((Cbarvec-eyevec),
     1  matmul(G,(Cbarvec-eyevec)))
       
       Stress = mu*exp(Q)*Pvec + K*J*(J-1)*Cinvvec          ! Stress is defined, finally :)
       sigma(1,1)=Stress(1)                                          ! storing it as a matrix
       sigma(1,2)=Stress(4)                                 !CHECK?
       sigma(1,3)=Stress(5)
       sigma(2,1)=sigma(1,2)
       sigma(2,2)=Stress(2)
       sigma(2,3)=Stress(6)
       sigma(3,1)=sigma(1,3)
       sigma(3,2)=sigma(2,3)
       sigma(3,3)=Stress(3)
       
       
       Omegaa(1,1)=Cinv(1,1)*Cinv(1,1)/2.d0
       Omegaa(1,2)=Cinv(1,2)*Cinv(1,2)
       Omegaa(1,3)=Cinv(1,3)*Cinv(1,3)
       Omegaa(1,4)=Cinv(1,1)*Cinv(1,2)
       Omegaa(1,5)=Cinv(1,1)*Cinv(1,3)
       Omegaa(1,6)=Cinv(1,2)*Cinv(1,3)       !1row complete
       Omegaa(2,2)=Cinv(2,2)*Cinv(2,2)/2.d0
       Omegaa(2,3)=Cinv(2,3)*Cinv(2,3)    
       Omegaa(2,4)=Cinv(2,1)*Cinv(2,2)
       Omegaa(2,5)=Cinv(2,1)*Cinv(2,3)
       Omegaa(2,6)=Cinv(2,2)*Cinv(2,3)       !2nd complete
       Omegaa(3,3)=Cinv(3,3)*Cinv(3,3)/2.d0
       Omegaa(3,4)=Cinv(3,1)*Cinv(3,2)
       Omegaa(3,5)=Cinv(3,1)*Cinv(3,3)
       Omegaa(3,6)=Cinv(3,2)*Cinv(3,3)       !3rd complete
       Omegaa(4,4)=1.d0/4*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(1,2))
       Omegaa(4,5)=1.d0/2*(Cinv(1,1)*Cinv(2,3)+Cinv(1,3)*Cinv(1,2))
       Omegaa(4,6)=1.d0/2*(Cinv(1,2)*Cinv(3,3)+Cinv(1,3)*Cinv(2,3))  !4th complete
       Omegaa(5,5)=1.d0/4*(Cinv(1,1)*Cinv(3,3)+Cinv(1,3)*Cinv(1,3))
       Omegaa(5,6)=1.d0/2*(Cinv(1,2)*Cinv(3,3)+Cinv(1,3)*Cinv(2,3))  !5th complete
       Omegaa(6,6)=1.d0/4*(Cinv(2,2)*Cinv(3,3)+Cinv(2,3)*Cinv(2,3))  !6th complete
       
       Omega = -((Omegaa+transpose(Omegaa)))           ! Omega is defined
       
       
       D=D+mu*exp(Q)*(ss**(2.d0)*(G-1.d0/3*(matmul(G,
     1  spread(Cstarvec,dim=2,ncopies=6)*spread(Cinvvec,
     2  dim=1,ncopies=6))+spread(Cinvvec,dim=2,ncopies=6)*spread(matmul
     3  (G,Cstarvec),dim=1,ncopies=6))-ss**(-1.d0)/3* dot_product( 
     4  Cstarvec,matmul(G,(Cbarstarvec-eyevec)))*Omega +1.d0/9*     !check if transpose needed?
     5  (dot_product(Cstarvec,matmul(G, Cstarvec)))*spread(Cinvvec,dim=2     !check again
     6  ,ncopies=6)*spread(Cinvvec,dim=1,ncopies=6)))           ! 1st line complete
     7  + mu*exp(Q)*(spread(2*Pvec,dim=2,ncopies=6)
     8  *spread((Pvec-1.d0/3*Cinvvec),dim=1,ncopies=6)- ss/3 * 
     9  spread(Cinvvec,dim=2,ncopies=6) *spread(matmul(G,       ! Check if (C-I) is OK? 
     1  (Cbarstarvec-eyevec)),dim=1, ncopies=6))                   !2nd term completed 
     2  +K*J*((2*J-1)*spread(Cinvvec,dim=2,ncopies=6)
     3  *spread(Cinvvec,dim=1,ncopies=6)+2*(J-1)* Omega)        ! 3rd term completed and D is defined, BC.  
     !!  !trB = sum(Bvec(1:3))/3.d0 + 2*(J-1)* Omega)
     !!  D=K*J*((2*J-1)*spread(Cinvvec,dim=2,ncopies=6)
     !!1  *spread(Cinvvec,dim=1,ncopies=6))
     !!  D = D + (ss*mu/3.d0)*( trB*spread(eyevec,dim=2,ncopies=6)*
     !!1                                   spread(Binvvec,dim=1,ncopies=6) 
     !!2              - spread(eyevec,dim=2,ncopies=6)*
     !!3                                    spread(eyevec,dim=1,ncopies=6)  
     !!4              - spread(Bvec,dim=2,ncopies=6)*
     !!5                                 spread(Binvvec,dim=1,ncopies=6) )
     !!
     !!  D = D + K*J*(J-0.5d0)*spread(eyevec,dim=2,ncopies=6)*
     !!1                                   spread(Binvvec,dim=1,ncopies=6)


       return

      end subroutine fung

      
      
      
      !Homework 5 
    ! Local Variables
    !!  integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
    !!  integer      :: face_node_list(3)                       ! List of nodes on an element face
    !!!
    !!  double precision  ::  xi(2,9)                          ! Area integration points
    !!  double precision  ::  w(9)                             ! Area integration weights
    !!  double precision  ::  N(9)                             ! 2D shape functions
    !!  double precision  ::  dNdxi(9,2),dNdxi0(9,2)                        ! 2D shape function derivatives
    !!  double precision  ::  dNdx(9,2)                        ! Spatial derivatives
    !!  double precision  ::  dxdxi(2,2), dxdxi0(2,2)                       ! Derivative of spatial coords wrt normalized coords
    !!
    !!!   Variables below are for computing integrals over element faces
    !!  double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
    !!  double precision  ::  xi1(6)                            ! 1D integration points
    !!  double precision  ::  w1(6)                              ! Integration weights
    !!  double precision  ::  N1(3)                             ! 1D shape functions
    !!  double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
    !!  double precision  ::  norm(2)                           ! Normal to an element face
    !!  double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !!!
    !!  double precision  ::  strain(4)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    !!  double precision  ::  stress(4)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    !!  double precision  ::  D(4,4)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
    !!  double precision  ::  B(4,22)                           ! strain = B*(dof_total)
    !!  double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
    !!  double precision  ::  rhs_temp(22)                      ! Temporary RHS vector (for incompatible mode elements)
    !!  double precision  ::  kuu(18,18)                        ! Upper diagonal stiffness
    !!  double precision  ::  kaa(4,4),kaainv(4,4)              ! Lower diagonal stiffness
    !!  double precision  ::  kau(4,18)                         ! Lower quadrant of stiffness
    !!  double precision  ::  kua(18,4)                         ! Upper quadrant of stiffness
    !!  double precision  ::  alpha(4)                          ! Internal DOF for incompatible mode element
    !!  double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    !!  double precision  ::  E, xnu, D44, D11, D12             ! Material properties
    !!  double precision  ::  det0, center(2,1)                 ! determinant at xi=0, center denotes {xi1=0,xi2=0}
    !!  double precision  ::  BigU(22,1)   
    !!  double precision  ::  ru(18,1), ra(4,1)
    !!  !
    !!!     Example ABAQUS UEL implementing 2D linear elastic elements
    !!!     Includes option for incompatible mode elements
    !!!     El props are:
    !!
    !!!     PROPS(1)         Young's modulus
    !!!     PROPS(2)         Poisson's ratio
    !!
    !!
    !!  if (NNODE == 3) n_points = 1              ! Linear triangle
    !!  if (NNODE == 4) n_points = 4               ! Linear rectangle
    !!  if (NNODE == 6) n_points = 4              ! Quadratic triangle
    !!  if (NNODE == 8) n_points = 9               ! Serendipity rectangle
    !!  if (NNODE == 9) n_points = 9             ! Quadratic rect
    !!
    !!! Write your code for a 2D element below
    !!!  write(6,*)'nnode'
    !!  call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)
    !!
    !!  if (MLVARX<2*NNODE) then                                    
    !!    write(6,*) ' Error in abaqus UEL '
    !!    write(6,*) ' Variable MLVARX must exceed 3*NNODE'
    !!    write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
    !!    stop
    !!  endif
    !!
    !!  RHS(1:MLVARX,1) = 0.d0
    !!  AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
    !!  ktemp(1:22,1:22)=0.d0
    !!  rhs_temp(1:22)=0.d0
    !!  alpha(1:4)=0.d0                                     !Initializing the local variables
    !!  kuu(1:18,1:18)=0.d0
    !!  kua(1:18,1:4)=0.d0
    !!  kau(1:4,1:18)=0.d0
    !!  kaa(1:4,1:4)=0.d0
    !!  ru(1:2*NNODE,1)=0.d0
    !!  ra(1:4,1)=0.d0
    !!  
    !!  D = 0.d0
    !!  E = PROPS(1)
    !!  xnu = PROPS(2)
    !!  d44 = 0.5D0*E/(1+xnu)
    !!  d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    !!  d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    !!  D(1:3,1:3) = d12
    !!  D(1,1) = d11
    !!  D(2,2) = d11
    !!  D(3,3) = d11
    !!  D(4,4) = d44
    !!  !D(5,5) = d44
    !!  !D(6,6) = d44
    !!  
    !!  ENERGY(1:8) = 0.d0
    !!  BigU (1:22,1) = 0.d0
    !!  det0=0.d0
    !!  center(1:2,1)=0.d0
    !!!     --  Loop over integration 
    !!  if(JTYPE==1) then              !Conventional elements
    !!  
    !!  do kint = 1, n_points
    !!    call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)      !changed
    !!    dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))          !changed
    !!    call abq_inverse_LU(dxdxi,dxidx,2)                              !changed
    !!    determinant = dxdxi(1,1)*dxdxi(2,2)- dxdxi(1,2)*dxdxi(2,1)      
    !!    dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)          
    !!    B = 0.d0
    !!    B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)                          
    !!    B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)                                                         
    !!    !B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)                           !zeros
    !!    B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)                            
    !!    B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)                              
    !!   ! B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
    !!   ! B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)                           
    !!    !B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
    !!    !B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)
    !!
    !!    strain = matmul(B(1:4,1:2*NNODE),U(1:2*NNODE))                  
    !! 
    !!    stress = matmul(D,strain)                                       
    !!    RHS(1:2*NNODE,1) = RHS(1:2*NNODE,1)
    !! 1   - matmul(transpose(B(1:4,1:2*NNODE)),stress(1:4))*           
    !! 2                                          w(kint)*determinant
    !! 
    !!    AMATRX(1:2*NNODE,1:2*NNODE) = AMATRX(1:2*NNODE,1:2*NNODE)       
    !! 1  + matmul(transpose(B(1:4,1:2*NNODE)),matmul(D,B(1:4,1:2*NNODE)))
    !! 2                                             *w(kint)*determinant
    !! 
    !!    ENERGY(2) = ENERGY(2)
    !! 1   + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy
    !! 
    !!    if (NSVARS>=n_points*4) then                                    ! Store stress at each integration point (if space was allocated to do so)
    !!        SVARS(4*kint-3:4*kint) = stress(1:4)                  
    !!    endif
    !!  end do
    !!  endif
    !!  
    !!  
    !!  
    !!  if (JTYPE==2) then                 ! Incompatible mode elements! 
    !!  do kint = 1, n_points
    !!  call abq_UEL_2D_shapefunctions(center(1:2,1),NNODE,N,dNdxi0)  
    !!  call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)      
    !!    dxdxi = matmul(coords(1:2,1:NNODE), dNdxi(1:NNODE,1:2))        
    !!    dxdxi0 = matmul(coords(1:2,1:NNODE),dNdxi0(1:NNODE,1:2))
    !!    
    !!    call abq_inverse_LU(dxdxi,dxidx,2)                              !changed
    !!    determinant = dxdxi(1,1)*dxdxi(2,2)- dxdxi(1,2)*dxdxi(2,1)      !a line added
    !!    det0 = dxdxi0(1,1)*dxdxi0(2,2)- dxdxi0(1,2)*dxdxi0(2,1) 
    !!    
    !!    dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)            !changed
    !!    B = 0.d0
    !!    B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)                            !changed
    !!    B(1, 2*NNODE+2:2*NNODE+4:2) = 0                                 ! Extra DOF (3 lines)
    !!    B(1, 2*NNODE+1) = det0/determinant *xi(1,kint)*dxidx(1,1)
    !!    B(1, 2*NNODE+3) = det0/determinant *xi(2,kint)*dxidx(2,1)
    !!    
    !!    B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)                              !changed                            
    !!    B(2, 2*NNODE+1:2*NNODE+3:2) = 0                                 ! Extra DOF (3 lines added)
    !!    B(2, 2*NNODE+2) = det0/determinant*xi(1,kint)*dxidx(1,2)
    !!    B(2, 2*NNODE+4) = det0/determinant*xi(2,kint)*dxidx(2,2)  
    !!   
    !!    !B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)                           !zeros
    !!    
    !!    B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)                            
    !!    B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)                              
    !!    
    !!     B(4, 2*NNODE+1) = det0/determinant *xi(1,kint)*dxidx(1,2)      ! Extra DOF (4 lines added)
    !!     B(4, 2*NNODE+2) = det0/determinant *xi(1,kint)*dxidx(1,1)
    !!     B(4, 2*NNODE+3) = det0/determinant *xi(2,kint)*dxidx(2,2)
    !!     B(4, 2*NNODE+4) = det0/determinant *xi(2,kint)*dxidx(2,1) 
    !!    
    !!    
    !!    ! B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
    !!   ! B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)                           
    !!    !B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
    !!    !B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)
    !!
    !!   !strain = matmul(B(1:4,1:2*NNODE),U(1:2*NNODE))                  
    !! 
    !!    !stress = matmul(D,strain)                                       
    !!    !RHS(1:2*NNODE,1) = RHS(1:2*NNODE,1)
    !! 1   !- matmul(transpose(B(1:4,1:2*NNODE)),stress(1:4))*             !changed
    !! 2    !                                      w(kint)*determinant
    !! 
    !!    ktemp(1:2*NNODE+4,1:2*NNODE+4) = ktemp(1:2*NNODE+4,1:2*NNODE+4)           !K_hat    Elemental stiffness 2n+4 cross 2n+4
    !! 1 + matmul(transpose(B(1:4,1:2*NNODE+4))
    !! 2 ,matmul(D,B(1:4,1:2*NNODE+4)))*w(kint)*determinant
    !! 
    !!    !ENERGY(2) = ENERGY(2)
    !! 1   !+ 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy
    !! 
    !!    !if (NSVARS>=n_points*4) then                                      !changed    ! Store stress at each integration point (if space was allocated to do so)
    !!     !   SVARS(4*kint-3:4*kint) = stress(1:4)                          !changed
    !!    !endif
    !!  end do
    !!  
    !!  
    !!  kuu=ktemp(1:2*NNODE,1:2*NNODE);      
    !!  kua=ktemp(1:2*NNODE, 2*NNODE+1:2*NNODE+4)
    !!  kau=ktemp(2*NNODE+1:2*NNODE+4, 1:2*NNODE)
    !!  kaa=ktemp(2*NNODE+1:2*NNODE+4, 2*NNODE+1:2*NNODE+4)
    !!  call abq_inverse_LU(kaa,kaainv,4) 
    !!  alpha=-matmul(kaainv,matmul(kau(1:4,1:2*NNODE),U))
    !!  
    !!  BigU(1:2*NNODE,1)=U(1:2*NNODE)
    !!  BigU(2*NNODE+1:2*NNODE+4,1) = alpha(1:4)
    !!  
    !!  
    !!  do kint = 1, n_points
    !!   call abq_UEL_2D_shapefunctions(center(1:2,1),NNODE,N,dNdxi0)  
    !!   call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)      !changed
    !!    dxdxi = matmul(coords(1:2,1:NNODE), dNdxi(1:NNODE,1:2))          !changed
    !!    dxdxi0 = matmul(coords(1:2,1:NNODE),dNdxi0(1:NNODE,1:2))
    !!    
    !!    call abq_inverse_LU(dxdxi,dxidx,2)                              !changed
    !!    determinant = dxdxi(1,1)*dxdxi(2,2)- dxdxi(1,2)*dxdxi(2,1)     !a line added
    !!    det0 = dxdxi0(1,1)*dxdxi0(2,2)- dxdxi0(1,2)*dxdxi0(2,1) 
    !!    
    !!    dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)            !changed
    !!    B = 0.d0
    !!    B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)                            !changed
    !!    B(1, 2*NNODE+2:2*NNODE+4:2) = 0                                ! Extra DOF (3 lines)
    !!    B(1, 2*NNODE+1) = det0/determinant *xi(1,kint)*dxidx(1,1)
    !!    B(1, 2*NNODE+3) = det0/determinant *xi(2,kint)*dxidx(2,1)
    !!    
    !!    B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)                              !changed                            
    !!    B(2, 2*NNODE+1:2*NNODE+3:2) = 0                                ! Extra DOF (3 lines added)
    !!    B(2, 2*NNODE+2) = det0/determinant*xi(1,kint)*dxidx(1,2)
    !!    B(2, 2*NNODE+4) = det0/determinant*xi(2,kint)*dxidx(2,2)  
    !!   
    !!    !B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)                           !zeros
    !!    
    !!    B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)                            !changed
    !!    B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)                              !changed
    !!    
    !!     B(4, 2*NNODE+1) = det0/determinant *xi(1,kint)*dxidx(1,2)             ! Extra DOF (4 lines added)
    !!     B(4, 2*NNODE+2) = det0/determinant *xi(1,kint)*dxidx(1,1)
    !!     B(4, 2*NNODE+3) = det0/determinant *xi(2,kint)*dxidx(2,2)
    !!     B(4, 2*NNODE+4) = det0/determinant *xi(2,kint)*dxidx(2,1) 
    !!  
    !!     strain = matmul(B(1:4,1:2*NNODE+4),BigU(1:2*NNODE+4,1))                  !changed
    !! 
    !!     stress = matmul(D,strain)     
    !!    ! 
    !!      if (NSVARS>=n_points*4) then                                  ! Store stress at each integration point (if space was allocated to do so)
    !!        SVARS(4*kint-3:4*kint) = stress(1:4)                          
    !!    endif
    !!     
    !!     rhs_temp(1:2*NNODE+4) = rhs_temp(1:2*NNODE+4)
    !! 1   - matmul(transpose(B(1:4,1:2*NNODE+4)),stress(1:4))*          !changed      1 cross 2n+4
    !! 2                                          w(kint)*determinant
    !!  
    !!  end do
    !!  
    !!  ru(1:2*NNODE,1) =  rhs_temp(1:2*NNODE)
    !!  ra(1:4,1) =  rhs_temp(2*NNODE+1:2*NNODE+4)
    !!  
    !!  ! Assembling
    !!  
    !!  AMATRX(1:2*NNODE, 1:2*NNODE)= kuu (1:2*NNODE, 1:2*NNODE) 
    !! 1- matmul(kua(1:2*NNODE, 1:4), matmul(kaainv, kau(1:4,1:2*NNODE)))
    !!  
    !!  RHS(1:2*NNODE,1) = ru(1:2*NNODE,1) 
    !! 1- matmul(kua(1:2*NNODE, 1:4),matmul(kaainv,ra(1:4 , 1)))
    !!  
    !!  PNEWDT = 1.d0 
    !!  endif
    !!  return 
    !!
    !!  END SUBROUTINE UEL_2D



      subroutine abq_UEL_2D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_UEL_2D_integrationpoints




      subroutine abq_UEL_2D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2

            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if

      end subroutine abq_UEL_2D_shapefunctions


      subroutine abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)


      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)



      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in subroutine abq_UEL_1D_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select







      end subroutine ABQ_UEL_1D_integrationpoints



      subroutine abq_facenodes_2D(nelnodes,face,list,nfacenodes)

      implicit none

      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes
    !
    !        Subroutine to return list of nodes on an element face for standard 2D solid elements
    !
      integer :: i3(3)
      integer :: i4(4)

      i3(1:3) = [2,3,1]
      i4(1:4) = [2,3,4,1]

      if (nelnodes == 3) then
        nfacenodes = 2
        list(1) = face
        list(2) = i3(face)
      else if (nelnodes == 4) then
        nfacenodes = 2
        list(1) = face
        list(2) = i4(face)
      else if (nelnodes == 6) then
        nfacenodes = 3
        list(1) = face
        list(2) = i3(face)
        list(3) = face+3
      else if (nelnodes == 8) then
        nfacenodes = 3
        list(1) = face
        list(2) = i4(face)
        list(3) = face+4
      else if (nelnodes == 9) then
        nfacenodes = 3
        if (face==1) list(1:3) = (/1,3,2/)
        if (face==2) list(1:3) = (/3,9,6/)
        if (face==3) list(1:3) = (/9,7,8/)
        if (face==4) list(1:3) = (/7,1,4/)
      endif

      end subroutine abq_facenodes_2d

      subroutine abq_inverse_LU(Ain,A_inverse,n)  ! Compute the inverse of an arbitrary matrix by LU decomposition

        implicit none

        integer, intent(in)  :: n

        double precision, intent(in)    :: Ain(n,n)
        double precision, intent(out)   :: A_inverse(n,n)

        double precision :: A(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
        double precision :: coeff
        integer :: i, j, k

        A(1:n,1:n) = Ain(1:n,1:n)
        L=0.d0
        U=0.d0
        b=0.d0

        do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
            end do
        end do

        forall (i=1:n)  L(i,i) = 1.d0
        forall (j=1:n) U(1:j,j) = A(1:j,j)
        do k=1,n
            b(k)=1.d0
            d(1) = b(1)
            do i=2,n
                d(i)=b(i)
                d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
            end do
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                x(i)=x(i)-dot_product(U(i,i+1:n),x(i+1:n))
                x(i) = x(i)/U(i,i)
            end do
            A_inverse(1:n,k) = x(1:n)
            b(k)=0.d0
        end do

      end subroutine abq_inverse_LU


