
program dg_1d
    implicit none 
    integer, parameter:: ngi=2,nloc=2,nele=100,nnod=nele+1
    real :: Mmat(nnod,nnod), Gmat(nnod,nnod), f(nnod), &
            phin(nnod), phin1(nnod) ! phin1 is phi^(n-1) - last step value
    real :: A(nnod,nnod), b(nnod), rhs(nnod) ! A x = rhs, rhs = b + explicit term
    real :: Lend(2) ! coordinates of domain ends
    real :: coor(nnod) ! coordinates of nodal points
    real :: elength ! length of element
    integer :: inod, jnod, & ! local node index
            globi, globj, & ! global node index
            ele, & ! element index
            gi ! gaussian points index
    real :: sf(nloc,ngi), sfdev(nloc,ngi) ! shape function and shape function derivatives
    real :: detwei(ngi) ! determinant of Jacobian x gaussian point weight
    real :: nn, nxnx
    integer :: ts, tstot ! timestep index, total timesteps
    real :: tend , tstep    ! end time and time step size

    ! geometry
    Lend = (/0., 1./)
    elength = (Lend(2) - Lend(1))/nele 
    do globi = 1, nnod
        coor(globi) = Lend(1) + real(globi-1) * elength
    enddo

    ! time
    tend = 2.
    tstep = 0.01
    tstot = ceiling(tend/tstep)

    ! shape function on reference element 
    ! use 1st order quadrature, 2 point, weight = 1
    sf(1,1)=(1.+1./sqrt(3.))/2.;    sf(1,2)=(1.-1./sqrt(3.))/2.;
    sf(2,1)=(1.-1./sqrt(3.))/2.;    sf(2,2)=(1.+1./sqrt(3.))/2.;
    ! derivatives of shape function on local elements
    ! because this is uniform mesh, sfdev is same for every element. we only compute once
    sfdev(1,1) = -0.5*2.0/elength;  sfdev(1,2) = -0.5*2.0/elength;   
    sfdev(2,1) = 0.5*2.0/elength;   sfdev(2,2) = 0.5*2.0/elength;
    detwei(1) = elength /2.0;       detwei(2) = elength / 2.0;

    ! assemble the matrix
    Mmat = 0.
    Gmat = 0.
    do ele = 1,nele 
        do inod = 1,nloc 
            globi = ele-1+inod 
            do jnod = 1,nloc 
                globj = ele-1+jnod 
                nn = .0
                nxnx = .0
                do gi = 1,ngi
                    ! print*, nxnx
                    nn = nn + sf(inod,gi) * sf(jnod,gi) * detwei(gi)
                    nxnx = nxnx + sfdev(inod, gi) * sfdev(jnod, gi) * detwei(gi)
                    ! print* ,sfdev(inod,gi), sfdev(jnod,gi), detwei(gi), nxnx
                enddo
                Mmat(globi, globj) = Mmat(globi, globj) + nn 
                Gmat(globi, globj) = Gmat(globi, globj) + nxnx
                ! print*, inod, jnod, nn, nxnx
            enddo
        enddo
    enddo

    ! rhs 
    f = 1.
    b = matmul(Mmat, f) 
    ! print*, b

    
    ! time loop
    phin1 = 0.  ! initial condition
    phin  = 0.  ! initiallizing vector
    
    open(unit=30, file='results.out', status='replace')
    do ts = 1, tstot
        A = Mmat/tstep + Gmat
        if (ts.eq.tstot) then 
            open(unit=10, file='matA.out', status='replace')
            open(unit=20, file='rhs.out', status='replace')
            do globi = 1,nnod 
                write(10,*) A(globi,:)
                write(20,*) rhs(globi)
            enddo
            close(10); close(20);
        endif
        rhs = b + matmul(Mmat, phin1) / tstep 
        ! Dirichlet boundary (assuming both ends are Diri. b.c.)
        A(1,:)=0.;       A(:,1)=0.;       A(1,1)=1.;
        A(nnod,:)=0.;    A(:,nnod)=0.;    A(nnod,nnod)=1.;
        rhs(1)=0.;            rhs(nnod)=0.;

        call linsol(A, rhs, phin, nnod)
        ! write out 
        write(30,*) ts*tstep, phin
        phin1 = phin 
    enddo
    close(30)

    ! ! let's try steady case

    ! open(unit=10, file='Gmat.out', status='replace')
    ! open(unit=20, file='rhs.out', status='replace')
    ! open(unit=30, file='Mmat.out', status='replace')
    ! do globi = 1,nnod 
    !     write(10,*) Gmat(globi,:)
    !     write(20,*) b(globi)
    !     write(30,*) Mmat(globi,:)
    ! enddo
    ! close(10); close(20); close(30)

    ! call linsol(Gmat, b, phin, nnod)
    ! print*, phin



end program


subroutine linsol(A, b, x, nnod)
    implicit none

    integer, intent(in):: nnod
    real, intent(inout)::A(nnod,nnod), b(nnod)
    real, intent(inout)::x(nnod)
    integer :: i 
    real :: w

    A(1,2) = A(1,2)/A(1,1)
    do i = 2, nnod
        w = A(i,i-1)/A(i-1,i-1)
        A(i,i) = A(i,i) - w * A(i-1,i)
        b(i) = b(i) - w * b(i-1)
    enddo
    x(nnod) = b(nnod) / A(nnod, nnod)
    do i = nnod-1, 1, -1
        x(i) = (b(i) - A(i,i+1) * x(i+1) ) / A(i,i)
    enddo
endsubroutine