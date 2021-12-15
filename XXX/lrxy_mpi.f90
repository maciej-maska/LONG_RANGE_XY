module configuration
    implicit none
    real(8), parameter ::pi = 3.14159265359_8, pi2 = 6.28318530718_8
    logical :: pbc
    integer :: N
    real(8), dimension (:), allocatable :: theta   
    real(8) :: temperature
end module configuration

program lr_xy
    use configuration
    !use lapack95
    implicit none
    include "mpif.h"
    integer :: i, j, k, n_warmup, mcs, isite, ierr, myrank, ss, len, status
    integer, dimension(:), allocatable :: seed
    real(8) :: beta, r, energy, e1, e2, de, x, theta_old, en, en2, nen, een, pars(60), ee(3), stiffness, stiff, stiff1
    character(len=10) :: arg1,arg2

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call RANDOM_SEED(SIZE=n)
    allocate(seed(n))
    call system_clock(ss)
    call random_seed(get = seed)
    seed = myrank * seed + abs( mod((ss*181), 104729))
    call random_seed(put = seed)


    call get_command_argument(1, arg1, len, status)
    read(arg1,'(i4)') N
    call get_command_argument(2, arg2, len, status)
    if (arg2 == 'P') then
        pbc = .true.
    else
        pbc = .false.
    endif
    allocate(theta(N))

    temperature = 0.2 + 0.05 * myrank ! 20 procesow => 0.2 <= T <= 1.2
    beta = 1._8/temperature
    n_warmup = 1000000
    een = 0._8
    en2 = 0._8
    stiff1 = 0._8
    do k = 1, 10
        if (pbc) then
            theta = 0
        else
            do i = 1, N
                 theta(i) = 0.5 * (i - 1)/(N - 1)
            enddo
        endif
        theta = pi2 * theta
        e1 = energy()
        en = 0._8
        nen = 0
        stiff = 0._8
        do mcs = 1, 2000000
            call random_number(x)
            isite = int((N-2)*x+2)     ! tylko wewnętrzne węzły sa losowane
            theta_old = theta(isite)
            call random_number(theta(isite))
            theta(isite) = pi2 * theta(isite)
            e2 = energy()
            if (e2 < e1) then
                e1 = e2
            else 
                call random_number(x)
                if (-beta*(e2-e1) > log(x)) then
                    e1 = e2
                else
                    theta(isite) = theta_old
                endif
            endif
            if (mcs > n_warmup .and. mod(mcs,10*N) == 0) then
                en = en + e1
                nen = nen + 1
                if (pbc) stiff = stiff + stiffness()  ! sztywność musi być liczona przy PBC
            endif
        enddo
        een = een + en/nen
        en2 = en2 + (en/nen)**2
        stiff1 = stiff1 + stiff/nen
enddo
    ee(1) = een/10                  ! energy averaged over 10 MC runs
    ee(2) = sqrt(en2/10-ee(1)**2)   ! energy fluctuations between different runs
    ee(3) = stiff1
    call MPI_Gather(ee, 3, MPI_DOUBLE, pars, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
        open(10,file='en_'//trim(arg2)//'_'//trim(arg1)//'.dat')
        do i = 1,20
            write(10,'(3f20.10)') 0.2 + 0.05*(i-1),pars(3*(i-1)+1),pars(3*(i-1)+2),pars(3*(i-1)+3)
        enddo
        open(11,file='spins_'//trim(arg2)//'.dat')
        do i=1,N
            write(11,*) theta(i)
        enddo
        close(11)
    endif
    call MPI_Finalize(ierr)
end program lr_xy

function energy()
    use configuration
    implicit none
    integer :: i,j
    real(8) :: energy, theta1
    energy = 0._8
    do j = 2, N
        theta1 = theta(j)
        do i = 1, j - 1
            energy = energy - (cos(theta1 - theta(i)))/(abs(j - i)**2)
        enddo
    enddo
end function energy

function stiffness()
    use configuration
    implicit none
    integer :: i,j
    real(8) :: h1, current, theta1, stiffness
    h1 = 0._8
    current = 0._8
    do j = 2, N
        theta1 = theta(j)
        do i = 1, j - 1

            h1 = h1 + (cos(theta1 - theta(i)))
            current = current + sin(theta1 - theta(i))/(j - i)
        enddo
    enddo
    stiffness = h1 / N**2 - current**2 / (temperature*N*N)
end function stiffness