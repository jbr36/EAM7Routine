
SUBROUTINE NIHEAM7(XC, DX, ENERGY) 
    ! *****************************************************************************************
    ! XC input coordinates
    ! DX output gradients
    ! Energy output
    ! *****************************************************************************************
    ! 
    ! This is the EAM7 Ni(100)/H for the interaction of a single H atom with a Ni(100)
    ! written by Judith B Rommel, if used cite:
    !
    ! M. S. Daw and M. I. Baskes, Phys. Rev. B 29, 6443 (1984).
    ! Wonchoba, S. E., Truhlar, D. G. (1996) Phys. Rev. B 53:11222-11241.
    ! REV B 51, 9985 (1995)
    ! James T. Kindt and John C. Tully, J. Chem. Phys. 111, 11060 (1999)
    !
    ! Input/output geometry in Angstrom, energies in eV, and gradients in eV/Angstrom!
    ! This version is for 1 hydrogen atom only!
    ! *****************************************************************************************

    USE MODHESS
    USE COMMONS, ONLY: NATOMS, ZSYM
    USE KEY, ONLY: FROZEN, NFREEZE

    IMPLICIT NONE
    DOUBLE PRECISION :: XC(3*NATOMS), DX(3*NATOMS),DDX(3*NATOMS,3*NATOMS) !Cartesian coords, 1st, 2nd derivatives
    DOUBLE PRECISION :: ENERGY
    integer          :: nclassic,nrigid,nb,i1,i2
    integer :: iflag ! -1 - coordinates to Angs, +1 - coordinates, energy and forces to a.u
    double precision :: boxlx,boxly,boxlz
    double precision :: qh(3,1,1),dvdqh(3,1,1)
    double precision :: qmclassic(3,NATOMS-NFREEZE-1,1),dvdqmclassic(3,NATOMS-NFREEZE-1,1)
    double precision :: qmrigid(3,NFREEZE,1),dvdqmrigid(3,NFREEZE,1)
    character(len=6)      :: boundary
    logical :: h_included,GTEST, STEST
    ! *****************************************************************************************
    !Description of the variables:
    !
    !nclassic - number of classical Ni atoms
    !nrigid - number of rigid Ni atoms
    !nb - number of the ring polymer beads
    !boxlx - length of the simulation box in x direction
    !boxly - length of the simulation box in y direction
    !boxlz - length of the simulation box in z direction
    !boundary - 'afixed', 'bfixed', or 'cfixed'. 'cfixed' or 'bfixed'- the lattices with 2D periodic boundary conditions
    !h_included - '.true.' or  '.false.'  for hydrogen (whether it is on the surface or not)
    !qmclassic  - array for classic Ni atoms coordinates
    !qmrigid  - array for rigid Ni atoms coordinates
    !qh - array for H atom coordinates
    !v - interaction potential
    !dvdqmclassic - array for classic Ni atoms forces
    !dvdqmrigid  - array for rigid Ni atoms forces
    !dvdqh  - array for H atom forces
    ! *****************************************************************************************


    nclassic=NATOMS-NFREEZE-1 !natoms-nfrozen !number of mobile atoms not including H
    nrigid=NFREEZE!nfrozen ! number of frozen atoms
    nb=1

    !boxlength (caluclated based on approx 2 times the atom radius 4.7 of NI in a.u.)
    boxlx = 24.64d0 ! 65.8d0
    boxly = 24.64d0 !65.8d0
    boxlz = 7.04d0 !14.1d0

    !boundary - 'afixed', 'bfixed', or 'cfixed'. 'cfixed' or 'bfixed'- the lattices with 2D periodic boundary conditions
    boundary='bfixed'
    !h_included - '.true.' or  '.false.'  for hydrogen (whether it is on the surface or not)
    h_included =.false.


    !initialise coords and forces
    qh(:,:,:) =0.d0
    dvdqh(:,:,:) =0.d0
    qmclassic(:,:,:)=0.D0
    dvdqmclassic(:,:,:)=0.D0
    qmrigid(:,:,:)=0.D0
    dvdqmrigid(:,:,:)=0.D0
    ENERGY = 0.0D0
    DX(:)=0
    DDX(:,:)=0


    ! The geometry of the mobile classical atoms
    DO  I1 = 1,nclassic+1
        IF ((ZSYM(I1).EQ.'NI').or.(ZSYM(I1).EQ.'Ni')) THEN
            qmclassic(1,I1-1,1)= XC(3*I1-2)!I1-1 because of H excluded
            qmclassic(2,I1-1,1)= XC(3*I1-1)
            qmclassic(3,I1-1,1)= XC(3*I1)
        ELSE IF (ZSYM(I1).EQ.'H') THEN
            qh(1,I1,1)= XC(3*I1-2)
            qh(2,I1,1)= XC(3*I1-1)
            qh(3,I1,1)= XC(3*I1)
            h_included =.true.
        END IF
    ENDDO


    ! The geometry of the rigid atoms
    IF (nclassic+2 .LE. NATOMS) THEN
        DO  I1 = nclassic+2,NATOMS
            IF (ZSYM(I1).EQ.'NI') THEN
                qmrigid(1,I1-nclassic-1,1)= XC(3*I1-2)!I1-nclassic-1 because of H excluded
                qmrigid(2,I1-nclassic-1,1)= XC(3*I1-1)
                qmrigid(3,I1-nclassic-1,1)= XC(3*I1)
            ELSE
                PRINT*, 'There are gas atoms in frozen zone'
                STOP
            END IF
        ENDDO
    ENDIF


    IF (GTEST) THEN
        ! PRINT*, 'Calculating energy and forces.  One moment please ...'

        !     Calculate the new forces and the interaction energy
        call EAM_EAMNIH_FORCES(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
            boundary,h_included,qmclassic,qmrigid,&
            qh,ENERGY,dvdqmclassic,dvdqmrigid,dvdqh)

        !print*,dvdqh,'dvdqh'
        !    print*,dvdqmclassic,'dvdqmclassic'
        !    print*,dvdqmclassic(1,1,1),dvdqmclassic(2,1,1),dvdqmclassic(3,1,1),'dvdqmclassic'

        !Gradient of Hydrogen
        DX(1)=dvdqh(1,1,1)
        DX(2)=dvdqh(2,1,1)
        DX(3)=dvdqh(3,1,1)
        !Gradient of mobile metal atoms
        DO I1 = 4,(nclassic+1)*3,3
            DX(I1) = dvdqmclassic(1,(I1-1)/3,1)
            DX(I1+1) = dvdqmclassic(2,(I1-1)/3,1)
            DX(I1+2) = dvdqmclassic(3,(I1-1)/3,1)
        ENDDO

        PRINT*, 'Gradient:'
        DO I1 = 1,15,3!1,(nclassic+1)*3,3
            PRINT*, DX(I1), DX(I1+1), DX(I1+2)
        ENDDO

    !           PRINT*, 'Energy = ', ENERGY, 'eV'
    !           PRINT*, 'Done with ENE'

    ELSE
        !PRINT*, 'Calculating energy.  One moment please ...'

        call EAM_EAMNIH_ENERGY(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
            boundary,h_included,qmclassic,qmrigid,qh,ENERGY)

    !            PRINT*, 'Energy = ', ENERGY,'eV'
    !            PRINT*, 'Done with ENE'

    ENDIF

    IF (STEST) THEN
        CALL EAM_EAMNIH_HESSIAN(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
            boundary,h_included,qmclassic,qmrigid,&
            qh,ENERGY,dvdqmclassic,dvdqmrigid,dvdqh,ddx)

        DO I1 = 1,3*natoms
            DO I2 = I1+1,3*natoms
                DDX(I2,I1) = DDX(I1,I2)
            ENDDO
        ENDDO
    ENDIF

    PRINT*, 'Energy = ', ENERGY,'eV'
    PRINT*, 'Done with ENE'


    call writespec_xyz(17,XC)

END SUBROUTINE NIHEAM7


subroutine EAM_EAMNIH_FORCES(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
    boundary,h_included,qmclassic,qmrigid,qh,v, &
    dvdqmclassic,dvdqmrigid,dvdqh)
    !
    Implicit None
    !
    !      EAM 6 Ni(100)/H. Input geometry in angstom Output energy, grad in eV
    !
    integer :: i,j,k,ilist
    integer, intent(in) :: nclassic,nrigid, nb
    double precision :: dx,dy,dz
    double precision :: vmm,vmh,fm,fh !interaction potene mm, mh, embedding ene fm, fh
    double precision, intent(in) :: boxlx,boxly,boxlz
    double precision :: onboxlx,onboxly,onboxlz
    character(len=6),intent(in) :: boundary
    double precision :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
    double precision :: ddvd2r,ddpad2r_m,ddpad2r_h
    double precision :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h,dpaqkdr_m(2)
    double precision :: dfidp,dfjdp,fjx,fjy,fjz,summ,fix,fiy,fiz,QK(2),PHQK(2)
    double precision :: dfidp1,dfidp2,fjx1,fjy1,fjz1,fjx2,fjy2,fjz2
    double precision,intent(out) :: v
    logical, intent(in) :: h_included
    !
    !
    double precision :: f_mclh(1:3,1:1,nclassic,1:nb), &
        f_mclhqk(1:3,1:1,nclassic,1:nb,1:2), &
        f_hmcl(1:3,1:1,nclassic,1:nb), &
        f_mrigh(1:3,1:1,nrigid,1:nb), &
        f_mrighqk(1:3,1:1,nrigid,1:nb,1:2), &
        f_hmrig(1:3,1:1,nrigid,1:nb), &
        f_mclcl(1:3,nclassic,nclassic,1:1), &
        f_mclrig(1:3,nrigid,nclassic,1:1), &
        f_mrigcl(1:3,nrigid,nclassic,1:1)
    !
    INTEGER :: nlist_mclh(nclassic,1:nb), &
        list_mclh(1,nclassic,1:nb),      &
        nlist_mrigh(nrigid,1:nb), &
        list_mrigh(1,nrigid,1:nb), &
        nlist_mclcl(nclassic,1:1), &
        list_mclcl(nclassic,nclassic,1:1), &
        nlist_mclrig(nclassic,1:1), &
        list_mclrig(nrigid,nclassic,1:1), &
        nlist_mrigrig(nrigid,1:1)
    !
    double precision :: dfdp_mclassic(nclassic,1:nb), &
        dfdp_mrigid(nrigid,1:nb), &
        dfdp_h(1:1,1:nb), &
        dfdpqk_h(1:1,1:nb,2)
    double precision :: ddfd2p_mclassic(nclassic,1:nb), &
        ddfd2p_mrigid(nrigid,1:nb), &
        ddfd2p_h(1:1,1:nb)
    !
    double precision,intent(in) :: qmclassic(1:3,nclassic,1:1), &
        qmrigid(1:3,nrigid,1:1)
    !
    double precision,intent(in) :: qh(1:3,1:1,1:nb)
    !
    double precision,intent(out) :: dvdqmclassic(1:3,nclassic,1:1), &
        dvdqmrigid(1:3,nrigid,1:1), &
        dvdqh(1:3,1:1,1:nb)
    !
    double precision :: eden_mrigid(nrigid,1:nb), &
        eden_mclassic(nclassic,1:nb), &
        eden_h(1:1,1:nb),&
        eden_hq1(1:1,1:nb), &
        eden_hq2(1:1,1:nb)

    !
    double precision :: smooth,dsmooth,ddsmooth

    !***************************************************************
    QK = (/ 2.D0 , 0.5D0  /) !EAM6 embedding H parameters for gradient correction


    !        nh = 1 ! this version is for 1 hydrogen atom only!
    vmm = 0.d0
    vmh = 0.d0
    v = 0.d0
    fm = 0.d0
    fh = 0.d0
    rc1 = 10.0d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
    rc1sq = rc1*rc1
    !
    !
    !
    !
    dvdqmclassic = 0.d0
    dvdqmrigid = 0.d0
    dvdqh = 0.d0
    dv = 0.d0
    dvdr = 0.d0

    !      Cl-H
    nlist_mclh = 0
    list_mclh = 0
    f_mclh = 0.d0
    f_hmcl = 0.d0
    f_mclhqk = 0.d0
    !
    !      RIG-H
    nlist_mrigh = 0
    list_mrigh = 0
    f_mrigh = 0.d0
    f_mrighqk = 0.d0
    f_hmrig = 0.d0
    !
    !      Cl-Cl
    nlist_mclcl = 0
    list_mclcl = 0
    f_mclcl = 0.d0
    !
    !     Cl-RIG
    nlist_mclrig = 0
    list_mclrig = 0
    f_mclrig = 0.d0
    f_mrigcl = 0.d0
    !
    !      RIG-RIG
    nlist_mrigrig = 0
    !       list_mrigrig = 0
    !       f_mrigrig = 0.d0
    !
    !
    dfdp_mclassic = 0.d0
    dfdp_mrigid = 0.d0
    dfdp_h = 0.d0
    !
    eden_mclassic = 0.d0
    eden_mrigid = 0.d0
    eden_h = 0.d0
    eden_hq1= 0.d0
    eden_hq2= 0.d0

    !
    !      ! Converting in eV/ANgs if input is not in Angstrom
    !       call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic,&
    !                    qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, -1)
    onboxlx = 1.d0/boxlx
    onboxly = 1.d0/boxly
    onboxlz = 1.d0/boxlz
    !
    IF(nclassic.ne.0) then
        do i = 1,nb
            do j = 1, nclassic
                !
                dx = qmclassic(1,j,1)-qh(1,1,i)
                dy = qmclassic(2,j,1)-qh(2,1,i)
                dz = qmclassic(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mclh(j,i) = nlist_mclh(j,i)+1
                    list_mclh(nlist_mclh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr
                    dvdx = dvdr*dxr
                    dvdy = dvdr*dyr
                    dvdz = dvdr*dzr
                    !
                    dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + dvdx
                    dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + dvdy
                    dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + dvdz
                    dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                    dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                    dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz
                    !
                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mclassic(j,i) = eden_mclassic(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                    eden_hq1(1,i) = eden_hq1(1,i) + pma**QK(1)
                    eden_hq2(1,i) = eden_hq2(1,i) + pma**QK(2)

                    !
                    !       !       DERIVATIVES HERE
                    !
                    dpaqkdr_m(1)=QK(1)*dpadr_m*pma**(QK(1)-1)
                    dpaqkdr_m(2)=QK(2)*dpadr_m*pma**(QK(2)-1)

                    ilist = nlist_mclh(j,i)
                    f_mclh(1,ilist,j,i) = dpadr_m*dxr
                    f_mclh(2,ilist,j,i) = dpadr_m*dyr
                    f_mclh(3,ilist,j,i) = dpadr_m*dzr
                    f_mclhqk(1,ilist,j,i,1) = dpaqkdr_m(1)*dxr
                    f_mclhqk(2,ilist,j,i,1) = dpaqkdr_m(1)*dyr
                    f_mclhqk(3,ilist,j,i,1) = dpaqkdr_m(1)*dzr
                    f_mclhqk(1,ilist,j,i,2) = dpaqkdr_m(2)*dxr
                    f_mclhqk(2,ilist,j,i,2) = dpaqkdr_m(2)*dyr
                    f_mclhqk(3,ilist,j,i,2) = dpaqkdr_m(2)*dzr
                    f_hmcl(1,ilist,j,i) = dpadr_h*dxr
                    f_hmcl(2,ilist,j,i) = dpadr_h*dyr
                    f_hmcl(3,ilist,j,i) = dpadr_h*dzr
                !
                !
                !
                end if
            end do
        end do
    END IF

    IF(nrigid.ne.0) then
        do i = 1,nb
            do j = 1, nrigid
                !
                dx = qmrigid(1,j,1)-qh(1,1,i)
                dy = qmrigid(2,j,1)-qh(2,1,i)
                dz = qmrigid(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mrigh(j,i) = nlist_mrigh(j,i)+1
                    list_mrigh(nlist_mrigh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr
                    dvdx = dvdr*dxr
                    dvdy = dvdr*dyr
                    dvdz = dvdr*dzr
                    !
                    dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                    dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                    dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz
                    !
                    !
                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mrigid(j,i) = eden_mrigid(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                    eden_hq1(1,i) = eden_hq1(1,i) + pma**QK(1)
                    eden_hq2(1,i) = eden_hq2(1,i) + pma**QK(2)

                    !
                    !       !       DERIVATIVES HERE
                    !
                    dpaqkdr_m(1)=QK(1)*dpadr_m*pma**(QK(1)-1)
                    dpaqkdr_m(2)=QK(2)*dpadr_m*pma**(QK(2)-1)

                    ilist = nlist_mrigh(j,i)
                    f_mrigh(1,ilist,j,i) = dpadr_m*dxr
                    f_mrigh(2,ilist,j,i) = dpadr_m*dyr
                    f_mrigh(3,ilist,j,i) = dpadr_m*dzr
                    f_mrighqk(1,ilist,j,i,1) = dpaqkdr_m(1)*dxr
                    f_mrighqk(2,ilist,j,i,1) = dpaqkdr_m(1)*dyr
                    f_mrighqk(3,ilist,j,i,1) = dpaqkdr_m(1)*dzr
                    f_mrighqk(1,ilist,j,i,2) = dpaqkdr_m(2)*dxr
                    f_mrighqk(2,ilist,j,i,2) = dpaqkdr_m(2)*dyr
                    f_mrighqk(3,ilist,j,i,2) = dpaqkdr_m(2)*dzr
                    f_hmrig(1,ilist,j,i) = dpadr_h*dxr
                    f_hmrig(2,ilist,j,i) = dpadr_h*dyr
                    f_hmrig(3,ilist,j,i) = dpadr_h*dzr
                !
                !
                !
                end if
            end do
        end do
    END IF

    !
    ! Classic-Classic

    do i=1,nclassic
        do j=i+1,nclassic
            dx = qmclassic(1,i,1) - qmclassic(1,j,1)
            dy = qmclassic(2,i,1) - qmclassic(2,j,1)
            dz = qmclassic(3,i,1) - qmclassic(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
            !            !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclcl(i,1) = nlist_mclcl(i,1)+1
                list_mclcl(nlist_mclcl(i,1),i,1) = j
                !
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm  + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr
                dvdx = dvdr*dxr
                dvdy = dvdr*dyr
                dvdz = dvdr*dzr
                dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz
                !
                dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) - dvdx
                dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) - dvdy
                dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) - dvdz
                !

                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do k = 1, nb
                    eden_mclassic(j,k) = eden_mclassic(j,k) + pma
                    eden_mclassic(i,k) = eden_mclassic(i,k) + pma
                end do
                !
                !        !      DERIVATIVES HERE
                !
                ilist = nlist_mclcl(i,1)
                f_mclcl(1,ilist,i,1) = dpadr_m*dxr
                f_mclcl(2,ilist,i,1) = dpadr_m*dyr
                f_mclcl(3,ilist,i,1) = dpadr_m*dzr
            !
            end if
        end do
    end do

    do i = 1,nclassic
        do k = 1,nrigid
            !
            dx = qmclassic(1,i,1) - qmrigid(1,k,1)
            dy = qmclassic(2,i,1) - qmrigid(2,k,1)
            dz = qmclassic(3,i,1) - qmrigid(3,k,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)

            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclrig(i,1) = nlist_mclrig(i,1)+1
                list_mclrig(nlist_mclrig(i,1),i,1) = k
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                !
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr
                dvdx = dvdr*dxr
                dvdy = dvdr*dyr
                dvdz = dvdr*dzr
                !
                dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz
                !
                !
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do  j = 1,nb
                    eden_mclassic(i,j) = eden_mclassic(i,j) + pma
                    eden_mrigid(k,j) = eden_mrigid(k,j) + pma
                end do
                !
                !
                !
                !     !     DERIVATIVES HERE
                !
                ilist = nlist_mclrig(i,1)
                f_mclrig(1,ilist,i,1) = dpadr_m*dxr
                f_mclrig(2,ilist,i,1) = dpadr_m*dyr
                f_mclrig(3,ilist,i,1) = dpadr_m*dzr
                f_mrigcl(1,ilist,i,1) = dpadr_m*dxr
                f_mrigcl(2,ilist,i,1) = dpadr_m*dyr
                f_mrigcl(3,ilist,i,1) = dpadr_m*dzr
            !
            !
            !
            end if
        end do
    end do


    do i = 1,nrigid
        do j = i+1,nrigid
            dx = qmrigid(1,i,1) - qmrigid(1,j,1)
            dy = qmrigid(2,i,1) - qmrigid(2,j,1)
            dz = qmrigid(3,i,1) - qmrigid(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
             !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mrigrig(i,1) = nlist_mrigrig(i,1)+1
                !
                dr = dsqrt(drsq)

                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !!
                do k = 1, nb
                    eden_mrigid(j,k) = eden_mrigid(j,k) + pma
                    eden_mrigid(i,k) = eden_mrigid(i,k) + pma
                end do
                !!      !       DERIVATIVES HERE
                !
                ilist = nlist_mrigrig(i,1)
            !                 f_mrigrig(1,ilist,i,1) = dpadr_m*dxr
            !                 f_mrigrig(2,ilist,i,1) = dpadr_m*dyr
            !                 f_mrigrig(3,ilist,i,1) = dpadr_m*dzr
            !
            end if
        end do
    end do
    !
    ! Embedding Energies
    !
    fm = 0.d0
    fh = 0.d0
    df = 0.d0
    !
    !!!

    if(nclassic.ne.0) then   !!!debug
        do k = 1,nb
            do j = 1,nclassic
                pm = eden_mclassic(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mclassic(j,k),ddfd2p_mclassic(j,k))
                fm = fm + df
            enddo
        enddo
    end if
    !
    !      do j = 1,15
    !        print*,j,eden_mrigid(j,1),'rhonickel'
    !      enddo
    if(nrigid.ne.0) then
        do k = 1,nb
            do j = 1,nrigid
                pm = eden_mrigid(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mclassic(j,k),ddfd2p_mclassic(j,k))
                fm = fm + df
            enddo
        enddo
    end if

    !       Embedding Forces for Ni-Ni
    !        Classic-Classic
    !
    if(nclassic.ne.0) then  !!! DEBUG
        do k = 1,1
            do j = 1,nclassic
                dfjdp = dfdp_mclassic(j,k)
                do ilist = 1,nlist_mclcl(j,k)
                    if(ilist.gt.nclassic) write(*,*) ' classic-classic listing error'
                    i = list_mclcl(ilist,j,k)
                    dfidp = dfdp_mclassic(i,k)
                    summ = dfidp + dfjdp
                    fjx = f_mclcl(1,ilist,j,k)
                    fjy = f_mclcl(2,ilist,j,k)
                    fjz = f_mclcl(3,ilist,j,k)
                    dvdqmclassic(1,j,k) = dvdqmclassic(1,j,k) + summ*fjx
                    dvdqmclassic(2,j,k) = dvdqmclassic(2,j,k) + summ*fjy
                    dvdqmclassic(3,j,k) = dvdqmclassic(3,j,k) + summ*fjz
                    dvdqmclassic(1,i,k) = dvdqmclassic(1,i,k) - summ*fjx
                    dvdqmclassic(2,i,k) = dvdqmclassic(2,i,k) - summ*fjy
                    dvdqmclassic(3,i,k) = dvdqmclassic(3,i,k) - summ*fjz
                enddo
            enddo
        enddo
    end if

    !        print*,dvdqmclassic(1,1,1),dvdqmclassic(2,1,1),dvdqmclassic(3,1,1),'there we are'
    !        STOP

    !
    !     Classic-RIGID

    if(nclassic.ne.0.and.nrigid.ne.0) then
        do k = 1, 1
            do j = 1, nclassic
                dfjdp = dfdp_mclassic(j,k)
                do ilist = 1, nlist_mclrig(j,k)
                    if(ilist.gt.nrigid) write(*,*) ' classic-rigid listing error'
                    i = list_mclrig(ilist,j,k)
                    dfidp = dfdp_mrigid(i,k)
                    fjx = f_mrigcl(1,ilist,j,k)
                    fjy = f_mrigcl(2,ilist,j,k)
                    fjz = f_mrigcl(3,ilist,j,k)
                    fix = f_mclrig(1,ilist,j,k)
                    fiy = f_mclrig(2,ilist,j,k)
                    fiz = f_mclrig(3,ilist,j,k)
                    dvdqmclassic(1,j,k)   =  dvdqmclassic(1,j,k) + (dfidp*fjx + dfjdp*fix)
                    dvdqmclassic(2,j,k)   =  dvdqmclassic(2,j,k) + (dfidp*fjy + dfjdp*fiy)
                    dvdqmclassic(3,j,k)   =  dvdqmclassic(3,j,k) + (dfidp*fjz + dfjdp*fiz)
                enddo
            enddo
        enddo
    end if

    if(h_included) then
        do k = 1,nb
            do j = 1,1
                ph = eden_h(j,k)
                phqk(1) = eden_hq1(j,k)
                phqk(2) = eden_hq2(j,k)

                !           print*,j,k,ph,'rhoh'
                call EAM_EMBEDING_HYDROGEN(ph,phqk,df,dfdp_h(j,k),ddfd2p_h(j,k),dfdpqk_h(j,k,1:2))
                fh = fh + df
            enddo
        enddo
    end if


    if(h_included) then
        !    !  Classic-H
        if(nclassic.ne.0) then
            do k = 1,nb
                do j = 1,nclassic
                    dfjdp = dfdp_mclassic(j,k)
                    do ilist = 1,nlist_mclh(j,k)
                        if(ilist.gt.nclassic) write(*,*) ' classic-H listing error'
                        i = list_mclh(ilist,j,k)
                        dfidp = dfdp_h(i,k)
                        dfidp1 = dfdpqk_h(i,k,1)
                        dfidp2 = dfdpqk_h(i,k,2)
                        fjx = f_mclh(1,ilist,j,k)
                        fjy = f_mclh(2,ilist,j,k)
                        fjz = f_mclh(3,ilist,j,k)
                        fjx1 = f_mclhqk(1,ilist,j,k,1)
                        fjy1 = f_mclhqk(2,ilist,j,k,1)
                        fjz1 = f_mclhqk(3,ilist,j,k,1)
                        fjx2 = f_mclhqk(1,ilist,j,k,2)
                        fjy2 = f_mclhqk(2,ilist,j,k,2)
                        fjz2 = f_mclhqk(3,ilist,j,k,2)
                        fix = f_hmcl(1,ilist,j,k)
                        fiy = f_hmcl(2,ilist,j,k)
                        fiz = f_hmcl(3,ilist,j,k)
                        dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + (dfidp*fjx + dfjdp*fix)
                        dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + (dfidp*fjy + dfjdp*fiy)
                        dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + (dfidp*fjz + dfjdp*fiz)
                        !              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
                        !              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
                        !              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
                        dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfidp1*fjx1 &
                            &  + dfidp2*fjx2 + dfjdp*fix)
                        dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfidp1*fjy1 &
                            & + dfidp2*fjy2 + dfjdp*fiy)
                        dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfidp1*fjz1 &
                            & + dfidp2*fjz2 + dfjdp*fiz)

                    enddo
                enddo
            enddo
        end if
        !

        !      !RIGID-H
        if(nrigid.ne.0) then
            do k = 1,nb
                do j = 1,nrigid
                    dfjdp = dfdp_mrigid(j,k)
                    do ilist = 1,nlist_mrigh(j,k)
                        if(ilist.gt.nrigid) write(*,*) ' rigid-H listing error'
                        i = list_mrigh(ilist,j,k)
                        dfidp = dfdp_h(i,k)
                        dfidp1 = dfdpqk_h(i,k,1)
                        dfidp2 = dfdpqk_h(i,k,2)
                        fjx = f_mrigh(1,ilist,j,k)
                        fjy = f_mrigh(2,ilist,j,k)
                        fjz = f_mrigh(3,ilist,j,k)
                        fjx1 = f_mrighqk(1,ilist,j,k,1)
                        fjy1 = f_mrighqk(2,ilist,j,k,1)
                        fjz1 = f_mrighqk(3,ilist,j,k,1)
                        fjx2 = f_mrighqk(1,ilist,j,k,2)
                        fjy2 = f_mrighqk(2,ilist,j,k,2)
                        fjz2 = f_mrighqk(3,ilist,j,k,2)
                        fix = f_hmrig(1,ilist,j,k)
                        fiy = f_hmrig(2,ilist,j,k)
                        fiz = f_hmrig(3,ilist,j,k)
                        !              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
                        !              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
                        !              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
                        dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfidp1*fjx1 + dfidp2*fjx2 + dfjdp*fix)
                        dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfidp1*fjy1 + dfidp2*fjy2 + dfjdp*fiy)
                        dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfidp1*fjz1 + dfidp2*fjz2 + dfjdp*fiz)

                    enddo
                enddo
            enddo
        end if
    !
    !
    end if  ! h_incl
    !
    !

    v = vmh + vmm + fh + fm
    !        print*,v,vmh, vmm, fh, fm,'v,vmh, vmm, fh, fm'
    !
    ! Converting forces, energy and distances in a.u.

    !        call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic, &
    !                       qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, 1)

    !       v=v/dble(nb)
    dvdqmrigid   = 0.d0
!
END subroutine EAM_EAMNIH_FORCES
!
!
!

subroutine EAM_EAMNIH_ENERGY(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
    boundary,h_included,qmclassic,qmrigid,qh,v)
    !
    Implicit None
    integer :: i,j,k,ilist
    integer, intent(in) :: nclassic,nrigid, nb
    double precision :: dx,dy,dz
    double precision :: vmm,vmh,fm,fh !interaction potene mm, mh, embedding ene fm, fh
    double precision, intent(in) :: boxlx,boxly,boxlz
    double precision :: onboxlx,onboxly,onboxlz
    character(len=6),intent(in) :: boundary
    double precision :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
    double precision :: ddvd2r,ddpad2r_m,ddpad2r_h
    double precision :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h
    double precision :: dfidp,dfjdp,fjx,fjy,fjz,summ,fix,fiy,fiz,QK(2),PHQK(2)
    double precision,intent(out) :: v
    logical, intent(in) :: h_included
    !
    !

    !
    INTEGER :: nlist_mclh(nclassic,1:nb), &
        list_mclh(1,nclassic,1:nb),      &
        nlist_mrigh(nrigid,1:nb), &
        list_mrigh(1,nrigid,1:nb), &
        nlist_mclcl(nclassic,1:1), &
        list_mclcl(nclassic,nclassic,1:1), &
        nlist_mclrig(nclassic,1:1), &
        list_mclrig(nrigid,nclassic,1:1), &
        nlist_mrigrig(nrigid,1:1)
    !
    double precision :: dfdp_mclassic(nclassic,1:nb), &
        dfdp_mrigid(nrigid,1:nb), &
        dfdp_h(1:1,1:nb), &
        dfdpqk_h(1:1,1:nb,2)
    double precision :: ddfd2p_mclassic(nclassic,1:nb), &
        ddfd2p_mrigid(nrigid,1:nb), &
        ddfd2p_h(1:1,1:nb)
    !
    double precision,intent(in) :: qmclassic(1:3,nclassic,1:1), &
        qmrigid(1:3,nrigid,1:1)
    !
    double precision,intent(in) :: qh(1:3,1:1,1:nb)
    !
    double precision :: eden_mrigid(nrigid,1:nb), &
        eden_mclassic(nclassic,1:nb), &
        eden_h(1:1,1:nb), &
        eden_hq1(1:1,1:nb), &
        eden_hq2(1:1,1:nb)
    !
    double precision :: smooth,dsmooth,ddsmooth

    !***************************************************************
    QK = (/ 2.D0 , 0.5D0  /) !EAM6 embedding H parameters

    !        nh = 1 ! this version is for 1 hydrogen atom only!
    vmm = 0.d0
    vmh = 0.d0
    v = 0.d0
    fm = 0.d0
    fh = 0.d0
    rc1 = 10.0d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
    rc1sq = rc1*rc1
    !
    !
    !
    !
    dv = 0.d0
    dvdr = 0.d0

    !      Cl-H
    nlist_mclh = 0
    list_mclh = 0
    !
    !      RIG-H
    nlist_mrigh = 0
    list_mrigh = 0
    !
    !      Cl-Cl
    nlist_mclcl = 0
    list_mclcl = 0
    !
    !     Cl-RIG
    nlist_mclrig = 0
    list_mclrig = 0
    !
    !      RIG-RIG
    nlist_mrigrig = 0
    !       list_mrigrig = 0

    !
    !
    dfdp_mclassic = 0.d0
    dfdp_mrigid = 0.d0
    dfdp_h = 0.d0
    dfdpqk_h = 0.d0
    !
    eden_mclassic = 0.d0
    eden_mrigid = 0.d0
    eden_h = 0.d0
    eden_hq1 = 0.d0
    eden_hq2 = 0.d0
    !
    !      ! Converting in eV/ANgs if input is not in Angstrom
    !       call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic,&
    !                    qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, -1)
    onboxlx = 1.d0/boxlx
    onboxly = 1.d0/boxly
    onboxlz = 1.d0/boxlz
    !
    IF(nclassic.ne.0) then
        do i = 1,nb
            do j = 1, nclassic
                !
                dx = qmclassic(1,j,1)-qh(1,1,i)
                dy = qmclassic(2,j,1)-qh(2,1,i)
                dz = qmclassic(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mclh(j,i) = nlist_mclh(j,i)+1
                    list_mclh(nlist_mclh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr

                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mclassic(j,i) = eden_mclassic(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                    eden_hq1(1,i) = eden_hq1(1,i) + pma**QK(1)
                    eden_hq2(1,i) = eden_hq2(1,i) + pma**QK(2)
                !

                end if
            end do
        end do
    END IF

    IF(nrigid.ne.0) then
        do i = 1,nb
            do j = 1, nrigid
                !
                dx = qmrigid(1,j,1)-qh(1,1,i)
                dy = qmrigid(2,j,1)-qh(2,1,i)
                dz = qmrigid(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mrigh(j,i) = nlist_mrigh(j,i)+1
                    list_mrigh(nlist_mrigh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr

                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mrigid(j,i) = eden_mrigid(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                !

                end if
            end do
        end do
    END IF

    !


    ! Classic-Classic

    do i=1,nclassic
        do j=i+1,nclassic
            dx = qmclassic(1,i,1) - qmclassic(1,j,1)
            dy = qmclassic(2,i,1) - qmclassic(2,j,1)
            dz = qmclassic(3,i,1) - qmclassic(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
            !            !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclcl(i,1) = nlist_mclcl(i,1)+1
                list_mclcl(nlist_mclcl(i,1),i,1) = j
                !
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm  + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr

                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do k = 1, nb
                    eden_mclassic(j,k) = eden_mclassic(j,k) + pma
                    eden_mclassic(i,k) = eden_mclassic(i,k) + pma
                end do
            !

            end if
        end do
    end do

    do i = 1,nclassic
        do k = 1,nrigid
            !
            dx = qmclassic(1,i,1) - qmrigid(1,k,1)
            dy = qmclassic(2,i,1) - qmrigid(2,k,1)
            dz = qmclassic(3,i,1) - qmrigid(3,k,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)

            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclrig(i,1) = nlist_mclrig(i,1)+1
                list_mclrig(nlist_mclrig(i,1),i,1) = k
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                !
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr

                !
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do  j = 1,nb
                    eden_mclassic(i,j) = eden_mclassic(i,j) + pma
                    eden_mrigid(k,j) = eden_mrigid(k,j) + pma
                end do
            !
            end if
        end do
    end do


    do i = 1,nrigid
        do j = i+1,nrigid
            dx = qmrigid(1,i,1) - qmrigid(1,j,1)
            dy = qmrigid(2,i,1) - qmrigid(2,j,1)
            dz = qmrigid(3,i,1) - qmrigid(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
             !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mrigrig(i,1) = nlist_mrigrig(i,1)+1
                !
                dr = dsqrt(drsq)

                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !!
                do k = 1, nb
                    eden_mrigid(j,k) = eden_mrigid(j,k) + pma
                    eden_mrigid(i,k) = eden_mrigid(i,k) + pma
                end do
            !
            end if
        end do
    end do
    !
    ! Embedding Energies
    !
    fm = 0.d0
    fh = 0.d0
    df = 0.d0
    !
    !!!

    if(nclassic.ne.0) then   !!!debug
        do k = 1,nb
            do j = 1,nclassic
                pm = eden_mclassic(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mclassic(j,k),ddfd2p_mclassic(j,k))

                fm = fm + df
            enddo
        enddo
    end if
    !
    !      do j = 1,15
    !        print*,j,eden_mrigid(j,1),'rhonickel'
    !      enddo
    if(nrigid.ne.0) then
        do k = 1,nb
            do j = 1,nrigid
                pm = eden_mrigid(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mclassic(j,k),ddfd2p_mclassic(j,k))
                fm = fm + df
            enddo
        enddo
    end if


    if(h_included) then
        do k = 1,nb
            do j = 1,1
                ph = eden_h(j,k)
                phqk(1) = eden_hq1(j,k)
                phqk(2) = eden_hq2(j,k)
                call EAM_EMBEDING_HYDROGEN(ph,phqk,df,dfdp_h(j,k),ddfd2p_h(j,k),dfdpqk_h(j,k,1:2))
                fh = fh + df
            enddo
        enddo
    end if


    !
    !
    v = vmh + vmm + fh + fm
    print*,v,vmh, vmm, fh, fm,'v,vmh, vmm, fh, fm'
!
      ! Converting forces, energy and distances in a.u.

!        call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic, &
!                       qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, 1)

!       v=v/dble(nb)

END subroutine EAM_EAMNIH_ENERGY


subroutine EAM_EAMNIH_HESSIAN(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
    boundary,h_included,qmclassic,qmrigid,qh,v, &
    dvdqmclassic,dvdqmrigid,dvdqh,ddx)
    !
    Implicit None
    !
    !      EAM 6 Ni(100)/H. Input geometry in angstom Output energy, grad in eV
    !
    integer :: i,j,k,ilist
    integer, intent(in) :: nclassic,nrigid, nb
    double precision :: dx,dy,dz
    double precision :: vmm,vmh,fm,fh !interaction potene mm, mh, embedding ene fm, fh
    double precision, intent(in) :: boxlx,boxly,boxlz
    double precision :: onboxlx,onboxly,onboxlz
    character(len=6),intent(in) :: boundary
    double precision :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
    double precision :: ddvd2r,ddpad2r_m,ddpad2r_h
    double precision :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h,dpaqkdr_m(2)
    double precision :: dfidp,dfjdp,fjx,fjy,fjz,summ,fix,fiy,fiz,QK(2),PHQK(2)
    double precision :: dfidp1,dfidp2,fjx1,fjy1,fjz1,fjx2,fjy2,fjz2
    double precision,intent(out) :: v
    logical, intent(in) :: h_included
    !
    !
    double precision :: f_mclh(1:3,1:1,nclassic,1:nb), &
        f_mclhqk(1:3,1:1,nclassic,1:nb,1:2), &
        f_hmcl(1:3,1:1,nclassic,1:nb), &
        f_mrigh(1:3,1:1,nrigid,1:nb), &
        f_mrighqk(1:3,1:1,nrigid,1:nb,1:2), &
        f_hmrig(1:3,1:1,nrigid,1:nb), &
        f_mclcl(1:3,nclassic,nclassic,1:1), &
        f_mclrig(1:3,nrigid,nclassic,1:1), &
        f_mrigcl(1:3,nrigid,nclassic,1:1)

    double precision :: ff2_mclh(1:3,1:3,1:1,nclassic,1:nb), &
        ff2_hmcl(1:3,1:3,1:1,nclassic,1:nb), &
        ff2_mrigh(1:3,1:3,1:1,nrigid,1:nb), &
        ff2_hmrig(1:3,1:3,1:1,nrigid,1:nb), &
        ff2_mclcl(1:3,1:3,nclassic,nclassic,1:1), &
        ff2_mclrig(1:3,1:3,nrigid,nclassic,1:1), &
        ff2_mrigcl(1:3,1:3,nrigid,nclassic,1:1)

    !
    INTEGER :: nlist_mclh(nclassic,1:nb), &
        list_mclh(1,nclassic,1:nb),      &
        nlist_mrigh(nrigid,1:nb), &
        list_mrigh(1,nrigid,1:nb), &
        nlist_mclcl(nclassic,1:1), &
        list_mclcl(nclassic,nclassic,1:1), &
        nlist_mclrig(nclassic,1:1), &
        list_mclrig(nrigid,nclassic,1:1), &
        nlist_mrigrig(nrigid,1:1)
    !
    double precision :: dfdp_mclassic(nclassic,1:nb), &
        dfdp_mrigid(nrigid,1:nb), &
        dfdp_h(1:1,1:nb), &
        dfdpqk_h(1:1,1:nb,2)
    double precision :: ddfd2p_mclassic(nclassic,1:nb), &
        ddfd2p_mrigid(nrigid,1:nb), &
        ddfd2p_h(1:1,1:nb)
    !
    double precision,intent(in) :: qmclassic(1:3,nclassic,1:1), &
        qmrigid(1:3,nrigid,1:1)
    !
    double precision,intent(in) :: qh(1:3,1:1,1:nb)
    !
    double precision,intent(out) :: dvdqmclassic(1:3,nclassic,1:1), &
        dvdqmrigid(1:3,nrigid,1:1), &
        dvdqh(1:3,1:1,1:nb)
    double precision,intent(out) :: ddx(3*(nclassic+nrigid+1),3*(nclassic+nrigid+1))
    !
    double precision :: eden_mrigid(nrigid,1:nb), &
        eden_mclassic(nclassic,1:nb), &
        eden_h(1:1,1:nb),&
        eden_hq1(1:1,1:nb), &
        eden_hq2(1:1,1:nb)

    !
    double precision :: smooth,dsmooth,ddsmooth,ddfid2p,ddfjd2p,summ2

    !***************************************************************
    QK = (/ 2.D0 , 0.5D0  /) !EAM6 embedding H parameters for gradient correction


    !        nh = 1 ! this version is for 1 hydrogen atom only!
    vmm = 0.d0
    vmh = 0.d0
    v = 0.d0
    fm = 0.d0
    fh = 0.d0
    rc1 = 10.D0 !10.0d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
    rc1sq = rc1*rc1
    !
    !
    !
    !
    dvdqmclassic = 0.d0
    dvdqmrigid = 0.d0
    dvdqh = 0.d0
    dv = 0.d0
    dvdr = 0.d0

    !      Cl-H
    nlist_mclh = 0
    list_mclh = 0
    f_mclh = 0.d0
    ff2_mclh = 0.d0
    ff2_hmcl = 0.d0
    f_hmcl = 0.d0
    f_mclhqk = 0.d0
    !
    !      RIG-H
    nlist_mrigh = 0
    list_mrigh = 0
    f_mrigh = 0.d0
    f_mrighqk = 0.d0
    f_hmrig = 0.d0
    ff2_mrigh = 0.d0
    ff2_hmrig = 0.d0
    !
    !      Cl-Cl
    nlist_mclcl = 0
    list_mclcl = 0
    f_mclcl = 0.d0
    ff2_mclcl = 0.d0
    !
    !     Cl-RIG
    nlist_mclrig = 0
    list_mclrig = 0
    f_mclrig = 0.d0
    f_mrigcl = 0.d0
    ff2_mclrig = 0.d0
    ff2_mrigcl = 0.d0
    !
    !      RIG-RIG
    nlist_mrigrig = 0
    !       list_mrigrig = 0
    !       f_mrigrig = 0.d0
    !
    !
    dfdp_mclassic = 0.d0
    dfdp_mrigid = 0.d0
    dfdp_h = 0.d0
    dfdpqk_h  = 0.d0

    ddfd2p_mclassic = 0.d0
    ddfd2p_mrigid = 0.d0
    ddfd2p_h = 0.d0
    !
    eden_mclassic = 0.d0
    eden_mrigid = 0.d0
    eden_h = 0.d0
    eden_hq1= 0.d0
    eden_hq2= 0.d0

    !
    !      ! Converting in eV/ANgs if input is not in Angstrom
    !       call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic,&
    !                    qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, -1)
    onboxlx = 1.d0/boxlx
    onboxly = 1.d0/boxly
    onboxlz = 1.d0/boxlz
    !
    IF(nclassic.ne.0) then
        do i = 1,nb
            do j = 1, nclassic
                dx = qmclassic(1,j,1)-qh(1,1,i)
                dy = qmclassic(2,j,1)-qh(2,1,i)
                dz = qmclassic(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mclh(j,i) = nlist_mclh(j,i)+1
                    list_mclh(nlist_mclh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr
                    dvdx = dvdr*dxr
                    dvdy = dvdr*dyr
                    dvdz = dvdr*dzr
                    !
                    dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + dvdx
                    dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + dvdy
                    dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + dvdz
                    dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                    dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                    dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz

                    !! H-H second deriv.
                    ddx(1,1)=ddx(1,1)+(1.D0/dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))

                    ddx(2,1)=ddx(2,1)+dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(2,2)=(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))

                    ddx(3,1)=ddx(3,1)+dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3,2)=ddx(3,2)+dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3,3)=ddx(3,3)+(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))

                    !!! H-Ni second deriv.
                    !!dxhdxni
                    ddx(3+3*j-2,1)=ddx(3+3*j-2,1)-(1.D0/dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))
                    ddx(3+3*j-2,2)=ddx(3+3*j-2,2)-dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3+3*j-2,3)=ddx(3+3*j-2,3)-dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                    !!dxhdyni
                    ddx(3+3*j-1,1)=ddx(3+3*j-1,1)-dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3+3*j-1,2)=ddx(3+3*j-1,2)-(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))
                    ddx(3+3*j-1,3)=ddx(3+3*j-1,3)-dy*dz/(dr**2)*(ddvd2r-dvdr/dr)
                    !!dxhdzni
                    ddx(3+3*j,1)=ddx(3+3*j,1)-dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3+3*j,2)=ddx(3+3*j,2)-dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3+3*j,3)=ddx(3+3*j,3)-(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))

                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mclassic(j,i) = eden_mclassic(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                    eden_hq1(1,i) = eden_hq1(1,i) + pma**QK(1)
                    eden_hq2(1,i) = eden_hq2(1,i) + pma**QK(2)

                    !
                    !       !       DERIVATIVES HERE
                    !
                    dpaqkdr_m(1)=QK(1)*dpadr_m*pma**(QK(1)-1)
                    dpaqkdr_m(2)=QK(2)*dpadr_m*pma**(QK(2)-1)

                    ilist = nlist_mclh(j,i)
                    f_mclh(1,ilist,j,i) = dpadr_m*dxr
                    f_mclh(2,ilist,j,i) = dpadr_m*dyr
                    f_mclh(3,ilist,j,i) = dpadr_m*dzr
                    f_mclhqk(1,ilist,j,i,1) = dpaqkdr_m(1)*dxr
                    f_mclhqk(2,ilist,j,i,1) = dpaqkdr_m(1)*dyr
                    f_mclhqk(3,ilist,j,i,1) = dpaqkdr_m(1)*dzr
                    f_mclhqk(1,ilist,j,i,2) = dpaqkdr_m(2)*dxr
                    f_mclhqk(2,ilist,j,i,2) = dpaqkdr_m(2)*dyr
                    f_mclhqk(3,ilist,j,i,2) = dpaqkdr_m(2)*dzr
                    f_hmcl(1,ilist,j,i) = dpadr_h*dxr
                    f_hmcl(2,ilist,j,i) = dpadr_h*dyr
                    f_hmcl(3,ilist,j,i) = dpadr_h*dzr


                    ff2_mclh(1,1,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dx**2+dpadr_m*(dr-dx**2/dr))
                    ff2_mclh(1,2,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(1,3,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(2,1,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(2,2,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dy**2+dpadr_m*(dr-dy**2/dr))
                    ff2_mclh(2,3,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(3,1,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(3,2,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mclh(3,3,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dz**2+dpadr_m*(dr-dz**2/dr))


                    ff2_hmcl(1,1,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dx**2+dpadr_h*(dr-dx**2/dr))
                    ff2_hmcl(1,2,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(1,3,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(2,1,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(2,2,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dy**2+dpadr_h*(dr-dy**2/dr))
                    ff2_hmcl(2,3,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(3,1,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(3,2,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmcl(3,3,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dz**2+dpadr_h*(dr-dz**2/dr))
                !
                !
                !
                end if
            end do
        end do
    END IF

    IF(nrigid.ne.0) then
        do i = 1,nb
            do j = 1, nrigid
                !
                dx = qmrigid(1,j,1)-qh(1,1,i)
                dy = qmrigid(2,j,1)-qh(2,1,i)
                dz = qmrigid(3,j,1)-qh(3,1,i)
                if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                    dx = dx-boxlx*nint(onboxlx*dx)
                    dy = dy-boxly*nint(onboxly*dy)
                !            !dz = dz-boxlz*nint(onboxlz*dz)
                end if
                drsq = dx*dx + dy*dy + dz*dz
                if (drsq.lt.rc1sq) then
                    ! Add to neighbour list
                    !
                    nlist_mrigh(j,i) = nlist_mrigh(j,i)+1
                    list_mrigh(nlist_mrigh(j,i),j,i) = 1
                    !
                    dr = dsqrt(drsq)
                    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                    call EAM_PAIR_NIH(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                    vmh = vmh + dv
                    dxr = dx/dr
                    dyr = dy/dr
                    dzr = dz/dr
                    dvdx = dvdr*dxr
                    dvdy = dvdr*dyr
                    dvdz = dvdr*dzr
                    !
                    dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                    dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                    dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz

                    !! H-H second deriv.
                    ddx(1,1)=ddx(1,1)+(1.D0/dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))

                    ddx(2,1)=ddx(2,1)+dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(2,2)=ddx(2,2)+(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))

                    ddx(3,1)=ddx(3,1)+dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3,2)=ddx(3,2)+dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                    ddx(3,3)=ddx(3,3)+(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))
                    !
                    !
                    call EAM_RHO_H(dr,pha,dpadr_h,ddpad2r_h,Smooth,DSmooth,DDSmooth)
                    call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                    !
                    eden_mrigid(j,i) = eden_mrigid(j,i) + pha
                    eden_h(1,i) = eden_h(1,i) + pma
                    eden_hq1(1,i) = eden_hq1(1,i) + pma**QK(1)
                    eden_hq2(1,i) = eden_hq2(1,i) + pma**QK(2)

                    !
                    !       !       DERIVATIVES HERE
                    !
                    dpaqkdr_m(1)=QK(1)*dpadr_m*pma**(QK(1)-1)
                    dpaqkdr_m(2)=QK(2)*dpadr_m*pma**(QK(2)-1)

                    ilist = nlist_mrigh(j,i)
                    f_mrigh(1,ilist,j,i) = dpadr_m*dxr
                    f_mrigh(2,ilist,j,i) = dpadr_m*dyr
                    f_mrigh(3,ilist,j,i) = dpadr_m*dzr
                    f_mrighqk(1,ilist,j,i,1) = dpaqkdr_m(1)*dxr
                    f_mrighqk(2,ilist,j,i,1) = dpaqkdr_m(1)*dyr
                    f_mrighqk(3,ilist,j,i,1) = dpaqkdr_m(1)*dzr
                    f_mrighqk(1,ilist,j,i,2) = dpaqkdr_m(2)*dxr
                    f_mrighqk(2,ilist,j,i,2) = dpaqkdr_m(2)*dyr
                    f_mrighqk(3,ilist,j,i,2) = dpaqkdr_m(2)*dzr
                    f_hmrig(1,ilist,j,i) = dpadr_h*dxr
                    f_hmrig(2,ilist,j,i) = dpadr_h*dyr
                    f_hmrig(3,ilist,j,i) = dpadr_h*dzr



                    ff2_mrigh(1,1,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dx**2+dpadr_m*(dr-dx**2/dr))
                    ff2_mrigh(1,2,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(1,3,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(2,1,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(2,2,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dy**2+dpadr_m*(dr-dy**2/dr))
                    ff2_mrigh(2,3,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(3,1,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(3,2,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                    ff2_mrigh(3,3,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_m*dz**2+dpadr_m*(dr-dz**2/dr))


                    ff2_hmrig(1,1,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dx**2+dpadr_h*(dr-dx**2/dr))
                    ff2_hmrig(1,2,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(1,3,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(2,1,ilist,j,i) = dx*dy/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(2,2,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dy**2+dpadr_h*(dr-dy**2/dr))
                    ff2_hmrig(2,3,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(3,1,ilist,j,i) = dx*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(3,2,ilist,j,i) = dy*dz/(dr**2)*(ddpad2r_h-dpadr_h/dr)
                    ff2_hmrig(3,3,ilist,j,i) = (1.D0/dr)**2*(ddpad2r_h*dz**2+dpadr_h*(dr-dz**2/dr))

                !
                !
                !
                end if
            end do
        end do
    END IF

    !
    ! Classic-Classic

    do i=1,nclassic
        do j=i+1,nclassic
            dx = qmclassic(1,i,1) - qmclassic(1,j,1)
            dy = qmclassic(2,i,1) - qmclassic(2,j,1)
            dz = qmclassic(3,i,1) - qmclassic(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
            !            !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclcl(i,1) = nlist_mclcl(i,1)+1
                list_mclcl(nlist_mclcl(i,1),i,1) = j
                !
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm  + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr
                dvdx = dvdr*dxr
                dvdy = dvdr*dyr
                dvdz = dvdr*dzr
                dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz
                !
                dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) - dvdx
                dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) - dvdy
                dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) - dvdz
                !

                !! Ni-Ni on diagonal second deriv.
                ddx(3+3*i-2,3+3*i-2)=ddx(3+3*i-2,3+3*i-2)+1.D0/(dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))

                ddx(3+3*i-1,3+3*i-2)=ddx(3+3*i-1,3+3*i-2)+dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i-1,3+3*i-1)=ddx(3+3*i-1,3+3*i-1)+(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))

                ddx(3+3*i,3+3*i-2)=ddx(3+3*i,3+3*i-2)+dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i,3+3*i-1)=ddx(3+3*i,3+3*i-1)+dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i,3+3*i)=ddx(3+3*i,3+3*i)+(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))



                !!! Ni-Ni second deriv. off diagonal
                ddx(3+3*j-2,3+3*i-2)=ddx(3+3*j-2,3+3*i-2)-(1.D0/dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))
                ddx(3+3*j-2,3+3*i-1)=ddx(3+3*j-2,3+3*i-1)-dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*j-2,3+3*i)=ddx(3+3*j-2,3+3*i)-dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                !!
                ddx(3+3*j-1,3+3*i-2)=ddx(3+3*j-1,3+3*i-2)-dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*j-1,3+3*i-1)=ddx(3+3*j-1,3+3*i-1)-(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))
                ddx(3+3*j-1,3+3*i)=ddx(3+3*j-1,3+3*i)-dy*dz/(dr**2)*(ddvd2r-dvdr/dr)
                !!
                ddx(3+3*j,3+3*i-2)=ddx(3+3*j,3+3*i-2)-dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*j,3+3*i-1)=ddx(3+3*j,3+3*i-1)-dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*j,3+3*i)=ddx(3+3*j,3+3*i)-(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))



                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do k = 1, nb
                    eden_mclassic(j,k) = eden_mclassic(j,k) + pma
                    eden_mclassic(i,k) = eden_mclassic(i,k) + pma
                end do
                !
                !        !      DERIVATIVES HERE
                !
                ilist = nlist_mclcl(i,1)
                f_mclcl(1,ilist,i,1) = dpadr_m*dxr
                f_mclcl(2,ilist,i,1) = dpadr_m*dyr
                f_mclcl(3,ilist,i,1) = dpadr_m*dzr

                ff2_mclcl(1,1,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dx**2+dpadr_m*(dr-dx**2/dr))
                ff2_mclcl(1,2,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclcl(1,3,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclcl(2,1,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclcl(2,2,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dy**2+dpadr_m*(dr-dy**2/dr))
                ff2_mclcl(2,3,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclcl(3,2,ilist,i,1) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclcl(3,3,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dz**2+dpadr_m*(dr-dz**2/dr))


            !
            end if
        end do
    end do

    do i = 1,nclassic
        do k = 1,nrigid
            !
            dx = qmclassic(1,i,1) - qmrigid(1,k,1)
            dy = qmclassic(2,i,1) - qmrigid(2,k,1)
            dz = qmclassic(3,i,1) - qmrigid(3,k,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)

            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mclrig(i,1) = nlist_mclrig(i,1)+1
                list_mclrig(nlist_mclrig(i,1),i,1) = k
                dr = dsqrt(drsq)
                !
                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                !
                call EAM_PAIR_NINI(dr,dv,dvdr,ddvd2r,Smooth,DSmooth,DDSmooth)
                vmm = vmm + dv
                dxr = dx/dr
                dyr = dy/dr
                dzr = dz/dr
                dvdx = dvdr*dxr
                dvdy = dvdr*dyr
                dvdz = dvdr*dzr
                !
                dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz


                !! Nicl-Nirig on diagonal second deriv.
                ddx(3+3*i-2,3+3*i-2)=ddx(3+3*i-2,3+3*i-2)+1.D0/(dr)**2*(ddvd2r*dx**2+dvdr*(dr-dx**2/dr))

                ddx(3+3*i-1,3+3*i-2)=ddx(3+3*i-1,3+3*i-2)+dx*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i-1,3+3*i-1)=ddx(3+3*i-1,3+3*i-1)+(1.D0/dr)**2*(ddvd2r*dy**2+dvdr*(dr-dy**2/dr))

                ddx(3+3*i,3+3*i-2)=ddx(3+3*i,3+3*i-2)+dx*dz/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i,3+3*i-1)=ddx(3+3*i,3+3*i-1)+dz*dy/(dr**2)*(ddvd2r-dvdr/dr)
                ddx(3+3*i,3+3*i)=ddx(3+3*i,3+3*i)+(1.D0/dr)**2*(ddvd2r*dz**2+dvdr*(dr-dz**2/dr))
                !
                !
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !
                do  j = 1,nb
                    eden_mclassic(i,j) = eden_mclassic(i,j) + pma
                    eden_mrigid(k,j) = eden_mrigid(k,j) + pma
                end do
                !
                !
                !
                !     !     DERIVATIVES HERE
                !
                ilist = nlist_mclrig(i,1)
                f_mclrig(1,ilist,i,1) = dpadr_m*dxr
                f_mclrig(2,ilist,i,1) = dpadr_m*dyr
                f_mclrig(3,ilist,i,1) = dpadr_m*dzr
                f_mrigcl(1,ilist,i,1) = dpadr_m*dxr
                f_mrigcl(2,ilist,i,1) = dpadr_m*dyr
                f_mrigcl(3,ilist,i,1) = dpadr_m*dzr

                ff2_mclrig(1,1,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dx**2+dpadr_m*(dr-dx**2/dr))
                ff2_mclrig(1,2,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(1,3,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(2,1,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(2,2,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dy**2+dpadr_m*(dr-dy**2/dr))
                ff2_mclrig(2,3,ilist,i,1) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(3,1,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(3,2,ilist,i,1) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mclrig(3,3,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dz**2+dpadr_m*(dr-dz**2/dr))

                ff2_mrigcl(1,1,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dx**2+dpadr_m*(dr-dx**2/dr))
                ff2_mrigcl(1,2,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(1,3,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(2,1,ilist,i,1) = dx*dy/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(2,2,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dy**2+dpadr_m*(dr-dy**2/dr))
                ff2_mrigcl(2,3,ilist,i,1) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(3,1,ilist,i,1) = dx*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(3,2,ilist,i,1) = dy*dz/(dr**2)*(ddpad2r_m-dpadr_m/dr)
                ff2_mrigcl(3,3,ilist,i,1) = (1.D0/dr)**2*(ddpad2r_m*dz**2+dpadr_m*(dr-dz**2/dr))
            !
            !
            !
            end if
        end do
    end do


    do i = 1,nrigid
        do j = i+1,nrigid
            dx = qmrigid(1,i,1) - qmrigid(1,j,1)
            dy = qmrigid(2,i,1) - qmrigid(2,j,1)
            dz = qmrigid(3,i,1) - qmrigid(3,j,1)
            if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
                dx = dx-boxlx*nint(onboxlx*dx)
                dy = dy-boxly*nint(onboxly*dy)
             !dz = dz-boxlz*nint(onboxlz*dz)
            end if
            drsq = dx*dx + dy*dy + dz*dz
            if (drsq.lt.rc1sq) then
                ! Add to neighbour list
                !
                nlist_mrigrig(i,1) = nlist_mrigrig(i,1)+1
                !
                dr = dsqrt(drsq)

                CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
                call EAM_RHO_NI(dr,pma,dpadr_m,ddpad2r_m,Smooth,DSmooth,DDSmooth)
                !!
                do k = 1, nb
                    eden_mrigid(j,k) = eden_mrigid(j,k) + pma
                    eden_mrigid(i,k) = eden_mrigid(i,k) + pma
                end do
                !!      !       DERIVATIVES HERE
                !
                ilist = nlist_mrigrig(i,1)
            !                 f_mrigrig(1,ilist,i,1) = dpadr_m*dxr
            !                 f_mrigrig(2,ilist,i,1) = dpadr_m*dyr
            !                 f_mrigrig(3,ilist,i,1) = dpadr_m*dzr
            !
            end if
        end do
    end do
    !
    ! Embedding Energies
    !
    fm = 0.d0
    fh = 0.d0
    df = 0.d0
    !
    !!!

    if(nclassic.ne.0) then   !!!debug
        do k = 1,nb
            do j = 1,nclassic
                pm = eden_mclassic(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mclassic(j,k),ddfd2p_mclassic(j,k))
                fm = fm + df
            enddo
        enddo
    end if
    !
    !      do j = 1,15
    !        print*,j,eden_mrigid(j,1),'rhonickel'
    !      enddo
    if(nrigid.ne.0) then
        do k = 1,nb
            do j = 1,nrigid
                pm = eden_mrigid(j,k)
                call EAM_EMBEDING_NICKEL(pm,df,dfdp_mrigid(j,k),ddfd2p_mrigid(j,k))
                fm = fm + df
            enddo
        enddo
    end if

    !       Embedding Forces for Ni-Ni
    !        Classic-Classic
    !
    if(nclassic.ne.0) then  !!! DEBUG
        do k = 1,1
            do j = 1,nclassic
                dfjdp = dfdp_mclassic(j,k)
                ddfjd2p = ddfd2p_mclassic(j,k)
                do ilist = 1,nlist_mclcl(j,k)
                    if(ilist.gt.nclassic) write(*,*) ' classic-classic listing error'
                    i = list_mclcl(ilist,j,k)
                    dfidp = dfdp_mclassic(i,k)
                    ddfid2p = ddfd2p_mclassic(i,k)
                    summ = dfidp + dfjdp
                    summ2= ddfjd2p + ddfid2p
                    fjx = f_mclcl(1,ilist,j,k)
                    fjy = f_mclcl(2,ilist,j,k)
                    fjz = f_mclcl(3,ilist,j,k)
                    dvdqmclassic(1,j,k) = dvdqmclassic(1,j,k) + summ*fjx
                    dvdqmclassic(2,j,k) = dvdqmclassic(2,j,k) + summ*fjy
                    dvdqmclassic(3,j,k) = dvdqmclassic(3,j,k) + summ*fjz
                    dvdqmclassic(1,i,k) = dvdqmclassic(1,i,k) - summ*fjx
                    dvdqmclassic(2,i,k) = dvdqmclassic(2,i,k) - summ*fjy
                    dvdqmclassic(3,i,k) = dvdqmclassic(3,i,k) - summ*fjz

                    !! Ni-Ni on diagonal second deriv.
                    ddx(3+3*j-2,3+3*j-2)=ddx(3+3*j-2,3+3*j-2)+summ2*fjx**2+summ*ff2_mclcl(1,1,ilist,j,k)

                    ddx(3+3*j-1,3+3*j-2)=ddx(3+3*j-1,3+3*j-2)+summ2*fjy*fjx+summ*ff2_mclcl(2,1,ilist,j,k)
                    ddx(3+3*j-1,3+3*j-1)=ddx(3+3*j-1,3+3*j-1)+summ2*fjy**2+summ*ff2_mclcl(2,2,ilist,j,k)

                    ddx(3+3*j,3+3*j-2)=ddx(3+3*j,3+3*j-2)+summ2*fjz*fjx+summ*ff2_mclcl(3,1,ilist,j,k)
                    ddx(3+3*j,3+3*j-1)=ddx(3+3*j,3+3*j-1)+summ2*fjz*fjy+summ*ff2_mclcl(3,2,ilist,j,k)
                    ddx(3+3*j,3+3*j)=ddx(3+3*j,3+3*j)+summ2*fjz**2+summ*ff2_mclcl(3,3,ilist,j,k)

                    !!! Ni-Ni second deriv. off diagonal, density no selfinteraction
                    ddx(3+3*i-2,3+3*j-2)=ddx(3+3*i-2,3+3*j-2)-(summ2*fjx**2+summ*ff2_mclcl(1,1,ilist,j,k))
                    ddx(3+3*i-2,3+3*j-1)=ddx(3+3*i-2,3+3*j-1)-(summ2*fjy*fjx+summ*ff2_mclcl(1,2,ilist,j,k))
                    ddx(3+3*i-2,3+3*j)=ddx(3+3*i-2,3+3*j)-(summ2*fjz*fjx+summ*ff2_mclcl(1,3,ilist,j,k))
                    !!
                    ddx(3+3*i-1,3+3*j-2)=ddx(3+3*i-1,3+3*j-2)-(summ2*fjy*fjx+summ*ff2_mclcl(2,1,ilist,j,k))
                    ddx(3+3*i-1,3+3*j-1)=ddx(3+3*i-1,3+3*j-1)-(summ2*fjy**2+summ*ff2_mclcl(2,2,ilist,j,k))
                    ddx(3+3*i-1,3+3*j)=ddx(3+3*i-1,3+3*j)-(summ2*fjy*fjz+summ*ff2_mclcl(2,3,ilist,j,k))
                    !!
                    ddx(3+3*i,3+3*j-2)=ddx(3+3*i,3+3*j-2)-(summ2*fjz*fjx+summ*ff2_mclcl(3,1,ilist,j,k))
                    ddx(3+3*i,3+3*j-1)=ddx(3+3*i,3+3*j-1)-(summ2*fjz*fjy+summ*ff2_mclcl(3,2,ilist,j,k))
                    ddx(3+3*i,3+3*j)=ddx(3+3*i,3+3*j)-(summ2*fjz**2+summ*ff2_mclcl(3,3,ilist,j,k))

                enddo
            enddo
        enddo
    end if

    !        print*,dvdqmclassic(1,1,1),dvdqmclassic(2,1,1),dvdqmclassic(3,1,1),'there we are'
    !        STOP

    !
    !     Classic-RIGID

    if(nclassic.ne.0.and.nrigid.ne.0) then
        do k = 1, 1
            do j = 1, nclassic
                dfjdp = dfdp_mclassic(j,k)
                ddfjd2p = ddfd2p_mclassic(j,k)
                do ilist = 1, nlist_mclrig(j,k)
                    if(ilist.gt.nrigid) write(*,*) ' classic-rigid listing error'
                    i = list_mclrig(ilist,j,k)
                    dfidp = dfdp_mrigid(i,k)
                    ddfid2p = ddfd2p_mrigid(i,k)
                    fjx = f_mrigcl(1,ilist,j,k)
                    fjy = f_mrigcl(2,ilist,j,k)
                    fjz = f_mrigcl(3,ilist,j,k)
                    fix = f_mclrig(1,ilist,j,k)
                    fiy = f_mclrig(2,ilist,j,k)
                    fiz = f_mclrig(3,ilist,j,k)
                    dvdqmclassic(1,j,k)   =  dvdqmclassic(1,j,k) + (dfidp*fjx + dfjdp*fix)
                    dvdqmclassic(2,j,k)   =  dvdqmclassic(2,j,k) + (dfidp*fjy + dfjdp*fiy)
                    dvdqmclassic(3,j,k)   =  dvdqmclassic(3,j,k) + (dfidp*fjz + dfjdp*fiz)

                    !! Ni-Ni on diagonal second deriv.
                    ddx(3+3*j-2,3+3*j-2)=ddx(3+3*j-2,3+3*j-2)+dfidp*fjx**2+dfjdp*fix**2+ &
                        & ddfjd2p*ff2_mclrig(1,1,ilist,j,k)+ddfid2p*ff2_mrigcl(1,1,ilist,j,k)

                    ddx(3+3*j-1,3+3*j-2)=ddx(3+3*j-1,3+3*j-2)+dfidp*fjy*fjx+dfjdp*fiy*fix+ &
                        & ddfjd2p*ff2_mclrig(2,1,ilist,j,k)+ddfid2p*ff2_mrigcl(2,1,ilist,j,k)
                    ddx(3+3*j-1,3+3*j-1)=ddx(3+3*j-1,3+3*j-1)+dfidp*fjy**2+dfjdp*fiy**2+ &
                        & ddfjd2p*ff2_mclrig(2,2,ilist,j,k)+ddfid2p*ff2_mrigcl(2,2,ilist,j,k)

                    ddx(3+3*j,3+3*j-2)=ddx(3+3*j,3+3*j-2)+dfidp*fjz*fjx+dfjdp*fiz*fix+ &
                        & ddfjd2p*ff2_mclrig(3,1,ilist,j,k)+ddfid2p*ff2_mrigcl(3,1,ilist,j,k)
                    ddx(3+3*j,3+3*j-1)=ddx(3+3*j,3+3*j-1)+dfidp*fjz*fjy+dfjdp*fiz*fiy+ &
                        & ddfjd2p*ff2_mclrig(3,2,ilist,j,k)+ddfid2p*ff2_mrigcl(3,2,ilist,j,k)
                    ddx(3+3*j,3+3*j)=ddx(3+3*j,3+3*j)+dfidp*fjz**2+dfjdp*fiz**2+ &
                        & ddfjd2p*ff2_mclrig(3,3,ilist,j,k)+ddfid2p*ff2_mrigcl(3,3,ilist,j,k)

                !!! NO rigid Ni- classic Ni second deriv. off diagonal, only selfinteraction

                enddo
            enddo
        enddo
    end if

    if(h_included) then
        do k = 1,nb
            do j = 1,1
                ph = eden_h(j,k)
                phqk(1) = eden_hq1(j,k)
                phqk(2) = eden_hq2(j,k)

                !           print*,j,k,ph,'rhoh'
                call EAM_EMBEDING_HYDROGEN(ph,phqk,df,dfdp_h(j,k),ddfd2p_h(j,k),dfdpqk_h(j,k,1:2))
                fh = fh + df
            enddo
        enddo
    end if


    if(h_included) then
        !    !  Classic-H
        if(nclassic.ne.0) then
            do k = 1,nb
                do j = 1,nclassic
                    dfjdp = dfdp_mclassic(j,k)
                    ddfjd2p = ddfd2p_mclassic(j,k)
                    do ilist = 1,nlist_mclh(j,k)
                        if(ilist.gt.nclassic) write(*,*) ' classic-H listing error'
                        i = list_mclh(ilist,j,k)
                        dfidp = dfdp_h(i,k)
                        dfidp1 = dfdpqk_h(i,k,1)
                        dfidp2 = dfdpqk_h(i,k,2)
                        ddfid2p = ddfd2p_h(i,k)
                        !              summ = dfidp + dfjdp
                        !              summ2= ddfjd2p + ddfid2p
                        fjx = f_mclh(1,ilist,j,k)
                        fjy = f_mclh(2,ilist,j,k)
                        fjz = f_mclh(3,ilist,j,k)
                        fjx1 = f_mclhqk(1,ilist,j,k,1)
                        fjy1 = f_mclhqk(2,ilist,j,k,1)
                        fjz1 = f_mclhqk(3,ilist,j,k,1)
                        fjx2 = f_mclhqk(1,ilist,j,k,2)
                        fjy2 = f_mclhqk(2,ilist,j,k,2)
                        fjz2 = f_mclhqk(3,ilist,j,k,2)
                        fix = f_hmcl(1,ilist,j,k)
                        fiy = f_hmcl(2,ilist,j,k)
                        fiz = f_hmcl(3,ilist,j,k)
                        dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + (dfidp*fjx + dfjdp*fix)
                        dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + (dfidp*fjy + dfjdp*fiy)
                        dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + (dfidp*fjz + dfjdp*fiz)
                        !              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
                        !              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
                        !              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
                        dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfidp1*fjx1 &
                            &  + dfidp2*fjx2 + dfjdp*fix)
                        dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfidp1*fjy1 &
                            & + dfidp2*fjy2 + dfjdp*fiy)
                        dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfidp1*fjz1 &
                            & + dfidp2*fjz2 + dfjdp*fiz)
                        !! H-H on diagonal second deriv.
                        ddx(1,1)=ddx(1,1)+dfidp*fjx**2+dfjdp*fix**2+ &
                            & ddfjd2p*ff2_mclh(1,1,ilist,j,k)+ddfid2p*ff2_hmcl(1,1,ilist,j,k)

                        ddx(2,1)=ddx(2,1)+dfidp*fjy*fjx+dfjdp*fiy*fix+ &
                            & ddfjd2p*ff2_mclh(2,1,ilist,j,k)+ddfid2p*ff2_hmcl(2,1,ilist,j,k)
                        ddx(2,2)=ddx(2,2)+dfidp*fjy**2+dfjdp*fiy**2+ &
                            & ddfjd2p*ff2_mclh(2,2,ilist,j,k)+ddfid2p*ff2_hmcl(2,2,ilist,j,k)

                        ddx(3,1)=ddx(3,1)+dfidp*fjz*fjx+dfjdp*fiz*fix+ &
                            & ddfjd2p*ff2_mclh(3,1,ilist,j,k)+ddfid2p*ff2_hmcl(3,1,ilist,j,k)
                        ddx(3,2)=ddx(3,2)+dfidp*fjz*fjy+dfjdp*fiz*fiy+ &
                            & ddfjd2p*ff2_mclh(3,2,ilist,j,k)+ddfid2p*ff2_hmcl(3,2,ilist,j,k)
                        ddx(3,3)=ddx(3,3)+dfidp*fjz**2+dfjdp*fiz**2+ &
                            & ddfjd2p*ff2_mclh(3,3,ilist,j,k)+ddfid2p*ff2_hmcl(3,3,ilist,j,k)
                        !!! H-Ni second deriv., off diagonal
                        ddx(3+3*i-2,1)=ddx(3+3*i-2,3+3*j-2)-(dfidp*fjx**2+dfjdp*fix**2+ &
                            & ddfjd2p*ff2_mclh(1,1,ilist,j,k)+ddfid2p*ff2_hmcl(1,1,ilist,j,k))
                        ddx(3+3*i-2,2)=ddx(3+3*i-2,3+3*j-1)-(dfidp*fjx*fjy+dfjdp*fix*fiy+ &
                            & ddfjd2p*ff2_mclh(1,2,ilist,j,k)+ddfid2p*ff2_hmcl(1,2,ilist,j,k))
                        ddx(3+3*i-2,3)=ddx(3+3*i-2,3+3*j)-(dfidp*fjx*fjz+dfjdp*fix*fiz+ &
                            & ddfjd2p*ff2_mclh(1,3,ilist,j,k)+ddfid2p*ff2_hmcl(1,3,ilist,j,k))
                        !!
                        ddx(3+3*i-1,1)=ddx(3+3*i-1,3+3*j-2)-(dfidp*fjy*fjx+dfjdp*fiy*fix+ &
                            & ddfjd2p*ff2_mclh(2,1,ilist,j,k)+ddfid2p*ff2_hmcl(2,1,ilist,j,k))
                        ddx(3+3*i-1,2)=ddx(3+3*i-1,3+3*j-1)-(dfjdp*fiy**2+ &
                            & ddfjd2p*ff2_mclh(2,2,ilist,j,k)+ddfid2p*ff2_hmcl(2,2,ilist,j,k))
                        ddx(3+3*i-1,3)=ddx(3+3*i-1,3+3*j)-(dfidp*fjy*fjz+dfjdp*fiy*fiz+ &
                            & ddfjd2p*ff2_mclh(2,3,ilist,j,k)+ddfid2p*ff2_hmcl(2,3,ilist,j,k))
                        !!
                        ddx(3+3*i,1)=ddx(3+3*i,3+3*j-2)-(dfidp*fjz*fjx+dfjdp*fiz*fix+ &
                            & ddfjd2p*ff2_mclh(3,1,ilist,j,k)+ddfid2p*ff2_hmcl(3,1,ilist,j,k))
                        ddx(3+3*i,2)=ddx(3+3*i,3+3*j-1)-(dfidp*fjz*fjy+dfjdp*fiz*fiy+ &
                            & ddfjd2p*ff2_mclh(3,2,ilist,j,k)+ddfid2p*ff2_hmcl(3,2,ilist,j,k))
                        ddx(3+3*i,3)=ddx(3+3*i,3+3*j)-(+dfidp*fjz**2+dfjdp*fiz**2+ &
                            & ddfjd2p*ff2_mclh(3,3,ilist,j,k)+ddfid2p*ff2_hmcl(3,3,ilist,j,k))



                    enddo
                enddo
            enddo
        end if
        !

        !      !RIGID-H
        if(nrigid.ne.0) then
            do k = 1,nb
                do j = 1,nrigid
                    dfjdp = dfdp_mrigid(j,k)
                    ddfjd2p = dfdp_mrigid(j,k)
                    do ilist = 1,nlist_mrigh(j,k)
                        if(ilist.gt.nrigid) write(*,*) ' rigid-H listing error'
                        i = list_mrigh(ilist,j,k)
                        dfidp = dfdp_h(i,k)
                        dfidp1 = dfdpqk_h(i,k,1)
                        dfidp2 = dfdpqk_h(i,k,2)
                        ddfid2p = ddfd2p_h(i,k)
                        fjx = f_mrigh(1,ilist,j,k)
                        fjy = f_mrigh(2,ilist,j,k)
                        fjz = f_mrigh(3,ilist,j,k)
                        fjx1 = f_mrighqk(1,ilist,j,k,1)
                        fjy1 = f_mrighqk(2,ilist,j,k,1)
                        fjz1 = f_mrighqk(3,ilist,j,k,1)
                        fjx2 = f_mrighqk(1,ilist,j,k,2)
                        fjy2 = f_mrighqk(2,ilist,j,k,2)
                        fjz2 = f_mrighqk(3,ilist,j,k,2)
                        fix = f_hmrig(1,ilist,j,k)
                        fiy = f_hmrig(2,ilist,j,k)
                        fiz = f_hmrig(3,ilist,j,k)
                        !              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
                        !              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
                        !              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
                        dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfidp1*fjx1 + dfidp2*fjx2 + dfjdp*fix)
                        dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfidp1*fjy1 + dfidp2*fjy2 + dfjdp*fiy)
                        dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfidp1*fjz1 + dfidp2*fjz2 + dfjdp*fiz)

                        !! H-H on diagonal second deriv.
                        ddx(1,1)=ddx(1,1)+dfidp*fjx**2+dfjdp*fix**2+ &
                            & ddfjd2p*ff2_mrigh(1,1,ilist,j,k)+ddfid2p*ff2_hmrig(1,1,ilist,j,k)

                        ddx(2,1)=ddx(2,1)+dfidp*fjy*fjx+dfjdp*fiy*fix+ &
                            & ddfjd2p*ff2_mrigh(2,1,ilist,j,k)+ddfid2p*ff2_hmrig(2,1,ilist,j,k)
                        ddx(2,2)=ddx(2,2)+dfidp*fjy**2+dfjdp*fiy**2+ &
                            & ddfjd2p*ff2_mrigh(2,2,ilist,j,k)+ddfid2p*ff2_hmrig(2,2,ilist,j,k)

                        ddx(3,1)=ddx(3,1)+dfidp*fjz*fjx+dfjdp*fiz*fix+ &
                            & ddfjd2p*ff2_mrigh(3,1,ilist,j,k)+ddfid2p*ff2_hmrig(3,1,ilist,j,k)
                        ddx(3,2)=ddx(3,2)+dfidp*fjz*fjy+dfjdp*fiz*fiy+ &
                            & ddfjd2p*ff2_mrigh(3,2,ilist,j,k)+ddfid2p*ff2_hmrig(3,2,ilist,j,k)
                        ddx(3,3)=ddx(3,3)+dfidp*fjz**2+dfjdp*fiz**2+ &
                            & ddfjd2p*ff2_mrigh(3,3,ilist,j,k)+ddfid2p*ff2_hmrig(3,3,ilist,j,k)

                    enddo
                enddo
            enddo
        end if
    !
    !
    end if  ! h_incl
    !
    !

    v = vmh + vmm + fh + fm
    !        print*,v,vmh, vmm, fh, fm,'v,vmh, vmm, fh, fm'
    !
    ! Converting forces, energy and distances in a.u.

    !        call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic, &
    !                       qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, 1)

    !       v=v/dble(nb)
    dvdqmrigid   = 0.d0
!
END subroutine EAM_EAMNIH_HESSIAN


















SUBROUTINE EAM_RHO_NI(r,p,dpdr,DDPD2R,S, DSDR,DDS2DR)
    IMPLICIT NONE
    !      DOUBLE PRECISION :: C(1:10)
    DOUBLE PRECISION :: EPS(1:10) !10, the total number of valence electrons in NI
    DOUBLE PRECISION :: NI(1:10)
    DOUBLE PRECISION :: R,P,DPDR,PS,PD,PPS,PPD,DPDRS,DPDRD
    INTEGER :: I
    DOUBLE PRECISION :: S, DSDR,PI,DDS2DR,DDPD2R,DDPD2RS,DDPD2RD
    !      DOUBLE PRECISION :: FACTORIAL(10)
    DOUBLE PRECISION :: COEF(10), RNI(1:10),DEXPEPS(1:10)
    !!
    PI = DACOS(-1.D0)!3.141592653589793115997963468544185162d0 !
    !!
    NI = 0.D0
    NI(1) = 1.D0
    NI(2) = 1.D0
    NI(3) = 2.D0
    NI(4) = 2.D0
    NI(5) = 3.D0
    NI(6) = 3.D0
    NI(7) = 4.D0
    NI(8) = 4.D0
    NI(9) = 3.D0
    NI(10) = 3.D0
    !!
    EPS = 0.D0
    EPS(1) = 54.88885D0
    EPS(2) = 38.48431D0
    EPS(3) = 27.42703D0
    EPS(4) = 20.88204D0
    EPS(5) = 10.95707D0
    EPS(6) = 7.31958D0
    EPS(7) = 3.92650D0
    EPS(8) = 2.15289D0
    EPS(9) = 12.67582D0
    EPS(10) = 5.43253D0
    !!
    !      C = 0.D0
    !      C(1) = -0.00389D0
    !      C(2) = -0.02991D0
    !      C(3) = -0.03189D0
    !      C(4) =  0.15289D0
    !      C(5) = -0.20048D0
    !      C(6) = -0.05423D0
    !      C(7) =  0.49292D0
    !      C(8) =  0.61875D0
    !      C(9) =  0.42120D0
    !      C(10) = 0.70658D0

    !      FACTORIAL(1) = 2.D0
    !      FACTORIAL(2) = 2.D0
    !      FACTORIAL(3) = 24.D0
    !      FACTORIAL(4) = 24.D0
    !      FACTORIAL(5) = 720.D0
    !      FACTORIAL(6) = 720.D0
    !      FACTORIAL(7) = 40320.D0
    !      FACTORIAL(8) = 40320.D0
    !      FACTORIAL(9) = 720.D0
    !      FACTORIAL(10) = 720.D0

    !     DO I = 1, 10
    !     COEF(I) = C(I)*(1.d0/DSQRT(FACTORIAL(I)))*(2.d0*EPS(I))**(NI(I)+0.5D0)
    !     END DO
    COEF(1) =      -3.163776491313133654159628349589d0
    COEF(2) =     -14.281438867475564791220676852390d0
    COEF(3) =    -145.067738725436015556624624878168d0
    COEF(4) =     351.787786416577375803171889856458d0
    COEF(5) =    -368.078386826022267541702603921294d0
    COEF(6) =     -24.259391731366324762575459317304d0
    COEF(7) =      26.162326823407976661428619991057d0
    COEF(8) =       2.197800222225523736341301628272d0
    COEF(9) =    1287.784878683165061374893411993980d0
    COEF(10) =     111.328804019789757262515195179731d0
    !!
    DO I = 1, 10
        RNI(I) = R**(NI(I)-1.D0)
        DEXPEPS(I) = DEXP(-EPS(I)*R)
    END DO
    !!

    !!
    PS = 0.D0
    !!     ATOMIC NI ELECTRONIC DENSITY for S-state
    DO I = 1, 8
        PS = PS + COEF(I)*RNI(I)*DEXPEPS(I)
    END DO
    !!
    PPS = PS
    PS = PS*PS/(4.d0*PI)
    !!
    PD = 0.D0
    !!     ATOMIC NI ELECTRONIC DENSITY for S-state
    DO I = 9,10
        PD = PD + COEF(I)*RNI(I)*DEXPEPS(I)
    END DO
    !!
    PPD = PD
    PD = PD*PD/(4.d0*PI)
    !!
    P = 2.D0*PS + 8.D0*PD
    !!
    !!     NOW DERIVATIVES
    DPDRS = 0.D0
    !!
    DO I = 1,2
        DPDRS = DPDRS + COEF(I)*(-EPS(I)*DEXPEPS(I))
    END DO
    !!     GOOD CHECK IS TO CYCLE 1,10 THIS ONE
    DO I = 3, 8
        DPDRS = DPDRS + COEF(I)*((NI(I)-1.D0)*R**(NI(I)-2.D0)-EPS(I)*R**(NI(I)-1.D0))*DEXPEPS(I)
    END DO
    !!
    DPDRD = 0.D0
    !!
    DO I = 9, 10
        DPDRD = DPDRD + COEF(I)*((NI(I)-1.D0)*R**(NI(I)-2.D0)-EPS(I)*RNI(I))*DEXPEPS(I)
    END DO
    !!
    !!    NOW SECOND DERIVATIVES
    DDPD2RS = 0.D0
    !!
    DO I = 1,2
        DDPD2RS = DDPD2RS + COEF(I)*(EPS(I)**2*DEXPEPS(I))
    END DO

    DO I = 3, 4
        DDPD2RS = DDPD2RS + COEF(I)*(-2.0D0*EPS(I)*(NI(I)-1.D0)*R**(NI(I)-2.D0)+EPS(I)**2*DEXPEPS(I))
    END DO
    !!     GOOD CHECK IS TO CYCLE 1,10 THIS ONE
    DO I = 5, 6
        DDPD2RS = DDPD2RS + COEF(I)*((NI(I)-1.D0)*(NI(I)-2.D0)*R**(NI(I)-3.D0)- &
            2.0D0*EPS(I)*(NI(I)-1.D0)*R**(NI(I)-2.D0)+EPS(I)**2*DEXPEPS(I))
    END DO
    !!
    DDPD2RD = 0.D0
    !!
    DO I = 9, 10
        DDPD2RD = DDPD2RD + COEF(I)*((NI(I)-1.D0)*(NI(I)-2.D0)*R**(NI(I)-3.D0)- &
            2.0D0*EPS(I)*(NI(I)-1.D0)*R**(NI(I)-2.D0)+EPS(I)**2*DEXPEPS(I))
    END DO

    !!    SECOND DERIVATIVES
    DDPD2RS = (DDPD2RS*PPS+DPDRS**2)/(2.D0*PI)
    DDPD2RD = (DDPD2RD*PPD+DPDRD**2)/(2.D0*PI)
    !!
    DDPD2R = 2.D0*DDPD2RS + 8.D0*DDPD2RD
    !!    FIRST DERIVATIVES
    DPDRS = DPDRS*PPS/(2.D0*PI)
    !!
    DPDRD = DPDRD*PPD/(2.D0*PI)
    !!
    DPDR = 2.D0*DPDRS + 8.D0*DPDRD

    !!   NOW ADDING SMOOTHING FUNCTION
    !     CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
    !!
    DDPD2R = DDPD2R*S+2.D0*DPDR*DSDR+P*DDS2DR
    DPDR = DPDR*S+P*DSDR
    P = P*S
    !!
    !!    CHECKED
    RETURN
END
!!
!
SUBROUTINE EAM_EMBEDING_HYDROGEN(P,RHOA,FH,DFHDP,DDFHD2P,DFHDPQK)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RHOA(2)
    DOUBLE PRECISION :: P, FH, DFHDP
    DOUBLE PRECISION :: XSECH,APPROXSECH,APPROXDSECH,DDFHD2P
    DOUBLE PRECISION :: DEXPSEC2,DEXPSEC,DFHDPQK(2),ALPHBETHGAMM(2,3)
    DOUBLE PRECISION :: AH1,BH1,CH1,DH1,EH1,FH1, ZETA(2,3), DEXPEPSHP !, QK(2)
    INTEGER :: IC
    DOUBLE PRECISION :: PCUT3,DELTA3
    DOUBLE PRECISION :: D5ROH,D5ROH2,D5ROH3,D5ROH4,D5ROH5
    DOUBLE PRECISION :: D6ROH,D6ROH2,D6ROH3,D6ROH4,D6ROH5

    !EAM6 (modified formula!)
    !       DOUBLE PRECISION :: FMEH1(2),EPSILONH1(2),DEXPEPSHPEAM6(2)
    !
    !      FMEH1 = (/ 476.1208, -545.3935 /)
    !      EPSILONH1 = (/ 5.0715, 5.2849 /)


    !EAM7 STUFF
    !        AH1 = 184100000000.D0
    AH1 = 184094523400.D0
    !        BH1 = 3276000000.D0
    BH1 = 3276017408.D0
    !        CH1 = 15100000.D0
    CH1 = 15102027.D0
    !        DH1 =-130400000000.D0
    DH1 =-130388582400.D0
    !        EH1 = 2376000000.D0
    EH1 = 2375889152.D0
    !        FH1 =-11950000D0
    FH1 =-11946214D0
    PCUT3 = 0.009D0
    DELTA3 = 0.007D0
    !        QK = (/ 2.D0 , 0.5D0  /) see earlier part in the programm where densities are calculated
    ALPHBETHGAMM = RESHAPE([-4.8029D0, -74.556D0, & !ALPHAh1, ALPHAH2
        &  34.462D0, 7.2784D0, &    !BETAH1, BETAH2
        &  0.38553D0, 1.0154D0],[2,3]) !GAMMAH1, GAMMAH2
    ! Values from Kindt Tully 1999
    ZETA = RESHAPE([-0.04084D0, -0.00325D0, &
        &  2712.68D0, 54.772D0, &
        &  0.001843D0, 0.812D0],[2,3])
    !Original EAM6 Values
    !        ZETA = RESHAPE([-0.04084D0, -0.00325D0, &
    !                      &  368600000D0, 13563D0, &
    !                      &  54.772D0, 0.812D0],[2,3])



    !        DEXPEPSHPEAM6(1) = FMEH1(1)*DEXP(EPSILONH1(1)*P)
    !        DEXPEPSHPEAM6(2) = FMEH1(2)*DEXP(EPSILONH1(2)*P)

    !Calculate FH
    FH=0.D0
    DFHDP=0.D0
    XSECH=0.D0
    APPROXSECH=0.D0
    APPROXDSECH=0.D0
    DFHDPQK=0.D0
    DDFHD2P=0.D0
    !    DO IC=1,2
    !        FH=FH+DEXPEPSHPEAM6(IC)*P
    !        DFHDP=DFHDP+DEXPEPSHPEAM6(IC)*(1.D0+FMEH1(IC)*P)
    !    ENDDO



    ! EAM7
    IF(0.D0 <= P .AND. P < PCUT3-DELTA3 ) THEN
        DEXPEPSHP=ALPHBETHGAMM(1,1)*DEXP(-ALPHBETHGAMM(1,2)*P)
        FH=DEXPEPSHP*P**ALPHBETHGAMM(1,3)
        DFHDP=DEXPEPSHP*P**(ALPHBETHGAMM(1,3)-1.D0)*(ALPHBETHGAMM(1,3)-ALPHBETHGAMM(1,2)*P)
        DDFHD2P=DEXPEPSHP*P**(ALPHBETHGAMM(1,3)-2.D0)*(ALPHBETHGAMM(1,3)*(ALPHBETHGAMM(1,3)-1.D0)- &
            & 2.D0*ALPHBETHGAMM(1,3)*ALPHBETHGAMM(1,2)*P-(ALPHBETHGAMM(1,2)*P)**2.D0)
    !        print*,FH,DFHDP,'DFHDP1'
    ELSEIF(PCUT3-DELTA3 <= P .AND. P <= PCUT3) THEN
        D5ROH = P - PCUT3
        D5ROH2= D5ROH**2
        D5ROH3= D5ROH**3
        D5ROH4= D5ROH2**2
        D5ROH5= D5ROH2*D5ROH3
        D6ROH = P - (PCUT3 - DELTA3)
        D6ROH2= D6ROH**2
        D6ROH3= D6ROH**3
        D6ROH4= D6ROH2**2
        D6ROH5= D6ROH2*D6ROH3

        FH=AH1*D5ROH5+BH1*D5ROH4+CH1*D5ROH3+DH1*D6ROH5+EH1*D6ROH4+FH1*D6ROH3

        DFHDP=AH1*5*D5ROH4+BH1*4*D5ROH3+CH1*3*D5ROH2+ &
            &  DH1*5*D6ROH4+EH1*4*D6ROH3+FH1*3*D6ROH2
        DDFHD2P=AH1*20*D5ROH3+BH1*12*D5ROH2+CH1*6*D5ROH+ &
            &  DH1*20*D6ROH3+EH1*12*D6ROH2+FH1*6*D6ROH

    !print*,FH,DFHDP,'DFHDP2'
    ELSEIF(P > PCUT3) THEN
        DEXPEPSHP=ALPHBETHGAMM(2,1)*DEXP(-ALPHBETHGAMM(2,2)*P)
        FH=DEXPEPSHP*P**ALPHBETHGAMM(2,3)
        DFHDP=DEXPEPSHP*P**(ALPHBETHGAMM(2,3)-1)*(ALPHBETHGAMM(2,3)-ALPHBETHGAMM(2,2)*P)
        DDFHD2P=DEXPEPSHP*P**(ALPHBETHGAMM(2,3)-2.D0)*(ALPHBETHGAMM(2,3)*(ALPHBETHGAMM(2,3)-1.D0)- &
            & 2.D0*ALPHBETHGAMM(2,3)*ALPHBETHGAMM(2,2)*P-(ALPHBETHGAMM(2,2)*P)**2.D0)
    ENDIF

    !     DO IC=1,2
    !         XSECH=ZETA(IC,2)*(RHOA(IC)-ZETA(IC,3))
    !         DEXPSEC=DEXP(XSECH)
    !         DEXPSEC2=DEXP(2*XSECH)
    !        !Approx sech(x)=1/cosh(x) for calcultions in F77
    !        !  APPROXSECH=1.D0/COSH(XSECH)
    !        !  APPROXDSECH= -APPROXSECH**2/SINH(XSECH)
    !         APPROXSECH=2*DEXPSEC/(DEXPSEC2+1)
    !         APPROXDSECH=-(DEXPSEC2-1)*2*DEXPSEC/(DEXPSEC2+1)**2
    !
    !        FH=FH+ZETA(IC,1)*APPROXSECH
    !        DFHDPQK(IC)= ZETA(IC,1)*ZETA(IC,2)*APPROXDSECH
    !     ENDDO


    RETURN
END
!
SUBROUTINE EAM_EMBEDING_NICKEL(P,FNI,DFNIDP,DDFNID2P)
    IMPLICIT NONE
    DOUBLE PRECISION :: P, FNI, DFNIDP,DDFNID2P
    DOUBLE PRECISION :: ANI,BNI,CNI,DNI,AS,BS,CS,DS,DELTA,PC
    DOUBLE PRECISION :: ALPHA, BETA, GAMMANI,DELTANI
    DOUBLE PRECISION :: ALPHAP,BETAP,GAMMAP,DELTANIP
    DOUBLE PRECISION :: DEXPALPHAP,DEXPBETAP,DEXPGAMMAP, DEXPDELTANIP
    DOUBLE PRECISION :: P2,P3,P4,P5,PmPC,PmPC3,PmPC4
    !
    ANI = -303.28911 !eV*Angstrom^3
    ALPHA = 7.64744D0 !Angstrom^3
    !      BNI = 87.987D0 ! 87.987D0**3 eV^1/3*Angstrom^3
    BNI = 87.98723D0**3 ! 87.987D0**3 eV^1/3*Angstrom^3 ^3
    BETA = 75.07472D0 !Angstrom^3
    CNI = -522.77364D0 !eV*Angstrom^3
    GAMMANI = 373.37826D0 !Angstrom^3
    !      DNI = 39.421D0 ! eV^1/5*Angstrom^3
    DNI = 39.42118D0**5 ! eV^1/5*Angstrom^3 ^5
    DELTANI = 56.34512D0 ! Angstrom^3


    PC = 0.11D0 !Angstrom^-3 !was 0.11 in EAM7
    DELTA = 0.01D0 !Angstrom^-3
    !      AS = -20150000000.d0 ! -2.015D0*10.D0**(10.D0)  !eV*Angstrom^15
    AS = -20145635330.d0 ! -2.015D0*10.D0**(10.D0)  !eV*Angstrom^15
    !      BS = -538400000.d0 !-5.384D0*10.D0**(8.D0) !eV*Angstrom^12
    BS = -538419456.d0 !-5.384D0*10.D0**(8.D0) !eV*Angstrom^12
    !      CS = -4060000.d0 !-4.060D0*10.D0**(6.D0) !eV*Angstrom^9
    CS = -4060125.d0 !-4.060D0*10.D0**(6.D0) !eV*Angstrom^9
    DS = -11.03098054D0 !eV
    !
    !
    IF(P.GT.PC) THEN
        FNI = DS
        DFNIDP = 0.D0
        DDFNID2P = 0.D0
    ELSEIF(P.GT.(PC-DELTA)) THEN
        PmPC = P-PC
        PmPC3 = PmPC**3.d0
        PmPC4 = PmPC**4.d0
        FNI = AS*PmPC**5.D0 + BS*PmPC4 + CS*PmPC3 + DS

        DFNIDP =5.D0*AS*PmPC4 + 4.D0*BS*PmPC3 + 3.D0*CS*PmPC**2.D0
        DDFNID2P = 20.D0*AS*PmPC3 + 12.D0*BS*PmPC**2 + 6.D0*CS*PmPC
    ELSEIF(P.GE.0.D0) THEN
        !
        !
        ALPHAP=alpha*p
        BETAP=beta*p
        GAMMAP=gammaNI*p
        DELTANIP=DELTANI*p
        !
        DEXPALPHAP=ANI*DEXP(-ALPHAP)
        DEXPBETAP=BNI*DEXP(-BETAP)
        DEXPGAMMAP=CNI*DEXP(-GAMMAP)
        DEXPDELTANIP=DNI*DEXP(-DELTANIP)
        !
        P2= P**2
        P3= P**3
        P4= P2**2
        P5= P2*P3

        !
        FNI = DEXPALPHAP*P+P3*DEXPBETAP+DEXPGAMMAP*P+P5*DEXPDELTANIP
        DFNIDP = DEXPALPHAP*(1.D0-ALPHAP) + DEXPBETAP*P2*(3.D0 - BETA*P) + &
            &  DEXPGAMMAP*(1.D0-GAMMAP) +DEXPDELTANIP*P4*(5.D0 - DELTANI*P)
        DDFNID2P = DEXPALPHAP*(1.D0-ALPHAP)*(-ALPHAP) + &
            &  DEXPBETAP*(6.D0*P-6.D0*BETA*P**2+BETA**2.D0*P**3.D0) + &
            &  DEXPGAMMAP*(1.D0-GAMMAP)*(-GAMMAP) + &
            &  DEXPDELTANIP*P3*(20.D0-10.D0*DELTANI*P+(DELTANI*P)**2)


    END IF
    RETURN
END
!
SUBROUTINE EAM_PAIR_NIH(R,PHINIH,DPHINIHDR,DDPHINIHD2R,S, DSDR,DDS2DR)
    !! Calculates the short range pairwise repulsive function phi and its gradient for NIH
    IMPLICIT NONE
    !!
    DOUBLE PRECISION :: C
    DOUBLE PRECISION :: Z0NI, BETANI, ALPHANI
    DOUBLE PRECISION :: R,PHINIH,DPHINIHDR,DDPHINIHD2R
    DOUBLE PRECISION :: DZHDR,DZNIDR,DDZNIDR2
    DOUBLE PRECISION :: ZH, ZNI,APPROXSECH,APPROXDSECH
    DOUBLE PRECISION :: S, DSDR,DDS2DR,DEXPSEC,DEXPSEC2
    DOUBLE PRECISION :: ZHCOEF,ZNICOEF, A1, A2
    DOUBLE PRECISION :: AH(3), BH(3), CH(3), DH(3)
    DOUBLE PRECISION :: ALPHAZ,BETAZ,GAMMAZ,RSW !EAM7 line
    DOUBLE PRECISION :: P2,P3,P4,P5,DDZHD2R
    INTEGER :: IC

    !!
    C = 14.3888D0 !eV*Angstrom
    AH = (/ 0.37473D0, 0.01012D0, 0.01073D0 /)    ! Angstrom ,...
    BH = (/ 0.93232D0, 0.54879D0, 0.37697D0 /)    ! Angstrom ,...
    CH = (/ 2.77511D0, 142.67761D0, 95.4552D0 /)  ! no unit
    DH = (/ 0.00752D0, 8.0036D0, 2.401D0 /) ! no unit, Angstrom^-1, Angstrom

    ALPHAZ = -5.957D0 !A^-5
    BETAZ = 9.804D0 !A^-4
    GAMMAZ = -3.705D0 !A^-3  so im Program so im Paper EAM7 GAMMAZ = -3.715D0 !A^-3
    RSW = 1.09D0 !A

    Z0NI = 10.0D0 !Z effective nuclear charge of Ni atom at distance R
    BETANI = 1.D0/1.116D0 ! Angstrom^-1
    ALPHANI = 1.D0/0.537D0 ! Angstrom^-1

    ZH=0.D0
    DZHDR=0.D0
    DEXPSEC=0.D0
    DEXPSEC2=0.D0
    DDZHD2R=0.D0

    !EAM7 Correction
    IF(R.LE.RSW) THEN
        P2= R**2
        P3= R**3
        P4= P2**2
        P5= P2*P3
        ZH=1+ALPHAZ*P5+BETAZ*P4+GAMMAZ*P3
        DZHDR=5*ALPHAZ*P4+4*BETAZ*P3+3*GAMMAZ*P2
        DDZHD2R=20*ALPHAZ*P3+12*BETAZ*P2+6*GAMMAZ*R
    ELSE

        DO IC=1,3
            ZHCOEF=((R/BH(IC))**(CH(IC)-1))*DEXP(-R/AH(IC))
            ZH=ZH+ZHCOEF*(R/BH(IC))
            DZHDR=DZHDR+ZHCOEF*(CH(IC)/BH(IC)-R/(BH(IC)*AH(IC)))
            DDZHD2R=DDZHD2R+DEXP(-R/AH(IC))/(BH(IC)**CH(IC))*R**(CH(IC)-2)* &
                & (CH(IC)*(CH(IC)-1)-2.D0*CH(IC)/AH(IC)*R+(R/AH(IC))**2)
        ENDDO
        A1=R-DH(3)
        DEXPSEC=DEXP(DH(2)*A1)
        DEXPSEC2=DEXP(2*DH(2)*A1)
        APPROXSECH=2*DEXPSEC/(DEXPSEC2+1)
        APPROXDSECH=-(DEXPSEC2-1)*2*DEXPSEC/(DEXPSEC2+1)**2
        ZH=ZH+DH(1)*APPROXSECH
        DZHDR=DZHDR+DH(1)*DH(2)*APPROXDSECH
        DDZHD2R = DDZHD2R+DH(1)*DH(2)**2*(APPROXSECH**2+ &
            & APPROXDSECH*(DEXPSEC2-1)/(DEXPSEC2+1))
    !!
    ENDIF

    ZNICOEF = Z0NI*DEXP(-ALPHANI*R)
    A2= 1.D0+BETANI*R
    ZNI = ZNICOEF*A2

    DZNIDR = ZNICOEF*(BETANI-ALPHANI*A2)
    DDZNIDR2=ZNICOEF*(-2.D0*ALPHANI*BETANI+ALPHANI**2.D0*A2)
    DDPHINIHD2R = C*((S**2*(DDZNIDR2*ZH+2.D0*DZNIDR*DZHDR+ZNI*DDZHD2R)/R- &
        & 2.D0*(DZNIDR*ZH+ZNI*DZHDR+ZNI*ZH/R)/R**2)+4.D0*(DZNIDR*ZH+&
        & ZNI*DZHDR-ZNI*ZH/R)/R+2.D0*ZNI*ZH/R*(S*DDS2DR+DSDR**2))

    !!    NOW ADDING SMOOTHING FUNCTION
    !!

    DZHDR = DZHDR*S + ZH*DSDR

    DZNIDR = DZNIDR*S + ZNI*DSDR


    !!
    ZH = ZH*S
    ZNI = ZNI*S
    !!
    !PHINIH = 0.D0
    !!
    PHINIH = C*ZH*ZNI/R
    !!
    DPHINIHDR = C*(DZHDR*ZNI + DZNIDR*ZH-ZNI*ZH/R)/R

    !!

    RETURN
END
!
!
!
SUBROUTINE EAM_PAIR_NINI(R,PHININI,DPHININIDR,DDPHININID2R,S,DSDR,DDS2DR)
    ! Calculates the short range pauirwise repulsive function phi and its gradient for NINI
    IMPLICIT NONE
    !
    DOUBLE PRECISION :: C, S, DSDR,DDS2DR,DDPHININID2R
    DOUBLE PRECISION :: Z0NI, BETANI, ALPHANI
    DOUBLE PRECISION :: DZNIDR1!, DZNIDR2
    DOUBLE PRECISION :: ZNI1!, ZNI2
    DOUBLE PRECISION :: R,PHININI,DPHININIDR,DDZNIDR2
    DOUBLE PRECISION :: ZNICOEF,A1,ZNI1S,ZNI1SR
    !
    C = 14.3888D0 !eV*Angstrom
    !
    Z0NI = 10.0D0 !Z effective nuclear charge of Ni atom at distance R
    ! EAM 4 Values
    !      BETANI = 0.8957D0 ! Angstrom^-1
    !      ALPHANI = 1.8633D0 ! Angstrom^-1
    ! EAM 6
    BETANI = 1.D0/1.116D0 ! Angstrom^-1
    ALPHANI = 1.D0/0.537D0 ! Angstrom^-1

    !
    !
    ZNICOEF = Z0NI*DEXP(-ALPHANI*R)
    A1=1.D0+BETANI*R
    ZNI1 = ZNICOEF*A1
    !
    DZNIDR1 = ZNICOEF*(BETANI-ALPHANI*A1)
    DDZNIDR2 = ZNICOEF*(-2.D0*ALPHANI*BETANI+ALPHANI**2.D0*A1)

    !
    !!    NOW ADDING SMOOTHING FUNCTION
    !     CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
    !!
    ! DZNIDR1 = DZNIDR1*S+DSDR*ZNI1
    !!
    !ZNI1 = ZNI1*S
    !
    !!
    !
    !PHININI = 0.D0
    !
    ZNI1S = ZNI1*S
    ZNI1SR=ZNI1S/R
    PHININI = C*ZNI1S*ZNI1SR

    DPHININIDR = (C*ZNI1SR)*(2.d0*(DZNIDR1*S+DSDR*ZNI1)-ZNI1SR)

    DDPHININID2R =2.D0*C*DZNIDR1**2/R*S**2+(C*ZNI1SR*S)*(2.D0*DDZNIDR2- &
        & 4.D0*DZNIDR1/R+2.D0*ZNI1/R**2)+4.D0*C*S*DSDR* &
        & (2.D0*ZNI1*DZNIDR1/R-(ZNI1/R)**2)+ &
        & 2.D0*C*ZNI1**2/R*(S*DDS2DR+DSDR**2)
    !
    RETURN
END
!!
subroutine EAM_RHO_H(r,p,dpdr,DDPD2R,S,DSDR,DDS2DR)
    implicit none
    ! ------------------------------------------------------------------
    ! H Electron Density and derivatives for NI(100)/H EAM4 parametrization: REV B 51, 9985 (1995)
    ! ------------------------------------------------------------------
    double precision :: r,delh,ch,dpdr,p,DDPD2R
    double precision :: s, dsdr,dds2dr
    !
    delh = 3.779452267842503321304548080661334097d0 !2.d0/0.5291772083D0                    !1.30927d0   Tom Markland's values
    ch =  2.148061616605061452389691112330183387d0 !1.d0/(PI*0.5291772083D0**3.d0)          !11.0025d0
    !
    p = ch*dexp(-delh*r)
    dpdr = -delh*p
    DDPD2R = delh**2*P
    !
    !    CALL EAM6_SMOOTHING(dR,Smooth,DSmooth,DDSmooth)
    !
    DDPD2R = DDPD2R*S+2.D0*DPDR*DSDR+P*DDS2DR
    DPDR = DPDR*S+P*DSDR
    P = P*S
    !
    RETURN
END

SUBROUTINE EAM6_SMOOTHING(R,S,DSDR,DDS2DR)
    IMPLICIT NONE
    DOUBLE PRECISION :: R, RC, DELTA,S,DSDR,DDS2DR
    DOUBLE PRECISION :: PI
    !
    PI = DACOS(-1.D0)
    RC = 5.D0
    DELTA = 5.D0
    !
    S = 0.D0
    DSDR = 0.D0
    IF(R.LT.RC) THEN
        S = 1.D0
        DSDR = 0.D0
        DDS2DR = 0.D0
    END IF
    IF(R.GE.RC.AND.R.LE.(RC+DELTA)) THEN
        S = (RC+DELTA-R)/DELTA - 1.D0/(2.D0*PI)*DSIN(PI*(2.D0*R-2.D0*RC-DELTA)/DELTA)
        DSDR = (-1.D0 - DCOS(PI*(2.D0*R-2.D0*RC-DELTA)/DELTA))/DELTA
        DDS2DR = DSIN(PI*(2.D0*R-2.D0*RC-DELTA)/DELTA)*PI*2.D0/(DELTA**2.D0)
    END IF
    IF(R.GT.(RC+DELTA)) THEN
        S = 0.D0
        DSDR = 0.D0
        DDS2DR = 0.D0
    END IF
!
END SUBROUTINE EAM6_SMOOTHING

