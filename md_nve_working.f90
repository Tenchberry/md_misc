program md_nve
  implicit none

  integer, parameter :: maxpart = 5000
  real, dimension (maxpart, 3) :: R
  real, dimension (maxpart) :: fX, fY, fZ, vx, vy, vz, pospx, pospy, pospz
  real :: en, delt, t, tmax, etot, pot, ekin
  real :: temp
  character (len=80) :: atomtype
  integer :: new, npart, i, j, l
  !real, dimension (maxpart, 3) :: dist, box
  real :: ecut, r2, r2inv, r6inv, ff, distx, disty, distz

  t = 0.0
  tmax = 40
  delt = 0.001
  temp = 0.6 !in KbT/epsilon

  en = 0
  etot = 0
  pot = 0
  ekin = 0

  do i = 1, npart
    pospx(i) = 0
    pospy(i) = 0
    pospz(i) = 0 !setting initial force on each particle to be zero
  end do

  call read()
  call init(pospx,pospy,pospz,temp)

  do i = 1, npart
    fX(i) = 0
    fY(i) = 0
    fZ(i) = 0 !setting initial force on each particle to be zero
  end do

  atomtype = "X"

  new = 0

  open(unit=600, file='potenergy.dat')
  open(unit=700, file='kinenergy.dat')
  open(unit=800, file='energytot.dat')
  open(unit=500, file='temperature.dat')

  do 50 while (t.lt.tmax)
    new = new + 1

    call calcforce (i, j, fX, fY, fZ, en, pot)
    !if (mod(new,1000) == 0) then
      !write(*,*) "Check forces for frame no.",t," "
      !do i = 1, npart
        !write(*,*)
        !write(*,*) "x - ",fX(i)," y - ",fY(i)," z - ",fZ(i)," "
        !write(*,*)
      !enddo
    !end if
    call integrate (pospx, pospy, pospz, fX, fY, fZ, en, temp, etot, ekin)
    write(600,*) new, pot
    write(700,*) new, ekin
    write(800,*) new, etot
    write(500,*) new, temp
    !if (mod(new,1000) == 0) then
      !write(*,*) "Check forces for frame no.",t," "
      !do i = 1, npart
        !write(*,*)
        !write(*,*) "x - ",fX(i)," y - ",fY(i)," z - ",fZ(i)," "
        !write(*,*)
      !enddo
      !write(*,*) "------Previous Initial Position Check no.",t,"-------"
      !do i = 1, npart
        !write(*,*)
        !write(*,*) pospx(i)
        !write(*,*) pospy(i)
        !write(*,*) pospz(i)
        !write(*,*)
      !enddo
    !end if
  if (mod(new,1000) == 0) then
    write(200,*) npart
    write(200,*)
    do i = 1, npart
      write(200,*) atomtype, R(i, 1), R(i, 2), R(i, 3), vx(i), vy(i), vz(i)
    end do
  end if
    t = t + delt
  50 end do

stop

  if (t.eq.tmax) then
    write(200,*) npart
    endfile (unit=200)
    close (unit=200)
    endfile (unit=500)
    close (unit=500)
    endfile (unit=600)
    close (unit=600)
    endfile (unit=700)
    close (unit=700)
    endfile (unit=800)
    close (unit=800)
  end if

contains
  subroutine read()

    !Local variables
    character (len=80) :: input, outfile, atomtype

    !300 write (6,*) 'Enter input filename: '
    !read (*,400) input
    !400 format (a)
    open(unit=400, file='LGcluster.xyz')

    !100 write (6,*) 'Enter output filename for trajectories: '
    !read (*,200) outfile
    !200 format (a)
    open (unit=200, file='jtest.xyz')

    read(400,*) npart
    read(400,*)

    write(200,*) npart
    write(200,*)

    !Reading coordinates of xyz file into array

    do i = 1, npart
      read(400,*) atomtype, R(i,1), R(i,2), R(i,3)
      write(200,*) atomtype, R(i,1), R(i,2), R(i,3), vx(i), vy(i), vz(i)
    end do

    !write(*,*) "Reading of file into array R complete"

  end subroutine read

  subroutine init(pospx, pospy, pospz, temp)

    !Initialisation of the trajectories - reading initial positions, calculating initial velocities

    !Local variables
    real, dimension (maxpart) :: vx, vy, vz, pospx, pospy, pospz, v
    real :: sumvx, sumvy, sumvz
    real :: sumv2
    real :: sf, temp

    !write(*,*) "Beginning Initialisation..."
    !write(*,*) npart

    !sumv = 0
    sumv2 = 0
    sumvx = 0
    sumvy = 0
    sumvz = 0
    !sumvx2 = 0
    !sumvy2 = 0
    !sumvz2 = 0

    do i = 1, npart
      !Reading initial positions and velocities (in 3 dimensions)

      !posx(i) = R(i, 1)
      !posy(i) = R(i, 2)
      !posz(i) = R(i, 3)
      !pos(i) = sqrt(R(i, 1)**2 + R(i, 2)**2 + R(i, 3)**2)

      vx(i) = rand() - 0.5
      vy(i) = rand() - 0.5
      vz(i) = rand() - 0.5
      !vx(i) = 0.0
      !vy(i) = 0.0
      !vz(i) = 0.0
      v(i) = vx(i)**2 + vy(i)**2 + vz(i)**2

      sumvx = sumvx + vx(i)
      sumvy = sumvy + vy(i)
      sumvz = sumvz + vz(i)

      !sumvx2 = sumvx2 + vx(i)**2
      !sumvy2 = sumvy2 + vy(i)**2
      !sumvz2 = sumvz2 + vz(i)**2

      !sumv = sumv + v(i)
      sumv2 = sumv2 + v(i) !summing velocities and square velocities in the x, y and z directions
    end do

      sumvx = sumvx/npart
      sumvy = sumvy/npart
      sumvz = sumvz/npart !setting velocity center of mass
      sumv2 = sumv2/npart !Mean-squared velocities

      !sumvx2 = sumvx2/npart
      !sumvy2 = sumvy2/npart
      !sumvz2 = sumvz2/npart !Mean-squared velocities

      sf = sqrt(3*temp/sumv2) !Scale factor for the velocities

    do i = 1, npart
      vx(i) = (vx(i) - sumvx)*sf
      vy(i) = (vy(i) - sumvy)*sf
      vz(i) = (vz(i) - sumvz)*sf !Scaling velocity vector magnitude such that the instantaneous and desired temperature are equal

      pospx(i) = R(i, 1) - vx(i)*delt
      pospy(i) = R(i, 2) - vy(i)*delt
      pospz(i) = R(i, 3) - vz(i)*delt !Estimating previous time step positions - for Verlet scheme

      !write(*,*) "------Previous Initial Position Check-------"
      !write(*,*)
      !write(*,*) pospx(i)
      !write(*,*) pospy(i)
      !write(*,*) pospz(i)
      !write(*,*)

    end do

    !write(*,*) "Initialisation complete"

  end subroutine init

  subroutine calcforce (val1, val2, fX, fY, fZ, en, pot)

    !Local variables
    integer :: val1, val2
    real, dimension(maxpart) :: fX, fY, fZ
    real :: en, ecut, pot
    !real :: boxr, rc, rc2, rc6

    !write(*,*) "Begin calculation of forces..."

    do i = 1, npart
      fX(i) = 0
      fY(i) = 0
      fZ(i) = 0 !setting initial force on each particle to be zero
    end do

    en = 0
    pot = 0
    !boxr = 0

    !do l = 1,3
      !box(l) = 2.0
      !boxr = boxr + box(l)**2
    !end do

    !rc2 = boxr/2.0
    !rc = 1/rc2
    !rc6 = rc**3
    !ecut = 4*rc6*(rc6-1) !Periodic boundary conditions on





    !write(*,*) "------Distances check no.",t,"-------"

    do val1 = 1, npart-1
      do val2 = val1+1, npart

        !r2 = dist(1)**2 + dist(2)**2 + dist(3)**2 !radial distance
        !if (r2.lt.rc2) then    !testing cutoff distances for pairwise interactions

        distx = R(val1, 1) - R(val2, 1)
        disty = R(val1, 2) - R(val2, 2)
        distz = R(val1, 3) - R(val2, 3)
        !distx = distx - boxX*nint(distx/boxX)
        !disty = disty - boxY*nint(disty/boxY)
        !distz = distz - boxZ*nint(distz/boxZ)
        r2 = distx**2 + disty**2 + distz**2
        !distx = distx - boxX*nint(distx/boxX)


        !r2 = dist(1)**2 + dist(2)**2 + dist(3)**2 !radial distance
        !if (r2.lt.rc2) then    !testing cutoff distances for pairwise interactions

        r2inv = 1.0/r2
        r6inv = r2inv**3
        ff = 48.0*r2inv*r6inv*(r6inv-0.5) !Compute Lennard-Jones potential derivative - forces

          !write(*,*) "------Distances check-------"
          !write(*,*)
          !write(*,*) "Radial = ",r2
          !write(*,*) "X component = ",dist(1)," Y component = ",dist(2)," Z component = ",dist(3)
          !write(*,*) "1/r6 = ",r6inv," 1/r2 = ",r2inv
          !write(*,*)

          !write(*,*) "------LJ Derivative check------"
          !write(*,*)
          !write(*,*) pot
          !write(*,*)

        !write(*,*)
        !write(*,*) "Radial = ",r2
        !write(*,*) "X component = ",distx," Y component = ",disty," Z component = ",distz
        !write(*,*) "1/r6 = ",r6inv," 1/r2 = ",r2inv
        !write(*,*) "Derivative of force = ",ff
        !write(*,*)

        fX(val1) = fX(val1) + ff*distx
        fY(val1) = fY(val1) + ff*disty
        fZ(val1) = fZ(val1) + ff*distz !update forces at position i

        fX(val2) = fX(val2) - ff*distx
        fY(val2) = fY(val2) - ff*disty
        fZ(val2) = fZ(val2) - ff*distz  !update forces at position j

          !write(*,*) "----Forces check----"
          !write(*,*)
          !write(*,*) "x - position i: ",fX(i)," position j: ",fX(j)
          !write(*,*) "y - position i: ",fY(i)," position j: ",fY(j)
          !write(*,*) "z - position i: ",fZ(i)," position j: ",fZ(j)
          !write(*,*)

          en = en + 4*r6inv*(r6inv-1)
          pot = en + 4*r6inv*(r6inv-1)

        !end if
      end do
    end do

  !write(*,*) "Forces successfully calculated"
  !write(*,*) "Total Potential = ",en," "

  end subroutine calcforce

  subroutine integrate (pospx, pospy, pospz, fX, fY, fZ, en, temp, etot, ekin)

    !Local variable
    real, dimension(maxpart) :: fX, fY, fZ, pospx, pospy, pospz, vx, vy, vz, v
    real :: posnx, posny, posnz
    real :: etot, en, sumv, sumv2, temp, ekin
    real :: sumvx, sumvy, sumvz, sumvxn, sumvyn, sumvzn, sf !vnew, vnewx, vnewy, vnewz, sumv2, sumv,
    !character (len=80) :: atomtype, outfile

    !write(*,*) "Finding new positions using Verlet intgeration scheme..."

    ekin = 0

    sumv = 0
    sumv2 = 0
    sumvx = 0
    sumvy = 0
    sumvz = 0
    temp = 0
    !sumvx2 = 0
    !sumvy2 = 0
    !sumvz2 = 0

    do i = 1, npart

      !write(*,*) "----Prev positions check----"
      !write(*,*)
      !write(*,*) pospx(i)
      !write(*,*) pospy(i)
      !write(*,*) pospz(i)
      !write(*,*)

      posnx = 2*R(i, 1) - pospx(i) + (delt**2)*fX(i)
      posny = 2*R(i, 2) - pospy(i) + (delt**2)*fY(i)
      posnz = 2*R(i, 3) - pospz(i) + (delt**2)*fZ(i) !New positions found using Verlet integration scheme

      vx(i) = (posnx - pospx(i))/(2*delt)
      vy(i) = (posny - pospy(i))/(2*delt)
      vz(i) = (posnz - pospz(i))/(2*delt) !New velocities

      !vnewx = 0
      !vnewy = 0
      !vnewz = 0

      sumvx = sumvx + vx(i)
      sumvy = sumvy + vy(i)
      sumvz = sumvz + vz(i)
      !sumvx2 = sumvx2 + vx(i)**2
      !sumvy2 = sumvy2 + vy(i)**2
      !sumvz2 = sumvz2 + vz(i)**2

      v(i) = vx(i)**2 + vy(i)**2 + vz(i)**2 !radial velocity
      !sumv = sumv + v(i)
      sumv2 = sumv2 + v(i) !v^2 computed for finding total kinetic energy

      pospx(i) = R(i, 1)
      pospy(i) = R(i, 2)
      pospz(i) = R(i, 3) !update positions of previous time (t - delt)
      R(i, 1) = posnx
      R(i, 2) = posny
      R(i, 3) = posnz !update positions of current time (t + delt)

      !write(*,*) "------New prev positions check # ",i," ------"
      !write(*,*)
      !write(*,*) pospx(i)
      !write(*,*) pospy(i)
      !write(*,*) pospz(i)
      !write(*,*)
      !write(*,*) R(i, 1)
      !write(*,*) R(i, 2)
      !write(*,*) R(i, 3)
      !write(*,*)
      !write(*,*) v(i)
      !write(*,*)

      ! Saving new coordinates into the array R

      !write(*,*) R(i, 1), R(i, 2), R(i, 3)

      !write(200,*) atomtype, R(i, 1), R(i, 2), R(i, 3) !writing new trajectories into outfile

    end do

    sumvxn = sumvx/npart
    sumvyn = sumvy/npart
    sumvzn = sumvz/npart

    !sumvx2n = sumvx2/npart
    !sumvy2n = sumvy2/npart
    !sumvz2n = sumvz2/npart

    !sf = sqrt(3*temp/sumv2)

    !do i = 1, npart
      !vx(i) = (vx(i) - sumvxn)*sf
      !vy(i) = (vy(i) - sumvyn)*sf
      !vz(i) = (vz(i) - sumvzn)*sf

      !v(i) = vx(i)**2 + vy(i)**2 + vz(i)**2 !radial velocity
    !end do

    temp = sumv2 / (3*npart) !updated temperature found using equipartition theorem
    ekin = 0.5*sumv2/npart
    etot = (en+0.5*sumv2)/npart !Total energy per particle



    !write(*,*) "Integration successful"
    !write(*,*) "Temperature = ",temp," "
    !write(*,*) sumv, sumv2
    !write(*,*) "Total Energy per particle = ",etot," "

  end subroutine integrate

end program md_nve
