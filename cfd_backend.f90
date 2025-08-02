!function that implements the jacobi poisson solver
subroutine jacobi(psi, omega, dx, dy, N, max_iter, tol)
    implicit none


    !defines the variables needed for input and output
    integer, intent(in) :: N, max_iter
    real, intent(in) :: dx, dy, tol
    real, intent(in) :: omega(N,N)
    real, intent(inout) :: psi(N,N)

    !defines the intermediate calculated variables needed
    real :: dx2, dy2, den, diff, max_diff
    real, allocatable :: psi_new(:,:)
    integer :: i, j, iter

    !calculates the constants needed
    dx2 = dx * dx
    dy2 = dy * dy
    den = 2.0 * (dx2 + dy2)

    !allocates the array
    allocate(psi_new(N, N))
    psi_new = psi

    !loop that handles the Jacobi iterations
    do iter = 1, max_iter
        max_diff = 0.0

        !updates the interior points
        do i = 2, N - 1
            do j = 2, N - 1
                !updates according to the Jacobi formula
                psi_new(i, j) = ((psi(i+1, j) + psi(i-1, j)) * dy2 + (psi(i, j+1) + psi(i, j-1)) * dx2 + (-omega(i, j) * dx2 * dy2)) / den

                !tracks the maximum difference
                diff = abs(psi_new(i, j) - psi(i, j))
                if (diff > max_diff) max_diff = diff
            end do 
        end do 
    
        !checks for convergence
        psi = psi_new
        if (max_diff < tol) exit
    end do
    deallocate(psi_new)
end subroutine jacobi

!function that implements the Gauss Seidel poisson solver
subroutine gauss_seidel(psi, omega, dx, dy, N, max_iter, tol)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N, max_iter
    real, intent(in) :: dx, dy, tol
    real, intent(in) :: omega(N,N)
    real, intent(inout) :: psi(N,N)

    !defines the intermediate calculated variables needed
    real :: dx2, dy2, den, diff, max_diff, new_val
    integer :: i, j, iter

    !calculates the constants needed
    dx2 = dx * dx
    dy2 = dy * dy
    den = 2.0 * (dx2 + dy2)

    !loop that handles the Jacobi iterations
    do iter = 1, max_iter
        max_diff = 0.0

        !updates the interior points
        do i = 2, N - 1
            do j = 2, N - 1
                !updates according to the Gauss Seidel formula
                new_val = ((psi(i+1, j) + psi(i-1, j)) * dy2 + (psi(i, j+1) + psi(i, j-1)) * dx2 + (-omega(i, j) * dx2 * dy2)) / den

                !tracks the maximum difference
                diff = abs(new_val - psi(i, j))
                if (diff > max_diff) max_diff = diff
                psi(i, j) = new_val
            end do 
        end do 
    
        !checks for convergence
        if (max_diff < tol) exit
    end do
end subroutine gauss_seidel

!function that implements the succesive over relaxation poisson solver
subroutine sor(psi, omega, dx, dy, N, max_iter, tol, omega_relax)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N, max_iter
    real, intent(in) :: dx, dy, tol, omega_relax
    real, intent(in) :: omega(N,N)
    real, intent(inout) :: psi(N,N)

    !defines the intermediate calculated variables needed
    real :: dx2, dy2, den, diff, max_diff, new_val
    integer :: i, j, iter

    !calculates the constants needed
    dx2 = dx * dx
    dy2 = dy * dy
    den = 2.0 * (dx2 + dy2)

    !loop that handles the Jacobi iterations
    do iter = 1, max_iter
        max_diff = 0.0

        !updates the interior points
        do i = 2, N - 1
            do j = 2, N - 1
                !updates according to the Gauss Seidel formula and adds relaxation factor
                new_val = ((psi(i+1, j) + psi(i-1, j)) * dy2 + (psi(i, j+1) + psi(i, j-1)) * dx2 + (-omega(i, j) * dx2 * dy2)) / den
                new_val = (1.0 - omega_relax) * psi(i,j) + omega_relax * new_val

                !tracks the maximum difference
                diff = abs(new_val - psi(i, j))
                if (diff > max_diff) max_diff = diff
                psi(i, j) = new_val
            end do 
        end do 
    
        !checks for convergence
        if (max_diff < tol) exit
    end do
end subroutine sor

!functin that handles the right hand side of the time solver
subroutine compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, rhs)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: rhs(N,N)   

    !defines the intermediate calculated variables needed
    integer :: i, j
    real :: dwdx, dwdy, d2wdx2, d2wdy2

    !computes the values for vorticity and streamfunction in the interior
    rhs = 0.0
    do i = 2, N - 1
        do j = 2, N - 1
            !computes the first derivative using the central difference formula
            dwdx = (omega(i+1,j) - omega(i-1,j)) / (2.0 * dx)
            dwdy = (omega(i,j+1) - omega(i,j-1)) / (2.0 * dy)

            !computes the second derivative using the central difference formula
            d2wdx2 = (omega(i+1,j) - 2.0*omega(i,j) + omega(i-1,j)) / (dx * dx)
            d2wdy2 = (omega(i,j+1) - 2.0*omega(i,j) + omega(i,j-1)) / (dy * dy)

            !updates vorticity transport with streamfunction
            rhs(i,j) = -x_vel(i,j) * dwdx - y_vel(i,j) * dwdy + nu * (d2wdx2 + d2wdy2)
        end do
    end do
end subroutine compute_rhs

!function that implements Euler's method (rk1)
subroutine rk1(omega, x_vel, y_vel, dx, dy, N, dt, nu, omega_new)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, dt, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: omega_new(N,N) 
    
    !defines the array to store the right hand side
    real, allocatable :: rhs(:, :)
    allocate(rhs(N, N))

    !calls the function the function to compute the rhs
    call compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, rhs)

    !updates according to Euler's method
    omega_new = omega + dt * rhs

    !deallocates arrays
    deallocate(rhs)
end subroutine rk1

!function that implements Heun's method (rk2)
subroutine rk2(omega, x_vel, y_vel, dx, dy, N, dt, nu, omega_new)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, dt, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: omega_new(N,N) 
    
    !defines the arrays to store the right hand side and the k values
    real, allocatable :: omega_predict(:, :), k1(:, :), k2(:, :)
    allocate(omega_predict(N, N), k1(N, N), k2(N, N))

    !k1
    call compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, k1)
    omega_predict = omega + dt * k1

    !k2
    call compute_rhs(omega_predict, x_vel, y_vel, dx, dy, N, nu, k2)

    !computes the final update
    omega_new = omega + 0.5 * dt * (k1 + k2)

    !deallocates arrays
    deallocate(omega_predict, k1, k2)
end subroutine rk2

!function that implements classical fourth order Runge-Kutta method (rk4)
subroutine rk4(omega, x_vel, y_vel, dx, dy, N, dt, nu, omega_new)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, dt, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: omega_new(N,N) 

    !defines the arrays to store the right hand side and the k values
    real, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
    real, allocatable :: omega_temp(:,:)
    allocate(k1(N,N), k2(N,N), k3(N,N), k4(N,N), omega_temp(N,N))

    !k1
    call compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, k1)

    !k2
    omega_temp = omega + 0.5 * dt * k1
    call compute_rhs(omega_temp, x_vel, y_vel, dx, dy, N, nu, k2)

    !k3
    omega_temp = omega + 0.5 * dt * k2
    call compute_rhs(omega_temp, x_vel, y_vel, dx, dy, N, nu, k3)

    !k4
    omega_temp = omega + dt * k3
    call compute_rhs(omega_temp, x_vel, y_vel, dx, dy, N, nu, k4)

    !computes the final update
    omega_new = omega + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)

    !deallocates arrays
    deallocate(k1, k2, k3, k4, omega_temp)
end subroutine rk4

!function that implements Verner's sixth order Runge-Kutta method (rk6)
subroutine rk6(omega, x_vel, y_vel, dx, dy, N, dt, nu, omega_new)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, dt, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: omega_new(N,N)

    !defines the arrays to store the right hand side and the k values
    real, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:), k5(:,:), k6(:,:), k7(:,:), k8(:,:)
    real, allocatable :: temp(:,:)

    !defines the coefficeints for Verners 6th order method
    real :: a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65
    real :: a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87
    real :: b1, b2, b3, b4, b5, b6, b7, b8
    a21 = 1.0/12.0; a31 = 0.0; a32 = 1.0/6.0
    a41 = 1.0/16.0; a42 = 0.0; a43 = 3.0/16.0
    a51 = 21.0/16.0; a52 = 0.0; a53 = -81.0/16.0; a54 = 9.0/2.0
    a61 = 1344688.0/250563.0; a62 = 0.0; a63 = -1709184.0/83521.0; a64 = 1365632.0/83521.0; a65 = -78208.0/250563.0
    a71 = -559.0/384.0; a72 = 0.0; a73 = 6.0; a74 = -204.0/47.0; a75 = 14.0/39.0; a76 = -4913.0/78208.0
    a81 = -625.0/224.0; a82 = 0.0; a83 = 12.0; a84 = -456.0/47.0; a85 = 48.0/91.0; a86 = 14739.0/136864.0; a87 = 6.0/7.0
    b1 = 7.0/90.0; b2 = 0.0; b3 = 0.0; b4 = 16.0/45.0; b5 = 16.0/45.0; b6 = 0.0; b7 = 2.0/15.0; b8 = 7.0/90.0

    !allocates memroy to the arrays
    allocate(k1(N,N), k2(N,N), k3(N,N), k4(N,N), k5(N,N), k6(N,N), k7(N,N), k8(N,N), temp(N,N))

    !k1
    call compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, k1)

    !k2
    temp = omega + dt * a21 * k1
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k2)

    !k3
    temp = omega + dt * (a31 * k1 + a32 * k2)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k3)

    !k4
    temp = omega + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k4)

    !k5
    temp = omega + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k5)

    !k6
    temp = omega + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k6)

    !k7
    temp = omega + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k7)

    !k8
    temp = omega + dt * (a81 * k1 + a82 * k2 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k8)

    !computes the final update
    omega_new = omega + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 + b8 * k8)

    !deallocates arrays
    deallocate(k1, k2, k3, k4, k5, k6, k7, k8, temp)
end subroutine rk6

!function that implements Verner's eighth order Runge-Kutta method (rk8)
subroutine rk8(omega, x_vel, y_vel, dx, dy, N, dt, nu, omega_new)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy, dt, nu
    real, intent(in) :: omega(N,N), x_vel(N,N), y_vel(N,N)
    real, intent(out) :: omega_new(N,N)

    !defines the arrays to store the right hand side and the k values
    real, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:), k5(:,:), k6(:,:), k7(:,:), k8(:,:), k9(:,:), k10(:,:), k11(:,:), k12(:,:), k13(:,:), k14(:,:), temp(:,:)

    !defines the coefficients for Verner’s 8th order method
    real :: s6
    real :: a21, a31, a32, a41, a42, a43, a51, a52, a53, a54
    real :: a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76
    real :: a81, a82, a83, a84, a85, a86, a87
    real :: a91, a92, a93, a94, a95, a96, a97, a98
    real :: a101, a102, a103, a104, a105, a106, a107, a108, a109
    real :: a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110
    real :: a121, a122, a123, a124, a125, a126, a127, a128, a129, a1210, a1211
    real :: a131, a132, a133, a134, a135, a136, a137, a138, a139, a1310, a1311, a1312
    real :: a141, a142, a143, a144, a145, a146, a147, a148, a149, a1410, a1411, a1412, a1413
    real :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14
    s6 = sqrt(6.0)
    a21 = 1.0/12.0
    a31 = 1.0/27.0; a32 = 2.0/27.0
    a41 = 1.0/24.0; a42 = 0.0; a43 = 1.0/8.0
    a51 = (4.0 + 94.0 * s6)/375.0; a52 = 0.0; a53 = (-94.0 - 84.0 * s6)/125.0; a54 = (328.0 + 208.0 * s6)/375.0
    a61 = (9.0 - s6)/150.0; a62 = 0.0; a63 = 0.0; a64 = (312.0 + 32.0 * s6)/1425.0; a65 = (69.0 + 29.0 * s6)/570.0
    a71 = (927.0 - 347.0 * s6)/1250.0; a72 = 0.0; a73 = 0.0; a74 = (-16248.0 + 7328.0 * s6)/9375.0; a75 = (-489.0 + 179.0 * s6)/3750.0; a76 = (14268.0 - 5798.0 * s6)/9375.0
    a81 = 2.0/27.0; a82 = 0.0; a83 = 0.0; a84 = 0.0; a85 = 0.0; a86 = (16.0 - s6)/54.0; a87 = (16.0 + s6)/54.0
    a91 = 19.0/256.0; a92 = 0.0; a93 = 0.0; a94 = 0.0; a95 = 0.0; a96 = (118.0 - 23.0 * s6)/512.0; a97 = (118.0 + 23.0 * s6)/512.0; a98 = -9.0/256.0
    a101 = 11.0/144.0; a102 = 0.0; a103 = 0.0; a104 = 0.0; a105 = 0.0; a106 = (266.0 - s6)/864.0; a107 = (266.0 + s6)/864.0; a108 = -1.0/16.0; a109 = -9.0/27.0
    a111 = (5034.0 - 271.0 * s6)/61440.0; a112 = 0.0; a113 = 0.0; a114 = 0.0; a115 = 0.0; a116 = 0.0; a117 = (7859.0 - 1626.0 * s6)/10240.0; a118 = (-2232.0 - 813.0 * s6)/20480.0
    a119 = (-594.0 + 271.0 * s6)/960.0; a1110 = 3065993473.0/597172653.0
    a121 = (5996.0 - 3795.0 * s6)/405.0; a122 = 0.0; a123 = 0.0; a124 = 0.0; a125 = 0.0
    a126 = (-4342.0 - 338.0 * s6)/9.0; a127 = (154922.0 - 40458.0 * s6)/135.0; a128 = (-4176.0 + 3794.0 * s6)/45.0
    a129 = (-340864.0 + 242816.0 * s6)/405.0; a1210 = (26304.0 - 15176.0 * s6)/45.0; a1211 = -26624.0/81.0
    a131 = (3793.0 + 2168.0 * s6)/103680.0; a132 = 0.0; a133 = 0.0; a134 = 0.0; a135 = 0.0
    a136 = (4042.0 + 2263.0 * s6)/13824.0; a137 = (-231278.0 + 40717.0 * s6)/69120.0
    a138 = (7947.0 - 2168.0 * s6)/11520.0; a139 = (1048.0 - 542.0 * s6)/405.0
    a1310 = (-1383.0 + 542.0 * s6)/720.0; a1311 = 2624.0/1053.0; a1312 = 3.0/1664.0
    a141 = -137.0/1296.0; a142 = 0.0; a143 = 0.0; a144 = 0.0; a145 = 0.0
    a146 = (5642.0 - 337.0 * s6)/864.0; a147 = (5642.0 + 337.0 * s6)/864.0
    a148 = -299.0/48.0; a149 = 184.0/81.0; a1410 = -44.0/9.0
    a1411 = -5120.0/1053.0; a1412 = -11.0/468.0; a1413 = 16.0/9.0
    b1 = 103.0/1680.0; b2 = 0.0; b3 = 0.0; b4 = 0.0; b5 = 0.0; b6 = 0.0; b7 = 0.0
    b8 = -27.0/140.0; b9 = 76.0/105.0; b10 = -201.0/280.0; b11 = 1024.0/1365.0
    b12 = 3.0/7280.0; b13 = 12.0/35.0; b14 = 9.0/280.0

    !allocates memory to the arrays
    allocate(k1(N,N), k2(N,N), k3(N,N), k4(N,N), k5(N,N), k6(N,N), k7(N,N), k8(N,N), k9(N,N), k10(N,N), k11(N,N), k12(N,N), k13(N,N), k14(N,N), temp(N,N))

    !k1
    call compute_rhs(omega, x_vel, y_vel, dx, dy, N, nu, k1)

    !k2
    temp = omega + dt * a21 * k1
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k2)

    !k3
    temp = omega + dt * (a31 * k1 + a32 * k2)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k3)

    !k4
    temp = omega + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k4)

    !k5
    temp = omega + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k5)

    !k6
    temp = omega + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k6)

    !k7
    temp = omega + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k7)

    !k8
    temp = omega + dt * (a81 * k1 + a82 * k2 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k8)

    !k9
    temp = omega + dt * (a91 * k1 + a92 * k2 + a93 * k3 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 + a98 * k8)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k9)

    !k10
    temp = omega + dt * (a101 * k1 + a102 * k2 + a103 * k3 + a104 * k4 + a105 * k5 + a106 * k6 + a107 * k7 + a108 * k8 + a109 * k9)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k10)

    !k11
    temp = omega + dt * (a111 * k1 + a112 * k2 + a113 * k3 + a114 * k4 + a115 * k5 + a116 * k6 + a117 * k7 + a118 * k8 + a119 * k9 + a1110 * k10)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k11)

    !k12
    temp = omega + dt * (a121 * k1 + a122 * k2 + a123 * k3 + a124 * k4 + a125 * k5 + a126 * k6 + a127 * k7 + a128 * k8 + a129 * k9 + a1210 * k10 + a1211 * k11)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k12)

    !k13
    temp = omega + dt * (a131 * k1 + a132 * k2 + a133 * k3 + a134 * k4 + a135 * k5 + a136 * k6 + a137 * k7 + a138 * k8 + a139 * k9 + a1310 * k10 + a1311 * k11 + a1312 * k12)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k13)

    !k14
    temp = omega + dt * (a141 * k1 + a142 * k2 + a143 * k3 + a144 * k4 + a145 * k5 + a146 * k6 + a147 * k7 + a148 * k8 + a149 * k9 + a1410 * k10 + a1411 * k11 + a1412 * k12 + a1413 * k13)
    call compute_rhs(temp, x_vel, y_vel, dx, dy, N, nu, k14)

    !computes the final update
    omega_new = omega + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 + b12 * k12 + b13 * k13 + b14 * k14)

    !deallocates arrays
    deallocate(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, temp)
end subroutine rk8

!function that comnputes the velocity from the streamfunction
subroutine compute_vel(psi, dx, dy, N, x_vel, y_vel)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy
    real, intent(in) :: psi(N,N)
    real, intent(out) :: x_vel(N,N), y_vel(N,N)
    integer :: i, j

    !calculates velocity using central difference formula
    do i = 2, N - 1
        do j = 2, N - 1
            x_vel(i,j) = (psi(i,j+1) - psi(i,j-1)) / (2.0 * dy)
            y_vel(i,j) = -(psi(i+1,j) - psi(i-1,j)) / (2.0 * dx)
        end do
    end do
end subroutine compute_vel

!function that computes the right hand side of the auxilary pressure function
subroutine compute_pressure_rhs(x_vel, y_vel, dx, dy, N, rhs_p)
    implicit none

    !defines the variables needed for input and output
    integer, intent(in) :: N
    real, intent(in) :: dx, dy
    real, intent(in) :: x_vel(N,N), y_vel(N,N)
    real, intent(out) :: rhs_p(N,N)

    !defines the intermediate calculated variables needed
    integer :: i, j
    real :: dudx, dudy, dvdx, dvdy

    !computes pressure using the central difference formula
    rhs_p = 0.0
    do i = 2, N-1
        do j = 2, N-1
            dudx = (x_vel(i+1,j) - x_vel(i-1,j)) / (2.0 * dx)
            dudy = (x_vel(i,j+1) - x_vel(i,j-1)) / (2.0 * dy)
            dvdx = (y_vel(i+1,j) - y_vel(i-1,j)) / (2.0 * dx)
            dvdy = (y_vel(i,j+1) - y_vel(i,j-1)) / (2.0 * dy)

            rhs_p(i,j) = - (dudx**2 + 2.0*dudy*dvdx + dvdy**2)
        end do
    end do
end subroutine compute_pressure_rhs

!function that handles the entire cfd calculation
program cfd_backend
    implicit none

    !defines the variables needed for input and output
    real :: L, U, v, T, cyl_x, cyl_y, cyl_r
    integer :: N, nt
    character(len=128) :: poisson_solver, time_solver

    !defines the simulation grid and arrays
    real :: dx, dy, dt, x, y, tol, omega_relax, time, cd, cl, drag_force, lift_force, pressure_contrib, area
    real, allocatable :: psi(:,:), omega(:,:), x_vel(:,:), y_vel(:,:)
    real, allocatable :: omega_new(:,:), pressure(:,:), rhs_p(:,:)
    logical, allocatable :: mask(:,:)
    logical :: cylinder_on

    !defines the loop counters
    integer :: step, i, j, max_iter

    !reads the parameters from file
    open(unit=10, file="parameters.txt")
    read(10, *) L, N, U, v, T, nt
    read(10, '(A)') poisson_solver
    read(10, '(A)') time_solver
    read(10, *) cyl_x, cyl_y, cyl_r
    close(10)

    !sets up domain variables
    dx = L / real(N)
    dy = dx
    dt = T / real(nt)
    cylinder_on = (cyl_r > 0.0)

    !allocates arrays
    allocate(psi(N,N), omega(N,N), x_vel(N,N), y_vel(N,N))
    allocate(omega_new(N,N), pressure(N,N), rhs_p(N,N), mask(N,N))

    !sets up arrays
    psi = 0.0
    omega = 0.0
    x_vel = 0.0
    y_vel = 0.0
    mask = .true.

    !initializes the inlet streamfunction profile
    do j = 1, N
        psi(1,j) = U * dy * real(j - 1)
    end do

    !seeds vorticity to initiate dynamics
    omega(N/2-2:N/2+2, N/2-2:N/2+2) = 5.0

    !sets up the cylinder mask if needed
    if (cylinder_on) then
        do i = 1, N
            do j = 1, N
                x = (i - 1) * dx
                y = (j - 1) * dy
                if ((x - cyl_x)**2 + (y - cyl_y)**2 < cyl_r**2) then
                    mask(i,j) = .false.
                end if
            end do
        end do
    end if

    !sets the solver parameters
    tol = 1e-6
    max_iter = 10000
    omega_relax = 1.5

    !solves the initial streamfunction
    select case (trim(poisson_solver))
        case ("Jacobi")
            call jacobi(psi, omega, dx, dy, N, max_iter, tol)
        case ("Gauss-Seidel")
            call gauss_seidel(psi, omega, dx, dy, N, max_iter, tol)
        case ("SOR")
            call sor(psi, omega, dx, dy, N, max_iter, tol, omega_relax)
    end select

    !computes the initial velocity
    call compute_vel(psi, dx, dy, N, x_vel, y_vel)

    !enforces the inlet velocity boundary condition
    x_vel(1,:) = U
    y_vel(1,:) = 0.0

    !opens the result files
    open(20, file="results/psi.csv", status="replace")
    open(21, file="results/omega.csv", status="replace")
    open(22, file="results/x_vel.csv", status="replace")
    open(23, file="results/y_vel.csv", status="replace")
    open(24, file="results/pressure.csv", status="replace")
    open(unit=99, file="progress.txt", status="replace")

    !main time loop
    do step = 1, nt
        !define the current time for saving to the file
        time = real(step-1) * dt

        !updates the current time for the frontend progress bar
        rewind(99)
        write(99,*) step
        flush(99)

        !writes the output to file
        write(20, '(A,F10.5)') '# time=', time
        do i = 1, N
            write(20,'(*(E20.8))') (psi(i,j), j=1,N)
        end do

        write(21, '(A,F10.5)') '# time=', time
        do i = 1, N
            write(21,'(*(E20.8))') (omega(i,j), j=1,N)
        end do

        write(22, '(A,F10.5)') '# time=', time
        do i = 1, N
            write(22,'(*(E20.8))') (x_vel(i,j), j=1,N)
        end do

        write(23, '(A,F10.5)') '# time=', time
        do i = 1, N
            write(23,'(*(E20.8))') (y_vel(i,j), j=1,N)
        end do

        write(24, '(A,F10.5)') '# time=', time
        do i = 1, N
            write(24,'(*(E20.8))') (pressure(i,j), j=1,N)
        end do

        !enforces the vorticity boundary conditions 
        do i = 2, N - 1
            omega(i,1) = -2.0 * psi(i,2) / (dy * dy)          
            omega(i,N) = -2.0 * psi(i,N-1) / (dy * dy)         
        end do
        do j = 2, N - 1
            omega(1,j) = -2.0 * (psi(2,j) - U * dy * (j - 1)) / (dx * dx)   
            omega(N,j) = -2.0 * psi(N-1,j) / (dx * dx)                    
        end do

        !solves the Poisson equation for the streamfunction
        select case (trim(poisson_solver))
            case ("Jacobi")
                call jacobi(psi, omega, dx, dy, N, max_iter, tol)
            case ("Gauss-Seidel")
                call gauss_seidel(psi, omega, dx, dy, N, max_iter, tol)
            case ("SOR")
                call sor(psi, omega, dx, dy, N, max_iter, tol, omega_relax)
        end select

        !enforces the cylinder mask for streamfunction
        if (cylinder_on) then
            do i = 1, N
                do j = 1, N
                    if (.not. mask(i,j)) psi(i,j) = 0.0
                end do
            end do
        end if

        !computes the velocity from the streamfunction
        call compute_vel(psi, dx, dy, N, x_vel, y_vel)

        !enforces the velocity boundary conditions
        x_vel(1,:) = U
        y_vel(1,:) = 0.0
        do i = 1, N
            x_vel(i,1) = 0.0
            x_vel(i,N) = 0.0
            y_vel(i,1) = 0.0
            y_vel(i,N) = 0.0
        end do
        do j = 2, N - 1
            x_vel(N,j) = x_vel(N-1,j)
            y_vel(N,j) = y_vel(N-1,j)
        end do

        !enforces the cylinder mask for velocity
        if (cylinder_on) then
            do i = 1, N
                do j = 1, N
                    if (.not. mask(i,j)) then
                        x_vel(i,j) = 0.0
                        y_vel(i,j) = 0.0
                    end if
                end do
            end do
        end if

        !computes the pressure right hand side and solves for pressure
        call compute_pressure_rhs(x_vel, y_vel, dx, dy, N, rhs_p)
        select case (trim(poisson_solver))
            case ("Jacobi")
                call jacobi(pressure, rhs_p, dx, dy, N, max_iter, tol)
            case ("Gauss-Seidel")
                call gauss_seidel(pressure, rhs_p, dx, dy, N, max_iter, tol)
            case ("SOR")
                call sor(pressure, rhs_p, dx, dy, N, max_iter, tol, omega_relax)
        end select

        !enforces the pressure boundary conditions
        do j = 2, N - 1
            pressure(1,j) = pressure(2,j)
            pressure(N,j) = pressure(N-1,j)
        end do
        do i = 2, N - 1
            pressure(i,1) = pressure(i,2)
            pressure(i,N) = pressure(i,N-1)
        end do
        pressure(1,1) = 0.0

        !enforces the cylinder mask for pressure
        if (cylinder_on) then
            do i = 1, N
                do j = 1, N
                    if (.not. mask(i,j)) pressure(i,j) = 0.0
                end do
            end do
        end if

        !performs time integration for vorticity
        select case (trim(time_solver))
            case ("Euler")
                call rk1(omega, x_vel, y_vel, dx, dy, N, dt, v, omega_new)
            case ("Heun")
                call rk2(omega, x_vel, y_vel, dx, dy, N, dt, v, omega_new)
            case ("RK4")
                call rk4(omega, x_vel, y_vel, dx, dy, N, dt, v, omega_new)
            case ("RK6")
                call rk6(omega, x_vel, y_vel, dx, dy, N, dt, v, omega_new)
            case ("RK8")
                call rk8(omega, x_vel, y_vel, dx, dy, N, dt, v, omega_new)
        end select

        !updates vorticity
        omega = omega_new

        !reapplies the vorticity boundary conditions after update
        do i = 2, N - 1
            omega(i,1) = -2.0 * psi(i,2) / (dy * dy)
            omega(i,N) = -2.0 * psi(i,N-1) / (dy * dy)
        end do
        do j = 2, N - 1
            omega(1,j) = -2.0 * (psi(2,j) - U * dy * (j - 1)) / (dx * dx)
            omega(N,j) = -2.0 * psi(N-1,j) / (dx * dx)
        end do

        !enforces the cylinder mask for vorticity
        if (cylinder_on) then
            do i = 1, N
                do j = 1, N
                    if (.not. mask(i,j)) omega(i,j) = 0.0
                end do
            end do
        end if

        !calcultes lift and drag forces
        drag_force = 0.0
        lift_force = 0.0
        do i = 2, N - 1
            do j = 2, N - 1
                if (.not. mask(i,j)) then
                    !computes the x normal forces
                    if (mask(i+1,j)) then
                        pressure_contrib = pressure(i+1,j)
                        drag_force = drag_force - pressure_contrib * dy
                    else if (mask(i-1,j)) then
                        pressure_contrib = pressure(i-1,j)
                        drag_force = drag_force + pressure_contrib * dy
                    end if

                    !computes the y normal forces
                    if (mask(i,j+1)) then
                        pressure_contrib = pressure(i,j+1)
                        lift_force = lift_force - pressure_contrib * dx
                    else if (mask(i,j-1)) then
                        pressure_contrib = pressure(i,j-1)
                        lift_force = lift_force + pressure_contrib * dx
                    end if

                end if
            end do
        end do

        !computes the forces as coefficients
        area = 2.0 * cyl_r
        cd = drag_force / (0.5 * U * U * area)
        cl = lift_force / (0.5 * U * U * area)

        !writes each quantity to its own file
        open(unit=91, file="results/time.csv", status="unknown", position="append")
        write(91,'(E20.8)') time
        close(91)

        open(unit=92, file="results/drag.csv", status="unknown", position="append")
        write(92,'(E20.8)') drag_force
        close(92)

        open(unit=93, file="results/lift.csv", status="unknown", position="append")
        write(93,'(E20.8)') lift_force
        close(93)

        open(unit=94, file="results/cd.csv", status="unknown", position="append")
        write(94,'(E20.8)') cd
        close(94)

        open(unit=95, file="results/cl.csv", status="unknown", position="append")
        write(95,'(E20.8)') cl
        close(95)
    end do

    !closes the result files
    close(20); close(21); close(22); close(23); close(24); close(99)

    !writes the cylinder mask to file
    open(25, file="results/mask.csv", status="replace")
    do i = 1, N
        write(25,*) (mask(i,j), j=1,N)
    end do
    close(25)
end program cfd_backend 
