close all;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%problem parameters

L = 1; %size of square domain
T = 1; %duration of time simulation

rho = 1.0; %density of fluid
mu = 1.7E-5; %viscosity of fluid
nu = .1; %kinematic viscosity

N_x = 20; %number of spatial grid points in each direction
N_t = 2E3; %number of steps in time

x = linspace(0,L,N_x); %spatial grid (row)
t = linspace(0,T,N_t); %time grid (row)

dx = x(2) - x(1); dt = t(2) - t(1); %step sizes

[X1,X2] = meshgrid(x,x); %2D interpolate error field

solution = zeros( 2 , N_x , N_x , N_t ); %store vector components of velocity solution in space and time
pressure_solution = zeros( N_x , N_x , N_t ); %store components of pressure solution in space and time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%build stuff for finite element problem

n_el = N_x-1; %number of elements per side
num_el = n_el^2; %total number of elements
n_dof = N_x^2; %number of nodal dofs

%compute and store nodal coordinates
Coord = zeros(N_x^2,2);
for i=1:N_x^2
    x1_idx = mod(i-1,N_x)+1;
    x2_idx = idivide( i-1,int16(N_x) )+1;
    Coord(i,1) = x(x1_idx);
    Coord(i,2) = x(x2_idx);
end

%assign DOFs to nodes
Dof = zeros(N_x^2,1);
for i=1:N_x^2
   Dof(i) = i;
end

%create elements in mesh
Edof = zeros(num_el,5);
count = 0;
for e=1:num_el
    Edof(e,1) = e;

    first = e + count;
    second = first + 1;
    third = first + N_x + 1;
    fourth = first + N_x;

    Edof(e,2:5) = [first,second,third,fourth];

    if mod(e,n_el) == 0
        count = count + 1;
    end
end

%elemental coordinates
[ex,ey] = coordxtr(Edof,Coord,Dof,4);

ep = [ 1 , 3 ]; %element thickness and quadrature points
D = [ 1 , 0 ; 0 , 1 ]; %constitutive relation

%precompute stiffness matrix
K = zeros(n_dof,n_dof);
for e=1:num_el
    edof = Edof(e,:);
    exe = ex(e,:); eye = ey(e,:);

    Ke = flw2i4e( exe , eye , ep , D );
    K = assem( edof , K , Ke );  
end

%boundary condition -- node 1 has zero displacement
bc = [ 1 , 0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time stepping
for i=1:(N_t-1)

    if 2*i == N_t
        disp('halfway!')
    end

    %rename each component solution at current time step for readability
    u1_t = solution(1,:,:,i); u1_t = squeeze(u1_t);
    u2_t = solution(2,:,:,i); u2_t = squeeze(u2_t);

    %compute two components of intermediate velocity with finite difference
    %only interior nodes are updated, boundary comes from next time step
    intermediate = zeros( 2 , N_x , N_x );

    convection1 = convection( u1_t , u2_t , 1 , N_x , dx );
    laplacian1 = laplacian( u1_t , N_x , dx );
    intermediate( 1 , 2:N_x-1 , 2:N_x-1 ) = dt * ( -convection1 + nu*laplacian1 ) + u1_t( 2:N_x-1 , 2:N_x-1 );

    convection2 = convection( u1_t , u2_t , 2 , N_x , dx );
    laplacian2 = laplacian( u2_t , N_x , dx );
    intermediate( 2 , 2:N_x-1 , 2:N_x-1 ) = dt * ( -convection2 + nu*laplacian2 ) + u2_t( 2:N_x-1 , 2:N_x-1 );

    %apply boundary condition for next time step to intermediate velocity
    intermediate( 1 , 1 , : ) = 1*boundary( 1 , x , t(i+1) , L );
    intermediate( 2 , : , 1 ) = 1*boundary( 2 , x , t(i+1) , L );
    intermediate( 1 , N_x , : ) = 1*boundary( 3 , x , t(i+1) , L );

    %use intermediate velocity to compute error from divergence free velocity field
    %how to deal with divergence at boundaries? set to zero?
    error = grad( squeeze( intermediate(1,:,:) ) , 1 , N_x , dx ) + grad( squeeze( intermediate(2,:,:) ) , 2 , N_x , dx );

    rhs = @(xx1,xx2) ( interp2( X1 , X2 , flipud( -(rho/dt) * error ) , xx1 , xx2 ) );

    %build for vector from interpolated rhs 
    F = zeros(n_dof,1);

    for e=1:num_el
        edof = Edof(e,:);
        exe = ex(e,:); eye = ey(e,:);

        %compute average source term over the element
        eq = 0;
        for ii=1:4
            eq = eq + 0.25*rhs( exe(ii) , eye(ii) );
        end

        [~,fe] = flw2i4e( exe , eye , ep , D , eq );
        F = insert( edof , F , fe );
    end

    %solve linear system for pressure dofs
    p = solveq( K , F , bc );

    %build pressure solution on grid
    pressure = zeros(N_x,N_x);

    %convert finite element solution to array
    for kk=1:N_x^2
        x2_idx = mod(kk-1,N_x)+1;
        x1_idx = idivide( kk-1,int16(N_x) )+1;

        pressure(x1_idx,x2_idx) = p(kk);
    end

    %store pressure array according to fd grid
    pressure = flipud( pressure );

    %update intermediate solution with pressure correction
    solution(1,:,:,i+1) = -( dt/rho ) * grad(pressure,1,N_x,dx) + squeeze( intermediate(1,:,:) );
    solution(2,:,:,i+1) = -( dt/rho ) * grad(pressure,2,N_x,dx) + squeeze( intermediate(2,:,:) );
    pressure_solution(:,:,i+1) = pressure;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plots and post-processing

iter = 1000; %number of iterations to plot
sf = 0.25; %scale factor for quiver plot

step = double( int16( N_t/iter ) ) - 1; %stepsize in solution array

figure(1)
for t=1:iter
    u1_t = flipud( squeeze( solution(1,:,:,step*t) ) );
    u2_t = flipud( squeeze( solution(2,:,:,step*t) ) );
    mag = ( u1_t.^2 + u2_t.^2 ).^0.5;
    contour( X1 , X2 , mag , 10 );
    hold on;
    quiver( X1 , X2 , sf*u1_t , sf*u2_t , 'AutoScale','off')
    hold off;
    title('Velocity Field')
    xlabel('x')
    ylabel('y')
    grid on
    pause(.01)
end

pmin = min(min(min(pressure_solution))); pmax = max(max(max(pressure_solution)));
figure(2)
for t=1:iter
    pp = flipud( squeeze( pressure_solution(:,:,step*t) ) );
    surf( X1 , X2 , pp )
    axis([ 0 L 0 L pmin pmax ])
    title('Pressure Field')
    %axis off
    xlabel('x')
    ylabel('y')
    grid on
    pause(.01)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%local functions

%function specifying applied velocity at boundary
function val = boundary( which , x , t , L )

    if which == 1
        val = 0.6*sin(2*pi*t) + 0.8*sin(pi*t);
    elseif which == 2
        val = 0.6*sin(2*pi*t) - 0.7*sin(pi*t);
    elseif which == 3
        val = sin(0.3*pi*t) + 0.3*sin(2*pi*t);
    end

    time_env = 0.4 + 0.6*exp( -5*(t-2).^2) + 0.6*exp( -5*(t-6).^2 );
    val = val * time_env * sin(pi*x/L);
   
end

%computes derivative of field on grid in given direction
%assumes dirichlet boundary conditions
function val = grad( field, idx , N_x , dx )
    
    if idx == 1
        fd = ( 1/(2*dx) ) * ( field( 2:N_x-1 , 3:N_x ) - field( 2:N_x-1 , 1:N_x-2 ) );
    elseif idx == 2
        fd = ( 1/(2*dx) ) * ( field( 1:N_x-2 , 2:N_x-1 ) - field( 3:N_x , 2:N_x-1 ) );
    end

    val = zeros(N_x,N_x);
    val( 2:N_x-1 , 2:N_x-1 ) = fd;

end

%computes one component of convective term
%assumes dirichlet boundary conditions
function val = convection( u1 , u2 , idx , N_x , dx )
    
    if idx == 1
        u = u1;
    elseif idx == 2
        u = u2;
    end
   
    x1_comp = ( 1/(2*dx) ) * ( u( 2:N_x-1 , 3:N_x ) - u( 2:N_x-1 , 1:N_x-2 ) );
    x2_comp = ( 1/(2*dx) ) * ( u( 1:N_x-2 , 2:N_x-1 ) - u( 3:N_x , 2:N_x-1 ) );

    update = u1( 2:N_x-1 , 2:N_x-1 ) .* x1_comp + u2( 2:N_x-1 , 2:N_x-1 ) .* x2_comp;

    val = update;
end

%computes laplacian
%assumes dirichlet boundary conditions
function val = laplacian( u , N_x , dx )
    
    x1_comp = (1/dx^2) * ( u( 2:N_x-1 , 3:N_x ) -2*u( 2:N_x-1 , 2:N_x-1 ) + u( 2:N_x-1 , 1:N_x-2 ) );
    x2_comp = (1/dx^2) * ( u( 1:N_x-2 , 2:N_x-1 ) -2*u( 2:N_x-1 , 2:N_x-1 ) + u( 3:N_x , 2:N_x-1 ) );
   
    update = x1_comp + x2_comp;

    val = update;
end

%gradient with boundaries computed
function val = p_grad( field, idx , N_x , dx )
    
    if idx == 1
        fd = ( 1/(2*dx) ) * ( field( 2:N_x-1 , 3:N_x ) - field( 2:N_x-1 , 1:N_x-2 ) );
    elseif idx == 2
        fd = ( 1/(2*dx) ) * ( field( 1:N_x-2 , 2:N_x-1 ) - field( 3:N_x , 2:N_x-1 ) );
    end

    val = zeros(N_x,N_x);
    val( 2:N_x-1 , 2:N_x-1 ) = fd;

    %edge derivatives
    val( 1 , : ) = (1/dx) * ( field( 1 , N_x-2 ) - field( 2 , N_x-2 ) );
    val( N_x , : ) = (1/dx) * ( field( N_x-1 , N_x-2 ) - field( N_x , N_x-2 ) );
    val( : , 1 ) = (1/dx) * ( field( N_x-2 , 2 ) - field( N_x-2 , 1 ) );
    val( : , N_x ) = (1/dx) * ( field( N_x-2 , N_x ) - field( N_x-1 ,N_x-1 ) );

end