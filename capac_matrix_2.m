% --------------------------------------------------------------
% Compute capacitance per unit length of
% a coaxial pair of rectangles
% --------------------------------------------------------------
function cap = capac_matrix_2(a, b, c, d, n,plot)

%%%%%%%%%%%%%%%%%%%%%%%%%% Copied from capactior.m
% Arguments:
%    a   =  width of inner conductor
%    b   = height of inner conductor
%    c   =  width of outer conductor
%    d   = height of outer conductor
%    n   = number of points in the x-direction (horizontal)
%          (optimum is 2-c/n, where c is about pi)
%   plot =  Set to true to see surface plot
% Returns:
%    cap = capacitance per unit length [pF/m]



% Make grids such that point land on inner and outer conductor
h  = 0.5*c/n;                % Grid size
r = ceil(0.5*a/h);
h = 0.5*a/r;            %assuming a = b
na = round(0.5*a/h);         % Number of segments on 'a'
qt = mod(c,h);
if qt ~= 0
    qm = abs(qt-h);
    if qm < qt
        c = c+qm;
        d = d+qm;
    else
        c = c-qt;
        d = d-qt;
    end
end
x  = linspace(0,0.5*c,n+1);  % Grid points along x-axis
n;
m = n;
% m  = round(0.5*d/h);       % Number of segments on 'd'
h;
mb = round(0.5*b/h);         % Number of segments on 'b'
y  = linspace(0,0.5*d,m+1);  % Grid points along y-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize potential and mask array
f = zeros(n+1,m+1);          % 2D-array with solution

A = zeros((n+1)*(m+1));
fnum = @(i,j) (j-1)*(n+1)+i; %Funciton telling us what gridpoint we are on

% Make the potential on the inner conductor one
for i = 1:na+1
    for j = 1:mb+1
        num = fnum(i,j);
        mask(i,j) = 0;
        f(i,j)    = 1;
        A(num,num) = 1;
    end
end
fvec = f(:);
% Fill Matrix A for potential that have no BC

for i = 2:n
    for j = 2:n
        num = fnum(i,j);
        if A(num,num) ~= 1 % Ignore the inner conductor
            A(num,num) = -4; % Set up the poisson eq
            A(num,num-1) = 1;
            A(num,num+1) = 1;
            A(num,num+n+1) = 1;
            A(num,num-n-1) = 1;
        end
    end
end

% Look at Vertical BC
i = 1;
for j = (mb+2):m
    num = (j-1)*(n+1)+i;
    A(num,num) = -4;
    %     A(num,num-1) = mask(num-1)/4;
    A(num,num+1) = 2;
    A(num,num+n+1) = 1;
    A(num,num-n-1) = 1;
end

%Look at horizontal BC
j = 1;
for i = (na+2):n
    num = (j-1)*(n+1)+i;
    A(num,num) = -4;
    A(num,num-1) = 1; % Using the Mask like this isn't correct
    A(num,num+1) = 1;
    A(num,num+n+1) = 2;
    %     A(num,num-n-1) = mask(num-n-1)/4;
end

% The outer conductor
i = n+1;
for j = 1:m+1
    num = (j-1)*(n+1)+i;
    A(num,num) = 1;
end

%The outer conductor
j = m+1;
for i = 1:n+1
    num = (j-1)*(n+1)+i;
    A(num,num) = 1;
end


A = sparse(A); %Sparsing the matrix

% Divide by grid spacing
A = A./h^2;
fvec = fvec./h^2;

fvec = A\fvec; % Solve the system of equations

% Convert potential to matrix form to plot and solve
ff = fvec;
ff = reshape(ff,[n+1,m+1]);

if plot == true
    surf(ff)
    title('\Phi')
    xlabel('x')
    ylabel('y')
end
cap = gauss(n,m,h,ff);

function cap = gauss(n,m,h,f)

% Arguments:
%    n    = number of points in the x-direction (horizontal)
%    m    = number of points in the y-direction (vertical)
%    h    = cell size
%    f    = 2D-array with solution
% Returns:
%    cap = capacitance per unit length [pF/m]

q = 0;

for i = 1:n
    q = q + (f(i,m)+f(i+1,m))*0.5; %integrate along upper boundary
end

for j = 1:m
    q = q + (f(n,j)+f(n,j+1))*0.5; %integrate along right boundary
end

cap = q*4;           % 4 quadrants
cap = cap*8.854187;  % epsilon0*1e12 gives answer in pF/m
