
% Zhao_4 SHPSG function

% VARIABLES TO MANIPULATE %%%%%%%%%%%%%%%
EI = 1;
FI = 1;
D2_8 = 0.4;
D9_15 = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff = SHPSG(EI,FI,D2_8,D9_15);

theta = linspace(0,pi,50);
phi = linspace(0,2*pi,50);
[theta,phi] = meshgrid(theta,phi);

figure;
view(3);
[x,y,z] = sph2car(coeff, phi, theta);
surf(x,y,z)
colormap('gray');
camlight;
axis equal;

% randomly generate SH coefficients given parameters
function [fvec] = SHPSG(EI,FI,D2_8,D9_15)

    Fvec = zeros(4,3);

    % determine C0 and C1 with a unit max principal dimension
    arr = [0 0 0;
        -1 1j 0;
        0 0 sqrt(2);
        1 1j 0];

    fvec_sphere = -sqrt(pi/6)*arr; % a sphere with a unit diameter
    Fvec(:,1) = fvec_sphere(:,1);
    Fvec(:,2) = EI*fvec_sphere(:,2);
    Fvec(:,3) = (FI*EI)*fvec_sphere(:,3);

    d1 = norm(Fvec,'fro');

    % Determine D2 and D9
    alpha = 1.387;
    beta = 1.426;

    D_2 = D2_8/(1+(2/3)^alpha+(2/4)^alpha+(2/5)^alpha+(2/6)^alpha+(2/7)^alpha+(2/8)^alpha)*d1;
    D_9 = D9_15/(1+(9/10)^beta+(9/11)^beta+(9/12)^beta+(9/13)^beta+(9/14)^beta+(9/15)^beta)*d1;

    % Determine D3-D8 and D10-D15
    % assume all descriptors have three identical decomposition at x, y, and z-axis
    I = zeros(16,3);
    I(2,:) = [1 1 1];
    I(3,:) = [D_2 D_2 D_2];

    for c = 3:8
        dn = D_2*((c-1)/2)^(-alpha)/sqrt(3);
        I(c+1,:) = [dn dn dn];
    end

    for c = 9:14
        dn = D_9*((c-1)/9)^(-beta)/sqrt(3); % or should it be c-1
        I(c+1,:) = [dn dn dn];
    end

    L = zeros(16^2, 3);
    N = zeros(16^2, 3);


    % Randomly generate P
    % Randomly generate P including C1'-C15' with c_n^(-m)=(-1)^m*c_n^m*

    for n = 1:15
        J = ones(n+1,3)-2*rand(n+1,3); % [-1,1]
        K = flipud(J(1:n,:));
        i = (-1)^n; A = [i i i];
        B = K.*A;
        L(n^2+1:(n+1)^2,:) = [J; B];
        M = ones(n,3)-2*rand(n,3);
        M1 = [M; 0 0 0];
        N(n^2+1:(n+1)^2,:) = [M1; flipud((-1)^(n+1)*M)];
    end

    P = L+N*1j;
    P(1,:) = ones(1,3)-2*rand(1,3);

    % Calculate D1-D16 with the SH coefficients of P
    Q = conj(P);
    R = zeros(16,3);

    for n = 0:15
        a = n^2+1;
        b = (n+1)^2; % matlab end is inclusive but python is not
        R(n+1,:) = sqrt(real(sum(P(a:b,:).*Q(a:b,:))));
    end

    % Determine C2-C15 by making the descriptors of P equal to D2-D15
    fvec = zeros(16^2,3);
    fvec(1:4,:) = Fvec;

    for d = 3:16
        a = (d-1)^2+1;
        b = d^2;

        for n = 1:3
            fvec(a:b,n) = P(a:b,n)/R(d,n)*I(d,n);
        end
    end
end

% convert spherical to cartesian coordinates with SH expansion
function [x,y,z] = sph2car(coeff, phi, theta)
    
    x = 0; y = 0; z = 0;
    index = 1;

    % compute spherical harmonic of degree N and order M
    for n = 0:8
        order = -n:n;

        for m = 0:2*n
            Y = harmonicY(n,order(m+1),theta,phi);
            x = x + coeff(index,1)*Y;
            y = y + coeff(index,2)*Y;
            z = z + coeff(index,3)*Y;
            index = index + 1;
        end
    end

    x = real(x);
    y = real(y);
    z = real(z);
end
