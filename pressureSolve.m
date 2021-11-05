clear;
clc
close all
tic

% Constant
rho = 1;
nu = 0.001;
U_in = 60;

% Domain/Discretisation
n = 5;
x = 8;
hx = x/2^n; nhx = x/hx + 2;
hy = 0.1/2^(n); nhy = 0.1/hy + 2;
dt = 0.00001;

[XX,YY] = meshgrid(0-hx/2:hx:4+hx/2,0-hy/2:hy:0.1+hy/2);

% Initial condition
U = zeros(nhy,nhx); % x-velocity
U(2:end-1,1:2) = U_in;
Unew = U;

V = zeros(nhy,nhx); % y-velocity
Vnew = V;

P = zeros(nhy,nhx); % Pressure
P(:,1) = 1;
trial = P;

% Residual initialisation
residual = zeros(nhy,nhx);

% Divergence initialisation
div = residual;

% Turbulence model initialisation
delY_u = residual;
delYY_u = delY_u;
delYX_u = delY_u;

% Iteration counter
n_count = 0;
w = 1.9;

% Error/Steady state criteria
U_change = 1;
U_change_max = 1e-4;  % Tested with residual independence
resid_pc = 1;
resid_pc_max = 1e-4;

% Mixing Length
distFromWall = YY(1:nhy/2,1);
distFromWall = [distFromWall;flip(distFromWall)];
del = 0.1/2;
l_mix = del*(0.14-0.08*(1-distFromWall/del).^2-0.06*(1-distFromWall/del).^4);



% For residual animation 
resid_show = [];


%% Solve
% while  U_change > U_change_max
for opop = 1:1
    n_count=n_count+1;
    
    for i = 2:nhy-1
        for j = 2:nhx - 1
            delY_u(i,j) = (Unew(i+1,j)-Unew(i-1,j))/(2*hy);
        end
    end
    
    delY_u(nhy,:) = 0;
    delY_u(1,2:nhx-1) = 0;
    delY_u(2:nhy,1) = -delY_u(2:nhy,2);
    delY_u(2:nhy,nhx) = -delY_u(2:nhy,nhx-1);
    
    l_mixM = repmat(l_mix,1,18);
    temp = l_mix.^2.*abs(delY_u).*delY_u;
    
    for i = 2:nhy - 1
        for j = 2:nhx - 1
            delYY_u(i,j) = (temp(i+1,j)-temp(i-1,j))/(2*hy);
        end
    end
  
    for i = 2:nhy-1
        for j = 2:nhx-1
            delYX_u(i,j) = (temp(i,j+1)-temp(i,j-1))/(2*hx);
        end
    end
    
    % U velocity 
    for i=2:nhy-1
        for j=2:nhx-1
            Unew(i,j)= (1-w)*U(i,j)+w*(U(i,j)-dt*(P(i+1,j)-P(i-1,j))/(2*hx)...
                +nu*dt*(1/(hx*hx)*(U(i+1,j)-2.*U(i,j)+U(i-1,j))...
                +1/(hy*hy)*(U(i,j+1)-2.*U(i,j)+U(i,j-1)))...
                -dt*U(i,j)/(hx)*(U(i,j)-U(i-1,j))...
                -dt*V(i,j)/(2*hy)*(U(i,j+1)-U(i,j-1))+dt*delYY_u(i,j));
        end
    end
    
    % Update Boundary Condition 
    Unew(nhy,:)=Unew(nhy-1,:);
    Unew(1,2:nhx-1)=-Unew(2,2:nhx-1)+2*U_in;
    Unew(2:nhy,1)=-Unew(2:nhy,2);
    Unew(2:nhy,nhx)=-Unew(2:nhy,nhx-1);
    U_change=max(max(abs((Unew-U)/dt)));
    
    % V velocity 
    for i=2:nhy-1
        for j=2:nhx-1
            Vnew(i,j)= (1-w)*V(i,j)+w*(V(i,j)-dt*(P(i,j+1)-P(i,j-1))/(2*hy)...
                +nu*dt*(1/(hx*hx)*(V(i+1,j)-2.*V(i,j)+V(i-1,j))...
                +1/(hy*hy)*(V(i,j+1)-2.*V(i,j)+V(i,j-1)))...
                -dt*U(i,j)/(hx)*(V(i,j)-V(i-1,j))...
                -dt*V(i,j)/(2*hy)*(V(i,j+1)-V(i,j-1))+dt*delYX_u(i,j));
        end
    end
    
    % Update Boundary Condition  
    Vnew(nhy,:)=Vnew(nhy-1,:);
    Vnew(1,:)=-Vnew(2,:);
    Vnew(:,1)=-Vnew(:,2);
    Vnew(:,nhx)=-Vnew(:,nhx-1);
    
    % Divergent term for continuity
    for i=2:nhy-1
        for j=2:nhx-1
            div(i,j)=(Unew(i+1,j)-Unew(i-1,j))/(2*hx)...
                +(Vnew(i,j+1)-Vnew(i,j-1))/(2*hy);
        end
    end
    
    div_sum=sum(sum(abs(div)))/((nhx-2)*(nhy-2));
    
    resid_pc=1;
    p_count=0;
    trial(:,:) = 0;
    trial_lag = trial;
    prs_coeff = -2/hx^2 - 2/hy^2;
    prs_coeff = 1/prs_coeff;
    
    relx_pc = 1.6;
    
    while resid_pc > resid_pc_max & p_count < 100
        p_count = p_count + 1;
        for i = 2:nhy - 1
            for j = 2:nhx - 1
                trial(i,j) = trial_lag(i,j)*(1-relx_pc)+relx_pc*0.25*(div(i,j)/dt - (trial_lag(i-1,j) + trial_lag(i+1,j))/hx^2 - ...
                        (trial_lag(i,j-1)+trial_lag(i,j+1))/hy^2)*prs_coeff;
            end
        end
        
        resid_pc=abs(sum(sum(abs(trial)-abs(trial_lag))))/((nhx-2)*(nhy-2));
        
        trial_lag = trial;

    end
    
    % Update Pressure 
    Pc_mine =  (1/(-2/(hx*hx)-2/(hy*hy)))*trial+trial;
    
    Pc_mine(1,:)=Pc_mine(2,:);
    Pc_mine(:,1)=Pc_mine(:,2);
    Pc_mine(:,nhx)=Pc_mine(:,nhx-1);
    Pc_mine(nhy,:)=0;
    P=P+Pc_mine;
    
    % New corrected velocities
    for i=2:nhy-1
        for j=2:nhx-1
            Unew(i,j)= Unew(i,j)-dt*(Pc_mine(i+1,j)-Pc_mine(i-1,j))/(2*hx);
            Vnew(i,j)= Vnew(i,j)-dt*(Pc_mine(i,j+1)-Pc_mine(i,j-1))/(2*hy);
        end
    end
    
    % Velocity boundary condition 
    Unew(nhy,:)=Unew(nhy-1,:);
    Unew(1,2:nhx-1)=-Unew(2,2:nhx-1)+2*U_in;
    Unew(2:nhy,1)=-Unew(2:nhy,2);
    Unew(2:nhy,nhx)=-Unew(2:nhy,nhx-1);
    
    Vnew(nhy,:)=Vnew(nhy-1,:);
    Vnew(1,:)=-Vnew(2,:);
    Vnew(:,1)=-Vnew(:,2);
    Vnew(:,nhx)=-Vnew(:,nhx-1);
    
    U=Unew;
    V=Vnew;
    
    
    resid_show(n_count) = U_change;
% 
%      semilogy(resid_show)
%      pause(0.01);
disp(U_change);
end

opop = 1;
% toc = performTime;



%% Sparse matrix construction 
function x = getSparse(nx,hx,hy)

ny = nx;

row = nx - 2;
col = ny - 2;

mainCoeff = 2/hx^2 + 2/hy^2;
offDiag1 = -1/hx^2;
offDiag2 = -1/hy^2;

J = length(1:row*col);
D = sparse(1:J,1:J,mainCoeff*ones(1,J),J,J);

offDiag = ones(1,J-1);
offDiag(1:row) = 0;
offDiag(row:row:end-row) = 0;
offDiag(end-row:end) = 0;

E = sparse(2:J,1:J-1,offDiag1*offDiag,J,J);
F = sparse(row+1:J,1:J-row,offDiag2*ones(1,J-row),J,J);

x = E + D + E' + F + F';

end






