% General parameters.
N = 100; R0 = 1e4; rho0 = 1e7; dt = 0.01;

% Parameter percentage.
pc = 10;

% Select pc = 0 for the simulation to run until a 1% of the initial radius
% is reached.

% Select some (positive) pc other than zero for the simulation to run
% until differences of pc% appear between the analytical and the numerical
% solution. Example: pc = 10 means that the simulation will run until
% discrepancies of 10% show up.

tf = zeros(5,1); Rfnum = zeros(5,1); Rfanal = zeros(5,1);
% Try different values of the integration parameter beta.
for ibeta = 0:4
    beta = 0.25*ibeta;
    [tvec,Rnum,Ranal,trho,rhomat] = funPW3(N,R0,rho0,beta,dt,pc);
    
    figure(ibeta+1)
    plot(tvec,Rnum,'o',tvec,Ranal); xlim([0 0.7]); ylim([0, R0]); grid
    xlabel('t [s]'); ylabel('R [cm]'); legend('Numerical','Analytical')
    title(['Radius of the star (\beta = ' num2str(beta) ')'])
    figure(ibeta+6)
    semilogy(1:N,rhomat)
    grid; xlabel('i'); ylabel('\rho_{i+1/2}');
    title(['Density profile at different times (\beta = '...
        num2str(beta) ')'])
    tf(ibeta+1) = tvec(end);
    Rfnum(ibeta+1) = Rnum(end);
    Rfanal(ibeta+1) = Ranal(end);
end

% Function which performs the free fall collapse simulation.
function [tvec,Rnum,Ranal,trho,rhomat] = funPW3(N,R0,rho0,beta,dt,pc)

% Constants.
G  = 6.67259*10^-8;       % Gravitational constant in [cm^3 g^-1 s^-2].
M = 4*pi/3*R0^3*rho0;     % Total mass in [g].
dm = M/N;                 % Mass difference between consecutive shells.
mvec = dm*(0:N);          % LENGTH N+1. Thus, m_1 = 0 and m_{N+1} = M.

% Initial conditions.
rhovec = rho0*ones(1,N);
uvec = zeros(1,N+1);
rvec = zeros(1,N+1);
for i=1:N
    rvec(i+1) = (3*dm/(4*pi*rhovec(i))+rvec(i)^3)^(1/3);
end

tvec(1) = 0;
Rnum(1) = R0;
Ranal(1) = R0;

% Time loop.
time = 0; it = 0; saverho = 0; delta = zeros(3,N);
while dt>1e-6
    time = time+dt;
    it = it+1;
    
    uvecn = uvec;  % So as to update uvec without modifying the previous time step velocities.
    rvecn = rvec;  % So as to update rvec without modifying the previous time step radii.
    rhovecn = rhovec;
    
    % While loop (epsilon) and update (x+dx).
    epsilonrho = 1; epsilonu = 1; epsilonr = 1;
    while epsilonrho > 1e-6 || epsilonu > 1e-6 || epsilonr > 1e-6
        
        % Internal shell
        C1 = 1/rhovec(1)-4/3*pi*rvec(2)^3/dm;
        C2 = (uvec(2)-uvecn(2))/dt+(1-beta)*G*mvec(2)/rvecn(2)^2+...
            beta*G*mvec(2)/rvec(2)^2;
        C3 = (rvec(2)-rvecn(2))/dt-(1-beta)*uvecn(2)-beta*uvec(2);
        A = [-1/rhovec(1)^2   0      -4*pi*rvec(2)^2/dm;...
             0                1/dt   -2*beta*G*mvec(2)/rvec(2)^3;...
             0                -beta  1/dt                           ];
        b = [-C1; -C2; -C3];
        dxvec = A\b;
        delta(1,1)=dxvec(1); delta(2,1)=dxvec(2); delta(3,1)=dxvec(3);
        
        % Intermediate shells and external shell
        for j=2:N
            F1 = 1/rhovec(j)-4/3*pi*(rvec(j+1)^3-rvec(j)^3)/dm;
            F2 = (uvec(j+1)-uvecn(j+1))/dt+(1-beta)*G*...
                mvec(j+1)/rvecn(j+1)^2+beta*G*mvec(j+1)/rvec(j+1)^2;
            F3 = (rvec(j+1)-rvecn(j+1))/dt-(1-beta)*uvecn(j+1)-...
                beta*uvec(j+1);
            A = [-1/rhovec(j)^2 0     -4*pi*rvec(j+1)^2/dm;...
                 0              1/dt  -2*beta*G*mvec(j+1)/rvec(j+1)^3;...
                 0              -beta 1/dt                           ];
            b = [-F1-4*pi/dm*rvec(j)^2*delta(3,j-1); -F2; -F3];
            dxvec = A\b;
        delta(1,j)=dxvec(1); delta(2,j)=dxvec(2); delta(3,j)=dxvec(3);
        end
        
        % UPDATE VARIABLES
        rhovec = rhovec+delta(1,:);
        uvec = [0 uvec(1,2:end)+delta(2,:)];
        rvec = [0 rvec(1,2:end)+delta(3,:)];
        
        % Check error
        epsilonrho = max(abs(delta(1,:)./rhovec));
        epsilonu = max(abs(delta(2,:)./uvec(1,2:end)));
        epsilonr = max(abs(delta(3,:)./rvec(1,2:end)));
    end
    
    % Calculate analytical solution.
    x0 = Ranal(end);
    Rana = newton(x0,time,1e-10,100);
    Rana = real(Rana);
    
    % Conditions for getting as close as possible to 1% of the initial
    % radius.
    if pc == 0 && rvec(end)<0.01*R0
        it = it-1;
        time = time-dt;
        dt = dt/2;
        rhovec = rhovecn;
        uvec = uvecn;
        rvec = rvecn;
    elseif pc == 0 && rvec(end)>0.01*R0
        % Save time.
        tvec(it+1) = time;
        
        % Save numerical solution.
        Rnum(it+1) = rvec(end);
        
        % Save analytical solution.
        Ranal(it+1) = Rana;
        
        % Save density profile once every ten time steps.
        if mod(it,10) == 0
            saverho = saverho+1;
            rhomat(saverho,:) = rhovec;
            trho(saverho) = time;
        end
    end
    
    % Conditions for getting as close as possible to a certain error pc.
    if pc ~= 0 && abs(rvec(end)-Rana)/Rana>pc/100
        it = it-1;
        time = time-dt;
        dt = dt/2;
        rhovec = rhovecn;
        uvec = uvecn;
        rvec = rvecn;
    elseif pc ~= 0 && abs(rvec(end)-Rana)/Rana<pc/100
        % Save time.
        tvec(it+1) = time;
        
        % Save numerical solution.
        Rnum(it+1) = rvec(end);
        
        % Save analytical solution.
        Ranal(it+1) = Rana;
        
        % Save density profile once every ten time steps.
        if mod(it,10) == 0
            saverho = saverho+1;
            rhomat(saverho,:) = rhovec;
            trho(saverho) = time;
        end
    end
end
end

function R = newton(x0,t,tol,kmax)
k=0; tolk=1; X=x0;
while tolk>tol && k<kmax
    m=length(x0);
    I=eye(m);h=sqrt(eps);
    for j=1:m
        f1=feval(@analytic,x0-I(:,j)*h,t);
        f2=feval(@analytic,x0+I(:,j)*h,t);
        DF(:,j)=(f2-f1)/(2*h);
    end
    Fk=feval(@analytic,X(:,k+1),t);
    deltaxk=DF\(-Fk);
    xk=X(:,k+1)+deltaxk;
    X=[X xk];
    tolk=max(abs(X(:,k+1)-X(:,k+2)));
    k=k+1;
end
R = X(:,end);
end

function F = analytic(Ranalyticalpc,t)
R0 = 1e4; rho0 = 1e7; G  = 6.67259*10^-8;
F = (sqrt(1-Ranalyticalpc/R0).*sqrt(Ranalyticalpc/R0)+...
    asin(sqrt(1-Ranalyticalpc/R0)))/sqrt(8*pi*G*rho0/3)-t;
end