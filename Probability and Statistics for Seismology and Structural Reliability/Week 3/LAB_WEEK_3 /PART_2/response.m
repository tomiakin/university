function SA=response(f,dt,T,xi)

% [SA]=RESPONSE(f,dt)
% Purpose: Computes elastic response spectra for SDOF systems with periods T
% and damping ratios xi subjected to a support acceleration f(t) sampled at
% a step dt. Zero initial conditions are assumed.
% ydd(t) + 2 xi w yd(t) + wA2 y(t) = -f(t)          y(0)=0, yd(0)=0

% Input:
%           f = excitation vector
%           dt = sample time in f;

% Output:
%           SA = Absolute-acceleration


% T=0.01:0.01:4;
% xi=0.05;

%         SD=zeros(length(T),length(xi));
%         PSV=zeros(length(T),length(xi));
%         PSA=zeros(length(T),length(xi));
%         SV=zeros(length(T),length(xi));
SA=zeros(length(T),length(xi));
%         ED=zeros(length(T),length(xi));

for j=1:length(xi)
    for i=1:length(T)
        
        if T(i)==0
            SA(i,j)=max(abs(f));
        else        
            w=2*pi/T(i); C=2*xi(j)*w; K=w*w; y(:,1)=[0; 0];
            A=[0 1; -K -C]; Ae=expm(A*dt); AeB=A\(Ae-eye(2))*[0; -1];
            for k=2:length(f)
                y(:,k)=Ae*y(:,k-1)+AeB*f(k);
            end

            SA(i,j)=max(abs(K*y(1, :)+C*y(2, :))) ;
        end

    end
end