%% User-provided partial derivatives of right-hand side |f|

function J=circadian_deri(xx,par,nx,np,v)

%% Parameter vector
tauB = par(1);
nB = par(2);
KeB = par(3);
KiB = par(4);
KdBP = par(5);
KdB = par(6);
nP = par(7);
KeP = par(8);
KaP = par(9);
KdP = par(10);
Ke2Tot = par(11);
coupleRatio = par(12);
B = xx(1,1);
P = xx(2,1);
BLag = xx(1,2);
J = [];

if length(nx)==1 && isempty(np) && isempty(v)
    %% First-order derivatives wrt to state nx+1
    if nx==0 % derivative wrt x(t)
        J(:,1) = [- KdB - KdBP*P; 0];
        J(:,2) = [-B*KdBP; -KdP];
    elseif nx==1 % derivative wrt x(t-tau1)
        J(:,1) = [-(KeB*nB*(BLag/KiB)^(nB - 1))/(KiB*((BLag/KiB)^nB + 1)^2);...
                  (KaP*KeP*nP*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^2)];
        J(:,2) = [0; 0]; % no PLags in system
    end
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% First order derivatives wrt parameters
    if np==1 % derivative wrt tauB
        J = [0; 0];
    elseif np==2 % derivative wrt nB
        J = [-(KeB*log(BLag/KiB)*(BLag/KiB)^nB)/((BLag/KiB)^nB + 1)^2; 0];
    elseif np==3 % derivative wrt KeB
        J = [1/((BLag/KiB)^nB + 1); 0];
    elseif np==4 % derivative wrt KiB
        J = [(BLag*KeB*nB*(BLag/KiB)^(nB - 1))/(KiB^2*((BLag/KiB)^nB + 1)^2); 0];
    elseif np==5 % derivative wrt KdBP
        J = [-B*P; 0];
    elseif np==6 % derivative wrt KdB
        J = [-B; 0];
    elseif np==7 % derivative wrt nP
        J = [0; -(KeP*log(KaP/BLag)*(KaP/BLag)^nP)/((KaP/BLag)^nP + 1)^2];
    elseif np==8 % derivative wrt KeP
        J = [0; 1/((KaP/BLag)^nP + 1)];
    elseif np==9 % derivative wrt KaP
        J = [0; -(KeP*nP*(KaP/BLag)^(nP - 1))/(BLag*((KaP/BLag)^nP + 1)^2)];
    elseif np==10 % derivative wrt KdP
        J = [0; -P];
    elseif np==11 % derivative wrt KeB2
        J = [(1-coupleRatio); coupleRatio];
    elseif np==12 % derivative wrt KeP2
        J = [-Ke2Tot; Ke2Tot];
    end
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% Mixed state, parameter derivatives
    if nx == 0 % derivative wrt x(t)
        if np == 5 % derivative wrt KdBP
            J(:,1) = [-P; 0];
            J(:,2) = [-B; 0];
        elseif np == 6 % der wrt KdB
            J(:,1) = [-1; 0];
            J(:,2) = [0; 0];
        elseif np == 10 % der wrt KdP
            J(:,1) = [0; 0];
            J(:,2) = [0; -1];
        else
            J = zeros(2,2);
        end
    elseif nx==1 % derivative wrt x(t-tau1)
        J(:,2) = [0; 0]; % no PLags
        if any(np == [1,5,6,10,11,12])% derivative wrt tauB, KdBP, KdB, KdP, KeB2, or KeP2
            J(:,1) = [0; 0];
        elseif np==2 % derivative wrt nB
            J(:,1) = [(2*KeB*nB*log(BLag/KiB)*(BLag/KiB)^nB*(BLag/KiB)^(nB - 1))/(KiB*((BLag/KiB)^nB + 1)^3)...
                - (KeB*nB*log(BLag/KiB)*(BLag/KiB)^(nB - 1))/(KiB*((BLag/KiB)^nB + 1)^2) -...
                (KeB*(BLag/KiB)^(nB - 1))/(KiB*((BLag/KiB)^nB + 1)^2);...
                  0];
        elseif np==3 % derivative wrt KeB
            J(:,1) = [-(nB*(BLag/KiB)^(nB - 1))/(KiB*((BLag/KiB)^nB + 1)^2); 0];
        elseif np==4 % derivative wrt KiB
            J = [(KeB*nB*(BLag/KiB)^(nB - 1))/(KiB^2*((BLag/KiB)^nB + 1)^2) -...
                (2*BLag*KeB*nB^2*(BLag/KiB)^(2*nB - 2))/(KiB^3*((BLag/KiB)^nB + 1)^3)...
                + (BLag*KeB*nB*(BLag/KiB)^(nB - 2)*(nB - 1))/(KiB^3*((BLag/KiB)^nB + 1)^2);...
                0];
        elseif np==7 % derivative wrt nP
            J = [0;...
                (KaP*KeP*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^2)...
                + (KaP*KeP*nP*log(KaP/BLag)*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^2)...
                - (2*KaP*KeP*nP*log(KaP/BLag)*(KaP/BLag)^nP*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^3)];
        elseif np==8 % derivative wrt KeP
            J = [0;...
                (KaP*nP*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^2)];
        elseif np==9 % derivative wrt KaP
            J = [0;...
                (KeP*nP*(KaP/BLag)^(nP - 1))/(BLag^2*((KaP/BLag)^nP + 1)^2)...
                - (2*KaP*KeP*nP^2*(KaP/BLag)^(2*nP - 2))/(BLag^3*((KaP/BLag)^nP + 1)^3)...
                + (KaP*KeP*nP*(KaP/BLag)^(nP - 2)*(nP - 1))/(BLag^3*((KaP/BLag)^nP + 1)^2)];
        end
    end
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %% Second order derivatives wrt state variables
    if nx(2)==0 % outer derivative wrt x(t)
        if nx(1)==0 % inner derivative wrt x(t)
            % compute [diff(diff(dy,B),B)*v(1) + diff(diff(dy,P),B)*v(2),
             % diff(diff(dy,B),P)*v(1) + diff(diff(dy,P),P)*v(2)];
            J(:,1) = [-KdBP*v(2); 0];
            J(:,2) = [-KdBP*v(1); 0];
        elseif nx(1)== 1 % inner derivative wrt x(t-tau)
            % compute [diff(diff(dy,BLag),B)*v(1) + 0,
             %         diff(diff(dy,BLag),P)*v(1) + 0];
             J = zeros(2,2);
        end
    elseif nx(2)==1 % outer derivative wrt x(t-tau)
        if nx(1)==0 % inner derivative wrt x(t)
            % compute [diff(diff(dy,B),BLag)*v(1) + diff(diff(dy,P),BLag)*v(2), [0;0]];
            J = zeros(2,2);
        elseif nx(1)== 1 % inner derivative wrt x(t-tau)
            % compute [diff(diff(dy,BLag),BLag)*v(1) + 0, [0;0]];
            J(:,1) = [v(1)*((2*KeB*nB^2*(BLag/KiB)^(2*nB - 2))/(KiB^2*((BLag/KiB)^nB + 1)^3)...
                - (KeB*nB*(BLag/KiB)^(nB - 2)*(nB - 1))/(KiB^2*((BLag/KiB)^nB + 1)^2));...
                     -v(1)*((2*KaP*KeP*nP*(KaP/BLag)^(nP - 1))/(BLag^3*((KaP/BLag)^nP + 1)^2)...
                     - (2*KaP^2*KeP*nP^2*(KaP/BLag)^(2*nP - 2))/(BLag^4*((KaP/BLag)^nP + 1)^3)...
                     + (KaP^2*KeP*nP*(KaP/BLag)^(nP - 2)*(nP - 1))/(BLag^4*((KaP/BLag)^nP + 1)^2))];
            J(:,2) = [0; 0];
        end
    end
end
%% Otherwise raise error
% Raise error if the requested derivative does not exist
if isempty(J)
    error(['SYS_DERI: requested derivative nx=%d, np=%d, size(v)=%d',...
        'could not be computed!'],nx,np,size(v));
end
end