function [A,PHI,P]=generateNRPSParameters(YRed,E,Pset,D,omega_nom)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function generates the coupling matrix A, 
    % the angle matrix PHI and the vector P, required for running the NRPS
    % model
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nRED=min(size(YRed)); %get system dimension
    %obtain absolute values and angles of reduced network
    YAbs=abs(YRed);
    PHI=atan2(-real(YRed),imag(YRed));
    PHI(isnan(PHI))=0;
    %calculate coupling matrix A and static Power terms
    P=zeros(nRED,1);
    A=zeros(nRED,nRED);
    for i=1:1:nRED
        for j=1:1:nRED
            if i==j
                P(i)=E(i)*E(i)*YAbs(i,i)*sin(PHI(i,i))+Pset(i)+D(i)*omega_nom;%22-02-27 it is +phi, but why needs to be checked.
            else
                A(i,j)=E(i)*E(j)*YAbs(i,j);
            end
        end
    end
end