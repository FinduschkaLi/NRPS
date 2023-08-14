function dXdt=NRPS(t,X,PHI,M,D,P,A,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equation of the standard NRPS model
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TET=X(1:n);
    OMG=X(n+1:n*2);
    %calculate angular difference matrix
    DELTATET=repmat(TET,1,n)-repmat(TET',n,1);
    %calculate model
    dOMG=M\(P-D*OMG-A.*sin(DELTATET-PHI)*ones(n,1));% attention, cannot make e infinitly small!!
    dTET=OMG;
    dXdt=[dTET;dOMG];
end
