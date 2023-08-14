function dXdt=FEPS(t,X,PHI,M,D,F,P,A,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equation of the Friction Enhanced Power System model
    % An extension of the NRPS model using Virtual Friction
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %retrieve state variables
    TET=X(1:n);
    OMG=X(n+1:n*2);
    %calculate the COI frequency
    omgC=sum(M*OMG)/sum(diag(M));
    %calculate angular difference matrix
    DELTATET=repmat(TET,1,n)-repmat(TET',n,1);
    %update state variables
    dOMG=M\(P-D*OMG-F*OMG+diag(F)*omgC-A.*sin(DELTATET-PHI)*ones(n,1));% attention, cannot make e infinitly small!!
    dTET=OMG;
    dXdt=[dTET;dOMG];
end
