function [pout]=getP(PHI,YRED,E0,TET)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function calculates the output power for the NRPS model
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=size(YRED,1);
    DELTATET=repmat(TET,1,n)-repmat(TET',n,1);
    EMAT=repmat(E0,1,n).*repmat(E0',n,1);
    YY=(abs(YRED));
    pout=(EMAT.*YY.*sin(DELTATET-PHI))*ones(n,1);