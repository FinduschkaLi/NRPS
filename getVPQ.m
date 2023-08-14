function [V,pout,qout]=getVPQ(Yv_a,PSI,Yev,E,TET)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function calculates the voltage at the capcitor,
    % the active and reactive output power for the RAPS or the NRPS model
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=size(PSI,1);
    v=Yev*(E.*exp(sqrt(-1)*TET));
    V=abs(v);
    VMAT=repmat(V,1,n).*repmat(V',n,1);
    ZET=angle(v);
    DELTAZET=repmat(ZET,1,n)-repmat(ZET',n,1);
    pout= VMAT.*Yv_a.*sin(DELTAZET-PSI)*ones(n,1);
    qout=-VMAT.*Yv_a.*cos(DELTAZET-PSI)*ones(n,1);
end
