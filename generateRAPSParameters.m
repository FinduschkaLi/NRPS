function [Yv_a,Yev,PSI]=generateRAPSParameters(YRED,Yg)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function calculates:
    % - the Matrix Yev, used to obtain the voltages v=Yev*e (phasor notation)
    % - the Matrix Yv, which is used to estimate the currents i=Yv*v
    % (phasor notation).
    % - Yv_a is the absolute value, PSI the angle-pi/2 of Yv
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=size(YRED,1);
    Yev=(eye(n)-Yg\YRED);
    Yv=YRED/Yev;
    PSI=angle(Yv)-pi/2;
    Yv_a=abs(Yv);
end