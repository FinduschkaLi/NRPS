function [y,dydt] = RungeK(fcthandle,t,y0,optVars)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function iterates through a time vector t and solves the
    % differential function given by fcthandle.
    % Initial conditions are provided by y0.
    % An optional time series optVars can be passed to the differential function
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %initialize result vectors****************************************
    y = zeros(size(y0,1),size(t,2));
    dydt=y;
    y(:,1)=y0;
    dydt(:,1)=0;
    % Loop through time vector****************************************
    % in case an additional time-series is passed to fcthandle, 
    % pass the value at each iteration
    if(nargin>3)
        for i=2:1:size(t,2)
            h=t(i)-t(i-1);
            k1=fcthandle(t(i-1),y(:,i-1),optVars(:,i-1));
            k2=fcthandle(t(i-1)+h/2,y(:,i-1)+h*k1/2,optVars(:,i-1));
            k3=fcthandle(t(i-1)+h/2,y(:,i-1)+h*k2/2,optVars(:,i-1));
            k4=fcthandle(t(i),y(:,i-1)+h*k3,optVars(:,i-1));
            y(:,i)=y(:,i-1)+1/6*h*(k1+2*k2+2*k3+k4);
            dydt(:,i)=fcthandle(t(i),y(:,i),optVars(:,i));
        end
    else %regular execution
        for i=2:1:size(t,2)
            h=t(i)-t(i-1);
            k1=fcthandle(t(i-1),y(:,i-1));
            k2=fcthandle(t(i-1)+h/2,y(:,i-1)+h*k1/2);
            k3=fcthandle(t(i-1)+h/2,y(:,i-1)+h*k2/2);
            k4=fcthandle(t(i),y(:,i-1)+h*k3);
            y(:,i)=y(:,i-1)+1/6*h*(k1+2*k2+2*k3+k4);
            dydt(:,i)=fcthandle(t(i),y(:,i));
        end
    end
end
