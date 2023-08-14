%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This model compares the NRPS and the RAPS model in a 4 generator grid
% where after 10 seconds, the Pset value is changed for one generator
% 2023 FREISSNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general definitions
i_=sqrt(-1);        %imaginary number
dt=0.004;           %time step of the model
nRED=4;             %number of generators
%nominal values
omega_n=pi*100;
V_n=110;

%VSM configurations*************************************************
Rg =[1.1437    1.3724    4.9014    2.6392];
Lg =[0.0286    0.0343    0.1225    0.0660];
Dp=[2.0264    1.6887    0.4728    0.8781];
J=[1.2523    0.6755    0.1064    0.1633];
Dq=[31.3636   26.1364    7.3182   13.5909];
Km=[4.1400    3.4500    0.9660    1.7940]*1e3;
%Filter parameters of the inverter (only the "CL" part of the LCL filter)
L_filter=7.5e-4;
R_filter=0.5;
C_filter=6.0000e-06;
R_damp = 5;

%Loads**************************************************************
RL=[25.8667,   30.0000,   45.0000,   45.0000];
LL=ones(1,nRED)*115.55;

%Tie-lines**********************************************************
Lline =[ 0.0100    0.0200    0.0300    0.0300];
Rline=[0.5 0.5 0.5 0.5];

% Power references**************************************************
PsetData = [
        1600        1600;
        1300        2000;
         600         600;
        1100        1100];

QsetData = zeros(4,2);
PowerRefTime=[10 20];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diagonal matrix Yg required for RAPS model:
Yg=diag(3./(Rg+Lg*i_*omega_n));

% Line impedances***************************************************
Ylines=[];
%Add the virtual stators of the VSMs
Ylines=[Ylines; 1 5 3/(Rg(1)+Lg(1)*i_*omega_n)];
Ylines=[Ylines; 2 6 3/(Rg(2)+Lg(2)*i_*omega_n)];
Ylines=[Ylines; 3 7 3/(Rg(3)+Lg(3)*i_*omega_n)];
Ylines=[Ylines; 4 8 3/(Rg(4)+Lg(4)*i_*omega_n)];
%Add the filters of the inverters
Ylines=[Ylines; 5 9  3/(R_filter(1)+L_filter(1)*i_*omega_n)];
Ylines=[Ylines; 6 10 3/(R_filter(1)+L_filter(1)*i_*omega_n)];
Ylines=[Ylines; 7 11 3/(R_filter(1)+L_filter(1)*i_*omega_n)];
Ylines=[Ylines; 8 12 3/(R_filter(1)+L_filter(1)*i_*omega_n)];
%Add the tie lines
TLNumbering=[1 12 9; 2 9 10; 3 10 11; 4 11 12; ];
for i_T=1:4
    Ylines=[Ylines;TLNumbering(i_T,3) TLNumbering(i_T,2) 3/(Rline(i_T)+Lline(i_T)*i_*omega_n)];
end

% Load impedances (impedances connected line-ground*****************
%add the capacitor of the LCL filter of the inverters
YLoads =       [ 5 3/(1/(C_filter*omega_n*i_)+R_damp);
                6 3/(1/(C_filter*omega_n*i_)+R_damp);
                7 3/(1/(C_filter*omega_n*i_)+R_damp);
                8 3/(1/(C_filter*omega_n*i_)+R_damp);
                ];
%add the actual loads in the grid
LoadNodes=[9 10 11 12];
for i_L=1:4
    YLoads= [YLoads;LoadNodes(i_L) 3/RL(i_L)+3/(i_*omega_n*LL(i_L));];
end

% translate VSM parameters to NRPS parameters***********************
M=diag(J)*omega_n;
D=diag(Dp)*omega_n;
R=diag(Dq)*sqrt(2);
K=diag(Km)/omega_n*sqrt(3);
F=zeros(nRED,nRED);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawNetwork(nRED,Ylines,50)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run both the NRPS and the RAPS model for each set of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial values for voltages, angles and frequencies****************
E0=[121 120 134 133]';
E0_3=ones(nRED,1)*V_n;
TET0=zeros(nRED,1);
OMG0=ones(nRED,1)*omega_n;
X0=[TET0; OMG0];
X0_3=[TET0;OMG0;E0_3];
%Init empty arrays**************************************************
%Time vector
Tspan=[];
%initialize values 3O
XNRPS_3=[];
V_3=[];
E_3=[];
POUT_3=[];
QOUT_3=[];
OMG_3=[];
GNDA_3=[];
E0_3=E0;
%initialize values NRPS
XNRPS=[];
V=[];
POUT=[];
QOUT=[];
OMG=[];
GNDA=[];

T_last=0;

%Run the actual model*************************************************
%For each set of configuration you must regenerate all parameters. In this 
%example, no impedance (load or line) changes, only Pset changes.

for i_T=1:2 %iterate through phases of constant configuration (i_T is used for Pset here)

    %generate reduced Admittance matrix (only used for NRPS model)
    YRED=reduceNetwork(nRED,Ylines,YLoads);
    [A,PHI,P]=generateNRPSParameters(YRED,E0,PsetData(:,i_T),Dp'*omega_n,omega_n);
    %calculate admittance matrices:
    %Yv_a and PSI are the absolute value and angle-pi/2 of the admittance
    %matrix, such that: i=Yv_a.*exp(i_*PSI)*v where
    %v: phasor of voltages at the filter capacitors
    %i: phasor of currents flowing into the grid
    [Yv_a,Yev,PSI]=generateRAPSParameters(YRED,Yg);
    %get the end time of this configuration (when does the configuration
    %change again?)
    T_next=PowerRefTime(i_T);
    %define a time vector for the current configuration
    tspan=T_last:dt:T_next;
    %run Simulations*****************************************************
    %t,X,PHI,M,R,D,Dq,F,Pset,Qset,Eset,Y,n
    [xNRPS_3, ~]=   RungeK(@(t,X)RAPS(t,X,PSI,M,K,D,R,F,PsetData(:,i_T),QsetData(:,i_T),V_n,Yv_a,Yev,omega_n,nRED),tspan,X0_3);
    [xNRPS, ~]=     RungeK(@(t,X)NRPS(t,X,PHI,M,D,P,A,nRED),tspan,X0);
    %get theta, omega, E, V, POUT, QOUT, grounded angles for each timestep
    %unpack values for the RAPS model (third order)
    tetNRPS_3=xNRPS_3(1:nRED,:);
    omgNRPS_3=xNRPS_3(nRED+1:nRED*2,:);
    e_3=xNRPS_3(nRED*2+1:nRED*3,:);
    omgCOI_3=J*omgNRPS_3./sum(J);
    tetCOI_3=angle(mean(exp(sqrt(-1)*tetNRPS_3)));
    gndNRPS_3=mod(tetNRPS_3(2:end,:)-tetNRPS_3(1,:)+pi,2*pi)-pi;
    %unpack values for the NRPS model (second order)
    tetNRPS=xNRPS(1:nRED,:);
    omgNRPS=xNRPS(nRED+1:nRED*2,:);
    omgCOI=J*omgNRPS./sum(J);
    tetCOI=angle(mean(exp(sqrt(-1)*tetNRPS)));
    gndNRPS=mod(tetNRPS(2:end,:)-tetNRPS(1,:)+pi,2*pi)-pi;
    %calculate the output powers****************************************
    pout_3=ones(nRED,length(tspan));
    qout_3=pout_3;
    pout=pout_3;
    qout=pout_3;
    v=pout_3;
    v_3=pout_3;
    for i=1:length(tspan)
        %3O[v,pout,qout]=getVPQ(Yv_a,PSI,Yev,E,TET)
        [v_3(:,i),pout_3(:,i),qout_3(:,i)]=getVPQ(Yv_a,PSI,Yev,e_3(:,i),tetNRPS_3(:,i));
        %NRPS
        [v(:,i),pout(:,i),qout(:,i)]=getVPQ(Yv_a,PSI,Yev,E0,tetNRPS(:,i));
    end

    %collect data from constant configurations into single matrices*****
    V_3=[V_3,v_3];
    E_3=[E_3,e_3];
    POUT_3=[POUT_3,pout_3];
    QOUT_3=[QOUT_3,qout_3];
    OMG_3=[OMG_3,omgNRPS_3];
    GNDA_3=[GNDA_3,gndNRPS_3];
    XNRPS_3=[XNRPS_3,xNRPS_3];
    Tspan=[Tspan,tspan];
    %collect all data to one vector NRPS
    V=[V,v];
    POUT=[POUT,pout];
    QOUT=[QOUT,qout];
    OMG=[OMG,omgNRPS];
    GNDA=[GNDA,gndNRPS];
    XNRPS=[XNRPS,xNRPS];
    %update last value for subsequent simulations
    X0=XNRPS(:,end);
    X0_3=XNRPS_3(:,end);
    %move along the scenario defined for the simulation
    T_last=T_next;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
figC=0;
figC=figC+1;
f=figure(figC);
f.Name="Omega";
plot(Tspan,OMG_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1)
plot(Tspan,OMG,':')
ylabel("\omega [rad/s]")
xlabel("time [s]")
xlim([Tspan(1) Tspan(end)])

figC=figC+1;
f=figure(figC);
f.Name="Active Power";
plot(Tspan,POUT_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1)
plot(Tspan,POUT',':')
ylim([0 4000])
ylabel("P [W]")
xlabel("time [s]")
xlim([Tspan(1) Tspan(end)])

figC=figC+1;
f=figure(figC);
f.Name="Reactive Power";
plot(Tspan,QOUT_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1)
plot(Tspan,QOUT',':')
ylabel("Q [VAr]")
xlabel("time [s]")
xlim([Tspan(1) Tspan(end)])

figC=figC+1;
f=figure(figC);
f.Name="Internal Voltages";
plot(Tspan,E_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1)
plot(Tspan,E0*ones(1,length(Tspan)),':')
ylabel("|E| [V]")
xlabel("time [s]")

figC=figC+1;
f=figure(figC);
f.Name="Grounded Angles";
plot(Tspan,GNDA_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1)
plot(Tspan,GNDA',':')
ylabel("\delta [rad]")
xlabel("time [s]")

figC=figC+1;
f=figure(figC);
f.Name="Voltage at PCC";
plot(Tspan,V_3','--')
hold all;
grid on;
set(gca,'ColorOrderIndex',1);
plot(Tspan,V',':')
ylabel("|V| [V]")
xlabel("time [s]")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% format plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LineStyle
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_axes=findobj(cfig,'type','axes');
    splots=length(h_axes);
    for p=1:splots
        h_axe=h_axes(p);
        if(isa(h_axe,'matlab.graphics.axis.Axes'))
            h_lines  = findall(h_axe,'Type','line');
            for i=1:1:length(h_lines)
                set(h_lines(i),'LineWidth',3);
            end
        end
    end
end
%% Font Size
sze=17;
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_lines = findall(cfig,'Type','axes');
    isax= strcmp(get(h_lines,'Tag'),'legend');
    h_lines(isax)=[];
    for i=1:1:length(h_lines)
        set(h_lines(i),'FontSize',sze)
        set(h_lines(i), 'FontName', 'Arial')
    end   
    tx=findall(cfig,'Type','text');
    for i=1:1:length(tx)
        set(tx(i),'FontSize',sze);
        set(tx(i),'Interpreter','Tex');
        set(tx(i), 'FontName', 'Arial');
    end
end
%% Sizing
sizeFig=[1000 600];
margins=[.1 .12 .02 .03];
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_axes=findobj(cfig,'type','axes');
    pos=get(cfig,"Position");
    pos(3)=sizeFig(1);
    pos(4)=sizeFig(2);
    set(cfig,"Position",pos);
    pos=get(h_axes(1),"Position");
    pos(1)=margins(1);
    pos(2)=margins(2);
    pos(3)=1-margins(1)-margins(3);
    pos(4)=(1-margins(2)-margins(4));
    
    set(h_axes(1),"Position",pos);
end
%% Legends
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    gcf=figHandles(k);
    set (0, 'currentfigure', gcf); 
    if(~isempty(strfind(gcf.Name,"Grounded")))
        hh = legend("$\!$","$\!$","$G_1$","$G_2$","$G_3$",'NumColumns',2);
    else
        hh = legend("$\!$","$\!$","$\!$","$\!$","$G_1$","$G_2$","$G_3$","$G_4$",'NumColumns',2);
    end
    hh.Interpreter = 'latex';
    title(hh,"RAPS\  NRPS$\ \ \ \ $");
    hh.ItemTokenSize = [35,18];
end