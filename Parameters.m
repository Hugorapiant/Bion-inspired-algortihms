%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                       %%%
%%%  TECHNICAL FEATURES SIX-PHASE I.M. W/ SYMMETRIC W.    %%%
%%%                                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%%
%%%%%%% FEEDING SYSTEM VALUES %%%%%
Fs = 100000;  Ts = 1/Fs;  dTs=Ts;
a = pi/3;                % Phase Shift
Us = sqrt(2)*180;        % [Volt] Voltage Supply VF - 180V
deg = 0;
f = 50;                  % [Hz] Frequency Supply
fcarrier = 5e3;          % [Hz] Frequency Of Carrier Signal   
L = 2;                   % [s] Length of signal         
npoint= Fs; 
A=exp(j*pi/3);
fd=(750*2)/60;

%%
%%%%%%% PARAMETERS E.C. %%%%%%%%%% 
             %Potencia   kW
Rcore = 811;           % [Ohm] Magnetic Losses Resistance 
Rs = 1.87;             % [Ohm] Stator Phase Resistance
Rr = 2.98;             % [Ohm] Rotor Phase Resistance
Lls = 0.0148;          % [H] Stator Phase Leakage Inductance
Llr = 0.0148;          % [H] Rotor Phase Leakage Inductance
Lm = 0.1861;           % [H] Mutual/Magnetization Inductance
J = 0.0243;            % [kgm^2] Moment of Inercia of Rotor
Tav = 0.23;
p = 2;                 % Number of Pole's Pairs

T = [[0.166667,  0.166667,  0.166667,     0.166667,  0.166667,  0.166667];...
    [0.333333,  0.166667, -0.166667,    -0.333333, -0.166667,  0.166667];...
    [       0,  0.288675,  0.288675,  1.29526e-16, -0.288675, -0.288675];...
    [0.333333, -0.166667, -0.166667,     0.333333, -0.166667, -0.166667];...
    [       0,  0.288675, -0.288675, -2.59052e-16,  0.288675, -0.288675];...
    [0.166667, -0.166667,  0.166667,    -0.166667,  0.166667, -0.166667]];

% MPC parameters
Lr=Llr+Lm;
Ls=Lls+Lm;
Tr=Lr/Rr;
c1=(Ls*Lr-(Lm^2));
c2=Lr/c1;
c3=Lm/c1;
c4=Ls/c1;
c5=1/Lls;
Tmax=50;

kp=2;   %2; 1 
ki=45;   %45; 70

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ABCDEF FRAME         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rss = Rs*eye(6,6);   %[Ohm] Matrix Stator Resistance
Rrr = Rr*eye(6,6);  %[Ohm] Matrix Rotor Resistance
Msr = (1/3)*Lm;     %[H] Maximum Value of the mutual inductance between stator and rotor phase
LA = (1/3)*Lm;      %[H] Stator n Rotor Phase Self-Inductance LA,LB,LC,LD,..

%%%% Matrix Indutancia propria do estator Lls+Lms
Lss = [((1/3*Lm)+Lls) 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm;...
       1/6*Lm ((1/3*Lm)+Lls) 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm;...
       -1/6*Lm 1/6*Lm ((1/3*Lm)+Lls) 1/6*Lm -1/6*Lm -1/3*Lm;...
       -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Lls) 1/6*Lm -1/6*Lm;...
       -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Lls) 1/6*Lm;...
       1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Lls)];
   
%%%% Matrix Indutancia propria do rotor Llr+Lms
Lrr = [((1/3*Lm)+Llr) 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm;...
       1/6*Lm ((1/3*Lm)+Llr) 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm;...
       -1/6*Lm 1/6*Lm ((1/3*Lm)+Llr) 1/6*Lm -1/6*Lm -1/3*Lm;...
       -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Llr) 1/6*Lm -1/6*Lm;...
       -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Llr) 1/6*Lm;...
       1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm ((1/3*Lm)+Llr)];

%%%%  Matrix Mutual Inductances Between STATOR n ROTOR
Lsr = [0 (pi/3) 2*(pi/3) 3*(pi/3) 4*(pi/3) 5*(pi/3);...
       -(pi/3) 0 (pi/3) 2*(pi/3) 3*(pi/3) 4*(pi/3);...
       -2*(pi/3) -(pi/3) 0 (pi/3) 2*(pi/3) 3*(pi/3);...
       -3*(pi/3) -2*(pi/3) -(pi/3) 0 (pi/3) 2*(pi/3);...
       -4*(pi/3) -3*(pi/3) -2*(pi/3) -(pi/3) 0 (pi/3);...
       -5*(pi/3) -4*(pi/3) -3*(pi/3) -2*(pi/3) -(pi/3) 0];
 
%%%%  Matrix Mutual Inductances Between ROTOR n STATOR
Lrs = transpose(Lsr);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC FAULT MATRIX & Factors   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Torque
ST = 0.5;          % Step Ttime
IV = 0;            % Initial Value
FV = 9.8;           % Final Value
Lto = 0;           % Load torque oscillation level
RBin= 7;
RU=4;

N = 138;           % Number of bobine turns
Ns = 12;           % Number of Short Circ turns
u = Ns/N;          % Factor
Rk = 2;            % Short Circuit Resistance

SC='phA';

if(SC=='phA')
    
    Rssc = [[Rss, [-u*Rs; 0; 0; 0; 0; 0]];[u*Rs 0 0 0 0 0 -(u*Rs+Rk)]];
    
    Llsc = [[Lls*eye(6,6), [-(u^2)*Lls; 0; 0; 0; 0; 0]];[u^2*Lls 0 0 0 0 0 -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  -(u/3)*Lm;...
        1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  -(u/6)*Lm;...
        -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm  (u/6)*Lm;...
        -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  (u/3)*Lm;...
        -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  (u/6)*Lm;...
        1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm  -(u/6)*Lm;...
        (u/3)*Lm (u/6)*Lm -(u/6)*Lm -(u/3)*Lm -(u/6)*Lm (u/6)*Lm -u^2*Lm];
    
    Lsc = (Llsc+Lssc);
    MM = [0 (pi/3) 2*(pi/3) (pi) 4*(pi/3) 5*(pi/3)];
end

if(SC=='phB')
    Rssc = [[Rss, [0; -u*Rs; 0; 0; 0; 0]];[0 u*Rs 0 0 0 0 -(u*Rs+Rk)]];
    
    Llsc = [[Lls*eye(6,6), [0; -(u^2)*Lls; 0; 0; 0; 0]];[0 u^2*Lls 0 0 0 0 -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  -(u/6)*Lm;...
        1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  -(u/3)*Lm;...
        -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm  -(u/6)*Lm;...
        -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  (u/6)*Lm;...
        -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  (u/3)*Lm;...
        1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm  (u/6)*Lm;...
        (u/6)*Lm (u/3)*Lm (u/6)*Lm -(u/6)*Lm -(u/3)*Lm -(u/6)*Lm -u^2*Lm];
    
    
    Lsc = (Llsc+Lssc);
    MM = [-(pi/3) 0  (pi/3) 2*(pi/3) 3*(pi/3) 4*(pi/3)];
end


if(SC=='phC')
    Rssc = [[Rss, [0; 0; -u*Rs; 0; 0; 0]];[0 0 u*Rs 0 0 0 -(u*Rs+Rk)]];
    
      Llsc = [[Lls*eye(6,6), [0; 0; -(u^2)*Lls; 0; 0; 0]];[0 0 u^2*Lls 0 0 0 -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  (u/6)*Lm;... 
            1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  -(u/6)*Lm;...
           -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm  -(u/3)*Lm;...
           -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  -(u/6)*Lm;... 
           -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  (u/6)*Lm;...
            1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm  (u/3)*Lm;...
      -(u/6)*Lm (u/6)*Lm (u/3)*Lm (u/6)*Lm -(u/6)*Lm -(u/3)*Lm -u^2*Lm];
  
   
    Lsc = (Llsc+Lssc);
    MM = [-(2*pi/3) -(pi/3) 0  (pi/3) 2*(pi/3) 3*(pi/3)];
end

if(SC=='phD')
    Rssc = [[Rss, [0; 0; 0; -u*Rs; 0; 0]];[0 0 0 u*Rs 0 0 -(u*Rs+Rk)]];
    
    Llsc = [[Lls*eye(6,6), [0; 0; 0; -(u^2)*Lls; 0; 0]];[0 0 0 u^2*Lls 0 0 -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  (u/3)*Lm;...
        1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  (u/6)*Lm;...
        -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm  -(u/6)*Lm;...
        -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  -(u/3)*Lm;...
        -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  -(u/6)*Lm;...
        1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm  (u/6)*Lm;...
        -(u/3)*Lm -(u/6)*Lm (u/6)*Lm (u/3)*Lm (u/6)*Lm -(u/6)*Lm -u^2*Lm];
    
    
    Lsc = (Llsc+Lssc);
    MM = [-3*(pi/3) -(2*pi/3) -(pi/3) 0  (pi/3) 2*(pi/3)];
end

if(SC=='phE')
    Rssc = [[Rss, [0; 0; 0; 0; -u*Rs; 0]];[0 0 0 0 u*Rs 0 -(u*Rs+Rk)]];
    
      Llsc = [[Lls*eye(6,6), [0; 0; 0; 0; -(u^2)*Lls; 0]];[0 0 0 0 u^2*Lls 0 -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  (u/6)*Lm;... 
            1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  (u/3)*Lm;...
           -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm   (u/6)*Lm;... 
           -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  -(u/6)*Lm;... 
           -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  -(u/3)*Lm;... 
            1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm -(u/6)*Lm;...
      -(u/6)*Lm -(u/3)*Lm -(u/6)*Lm (u/6)*Lm (u/3)*Lm (u/6)*Lm -u^2*Lm];
  
   
    Lsc = (Llsc+Lssc);
    MM = [-(4*pi/3) -(3*pi/3) -(2*pi/3) -(pi/3) 0  (pi/3)];
end

if(SC=='phF')
    Rssc = [[Rss, [0; 0; 0; 0; 0; -u*Rs]];[0 0 0 0 0 u*Rs -(u*Rs+Rk)]];
    
      Llsc = [[Lls*eye(6,6), [0; 0; 0; 0; 0; -(u^2)*Lls]];[0 0 0 0 0 u^2*Lls -(u^2)*Lls]];
    
    Lssc = [1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm  -(u/6)*Lm;...
            1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm  (u/6)*Lm;...
           -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm -1/3*Lm   (u/3)*Lm;... 
           -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm -1/6*Lm  (u/6)*Lm;...
           -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm 1/6*Lm  -(u/6)*Lm;... 
            1/6*Lm -1/6*Lm -1/3*Lm -1/6*Lm 1/6*Lm 1/3*Lm -(u/3)*Lm;...
      (u/6)*Lm -(u/6)*Lm -(u/3)*Lm -(u/6)*Lm (u/6)*Lm (u/3)*Lm -u^2*Lm];
  
   
    Lsc = (Llsc+Lssc);
    MM = [-(5*pi/3) -(4*pi/3) -(3*pi/3) -(2*pi/3) -(pi/3) 0];
end
%%
disp('----- READY -----')
