classdef u<uBase
%{
   Description: unit handling class

   Copyright (c) 2023, Karl Erik Thoresen 
   All rights reserved.

%}
    properties(Constant)
        %============ START THE ACTUAL CODE TO DEFINE THE UNITS STRUCT =========
        %------- SI fundamental -------


        %--------Named units derived from SI base units
        % Expressed in terms of other SI units
       rad  = u.m/u.m;           % Angle Radian
       sr   = u.m^2/u.m^2;       % Solid angle Steradian
       Hz   = 1/u.s;             % Frequency hertz
       N    = u.kg*u.m/(u.s^2);  % Force newton
       Pa   = u.N/u.m^2;         % Pressure pascal
       J    = u.N*u.m;           % Energy joule
       W    = u.J/u.s;           % Power watt
       C    = u.s*u.A;           % Electric charge coulomb
       V    = u.W/u.A;           % Voltage volt
       F    = u.C/u.V;           %  Electric capacitance farad
       ohm  = u.V/u.A;           % Electric resistance ohm
       S    = u.A/u.V;       	% Electrical conductance siemens
       Wb   = u.V*u.s;           % Magnetic flux weber
       T	   = u.Wb/u.m^2;        % Magnetic field strength tesla
       H    = u.Wb/u.A;          % Inductance henry
       lm   = u.cd*u.sr;         % Luminous flux lumen
       lx   = u.lm/u.m^2;        % Illuminance lux
       Bq   = 1/u.s;             % Radioactivity (decays per unit time) becquerel
       Gy   = u.J/u.kg;          % Absorbed dose (of ionizing radiation)gray
       Sv   = u.J/u.kg;          % Equivalent dose (of ionizing radiation)sievert
       kat  = u.mol/u.s;         % Catalytic activity katal
        
        %------- length ----
       km = 1e3*u.m;
       cm = 1e-2*u.m;
       dm=  1e-1*u.m;
       mm = 1e-3*u.m;
       um = 1e-6*u.m;
       nm = 1e-9*u.m;
       ang = 1e-10*u.m;
       in = 2.54*u.cm;
       mil = 1e-3*u.in;
       ft = 12*u.in;
       yd = 3*u.ft;
       mi = 5280*u.ft;
       a0 = .529e-10*u.m;
        
        %------- Volume -------
       cc = (u.cm)^3;
       L = 1000*u.cc;
       mL = u.cc;
       gal  = 231*u.in^3;        % US gallon
       quart = u.gal/4;
       floz = u.gal/128;
       pint = 16*u.floz;         % pint (USA, wet) = 16 US fluid ounces
       bbl  = 42*u.gal;          % US oil barrel
        
        %----- mass ---------
       gm = 1e-3*u.kg;           % gram
       mg = 1e-3*u.gm;           % milligram
       lb = 0.45359237*u.kg;     % pound mass
       lbm = 0.45359237*u.kg;     % pound mass
       oz = (1/16)*u.lb;
       amu = 1.66e-27*u.kg;
      % ton = 2000*u.lbm % ton is often defined differently
      %  tonn = 1000*u.kg;
       tonne = 1000*u.kg; % t or tonne is the correct unit for metric ton
       tn = 2000*u.lbm; % Short ton, 2000 lbm, 907.18 kg
        % u.lbm = 0.45359237;         % lbm -> kg
        % u.lbs = 0.4536;             % lbs -> kg
 
        % densities
        ppg = u.lb/u.gal;         % pound pr gallon
        sg  = 1000*u.kg/u.m^3    % specific gravity       
  

        %---- time -------
       ms = 1e-3*u.s;
       us = 1e-6*u.s;
       ns = 1e-9*u.s;
       ps = 1e-12*u.s;
       min = 60*u.s;
       hr  = 60*u.min;
       day = 24*u.hr;
       yr  = 365.242199*u.day;
        
        %---- angular units --------------
       deg = 2*pi/360;             % deg -> rad
       rev = 2*pi;
        
        %---- frequency ----
       kHz = 1e3 *u.Hz;
       MHz = 1e6 *u.Hz;
       GHz = 1e9 *u.Hz;
       rpm = u.rev/u.min;        % rpm -> rad/s
       spm = u.rpm;              % strokes pr minute (for mud pumps)
        
        %---- force -------
       dyne = 1e-5*u.N;
       lbf  = 4.44822*u.N;
       klbf = 1000*u.lbf;
       tf   = 1000*u.kg*c.g;     % Metric ton - force
       tnf  = 2000*u.lbf;      % Short ton - force
       kN   = 1000*u.N; 
       %---- Speed -------
       kmh = u.km/u.hr;
        %---- torque -------
       Nm  = u.N*u.m;
       kNm = 1000*u.Nm;
        
        %---- flow ---------
       lpm = u.L/u.min;            % lpm -> m3/s
       gpm = u.gal/u.min;          % gpm -> m3/s
        
        %---- fluid viscosity data ---------
        % Based on   %g300 =  511 %https://www.fann.com/content/dam/fann/Manuals/Model%2035%20Viscometer.pdf
        % The g300 is probably 510.5 based on the values of the g600 of
        % 1021  something that is comfirmed in:
        % http://mycommittees.api.org/standards/ecs/sc10/Meeting%20Materials/2014/January%20-%20Davison%20RP%2010B2%20Review/Bernard/Fann%20coef%20calculation-2.pdf
        % Further the Fann viscometer documentation gives
        % ktau=1.065*u.lbf/(100*u.ft^2)  while the api says 1.067
        % The paper  at api.org also says 1.067 or actually 1.0666
        %kTau = 1.0678*u.lbf/(100*u.ft^2);% shear stress factor(Pa/deg) of fann viscometer frp, API->19
        fann = 1.065*u.lbf/(100*u.ft^2); % from fann manual
        %   fann = 1.0678*u.lbf/(100*u.ft^2);% shear stress factor(Pa/deg) of fann viscometer frp, API->19

        g300 = 511;   % shear rate (1/s) at fann viscometer speed = 300 rpm

      %g300 =510.69
        kGam = u.g300/300;            % shear rate coefficient (1/s)/rpm
        cP   = 1e-3*u.Pa*u.s;      % centiPose - dynamic viscosity unit
        cSt= 1 *u.mm^2/u.s;        % centiStokes- kinematic viscosity;
        %----- energy -----
       MJ = 1e6*u.J;
       kJ = 1e3*u.J;
       mJ = 1e-3*u.J;
       uJ = 1e-6*u.J;
       nJ = 1e-9*u.J;
       eV = 1.6022e-19*u.J;
       BTU = 1.0550559e3*u.J;
       kWh = 3.6e6*u.J;
       cal = 4.1868*u.J;
       kCal = 1e3*u.cal;
        
        %---- temperature ---
       mK   = 1e-3*u.K;
       uK   = 1e-6*u.K;
       nK   = 1e-9*u.K;
       
       %%-----thermal
       k_cond= u.W/(u.m*u.K);



      % degC = 273.15;% NB uses adittion
       degF = u_degF;
       degC  = u_degC;
        %---- pressure -----
       torr = 133.322*u.Pa;
       mtorr = 1e-3*u.torr;
       bar = 1e5*u.Pa;
       mbar = 1e-3*u.bar;
       atm = 1.013e5*u.Pa;
       psi = 6.8947548e3*u.Pa;
        
        %----- power --- ---
       MW = 1e6*u.W;
       kW = 1e3*u.W;
       mW = 1e-3*u.W;
       uW = 1e-6*u.W;
       nW = 1e-9*u.W;
       pW = 1e-12*u.W;
       hp = 745.69987*u.W;
        
        %------ Voltage -----
       kV = 1e3*u.V;
       mV = 1e-3*u.V;
       uV = 1e-6*u.V;
        
        %----- Current ------
       mA = 1e-3*u.A;
       uA = 1e-6*u.A;
       nA = 1e-9*u.A;
        
        %----magnetic field -----
       gauss = 1e-4*u.T;
        
        %% Choke data
       rho_ref = 998.203*u.kg/u.m^3;    % density of water @20 degC
       CvUS =  real(sqrt(u.rho_ref*u.gal/u.lb)*u.gal/u.min*sqrt(u.lb/u.gal/u.psi)); % conversion from Cv given in imperial units to SI
        
    %% unit handle
       P    = UnitType.Pressure;
       Ecd  = UnitType.ECD;
       Flw  = UnitType.Flowrate;
       Frc  = UnitType.Force;% kN
       Trq  = UnitType.Torque;
       Len  = UnitType.Length;
       Dia  = UnitType.Diameter;
       Rel  = UnitType.Ratio;
       Vol  = UnitType.Volume;
       Vel  = UnitType.Velocity;
       Ang  = UnitType.Angle;
       Dens = UnitType.Density;
       None = UnitType.None;
       Rot =  UnitType.Rotation;
       WOB  = UnitType.WOB;
       ROP   =UnitType.ROP;
       Area  =  UnitType.Area
       Nozzle =UnitType.Nozzle;
       Temp   = UnitType.Temperature;
       TempGrad = UnitType.TemperatureGradient;
       Time     = UnitType.Time;
       Stress =   UnitType.Stress;
       ThrmC   = UnitType.ThermalConductivity; 
       HeatC   = UnitType.SpecHeatCapacity;
       Fann    =UnitType.Fann;
       Visc     = UnitType.Viscosity;
       dP      =UnitType.DeltaPressure;
       MassFlwTon =UnitType.MassFlowRateTon;
       PowerMega = UnitType.PowerMega;
       Hookload = UnitType.Hookload;
       Mass = UnitType.Mass;
    end
    
    methods (Static)
       % function a=ton()
       %     a=nan; %#ok<NASGU>
       %     error('disambigous errror, ton was defined as tonn or tonne, which was wrong.Please check what to yse');
       % end
        function ret=lpmd(~),   ret= 51 ;        end
        function Kelvin=FtoK(F)
            %% Kelvin=FtoK(F)
            Kelvin =(F+459.67)/1.8;
        end
        function Deg=FtoDeg(F)
            %% Kelvin=FtoK(F)
            
            Kelvin =(F+459.67)/1.8;
            Deg=Kelvin-273.15;
        end
        function F = KtoF(Kelvin)
            Deg = Kelvin-273.15;
            F = Deg * (9/5) + 32;
        end

        function str=str(val,UnitType,inspec)
            %% str=prt(val,UnitType,spec)
            % automatic create string of umber with units
            %
            if nargin ==3
                spec=inspec;
            else
                switch UnitType
                    % Try to guess the most common digits representaion to
                    % use
                    case u.Len
                        spec='%.0f';
                    case u.Fann
                        spec='%.1f';
                    case u.Vel
                        spec='%.2f';
                    otherwise
                        spec=' %.1f';
                end

            end
            if  length(val)>1
                spec=strcat(spec,',');
                str=""+ sprintf(spec,val/UnitType)+" " +UnitType.Str;
            else
                str=""+ sprintf(spec,val/UnitType)+" " +UnitType.Str;
            end
        end

        function str=st(val,UnitType,inspec)
            name= inputname(1);
             if nargin ==3,str=u.str(val,UnitType,inspec);
             else,         str=u.str(val,UnitType);end
             str=sprintf('%s =%s',name,str);
        end
        function setSymBase()
            sympref('FloatingPointOutput',true)
            p = mfilename('fullpath');
            [curDir, ~, ~] = fileparts(p);
            copyfile(curDir+"\uSymBase.m", curDir+"\uBase.m");
        end

        function setNumBase()
            p = mfilename('fullpath');
            [curDir, ~, ~] = fileparts(p);
            copyfile(curDir+"\uNumBase.m", curDir+"\uBase.m");

        end
    end

end