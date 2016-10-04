function varargout = SolarSUAV_GUI(varargin)
% SOLARSUAV_GUI MATLAB code for SolarSUAV_GUI.fig
%      SOLARSUAV_GUI, by itself, creates a new SOLARSUAV_GUI or raises the existing
%      singleton*.
%
%      H = SOLARSUAV_GUI returns the handle to a new SOLARSUAV_GUI or the handle to
%      the existing singleton*.
%
%      SOLARSUAV_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOLARSUAV_GUI.M with the given input arguments.
%
%      SOLARSUAV_GUI('Property','Value',...) creates a new SOLARSUAV_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SolarSUAV_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SolarSUAV_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SolarSUAV_GUI

% Last Modified by GUIDE v2.5 25-May-2015 14:38:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SolarSUAV_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SolarSUAV_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SolarSUAV_GUI is made visible.
function SolarSUAV_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SolarSUAV_GUI (see VARARGIN)

% Choose default command line output for SolarSUAV_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SolarSUAV_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SolarSUAV_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnRun.
function btnRun_Callback(hObject, eventdata, handles)
% hObject    handle to btnRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% CUT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clear Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cla(handles.axes1);
        cla(handles.axes2);
        cla(handles.axes3);
        cla(handles.axes4);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Clear Outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        set(handles.txtWingArea,'String','');
        set(handles.txtTakeoffWeight,'String','');
        set(handles.txtPowerReq,'String','');
        set(handles.txtBattWh,'String','');
        set(handles.wingSpan,'String','');
        set(handles.airframe_mass,'String','');
        set(handles.battery_mass,'String','');
        set(handles.solarcells_mass,'String','');
        set(handles.propgroup_mass,'String','');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Geographic Position - Input Software
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        geo_Latitude = str2double(get(handles.geo_Latitude,'String'));
        geo_Altitude = str2double(get(handles.geo_Altitude,'String'));
        [geo_rho, geo_press]=stdatmo(geo_Altitude);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Software Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pv_eff = str2double(get(handles.pv_eff,'String'))/100;
        mppt_eff = str2double(get(handles.mppt_eff,'String'))/100;
        Encapsulation = str2double(get(handles.encapsulation,'String')); %kg/m²
        prop_eff = str2double(get(handles.prop_eff,'String'))/100;
        MotorEff = str2double(get(handles.motor_eff,'String'))/100;
        MotorControlEff = str2double(get(handles.control_eff,'String'))/100;
        BEC_eff = str2double(get(handles.BEC_eff,'String'))/100;
        gearbox_eff = str2double(get(handles.gearbox_eff,'String'))/100;
        BattEff = str2double(get(handles.BattEff,'String'))/100;
        Kbatt = str2double(get(handles.battSpecEnergy,'String'));
        e = str2double(get(handles.OswaldCte,'String'));
        AR = str2double(get(handles.AspectRatio,'String'));
        Payload = str2double(get(handles.payloadWeight,'String'));
        PayloadPower = str2double(get(handles.payloadPower,'String'));
        cruiseSpeed = str2double(get(handles.cruiseSpeed,'String'));
        SolarCellDensity = str2double(get(handles.solarCellDensity,'String'));
        CD0 = str2double(get(handles.CD0,'String'));
        Rmax  = str2double(get(handles.Rmax,'String'));
        RdesignInput = str2double(get(handles.txtRDesign,'String'));
        W_SInput = str2double(get(handles.txtW_S,'String'));
        
        eff=MotorControlEff*MotorEff*gearbox_eff*prop_eff*pv_eff*mppt_eff; % From Brandt => prop*gear*motor*cells*(0.45+0.55*batt)
        g = 9.81;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solar Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(handles.axes1)
        [MonthlyIrradiance, HourlyIrradiance]=SolarModel(geo_Latitude,geo_press);
        Irr=max(MonthlyIrradiance); %[W/m²] Imax to constraint analysis 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Day Length
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(handles.axes2)
        [mindaylength,maxdaylength]=daylength(geo_Latitude);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraint Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h3=subplot(handles.axes3);
            cla(h3)
            title('Constraint Diagram','FontWeight','Bold')
            WS_x=0:100;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constant Altitude/Speed Cruise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                V=cruiseSpeed;
                q=0.5*geo_rho*V^2;
                n=1;
                climbrate=0;
                i=1;
                for WS=0:100
                    R(i) = ((q*V)/(eff*Irr))*((1/(pi*e*AR))*((n/q)*(WS))^2+CD0) + (1/(eff*Irr))*WS*climbrate;
                    i=i+1;
                end
                plot(WS_x,R,'b')
                hold on

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constant Altitude/Speed Turn
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                V=cruiseSpeed;
                q=0.5*geo_rho*V^2;
                n = str2double(get(handles.turnLoadFactor,'String'));       
                climbrate=0;
                i=1;
                for WS=0:100
                    R(i) = ((q*V)/(eff*Irr))*((1/(pi*e*AR))*((n/q)*(WS))^2+CD0) + (1/(eff*Irr))*WS*climbrate;
                    i=i+1;
                end
                plot(WS_x,R,'g')
                hold on

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constant Speed Climb
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                V=cruiseSpeed;
                q = 0.5*geo_rho*V^2;
                n = 1; %str2double(get(handles.climbLoadFactor,'String')); 
                climbrate = str2double(get(handles.climbRate,'String'));
                i=1;
                for WS=0:100
                    R(i) = ((q*V)/(eff*Irr))*((1/(pi*e*AR))*((n/q)*(WS))^2+CD0) + (1/(eff*Irr))*WS*climbrate;
                    i=i+1;
                end
                plot(WS_x,R,'r')
                hold on

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rmax
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hline=refline([0 Rmax]);
                set(hline,'Color','y');
                
            legend('Cruise', 'Turn', 'Climb', 'Rmax')
            grid on
            ylim([0 1.5])
            xlim([0 100])
            xlabel('Wing Loading (W/S) [N/m²]')
            ylabel('R (Solar Cell Area/Wing Area)')
            title('Constraint Diagram','FontWeight','Bold')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sizing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Tnight = 24 - mindaylength;
         q=.5*1.1655*cruiseSpeed^2; %0.5*geo_rho*cruiseSpeed^2;
         L_Dmax = 1/2*sqrt(pi*AR*e/CD0); %max endurance
         
         if not (RdesignInput > 0 && W_SInput >0)
            [Wto_S, Rdesign] = ginput(1)
            set(handles.txtW_S,'String',num2str(Wto_S,3));
            set(handles.txtRDesign,'String',num2str(Rdesign,2));
         else
             Wto_S = W_SInput;
             Rdesign = RdesignInput;
         end
         
         % First Iteration
         Wairframe_S = 0.98*g;
         Wsolarcell_S = Rdesign*(SolarCellDensity+Encapsulation)*g;
         CL = sqrt(CD0*pi*e*AR); %SAE Paper - Drag Polar
         CD = CD0 + (CL^2)/(pi*AR*e);
         PlevelFlight = q*CD*cruiseSpeed/(prop_eff*MotorEff*MotorControlEff*gearbox_eff);
         Energy_levelflight = (PlevelFlight/BattEff)*Tnight; %Wh/m²
         Wbattery_S = (Energy_levelflight/(Kbatt))*g; %N/m²
         Wmotor_S = 0.008*PlevelFlight*g; %
         S = (Payload*g) / (Wto_S - Wairframe_S - 1.1*(Wbattery_S + Wmotor_S + Wsolarcell_S));
         Wto = S * Wto_S;
         
         i=0;
         Wto0 = Wto;
         R=1;
         % Battery Resize due to Payload Power not considered in First
         % Iteration
         while (abs(R)>0.001)
            if (S>0)
                PlevelFlight = (q*CD*cruiseSpeed/(prop_eff*MotorEff*MotorControlEff*gearbox_eff))+PayloadPower/S;
            else
                disp('[Error] Increase Wing Loading');
                return;
            end
            Energy_levelflight = (PlevelFlight/BattEff)*Tnight; %Wh/m²
            Wbattery_S = (Energy_levelflight/(Kbatt))*g; %N/m² 
            S = (Payload*g) / (Wto_S - Wairframe_S - 1.1*(Wbattery_S + Wmotor_S + Wsolarcell_S));
            Wto = S * Wto_S;
            R=(Wto-Wto0)/Wto;
            Wto0 = Wto;
            i=i+1;
         end
         msg = ['-> Sizing Iterations: ', num2str(i)];
         disp(msg);
         msg = ['-> Wing Area [m²]: ', num2str(S)];
         disp(msg);
         msg = ['-> Take-off Weight [N]: ', num2str(Wto)];
         disp(msg);
         msg = ['-> Power Required [W]: ~', num2str(PlevelFlight*S)];
         disp(msg);
         msg = ['-> Total mass: ', num2str(Wto/g)];
         disp(msg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Energy Balance Profile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(handles.axes4)
        Hourly = (S*pv_eff*Rdesign).*[HourlyIrradiance];
        
        BattCapacity = Kbatt * (Wbattery_S/g) * S; % [Wh]
        PowerLevelFlight = 0.5*geo_rho*CD*S*cruiseSpeed^3; % P = T*V [W]
        
        % Power Required only Solar supplying
        PowerRequired = PowerLevelFlight/(prop_eff*MotorEff*MotorControlEff*gearbox_eff*mppt_eff) + PayloadPower/(BEC_eff*mppt_eff);
        % Power Required only Battery supplying
        PowerRequiredOnlyBatt = PowerLevelFlight/(prop_eff*MotorEff*MotorControlEff*gearbox_eff*BattEff) + PayloadPower/(BEC_eff*BattEff);
        
        % This piece of code analyzes the transition between only-battery
        % to only-solar
            first=find(Hourly>PowerRequiredOnlyBatt, 1, 'first');
            last=find(Hourly>PowerRequiredOnlyBatt, 1, 'last');
            PwrReq1=ones(1,first).*PowerRequiredOnlyBatt;
            PwrReq2=ones(1,last-first).*PowerRequired;
            PwrReq3=ones(1,24-last).*PowerRequiredOnlyBatt;
            PwrReq=[PwrReq1 PwrReq2 PwrReq3];

        plot(1:24,PwrReq,1:24,Hourly,'LineWidth',2)
        legend('Required','Solar')
        title('Steady Flight Power Profile','FontWeight','Bold')
        xlabel('Hour')
        ylabel('Power (W)')
        hold on
        
        
        TotalEnergyRequired = trapz(PwrReq); 
        SolarEnergyGenerated = trapz(Hourly);
        ExcessSolarPower = trapz(Hourly(first:last))-trapz(PwrReq(first:last));
        BatteryOnlyEnergyMorning = trapz(PwrReq(1:first))-trapz(Hourly(1:first));
        BatteryOnlyEnergyAfternoon = trapz(PwrReq(last:24))-trapz(Hourly(last:24));
        text(2,mean(PowerRequired)/1.5,'Batt Only','FontWeight','Bold','FontSize',8)
        text(2,mean(PowerRequired)/2.8,strcat(num2str(BatteryOnlyEnergyMorning,4),' Wh'),'FontSize',8)
        text(19,mean(PowerRequired)/1.5,'Batt Only','FontWeight','Bold','FontSize',8)
        text(19,mean(PowerRequired)/2.8,strcat(num2str(BatteryOnlyEnergyAfternoon,4),' Wh'),'FontSize',8)
        text(9,(max(Hourly)/2)*1.1,'Excess Energy')
        text(10,max(Hourly)/2,strcat(num2str(ExcessSolarPower,4),' Wh'),'FontSize',8)
        
        if (ExcessSolarPower < (BatteryOnlyEnergyMorning+BatteryOnlyEnergyAfternoon))
            text(4,80,'Excess Solar Power does not recharge battery fully!','FontWeight','Bold','FontSize',7,'Color','r')
        end
        
        xlim([1 24])
        ylim([0 max(Hourly)])
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Software Outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        wingSpan = sqrt(AR*S);
        set(handles.txtWingArea,'String',num2str(S,5));
        set(handles.txtTakeoffWeight,'String',num2str(Wto,5));
        set(handles.wingSpan,'String',num2str(wingSpan,5));
        set(handles.txtPowerReq,'String',num2str(PowerRequired,5));
        set(handles.txtBattWh,'String',num2str(BattCapacity,5));
        set(handles.airframe_mass,'String',num2str(Wairframe_S*S/g,4));
        set(handles.battery_mass,'String',num2str(Wbattery_S*S/g,4));
        set(handles.solarcells_mass,'String',num2str(Wsolarcell_S*S/g,4));
        set(handles.propgroup_mass,'String',num2str(Wmotor_S*S/g,4));
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% CUT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
