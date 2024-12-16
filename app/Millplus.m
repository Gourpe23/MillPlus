classdef Millplus < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        TabGroup                        matlab.ui.container.TabGroup
        MILLTab                         matlab.ui.container.Tab
        GridLayout2                     matlab.ui.container.GridLayout
        CraftedbyDptofMechEngineeringUPVEHULabel  matlab.ui.control.Label
        Label_2                         matlab.ui.control.Label
        Label_3                         matlab.ui.control.Label
        Label_4                         matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        Image5                          matlab.ui.control.Image
        ToolmaterialandmodalparametersTab  matlab.ui.container.Tab
        ModesinYPanel                   matlab.ui.container.Panel
        Image5_4                        matlab.ui.control.Image
        DampingEditField_4              matlab.ui.control.NumericEditField
        DampingEditField_4Label         matlab.ui.control.Label
        StiffnesskNmEditField_4         matlab.ui.control.NumericEditField
        StiffnesskNmEditField_4Label    matlab.ui.control.Label
        FrequencyfHzEditField_4         matlab.ui.control.NumericEditField
        FrequencyfHzEditField_4Label    matlab.ui.control.Label
        ModalparameteresPanel           matlab.ui.container.Panel
        ModesinXPanel                   matlab.ui.container.Panel
        Image5_3                        matlab.ui.control.Image
        DampingEditField_3              matlab.ui.control.NumericEditField
        DampingEditField_3Label         matlab.ui.control.Label
        StiffnesskNmEditField_3         matlab.ui.control.NumericEditField
        StiffnesskNmEditField_3Label    matlab.ui.control.Label
        FrequencyfHzEditField_3         matlab.ui.control.NumericEditField
        FrequencyfHzEditField_3Label    matlab.ui.control.Label
        CuttingcoefficientsPanel        matlab.ui.container.Panel
        KnMPaEditField                  matlab.ui.control.NumericEditField
        KnMPaEditFieldLabel             matlab.ui.control.Label
        MaterialslibraryLabel           matlab.ui.control.Label
        KtMPaEditField                  matlab.ui.control.NumericEditField
        KtMPaEditFieldLabel             matlab.ui.control.Label
        ToolparametersPanel             matlab.ui.container.Panel
        Helixangle60EditField           matlab.ui.control.NumericEditField
        Helixangle60EditFieldLabel      matlab.ui.control.Label
        PLOTMILLButton                  matlab.ui.control.Button
        NofflutesZ10EditField           matlab.ui.control.NumericEditField
        NofflutesZ10EditFieldLabel      matlab.ui.control.Label
        TooldiameterDmmEditField        matlab.ui.control.NumericEditField
        TooldiameterDmmEditFieldLabel   matlab.ui.control.Label
        UIAxes4                         matlab.ui.control.UIAxes
        MillingoperationTab             matlab.ui.container.Tab
        SimulationLabel                 matlab.ui.control.Label
        SimulationprogressGauge         matlab.ui.control.Gauge
        SimulationparametersPanel       matlab.ui.container.Panel
        Periodstosimulate10EditField    matlab.ui.control.NumericEditField
        Periodstosimulate10EditFieldLabel  matlab.ui.control.Label
        IsolinesEditField               matlab.ui.control.NumericEditField
        IsolinesEditFieldLabel          matlab.ui.control.Label
        OthersimulationparametersLabel  matlab.ui.control.Label
        nyEditFieldLabel                matlab.ui.control.Label
        nxEditFieldLabel                matlab.ui.control.Label
        ap_maxEditField                 matlab.ui.control.NumericEditField
        ap_maxEditFieldLabel            matlab.ui.control.Label
        ap_minEditField                 matlab.ui.control.NumericEditField
        ap_minEditFieldLabel            matlab.ui.control.Label
        nyEditField                     matlab.ui.control.NumericEditField
        XLabel                          matlab.ui.control.Label
        nxEditField                     matlab.ui.control.NumericEditField
        AxialdepthofcutmmLabel          matlab.ui.control.Label
        MESHLabel                       matlab.ui.control.Label
        S_maxEditField                  matlab.ui.control.NumericEditField
        S_maxEditFieldLabel             matlab.ui.control.Label
        S_minEditField                  matlab.ui.control.NumericEditField
        S_minEditFieldLabel             matlab.ui.control.Label
        SpindlespeedrangerpmLabel       matlab.ui.control.Label
        SaveoptionsinAStabilitychartsPanel  matlab.ui.container.Panel
        GeneratereportButton            matlab.ui.control.Button
        ClearandreleaseButton_3         matlab.ui.control.Button
        SavesimulationdataButton        matlab.ui.control.Button
        CuttingparametersButtonGroup    matlab.ui.container.ButtonGroup
        SpindlespeednrpmEditFieldLabel  matlab.ui.control.Label
        vcDn1000Label                   matlab.ui.control.Label
        SpindlespeednrpmEditField       matlab.ui.control.EditField
        AxialdepthofcutapmmEditField    matlab.ui.control.EditField
        AxialdepthofcutapmmEditFieldLabel  matlab.ui.control.Label
        RadialdepthofcutaemmEditField   matlab.ui.control.NumericEditField
        RadialdepthofcutaemmEditFieldLabel  matlab.ui.control.Label
        FeedpertoothfzmmZrevEditField   matlab.ui.control.NumericEditField
        FeedpertoothfzmmZrevEditFieldLabel  matlab.ui.control.Label
        TypeofmillingButtonGroup        matlab.ui.container.ButtonGroup
        Image_2                         matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        DownmillingButton               matlab.ui.control.RadioButton
        UpmillingButton_2               matlab.ui.control.RadioButton
        Panel                           matlab.ui.container.Panel
        RUNButton                       matlab.ui.control.Button
        Lamp                            matlab.ui.control.Lamp
        AStabilitychartsTab             matlab.ui.container.Tab
        Panel_3                         matlab.ui.container.Panel
        UIAxes2                         matlab.ui.control.UIAxes
        RaButton                        matlab.ui.control.Button
        p2pXYButton                     matlab.ui.control.Button
        p2pYButton                      matlab.ui.control.Button
        p2pXButton                      matlab.ui.control.Button
        Panel_2                         matlab.ui.container.Panel
        UIAxes1                         matlab.ui.control.UIAxes
        BDynamicforcesTab               matlab.ui.container.Tab
        CuttingandsimulationparametersPanel  matlab.ui.container.Panel
        ClearandreleaseButton           matlab.ui.control.Button
        ChatterButton                   matlab.ui.control.Button
        RUNButton_2                     matlab.ui.control.Button
        PeriodstosimulateT10EditField   matlab.ui.control.NumericEditField
        PeriodstosimulateT10EditFieldLabel  matlab.ui.control.Label
        SpindlespeednrpmEditField_2     matlab.ui.control.NumericEditField
        SpindlespeednrpmEditField_2Label  matlab.ui.control.Label
        AxialdepthofcutapmmEditField_2  matlab.ui.control.NumericEditField
        AxialdepthofcutapmmEditField_2Label  matlab.ui.control.Label
        RadialdepthofcutaemmEditField_2  matlab.ui.control.NumericEditField
        RadialdepthofcutaemmEditField_2Label  matlab.ui.control.Label
        FeedpertoothfzmmZrevEditField_2  matlab.ui.control.NumericEditField
        FeedpertoothfzmmZrevEditField_2Label  matlab.ui.control.Label
        UIAxes3                         matlab.ui.control.UIAxes
        UIAxes3_3                       matlab.ui.control.UIAxes
        UIAxes3_2                       matlab.ui.control.UIAxes
        COptimizationtoolTab            matlab.ui.container.Tab
        Panel_4                         matlab.ui.container.Panel
        Panel_7                         matlab.ui.container.Panel
        cmminLabel_2                    matlab.ui.control.Label
        ProductivityimprovedbyEditField_9  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_9  matlab.ui.control.Label
        ProductivityimprovedbyEditField_8  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_8  matlab.ui.control.Label
        mLabel_2                        matlab.ui.control.Label
        ProductivityimprovedbyEditField_7  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_7  matlab.ui.control.Label
        ProductivityimprovedbyEditField_6  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_6  matlab.ui.control.Label
        Panel_6                         matlab.ui.container.Panel
        SurfroughnessimprovedbyEditField  matlab.ui.control.NumericEditField
        SurfroughnessimprovedbyEditFieldLabel  matlab.ui.control.Label
        mLabel                          matlab.ui.control.Label
        ProductivityimprovedbyEditField_5  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_5  matlab.ui.control.Label
        ProductivityimprovedbyEditField_4  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_4  matlab.ui.control.Label
        Label_6                         matlab.ui.control.Label
        Panel_5                         matlab.ui.container.Panel
        ProductivityimprovedbyEditField  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel  matlab.ui.control.Label
        cmminLabel                      matlab.ui.control.Label
        ProductivityimprovedbyEditField_3  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_3  matlab.ui.control.Label
        ProductivityimprovedbyEditField_2  matlab.ui.control.NumericEditField
        ProductivityimprovedbyEditFieldLabel_2  matlab.ui.control.Label
        Label_5                         matlab.ui.control.Label
        UIAxes7                         matlab.ui.control.UIAxes
        OptimizationobjectivePanel      matlab.ui.container.Panel
        ChoosecriterionLabel            matlab.ui.control.Label
        ClearfigureButton               matlab.ui.control.Button
        MaximumRamSlider                matlab.ui.control.Slider
        MaximumRamSliderLabel           matlab.ui.control.Label
        MinimumMRRmmminSlider           matlab.ui.control.Slider
        MRRmm3minLabel                  matlab.ui.control.Label
        EditField_4                     matlab.ui.control.NumericEditField
        EditField_3                     matlab.ui.control.NumericEditField
        EditField_2                     matlab.ui.control.NumericEditField
        EditField                       matlab.ui.control.NumericEditField
        MaximumRamSlider_2              matlab.ui.control.Slider
        MaximumalowedRamumLabel         matlab.ui.control.Label
        FxymeanFxySlider                matlab.ui.control.Slider
        FxymeanFxySliderLabel           matlab.ui.control.Label
        StartingpointPanel              matlab.ui.container.Panel
        SpindlespeednrpmEditField_6     matlab.ui.control.NumericEditField
        SpindlespeednrpmEditField_6Label  matlab.ui.control.Label
        AxialdepthofcutapmmEditField_6  matlab.ui.control.NumericEditField
        AxialdepthofcutapmmEditField_6Label  matlab.ui.control.Label
        RadialdepthofcutaemmEditField_4  matlab.ui.control.NumericEditField
        RadialdepthofcutaemmEditField_4Label  matlab.ui.control.Label
        FeedpertoothfzmmZrevEditField_4  matlab.ui.control.NumericEditField
        FeedpertoothfzmmZrevEditField_4Label  matlab.ui.control.Label
        OptimumpointPanel               matlab.ui.container.Panel
        SpindlespeednrpmEditField_7     matlab.ui.control.NumericEditField
        SpindlespeednrpmEditField_7Label  matlab.ui.control.Label
        AxialdepthofcutapmmEditField_7  matlab.ui.control.NumericEditField
        AxialdepthofcutapmmEditField_7Label  matlab.ui.control.Label
        RadialdepthofcutaemmEditField_5  matlab.ui.control.NumericEditField
        RadialdepthofcutaemmEditField_5Label  matlab.ui.control.Label
        FeedpertoothfzmmZrevEditField_5  matlab.ui.control.NumericEditField
        FeedpertoothfzmmZrevEditField_5Label  matlab.ui.control.Label
        OPTIMUMButton                   matlab.ui.control.Button
        Label_7                         matlab.ui.control.Label
        StartingpointDropDown           matlab.ui.control.DropDown
    end

properties (Access = public)
    Diameter          
    Z                  
    Helixangle 
    SwitchTypeMilling
    ptpx               
    ptpy
    Pot
    Ra
    Qc
    Fx
    Fy
    Time   
    omg_x 
    k_x 
    zta_x
    omg_y
    k_y 
    zta_y  
    Fx_mean
    Fy_mean
    optimal_n_index
    optimal_ap_index
    optimal_n
    optimal_ap
    SimulationResults   
    SpindleSpeedMatrix
    AxialDepthOfCutMatrix
    FeedPerTooth
    StartingPointMarker
    StartingPointText
    OptimalPointMarker
    OptimalPointText
    ForceRatioMarkers
    SurfaceRoughnessMarkers
    RaMarkers
    OptimalPoint
end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: FrequencyfHzEditField_3
        function FrequencyfHzEditFieldValueChanged(app, event)
            
        end

        % Value changed function: StiffnesskNmEditField_3
        function StiffnesskNmEditFieldValueChanged(app, event)
            
        end

        % Value changed function: DampingEditField_3
        function DampingEditFieldValueChanged(app, event)
            
        end

        % Value changed function: FrequencyfHzEditField_4
        function FrequencyfHzEditField_2ValueChanged(app, event)
           
        end

        % Value changed function: StiffnesskNmEditField_4
        function StiffnesskNmEditField_2ValueChanged(app, event)
           
        end

        % Value changed function: DampingEditField_4
        function DampingEditField_2ValueChanged(app, event)
         
        end

        % Value changed function: TooldiameterDmmEditField
        function TooldiameterDmmEditFieldValueChanged(app, event)
            %D = app.TooldiameterDmmEditField.Value;
        end

        % Value changed function: NofflutesZ10EditField
        function NofflutesZ10EditFieldValueChanged(app, event)
            %Z = app.NofflutesZ10EditField.Value;
        end

        % Value changed function: Helixangle60EditField
        function Helixangle60EditFieldValueChanged(app, event)
            %Helixangle = app.Helixangle60EditField.Value;
        end

        % Button pushed function: PLOTMILLButton
        function PLOTMILLButtonPushed(app, event)
            D = app.TooldiameterDmmEditField.Value;
            Z = app.NofflutesZ10EditField.Value;
            beta = app.Helixangle60EditField.Value;

            H = 2*D;         
            R = D/2;          
            n_points = 100;            
            cla(app.UIAxes4);
            hold(app.UIAxes4, 'on');
            theta = linspace(0, 2*pi, n_points);  
            z_cyl = linspace(0, H, n_points);   
            [theta_grid, Z_cyl] = meshgrid(theta, z_cyl);

            X_cyl = R * cos(theta_grid);  
            Y_cyl = R * sin(theta_grid);  
            cylinder_surface = surf(app.UIAxes4, X_cyl, Y_cyl, Z_cyl, 'EdgeColor', 'none');
            material(cylinder_surface, 'metal');
            set(cylinder_surface, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1);
            light(app.UIAxes4, 'Position', [1 0 1], 'Style', 'infinite');
            light(app.UIAxes4, 'Position', [-1 0 1], 'Style', 'infinite');
            lighting(app.UIAxes4, 'gouraud');
            camlight(app.UIAxes4, 'headlight');
            axis(app.UIAxes4, 'equal');
 
            for i = 1:Z
                theta0 = (i - 1) * 2 * pi / Z;
                t = linspace(0, H, n_points);  % Altura
                x_helix = R * cos(theta0 + t * tand(beta) / R);
                y_helix = R * sin(theta0 + t * tand(beta) / R);
                plot3(app.UIAxes4, x_helix, y_helix, t, 'r', 'LineWidth', 3);
            end

            xlabel(app.UIAxes4, 'x (mm)');
            ylabel(app.UIAxes4, 'y (mm)');
            zlabel(app.UIAxes4, 'z (mm)');
            title(app.UIAxes4, 'Solid mill');
            grid(app.UIAxes4, 'on');
            view(app.UIAxes4, 3);
            light(app.UIAxes4, 'Position', [0 1 0], 'Style', 'local');
            light(app.UIAxes4, 'Position', [0 -1 1], 'Style', 'local');
        end

        % Value changed function: FeedpertoothfzmmZrevEditField
        function FeedpertoothfzmmZrevEditFieldValueChanged(app, event)
       
        end

        % Value changed function: RadialdepthofcutaemmEditField
        function RadialdepthofcutaemmEditFieldValueChanged(app, event)
        
        end

        % Callback function
        function InmersionlevelaeDEditFieldValueChanged(app, event)

        end

        % Selection changed function: CuttingparametersButtonGroup
        function updateInmersionLevel(app, event)
       
        end

        % Value changed function: KtMPaEditField
        function KtMPaEditFieldValueChanged(app, event)
            Ktc = app.KtMPaEditField.Value; 
        end

        % Value changed function: KnMPaEditField
        function KnMPaEditFieldValueChanged(app, event)
            Knc = app.KnMPaEditField.Value;
        end

        % Value changed function: S_minEditField
        function S_minEditFieldValueChanged(app, event)
            Smin = app.S_minEditField.Value;
        end

        % Value changed function: S_maxEditField
        function S_maxEditFieldValueChanged(app, event)
            Smax = app.S_maxEditField.Value;
        end

        % Value changed function: nxEditField
        function nxEditFieldValueChanged(app, event)
            nx = app.nxEditField.Value;
        end

        % Value changed function: ap_minEditField
        function ap_minEditFieldValueChanged(app, event)
            apmin = app.ap_minEditField.Value;
        end

        % Value changed function: ap_maxEditField
        function ap_maxEditFieldValueChanged(app, event)
            apmax = 1e-3*app.ap_maxEditField.Value;
        end

        % Value changed function: nyEditField
        function nyEditFieldValueChanged(app, event)
            ny = app.nyEditField.Value;
        end

        % Value changed function: IsolinesEditField
        function IsolinesEditFieldValueChanged(app, event)
            iso = app.IsolinesEditField.Value;
        end

        % Value changed function: Periodstosimulate10EditField
        function Periodstosimulate10EditFieldValueChanged(app, event)
            revs = app.Periodstosimulate10EditField.Value;
        end

        % Button pushed function: RUNButton
        function RUNButtonPushed(app, event)
            
            app.Panel_5.Visible = 'off';
            app.Panel_6.Visible = 'off';
            app.Panel_7.Visible = 'off';

            app.SimulationprogressGauge.Value = 0; 
            app.Lamp.Color = [1 0 0]; 
            drawnow;

            app.FeedpertoothfzmmZrevEditField.Enable = 'off';       
            app.RadialdepthofcutaemmEditField.Enable = 'off';       
            app.S_minEditField.Enable = 'off';                      
            app.S_maxEditField.Enable = 'off';                      
            app.ap_minEditField.Enable = 'off';                     
            app.ap_maxEditField.Enable = 'off';                     
            app.nxEditField.Enable = 'off';                         
            app.nyEditField.Enable = 'off';                         
            app.IsolinesEditField.Enable = 'off';  
            app.Periodstosimulate10EditField.Enable = 'off';  
            
            % Bloquear campos en Dynamic forces
            app.PeriodstosimulateT10EditField.Value = app.Periodstosimulate10EditField.Value;
            app.RadialdepthofcutaemmEditField_2.Value = app.RadialdepthofcutaemmEditField.Value;
            app.FeedpertoothfzmmZrevEditField_2.Value = app.FeedpertoothfzmmZrevEditField.Value;
            
            app.PeriodstosimulateT10EditField.Enable = 'off';       
            app.RadialdepthofcutaemmEditField_2.Enable = 'off';     
            app.FeedpertoothfzmmZrevEditField_2.Enable = 'off';    
    
            %% Modal parameters
            omg_x = app.FrequencyfHzEditField_3.Value * 2 * pi;  
            k_x = app.StiffnesskNmEditField_3.Value;           
            zta_x = app.DampingEditField_3.Value;           
            m_x = k_x / omg_x^2;                           
            omg_y = app.FrequencyfHzEditField_4.Value * 2 * pi;  
            k_y = app.StiffnesskNmEditField_4.Value;           
            zta_y = app.DampingEditField_4.Value;          
            m_y = k_y / omg_y^2;                        

            %% Tool parameters
            Ktc = 1e6*app.KtMPaEditField.Value;
            Knc = 1e6*app.KnMPaEditField.Value;
            D = 1e-3*app.TooldiameterDmmEditField.Value;
            gamma = (pi/180)*app.Helixangle60EditField.Value;
            fz = 1e-3*app.FeedpertoothfzmmZrevEditField.Value;
            zn = app.NofflutesZ10EditField.Value; 
            ae = 1e-3*app.RadialdepthofcutaemmEditField.Value;
            aD = ae/D;
    
            typeMilling = app.TypeofmillingButtonGroup.SelectedObject;
            if typeMilling == app.UpmillingButton_2
                updw = 1;  % down milling
            else
                updw = -1; % up milling
            end
    
            if updw == 1, phist = 0; phiex = acos(1-2*aD);
            elseif updw == -1, phist = acos(2*aD-1);  phiex = pi;
            end
            
            %% Simulations parameters
            iso = app.IsolinesEditField.Value;
            drev = 1/360; 
            revs = app.Periodstosimulate10EditField.Value;
            steps_rev = round(1/zn/drev)*zn;
            Smin = app.S_minEditField.Value;
            Smax = app.S_maxEditField.Value;
            apmin = 1e-3*app.ap_minEditField.Value;
            apmax = 1e-3*app.ap_maxEditField.Value;
            nx = app.nxEditField.Value;
            ny = app.nyEditField.Value;
            n_s = linspace(Smin,Smax,nx); 
            ap_s = linspace(apmin,apmax,ny);
            varx = zeros(length(n_s),length(ap_s)); 
            vary = varx; ptpx = varx; ptpy = varx; Ra = varx;
            
            %% Solucion por euler
            for in = 1:length(n_s)
              app.SimulationprogressGauge.Value = (in / length(n_s)) * 100;
              n = n_s(in);
              for iap = 1:length(ap_s)
                [F, z] = app.simulation(n,ap_s(iap),revs,drev,fz,D,zn,gamma,phist,phiex,omg_x,omg_y,zta_x,zta_y,m_x,m_y,Ktc,Knc);
                Fx = F(:,1); 
                Fy = F(:,2);
                x = z(:,1); 
                y = z(:,2);
    
                ind_rev_ra = steps_rev*10 + 1:steps_rev/zn:length(x);
                Ra(in, iap) = app.roughness(x(ind_rev_ra), y(ind_rev_ra), fz, D, revs - 10, zn); % roughness Ra [um]
                Qc(in, iap) = (1000*ae) * (1000 * ap_s(iap)) * (1000*fz) * zn * n;  % en mm³/min
            
                varx(in, iap) = var(x);
                vary(in, iap) = var(y);
                ptpx(in, iap) = max(Fx) - min(Fx);
                ptpy(in, iap) = max(Fy) - min(Fy);
                
                Fx_mean(in, iap) = mean(Fx);
                Fy_mean(in, iap) = mean(Fy);
                
                force_xy(in, iap) = sqrt(ptpx(in, iap)^2 + ptpy(in, iap)^2);
                force_mean(in, iap) = sqrt((Fx_mean(in, iap))^2+(Fy_mean(in, iap))^2);
               
                Pot(in, iap) = (1e-3*pi*n*D/60)*sqrt((mean(Fx))^2+(mean(Fy))^2);  % power P [kW]
              end
              pause(0.01);
            end
    
            Time = linspace (0,steps_rev*revs,length(F(:,1)));
    
            app.ptpx = ptpx;
            app.ptpy = ptpy;
            app.Pot = Pot;
            app.Ra = Ra;
            app.Qc = Qc;
            app.Fx = Fx;
            app.Fy = Fy;
            app.Time = Time;
    
            % Fxy-left figure
            contourf(app.UIAxes1, n_s, 1000*ap_s, sqrt(ptpx'.^2 + ptpy'.^2), iso,'LineColor', 'none');
            colormap(app.UIAxes1, 'jet');
            colorbar(app.UIAxes1);
            xlabel(app.UIAxes1, 'Spindle speed [rpm]'); ylabel(app.UIAxes1, 'Axial depth of cut [mm]');
            title(app.UIAxes1, 'Peak-to-peak Force XY [N]');
            hold(app.UIAxes1, 'on');
            [C, h] = contour(app.UIAxes1, n_s, 1000*ap_s, sqrt(ptpx'.^2 + ptpy'.^2), 10, 'LineColor', 'w');
            clabel(C, h, 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
            hold(app.UIAxes1, 'off');
            
            % Multi-right figure
            colorbar(app.UIAxes2, 'off');    
            [C1, h1] = contour(app.UIAxes2, n_s, 1000*ap_s, sqrt(app.ptpx'.^2 + app.ptpy'.^2), 20, 'LineColor', [0.5, 0.5, 0.5]);
            clabel(C1, h1, 'Color', [0.5, 0.5, 0.5], 'FontSize', 8, 'FontWeight', 'bold');  
    
            hold(app.UIAxes2, 'on');
            [C2, h2] = contour(app.UIAxes2, n_s, 1000*ap_s, app.Pot', 10, 'LineColor', [0.5, 0, 0]);
            clabel(C2, h2, 'Color', [0.5, 0, 0], 'FontSize', 8, 'FontWeight', 'bold');
    
            [C3, h3] = contour(app.UIAxes2, n_s, 1000*ap_s, 1e-3*app.Qc', 10, 'LineColor', 'r');
            clabel(C3, h3, 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
    
            [C4, h4] = contour(app.UIAxes2, n_s, 1000*ap_s, Ra', 10, 'LineColor', 'g');
            clabel(C4, h4, 'Color', 'g', 'FontSize', 8, 'FontWeight', 'bold');
    
            xlabel(app.UIAxes2, 'Spindle speed [rpm]');
            ylabel(app.UIAxes2, 'Axial depth of cut [mm]');
            legend([h1, h2, h3, h4], {'Force XY [N]', 'Power [kW]', 'MRR (Qc = vf x ap x ae) [cm^3/min]', 'Ra [\mum]'}, 'Location', 'best');
            grid(app.UIAxes2, 'on');
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.XMinorGrid = 'on';
            app.UIAxes2.YMinorGrid = 'on';
            hold(app.UIAxes2, 'off');
   
            cla(app.UIAxes7);
            grid(app.UIAxes7, 'on');  
            app.UIAxes7.XGrid = 'on';  
            app.UIAxes7.YGrid = 'on';  
            app.UIAxes7.XMinorGrid = 'on';  
            app.UIAxes7.YMinorGrid = 'on';
                                
            [C1, h1] = contour(app.UIAxes7, n_s, 1e3*ap_s, sqrt(app.ptpx'.^2 + app.ptpy'.^2), 20, 'LineColor', 'k');
            clabel(C1, h1, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');
            xlabel(app.UIAxes7, 'Spindle speed [rpm]');
            ylabel(app.UIAxes7, 'Axial depth of cut [mm]');
    
            app.SimulationprogressGauge.Value = 100;
            app.Lamp.Color = [0 1 0];
            app.FeedpertoothfzmmZrevEditField_4.Value = app.FeedpertoothfzmmZrevEditField.Value;
            app.RadialdepthofcutaemmEditField_4.Value = app.RadialdepthofcutaemmEditField.Value;
        
            ap_min = app.ap_minEditField.Value; % mm
            ap_max = app.ap_maxEditField.Value; % mm
            s_min = app.S_minEditField.Value;
            s_max = app.S_maxEditField.Value;
            app.AxialdepthofcutapmmEditField_6.Value = ((ap_min + ap_max) / 2);
            app.SpindlespeednrpmEditField_6.Value = ((s_min + s_max) / 2);

            app.Panel_5.Visible = 'off';
            app.Panel_6.Visible = 'off';
            app.Panel_7.Visible = 'off';

            % Sliders initialization by default
            app.FxymeanFxySlider.Visible = 'on';          
            app.FxymeanFxySlider.Enable = 'on';   
            app.EditField.Visible = 'on';         
            app.EditField.Enable = 'off';              
            app.MaximumRamSlider_2.Visible = 'off';  
            app.EditField_2.Visible = 'off'; 
            app.MinimumMRRmmminSlider.Visible = 'off';     
            app.EditField_3.Visible = 'off';         
            app.MaximumRamSlider.Visible = 'off';
            app.EditField_4.Visible = 'off';

            choice = uiconfirm(app.UIFigure,'¿Do you want to start the Optimization tool tab?', 'Confirmation', 'Options', {'Yes', 'No'}, 'DefaultOption', 1,'CancelOption', 2);
            % Check the user's choice
            if strcmp(choice, 'Yes')
                app.TabGroup.SelectedTab = app.COptimizationtoolTab;
            end
            app.ProductivityimprovedbyEditField.Visible = false;
            app.ProductivityimprovedbyEditField_2.Visible = false;
            app.ProductivityimprovedbyEditField_3.Visible = false;
            app.ptpx = ptpx;        
            app.ptpy = ptpy;       
            app.Fx_mean = Fx_mean;
            app.Fy_mean = Fy_mean;
        end
        
        function [F, z] = simulation(app,sps,ap,revs,drev,fz,D,zn,gamma,phist,phiex,omg_x,omg_y,zta_x,zta_y,m_x,m_y,Ktc,Knc)
            nmd = length(omg_x);
            steps_rev = round(1/zn/drev)*zn;
            dphi = 2*pi/steps_rev; 
            dz_phi = D*dphi/(2*tan(gamma));
            if ap<dz_phi
              dz = ap;
            else
              dphi = 2*ap*tan(gamma)/D/ceil(ap/dz_phi);
              steps_rev = round(2*pi/zn/dphi)*zn;
              dz = D*dphi/(2*tan(gamma));
            end
            steps_axial = round(ap/dz);
            teeth = 1+(0:zn-1)*steps_rev/zn;
            phi = dphi*(0:steps_rev*revs-1);
            dt = 60/(steps_rev*sps);
            
            % inicialización
            sup_pre = zeros(steps_axial,length(phi)); sup_act = sup_pre; Fx = zeros(1,length(phi)); Fy = Fx;
            h = zeros(steps_axial, length(phi) - 1);
            ux = zeros(length(phi), 2, nmd); uy = ux; x = zeros(length(phi), 2); y = x;
            for i_phi = 1:length(phi)-1
              for iN = 1:zn
                for ia = 1:steps_axial
                  phi_tth_dsk = mod(phi(teeth(iN)) + dphi*(i_phi-ia+1), 2*pi);
                  i_tth_dsk = round(phi_tth_dsk/dphi);
                  if ~i_tth_dsk, i_tth_dsk = steps_rev; end
                  if (phi_tth_dsk>phist)*(phi_tth_dsk<phiex)
                    sup_act(ia,i_phi) = x(i_phi,1)*sin(phi_tth_dsk) - y(i_phi,1)*cos(phi_tth_dsk);
                    h(ia,i_phi) = fz*sin(phi_tth_dsk) - sup_pre(ia,i_tth_dsk) + sup_act(ia,i_phi);
                    if h(ia,i_phi) > 0
                      sup_pre(ia,i_tth_dsk) = sup_act(ia,i_phi);
                    else
                      h(ia,i_phi) = 0;
                      sup_pre(ia,i_tth_dsk) = sup_pre(ia,i_tth_dsk) - fz*sin(phi_tth_dsk);
                    end
                    Fn = Knc*dz*h(ia,i_phi);
                    Ft = Ktc*dz*h(ia,i_phi);
                    Fx(i_phi) = Fx(i_phi) + Ft*cos(phi_tth_dsk) + Fn*sin(phi_tth_dsk); %sobre herramienta
                    Fy(i_phi) = Fy(i_phi) + Ft*sin(phi_tth_dsk) - Fn*cos(phi_tth_dsk);
                  end
                end
              end
              for im = 1:nmd
                z_eul = app.euler_milling2([ux(i_phi,:,im)', uy(i_phi,:,im)'],[-Fx(i_phi), -Fy(i_phi)],[m_x(im), m_y(im)],[zta_x(im), zta_y(im)],[omg_x(im), omg_y(im)], dt);
                ux(i_phi+1,:,im) = z_eul(:,1)';
                uy(i_phi+1,:,im) = z_eul(:,2)';
                x(i_phi+1,:) = x(i_phi+1,:) + ux(i_phi+1,:,im);
                y(i_phi+1,:) = y(i_phi+1,:) + uy(i_phi+1,:,im);
              end
            end
            F = [Fx;Fy]';
            z = [x(:,1), y(:,1)];
        end
               
        %% Integración por Euler
        function z_out = euler_milling2(~,z,F,m,zeta,omega_n,dt)
            z0 = z(1,:); z0p = z(2,:);
            zpp = F./m-2*zeta.*omega_n.*z0p - omega_n.^2.*z0;
            zp = z0p + zpp*dt;
            z = z0 + zp*dt;
            z_out = [z; zp];
        end
        
        function Ra_ = roughness(~,x,y,fz,D,revs,zn)
            R = D/2;
            theta_fz = asin(fz/2/R); samples_prf = 11;
            phi_sur = fliplr(linspace(pi-theta_fz,pi+theta_fz,samples_prf));
            x_prf = R*sin(phi_sur); y_prf = -sqrt(R^2 - (x_prf).^2) + R;
            x_srf = zeros(1, samples_prf*revs*zn - samples_prf); y_srf = x_srf;
            for ii = 1:(revs*zn-1)
            y_srf(1,1+samples_prf*(ii-1):samples_prf*ii) = linspace(y(ii),y(ii+1),samples_prf) + y_prf;
            x_srf(1,1+samples_prf*(ii-1):samples_prf*ii) = linspace(x(ii),x(ii+1),samples_prf) + x_prf + fz*(ii-1);
            end
            Ra_ = sum(abs(y_srf - mean(y_srf)))/length(y_srf)*1E6;
        end

        % Button pushed function: p2pXButton
        function p2pXButtonPushed(app, event)
                % Clear the current plot in UIAxes1
                cla(app.UIAxes1);
                legend(app.UIAxes1, 'off');
                Smin = app.S_minEditField.Value;
                Smax = app.S_maxEditField.Value;
                apmin = 1e-3*app.ap_minEditField.Value;
                apmax = 1e-3*app.ap_maxEditField.Value;
                nx = app.nxEditField.Value;
                ny = app.nyEditField.Value;
                iso = app.IsolinesEditField.Value;
                n_s = linspace(Smin,Smax,nx); 
                ap_s = linspace(apmin,apmax,ny);
               
                contourf(app.UIAxes1, n_s, 1000*ap_s, app.ptpx', iso,'LineColor', 'none');
                colormap(app.UIAxes1, 'jet');
                colorbar(app.UIAxes1);
                xlabel(app.UIAxes1, 'Spindle speed [rpm]'); ylabel(app.UIAxes1, 'Axial depth of cut [mm]');
                title(app.UIAxes1, 'Peak-to-peak Force X');
                hold(app.UIAxes1, 'on');
                [C, h] = contour(app.UIAxes1, n_s, 1000*ap_s, app.ptpx', 10, 'LineColor', 'w');
                clabel(C, h, 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold'); 
                hold(app.UIAxes1, 'off');
        end

        % Button pushed function: p2pYButton
        function p2pYButtonPushed(app, event)
                % Clear the current plot in UIAxes1
                cla(app.UIAxes1);
                legend(app.UIAxes1, 'off');
                Smin = app.S_minEditField.Value;
                Smax = app.S_maxEditField.Value;
                apmin = 1e-3*app.ap_minEditField.Value;
                apmax = 1e-3*app.ap_maxEditField.Value;
                nx = app.nxEditField.Value;
                ny = app.nyEditField.Value;
                iso = app.IsolinesEditField.Value;
                n_s = linspace(Smin,Smax,nx); 
                ap_s = linspace(apmin,apmax,ny);

                contourf(app.UIAxes1, n_s, 1000*ap_s, app.ptpy', iso,'LineColor', 'none');
                colormap(app.UIAxes1, 'jet');
                colorbar(app.UIAxes1);
                xlabel(app.UIAxes1, 'Spindle speed [rpm]'); ylabel(app.UIAxes1, 'Axial depth of cut [mm]');
                title(app.UIAxes1, 'Peak-to-peak Force Y');
                hold(app.UIAxes1, 'on');
                [C, h] = contour(app.UIAxes1, n_s, 1000*ap_s, app.ptpy', 10, 'LineColor', 'w');
                clabel(C, h, 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
                hold(app.UIAxes1, 'off');
        end

        % Button pushed function: p2pXYButton
        function p2pXYButtonPushed(app, event)
                % Clear the current plot in UIAxes1
                cla(app.UIAxes1);
                legend(app.UIAxes1, 'off');
                Smin = app.S_minEditField.Value;
                Smax = app.S_maxEditField.Value;
                apmin = 1e-3*app.ap_minEditField.Value;
                apmax = 1e-3*app.ap_maxEditField.Value;
                nx = app.nxEditField.Value;
                ny = app.nyEditField.Value;
                iso = app.IsolinesEditField.Value;
                n_s = linspace(Smin,Smax,nx); 
                ap_s = linspace(apmin,apmax,ny);

                contourf(app.UIAxes1, n_s, 1000*ap_s, sqrt(app.ptpx'.^2 + app.ptpy'.^2), iso,'LineColor', 'none');
                colormap(app.UIAxes1, 'jet');
                colorbar(app.UIAxes1);
                xlabel(app.UIAxes1, 'Spindle speed [rpm]'); ylabel(app.UIAxes1, 'Axial depth of cut [mm]');
                title(app.UIAxes1, 'Peak-to-peak Force XY');
                hold(app.UIAxes1, 'on');
                [C, h] = contour(app.UIAxes1, n_s, 1000*ap_s, sqrt(app.ptpx'.^2 + app.ptpy'.^2), 10, 'LineColor', 'w');
                clabel(C, h, 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
        end

        % Callback function
        function IsolinesSwitchValueChanged(app, event)

        end

        % Value changed function: FeedpertoothfzmmZrevEditField_2
        function FeedpertoothfzmmZrevEditField_2ValueChanged(app, event)
            value = app.FeedpertoothfzmmZrevEditField_2.Value;       
        end

        % Value changed function: RadialdepthofcutaemmEditField_2
        function RadialdepthofcutaemmEditField_2ValueChanged(app, event)
            value = app.RadialdepthofcutaemmEditField_2.Value; 
        end

        % Value changed function: AxialdepthofcutapmmEditField_2
        function AxialdepthofcutapmmEditField_2ValueChanged(app, event)
            value = app.AxialdepthofcutapmmEditField_2.Value; 
        end

        % Value changed function: SpindlespeednrpmEditField_2
        function SpindlespeednrpmEditField_2ValueChanged(app, event)
            value = app.SpindlespeednrpmEditField_2.Value;
        end

        % Button pushed function: RUNButton_2
        function RUNButton_2Pushed(app, event)
        app.FeedpertoothfzmmZrevEditField_2.Enable = 'off';      
        app.RadialdepthofcutaemmEditField_2.Enable = 'off';
        app.AxialdepthofcutapmmEditField_2.Enable = 'off';
        app.SpindlespeednrpmEditField_2.Enable = 'off';
        app.PeriodstosimulateT10EditField.Enable = 'off';

        %% Modal parameters
        omg_x = app.FrequencyfHzEditField_3.Value * 2 * pi;
        k_x = app.StiffnesskNmEditField_3.Value;          
        zta_x = app.DampingEditField_3.Value;        
        m_x = k_x / omg_x^2;                         
        omg_y = app.FrequencyfHzEditField_4.Value * 2 * pi;
        k_y = app.StiffnesskNmEditField_4.Value;  
        zta_y = app.DampingEditField_4.Value;            
        m_y = k_y / omg_y^2;

        %% Tool parameters
        Ktc = 1e6*app.KtMPaEditField.Value;
        Knc = 1e6*app.KnMPaEditField.Value;
        D = 1e-3*app.TooldiameterDmmEditField.Value;
        gamma = (pi/180)*app.Helixangle60EditField.Value;
        fz = 1e-3*app.FeedpertoothfzmmZrevEditField.Value;
        zn = app.NofflutesZ10EditField.Value; 
        ae = 1e-3*app.RadialdepthofcutaemmEditField.Value;
        aD = ae/D;

        typeMilling = app.TypeofmillingButtonGroup.SelectedObject;
    
        if typeMilling == app.UpmillingButton_2
            updw = 1;  % down milling
        else
            updw = -1; % up milling
        end

        if updw == 1, phist = 0; phiex = acos(1-2*aD);
        elseif updw == -1, phist = acos(2*aD-1);  phiex = pi;
        end
        
        %% Simulations parameters
        iso = app.IsolinesEditField.Value;
        drev = 1/360;
        revs = app.Periodstosimulate10EditField.Value;
            steps_rev = round(1/zn/drev)*zn;
            n_sp = app.SpindlespeednrpmEditField_2.Value;
            ap_sp = 1e-3*app.AxialdepthofcutapmmEditField_2.Value;          
            [F, z] = app.simulation(n_sp,ap_sp,revs,drev,fz,D,zn,gamma,phist,phiex,omg_x,omg_y,zta_x,zta_y,m_x,m_y,Ktc,Knc);
            Fx = F(:,1); Fy = F(:,2);
            x = z(:,1); y = z(:,2);
            Time = linspace (0,revs*60/n_sp,length(F(:,1)));

            app.Fx = Fx;
            app.Fy = Fy;
            app.Time = Time;

            function [F, z] = simulation(app, sps, ap, revs, drev, fz, D, zn, gamma, phist, phiex, omg_x, omg_y, zta_x, zta_y, m_x, m_y, Ktc, Knc)
            nmd = length(omg_x);            
            dphi = 2*pi/steps_rev;
            dz_phi = D*dphi/(2*tan(gamma));
            if ap<dz_phi
            dz = ap;
            else
            dphi = 2*ap*tan(gamma)/D/ceil(ap/dz_phi);
            steps_rev = round(2*pi/zn/dphi)*zn;
            dz = D*dphi/(2*tan(gamma));
            end
            steps_axial = round(ap/dz);
            teeth = 1+(0:zn-1)*steps_rev/zn; 
            phi = dphi*(0:steps_rev*revs-1);
            dt = 60/(steps_rev*sps);
            sup_pre = zeros(steps_axial,length(phi)); sup_act = sup_pre; Fx = zeros(1,length(phi)); Fy = Fx;
            h = zeros(steps_axial, length(phi) - 1);
            ux = zeros(length(phi), 2, nmd); uy = ux; x = zeros(length(phi), 2); y = x;
            for i_phi = 1:length(phi)-1
            for iN = 1:zn
                for ia = 1:steps_axial
                phi_tth_dsk = mod(phi(teeth(iN)) + dphi*(i_phi-ia+1), 2*pi);
                i_tth_dsk = round(phi_tth_dsk/dphi);
                if ~i_tth_dsk, i_tth_dsk = steps_rev; end
                if (phi_tth_dsk>phist)*(phi_tth_dsk<phiex)
                    sup_act(ia,i_phi) = x(i_phi,1)*sin(phi_tth_dsk) - y(i_phi,1)*cos(phi_tth_dsk);
                    h(ia,i_phi) = fz*sin(phi_tth_dsk) - sup_pre(ia,i_tth_dsk) + sup_act(ia,i_phi);
                    if h(ia,i_phi) > 0
                    sup_pre(ia,i_tth_dsk) = sup_act(ia,i_phi);
                    else
                    h(ia,i_phi) = 0;
                    sup_pre(ia,i_tth_dsk) = sup_pre(ia,i_tth_dsk) - fz*sin(phi_tth_dsk);
                    end
                    Fn = Knc*dz*h(ia,i_phi);
                    Ft = Ktc*dz*h(ia,i_phi);
                    Fx(i_phi) = Fx(i_phi) + Ft*cos(phi_tth_dsk) + Fn*sin(phi_tth_dsk);
                    Fy(i_phi) = Fy(i_phi) + Ft*sin(phi_tth_dsk) - Fn*cos(phi_tth_dsk);
                end
                end
            end
            for im = 1:nmd
                z_eul = app.euler_milling2([ux(i_phi,:,im)', uy(i_phi,:,im)'],[-Fx(i_phi), -Fy(i_phi)],[m_x(im), m_y(im)],[zta_x(im), zta_y(im)],[omg_x(im), omg_y(im)], dt);
                ux(i_phi+1,:,im) = z_eul(:,1)';
                uy(i_phi+1,:,im) = z_eul(:,2)';
                x(i_phi+1,:) = x(i_phi+1,:) + ux(i_phi+1,:,im);
                y(i_phi+1,:) = y(i_phi+1,:) + uy(i_phi+1,:,im);
            end
            end
            F = [Fx;Fy]';
            z = [x(:,1), y(:,1)];
            end

            % Crear el gráfico en UIAxes3
            plot(app.UIAxes3, Time, Fx, 'Color', [0 0.4470 0.7410]);
            hold(app.UIAxes3, 'on');
            plot(app.UIAxes3, Time, Fy, 'Color', [0.8500 0.3250 0.0980]);
            hold(app.UIAxes3, 'off');    
            xlabel(app.UIAxes3, 'Time [s]'); ylabel(app.UIAxes3, 'Forces [N]'); legend(app.UIAxes3, 'Fx', 'Fy');
            grid(app.UIAxes3, 'on'); 
            app.UIAxes3.XGrid = 'on'; 
            app.UIAxes3.YGrid = 'on'; 
            app.UIAxes3.XMinorGrid = 'on'; 
            app.UIAxes3.YMinorGrid = 'on'; 
            axis(app.UIAxes3, 'tight');   
        end

        % Value changed function: PeriodstosimulateT10EditField
        function PeriodstosimulateT10EditFieldValueChanged(app, event)
            value = app.PeriodstosimulateT10EditField.Value;
            
        end

        % Button pushed function: SavesimulationdataButton
        function SavesimulationdataButtonPushed(app, event)
                simulationData = app.SimulationResults;
                [file, path] = uiputfile('*.mat', 'Guardar datos de simulación');
                if ischar(file)
                    fullFileName = fullfile(path, file);
                    save(fullFileName, 'simulationData');
                    uialert(app.UIFigure, 'Datos guardados correctamente', 'Éxito');
                else
                    uialert(app.UIFigure, 'La operación fue cancelada', 'Cancelado');
                end
        end

        % Button down function: UIAxes4
        function UIAxes4ButtonDown(app, event)
            
        end

        % Callback function
        function UIAxes5ButtonDown(app, event)

        end

        % Button pushed function: RaButton
        function RaButtonPushed(app, event)
            cla(app.UIAxes1);
            Smin = app.S_minEditField.Value;
            Smax = app.S_maxEditField.Value;
            apmin = 1e-3*app.ap_minEditField.Value;
            apmax = 1e-3*app.ap_maxEditField.Value;
            nx = app.nxEditField.Value;
            ny = app.nyEditField.Value;
            iso = app.IsolinesEditField.Value;
            n_s = linspace(Smin,Smax,nx); 
            ap_s = linspace(apmin,apmax,ny);
 
            contourf(app.UIAxes1, n_s, 1000*ap_s, app.Ra', iso,'LineColor', 'none');
            colormap(app.UIAxes1, 'jet');
            colorbar(app.UIAxes1);
            xlabel(app.UIAxes1, 'Spindle speed [rpm]'); ylabel(app.UIAxes1, 'Axial depth [mm]');
            title(app.UIAxes1, 'Surface roughness Ra [\mum]'); 
        
            hold(app.UIAxes1, 'on');
            [C1, h1] = contour(app.UIAxes1, n_s, 1000*ap_s, app.Ra', 10, 'LineColor', 'w');
            clabel(C1, h1, 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
        
        end

        % Button pushed function: ChatterButton
        function ChatterButtonPushed(app, event)
            cla(app.UIAxes3_2);
            cla(app.UIAxes3_3);
            L = length(app.Fx); 
            Fs = 1 / (app.Time(2) - app.Time(1));

            %% FFT-X
            Y_Fx = fft(app.Fx);
            P2_Fx = abs(Y_Fx / L);
            P1_Fx = P2_Fx(1:L/2+1);
            P1_Fx(2:end-1) = 2*P1_Fx(2:end-1);
            f = Fs * (0:(L/2)) / L;

            %% FFT-Y
            Y_Fy = fft(app.Fy);
            P2_Fy = abs(Y_Fy / L);
            P1_Fy = P2_Fy(1:L/2+1);
            P1_Fy(2:end-1) = 2*P1_Fy(2:end-1);

            plot(app.UIAxes3_2, f, P1_Fx, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
            hold(app.UIAxes3_2, 'on');
            xlabel(app.UIAxes3_2, 'Frequency [Hz]');
            ylabel(app.UIAxes3_2, 'FFT - Fx');
            grid(app.UIAxes3_2, 'on');
            plot(app.UIAxes3_3, f, P1_Fy, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
            hold(app.UIAxes3_3, 'on');
            xlabel(app.UIAxes3_3, 'Frequency [Hz]');
            ylabel(app.UIAxes3_3, 'FFT - Fy');
            grid(app.UIAxes3_3, 'on');

            spindle_speed_rpm = app.SpindlespeednrpmEditField_2.Value;
            num_teeth = app.NofflutesZ10EditField.Value;
            spindle_frequency_Hz = spindle_speed_rpm / (60 * num_teeth);

            num_multiples = floor(max(f) / spindle_frequency_Hz);
            multiples_of_spindle_frequency = spindle_frequency_Hz * (1:num_multiples);

            y_limits_Fx = ylim(app.UIAxes3_2);
            y_limits_Fy = ylim(app.UIAxes3_3);

            for i = 1:num_multiples
                plot(app.UIAxes3_2, multiples_of_spindle_frequency(i), y_limits_Fx(1), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
            end

            for i = 1:num_multiples
                plot(app.UIAxes3_3, multiples_of_spindle_frequency(i), y_limits_Fy(1), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
            end
            
            xlim(app.UIAxes3_2, [0 4000]);
            xlim(app.UIAxes3_3, [0 4000]);            
            ylim(app.UIAxes3_2, 'auto');
            ylim(app.UIAxes3_3, 'auto');
            
            [sorted_P1_Fx, sorted_indices_Fx] = sort(P1_Fx, 'descend');
            peak_frequencies_Fx = f(sorted_indices_Fx(1:3));
            [sorted_P1_Fy, sorted_indices_Fy] = sort(P1_Fy, 'descend');
            peak_frequencies_Fy = f(sorted_indices_Fy(1:3)); 
        
            peak_frequencies_Fx = peak_frequencies_Fx(peak_frequencies_Fx > 10);
            peak_frequencies_Fy = peak_frequencies_Fy(peak_frequencies_Fy > 10);
        
            if length(peak_frequencies_Fx) >= 2
                peak_frequencies_Fx = peak_frequencies_Fx(1:2);
            end
            if length(peak_frequencies_Fy) >= 2
                peak_frequencies_Fy = peak_frequencies_Fy(1:2);
            end
        
            tolerance = 2;
            function is_chatter = detect_chatter(peak_frequencies, multiples_of_spindle_frequency, tolerance)
                is_chatter = false;  % Asumimos que no hay chatter inicialmente
                for i = 1:length(peak_frequencies)
                    % Para cada pico, encontrar el armónico más cercano
                    [min_diff, ~] = min(abs(peak_frequencies(i) - multiples_of_spindle_frequency));
                    
                    % Comprobar si el pico está dentro de la tolerancia del armónico más cercano
                    if min_diff > tolerance
                        is_chatter = true;  % Si está fuera de la tolerancia, hay chatter
                        break;  % No es necesario seguir si encontramos chatter
                    end
                end
            end
        
            % Detectar chatter en Fx y Fy
            chatter_Fx = detect_chatter(peak_frequencies_Fx, multiples_of_spindle_frequency, tolerance);
            chatter_Fy = detect_chatter(peak_frequencies_Fy, multiples_of_spindle_frequency, tolerance);
        
            % Decisión sobre chatter
            if chatter_Fx || chatter_Fy
                uialert(app.UIFigure, 'Chatter detected!', 'Chatter Alert', 'Icon', 'warning');
            else
                uialert(app.UIFigure, 'No chatter detected (forced vibration).', 'Stable Cut', 'Icon', 'success');
            end
        
            hold(app.UIAxes3_2, 'off');
            hold(app.UIAxes3_3, 'off');
        end

        % Button pushed function: ClearandreleaseButton
        function ClearandreleaseButtonPushed(app, event)

                % Limpiar completamente los ejes, incluyendo las gráficas, leyendas, colorbars, etc.
                cla(app.UIAxes3, 'reset');   % Limpiar y restablecer el eje UIAxes3
                cla(app.UIAxes3_2, 'reset'); % Limpiar y restablecer el eje UIAxes3_2
                cla(app.UIAxes3_3, 'reset'); % Limpiar y restablecer el eje UIAxes3_3
            
                % Opcional: Si tienes barras de color o leyendas que quieras eliminar manualmente
                colorbar(app.UIAxes3, 'off');   % Quitar la barra de colores de UIAxes3
                colorbar(app.UIAxes3_2, 'off'); % Quitar la barra de colores de UIAxes3_2
                colorbar(app.UIAxes3_3, 'off'); % Quitar la barra de colores de UIAxes3_3
            
                % Asegurar que no quede ninguna leyenda visible
                legend(app.UIAxes3, 'off');   % Quitar leyenda de UIAxes3
                legend(app.UIAxes3_2, 'off'); % Quitar leyenda de UIAxes3_2
                legend(app.UIAxes3_3, 'off'); % Quitar leyenda de UIAxes3_3

                % Desbloquear los campos que se habían bloqueado al pulsar Run
                app.FeedpertoothfzmmZrevEditField_2.Enable = 'on';      
                app.RadialdepthofcutaemmEditField_2.Enable = 'on';
                app.AxialdepthofcutapmmEditField_2.Enable = 'on';
                app.SpindlespeednrpmEditField_2.Enable = 'on';
                app.PeriodstosimulateT10EditField.Enable = 'on';

        end

        % Callback function
        function RUNButton_3Pushed(app, event)
             
        end

        % Callback function
        function SpindlespeednrpmEditField_3ValueChanged(app, event)
            SS = app.SpindlespeednrpmEditField_3.Value;
        end

        % Selection changed function: TypeofmillingButtonGroup
        function TypeofmillingButtonGroupSelectionChanged(app, event)
            selectedButton = app.TypeofmillingButtonGroup.SelectedObject;
        end

        % Callback function
        function ClearandreleaseButton_2Pushed(app, event)
                app.FeedpertoothfzmmZrevEditField_3.Enable = 'on';
                app.AxialdepthofcutapmmEditField_3.Enable = 'on';
                app.RadialdepthofcutaemmEditField_3.Enable = 'on';
                app.SpindlespeednrpmEditField_3.Enable = 'on';
                app.PeriodstosimulateT10EditField_2.Enable = 'on';
                app.Axialzdiscretizationsteps20EditField.Enable = 'on';
                app.Angulardiscretizationsteps8000EditField.Enable = 'on';
            
                cla(app.UIAxes3_4);
                cla(app.UIAxes3_5);
                cla(app.UIAxes3_6);
                colorbar(app.UIAxes3_4, 'off');
                colorbar(app.UIAxes3_5, 'off');
                colorbar(app.UIAxes3_6, 'off');
        end

        % Callback function
        function SelectintervalButtonPushed(app, event)

        end

        % Button pushed function: ClearandreleaseButton_3
        function ClearandreleaseButton_3Pushed(app, event)
            app.FeedpertoothfzmmZrevEditField.Enable = 'on';       
            app.RadialdepthofcutaemmEditField.Enable = 'on';       
            app.S_minEditField.Enable = 'on';                      
            app.S_maxEditField.Enable = 'on';            
            app.ap_minEditField.Enable = 'on';          
            app.ap_maxEditField.Enable = 'on';       
            app.nxEditField.Enable = 'on';
            app.nyEditField.Enable = 'on';
            app.IsolinesEditField.Enable = 'on';  
            app.Periodstosimulate10EditField.Enable = 'on';  
            app.PeriodstosimulateT10EditField.Enable = 'on';       
            app.RadialdepthofcutaemmEditField_2.Enable = 'on';     
            app.FeedpertoothfzmmZrevEditField_2.Enable = 'on';    

            cla(app.UIAxes1);
            title(app.UIAxes1, ''); 
            legend(app.UIAxes1, 'off');
            if isfield(app, 'ColormapTab')
                delete(app.ColormapTab);
            end
            cla(app.UIAxes1, 'reset');
            title(app.UIAxes2, '');
            legend(app.UIAxes2, 'off');
            cla(app.UIAxes2, 'reset');

            cla(app.UIAxes7);
            app.Panel_5.Visible = 'off';
            app.Panel_6.Visible = 'off';
            app.Panel_7.Visible = 'off';
        end

        % Value changed function: StartingpointDropDown
        function StartingpointDropDownValueChanged(app, event)
            selectedCriterion = app.StartingpointDropDown.Value;    
            Smin = app.S_minEditField.Value;
            Smax = app.S_maxEditField.Value;
            apmin = app.ap_minEditField.Value; % mm
            apmax = app.ap_maxEditField.Value; % mm
            nx = app.nxEditField.Value;
            ny = app.nyEditField.Value;
            n_s = linspace(Smin, Smax, nx); 
            ap_s = linspace(apmin, apmax, ny);

            app.FxymeanFxySlider.Enable = 'off';
            app.MaximumRamSlider_2.Enable = 'off';

            switch selectedCriterion
                case 'Maximize productivity'
                    app.FxymeanFxySlider.Visible = 'on';             
                    app.FxymeanFxySlider.Enable = 'on';             
                    app.EditField.Visible= 'on';                    
                    app.EditField.Enable = 'off';                    
                    app.MaximumRamSlider_2.Visible = 'off';     
                    app.EditField_2.Visible = 'off';             
                    
                    app.MinimumMRRmmminSlider.Visible = 'off';       
                    app.EditField_3.Visible = 'off';                
                    app.MaximumRamSlider.Visible = 'off';  
                    app.EditField_4.Visible = 'off';          
                    
                    cla(app.UIAxes7);
                    [C1, h1] = contour(app.UIAxes7, n_s, ap_s, sqrt(app.ptpx'.^2 + app.ptpy'.^2), 20, 'LineColor', 'k');
                    clabel(C1, h1, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');
                    xlabel(app.UIAxes7, 'Spindle speed [rpm]');
                    ylabel(app.UIAxes7, 'Axial depth of cut [mm]');

                case 'Minimize surface roughness'
                    app.FxymeanFxySlider.Visible = 'off';           
                    app.EditField.Visible = 'off';                  
                    app.MaximumRamSlider_2.Visible = 'on';   
                    app.MaximumRamSlider_2.Enable = 'on';  
                    app.EditField_2.Visible = 'on';             
                    app.EditField_2.Enable = 'off';     

                    app.MinimumMRRmmminSlider.Visible = 'off'; 
                    app.EditField_3.Visible = 'off';              
                    app.MaximumRamSlider.Visible = 'off';  
                    app.EditField_4.Visible = 'off';  

                    cla(app.UIAxes7);
                    [C4, h4] = contour(app.UIAxes7, n_s, ap_s, app.Ra', 20, 'LineColor', 'k');
                    clabel(C4, h4, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');  % Etiquetas en negro
                    xlabel(app.UIAxes7, 'Spindle speed [rpm]');
                    ylabel(app.UIAxes7, 'Axial depth of cut [mm]');
                
                case 'Productivity and surface accuracy'
                    app.FxymeanFxySlider.Visible = 'off';      
                    app.EditField.Visible = 'off';               
                    app.MaximumRamSlider_2.Visible = 'off';
                    app.EditField_2.Visible = 'off';                 
                  
                    app.MinimumMRRmmminSlider.Visible = 'on';      
                    app.MinimumMRRmmminSlider.Enable = 'on';     
                    app.EditField_3.Visible = 'on';                  
                    app.EditField_3.Enable = 'off';              
                    app.MaximumRamSlider.Visible = 'on';    
                    app.MaximumRamSlider.Enable = 'on';     
                    app.EditField_4.Visible = 'on';              
                    app.EditField_4.Enable = 'off';    

                    cla(app.UIAxes7);
                    [C4, h4] = contour(app.UIAxes7, n_s, ap_s, app.Ra', 20, 'LineColor', 'k');
                    clabel(C4, h4, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');
                    xlabel(app.UIAxes7, 'Spindle speed [rpm]');
                    ylabel(app.UIAxes7, 'Axial depth of cut [mm]');
            end            
        end

        % Size changed function: COptimizationtoolTab
        function COptimizationtoolTabSizeChanged(app, event)
            app.FeedpertoothfzmmZrevEditField_4.Value = app.FeedpertoothfzmmZrevEditField.Value;
            app.RadialdepthofcutaemmEditField_4.Value = app.RadialdepthofcutaemmEditField.Value;

            ap_min = app.ap_minEditField.Value;  
            ap_max = app.ap_maxEditField.Value;  
            n_min = app.S_minEditField.Value;    
            n_max = app.S_maxEditField.Value;

            average_ap = (ap_min + ap_max) / 2;
            average_n = (n_min + n_max) / 2;

            app.AxialdepthofcutapmmEditField_4.Value = average_ap;
            app.SpindlespeednrpmEditField_4.Value = average_n;

            app.AxialdepthofcutapmmEditField_4.Enable = 'off';
            app.SpindlespeednrpmEditField_4.Enable = 'off';
            app.FeedpertoothfzmmZrevEditField_4.Enable = 'off';
            app.RadialdepthofcutaemmEditField_4.Enable = 'off';

            app.ProductivityimprovedbyEditField.Visible = 'off';
            app.ProductivityimprovedbyEditField_2.Visible = 'off';
            app.ProductivityimprovedbyEditField_3.Visible = 'off';

            app.Panel_5.Visible = 'off';
            app.Panel_6.Visible = 'off';
            app.Panel_7.Visible = 'off';
        end

        % Button down function: UIAxes7
        function UIAxes7ButtonDown(app, event)
            
        end

        % Button pushed function: OPTIMUMButton
        function OPTIMUMButtonPushed(app, event)
            criterion = app.StartingpointDropDown.Value;
            zn = app.NofflutesZ10EditField.Value;
            Smin = app.S_minEditField.Value;
            Smax = app.S_maxEditField.Value;
            apmin = app.ap_minEditField.Value; % mm
            apmax = app.ap_maxEditField.Value; % mm
            nx = app.nxEditField.Value;
            ny = app.nyEditField.Value;
            n_s = linspace(Smin, Smax, nx); 
            ap_s = linspace(apmin, apmax, ny);

            if isprop(app, 'ForceRatioMarkers') && ~isempty(app.ForceRatioMarkers)
                delete(app.ForceRatioMarkers);
                app.ForceRatioMarkers = [];
            end
            if isprop(app, 'RaMarkers') && ~isempty(app.RaMarkers)
                delete(app.RaMarkers);
                app.RaMarkers = [];
            end
            if isprop(app, 'SurfaceRoughnessMarkers') && ~isempty(app.RaMarkers)
                delete(app.SurfaceRoughnessMarkers);
                app.SurfaceRoughnessMarkers = [];
            end
            if isprop(app, 'BothConditionMarkers') && ~isempty(app.BothConditionMarkers)
                delete(app.BothConditionMarkers);
                app.BothConditionMarkers = [];
            end

            % Delete starting and optimal points
            if isprop(app, 'StartingPointMarker') && ~isempty(app.StartingPointMarker)
                delete(app.StartingPointMarker);
                app.StartingPointMarker = [];
            end
            if isprop(app, 'StartingPointText') && ~isempty(app.StartingPointText)
                delete(app.StartingPointText);
                app.StartingPointText = [];
            end
            if isprop(app, 'OptimalPointMarker') && ~isempty(app.OptimalPointMarker)
                delete(app.OptimalPointMarker);
                app.OptimalPointMarker = [];
            end
            if isprop(app, 'OptimalPointText') && ~isempty(app.OptimalPointText)
                delete(app.OptimalPointText);
                app.OptimalPointText = [];
            end
        
            % Switch based on criterion
            switch criterion
                case 'Maximize productivity'
                    hold(app.UIAxes7, 'on');
                    slider_value = app.FxymeanFxySlider.Value / 100;
                    app.EditField.Enable = 'off';
                    n_force_points = [];
                    ap_force_points = [];
                    force_markers = [];
                    max_Qc = 0;
                    optimal_in = NaN;
                    optimal_iap = NaN;
        
                    for in = 1:nx
                        for iap = 1:ny
                            force_xy = sqrt(app.ptpx(in, iap)^2 + app.ptpy(in, iap)^2);
                            force_mean = sqrt(app.Fx_mean(in, iap)^2 + app.Fy_mean(in, iap)^2);
                            force_ratio = force_xy / force_mean;
                            Qcc = app.RadialdepthofcutaemmEditField_4.Value * ap_s(iap) * ...
                                  app.FeedpertoothfzmmZrevEditField_4.Value * zn * n_s(in);
        
                            if force_ratio <= slider_value
                                n_force_points = [n_force_points; n_s(in)];
                                ap_force_points = [ap_force_points; ap_s(iap)];
                            end

                            if force_ratio <= slider_value && Qcc > max_Qc
                                max_Qc = Qcc;
                                optimal_in = n_s(in);
                                optimal_iap = ap_s(iap);
                            end
                        end
                    end

                    if ~isempty(n_force_points) && ~isempty(ap_force_points)
                        app.ForceRatioMarkers = plot(app.UIAxes7, n_force_points, ap_force_points, 'bo', 'MarkerSize', 3, 'MarkerFaceColor', 'b');
                    else
                        uialert(app.UIFigure, 'No points found with force ratio below the selected value', 'Optimization Alert');
                    end
        
                    if ~isnan(optimal_in) && ~isnan(optimal_iap)
                        if isprop(app, 'FeedpertoothfzmmZrevEditField_5')
                            app.FeedpertoothfzmmZrevEditField_5.Value = app.FeedpertoothfzmmZrevEditField_4.Value;
                        end
                        if isprop(app, 'RadialdepthofcutaemmEditField_5')
                            app.RadialdepthofcutaemmEditField_5.Value = app.RadialdepthofcutaemmEditField_4.Value;
                        end
                        if isprop(app, 'AxialdepthofcutapmmEditField_7')
                            app.AxialdepthofcutapmmEditField_7.Value = optimal_iap;
                            app.AxialdepthofcutapmmEditField_7.Visible = 'on';
                        end
                        if isprop(app, 'SpindlespeednrpmEditField_7')
                            app.SpindlespeednrpmEditField_7.Value = optimal_in;
                            app.SpindlespeednrpmEditField_7.Visible = 'on';
                        end
                        app.SpindlespeednrpmEditField_7.Visible = true;
                        app.AxialdepthofcutapmmEditField_7.Visible = true;
        
                        Qc_st = app.RadialdepthofcutaemmEditField_4.Value * app.AxialdepthofcutapmmEditField_6.Value * app.FeedpertoothfzmmZrevEditField_4.Value * zn * app.SpindlespeednrpmEditField_6.Value;
                        impro_produc = ((max_Qc - Qc_st) / Qc_st) * 100;
        
                        if isprop(app, 'ProductivityimprovedbyEditField')
                            app.ProductivityimprovedbyEditField.Value = round(impro_produc);
                            app.ProductivityimprovedbyEditField.Visible = 'on';
                        end
                        if isprop(app, 'ProductivityimprovedbyEditField_2')
                            app.ProductivityimprovedbyEditField_2.Value = round(Qc_st)*1e-3;
                            app.ProductivityimprovedbyEditField_2.Visible = 'on';
                        end
                        if isprop(app, 'ProductivityimprovedbyEditField_3')
                            app.ProductivityimprovedbyEditField_3.Value = round(max_Qc)*1e-3;
                            app.ProductivityimprovedbyEditField_3.Visible = 'on';
                        end
                        if ~isnan(optimal_in) && ~isnan(optimal_iap)
                            app.ProductivityimprovedbyEditField.Visible = true;
                            app.ProductivityimprovedbyEditField_2.Visible = true;
                            app.ProductivityimprovedbyEditField_3.Visible = true;
                        end

                        app.StartingPointMarker = plot(app.UIAxes7, app.SpindlespeednrpmEditField_6.Value, app.AxialdepthofcutapmmEditField_6.Value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                        app.StartingPointText = text(app.UIAxes7, app.SpindlespeednrpmEditField_6.Value, app.AxialdepthofcutapmmEditField_6.Value, 'Starting Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');
                        app.OptimalPointMarker = plot(app.UIAxes7, optimal_in, optimal_iap, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
                        app.OptimalPointText = text(app.UIAxes7, optimal_in, optimal_iap, 'Optimal Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');

                        app.Panel_5.Visible = 'on';
                        app.Panel_6.Visible = 'off';
                        app.Panel_7.Visible = 'off';
                    else
                        uialert(app.UIFigure, 'No optimum satisfying vibration constraint. Revise the axial depth of cut or productivity slider', 'Optimization Alert');
                    end

                    case 'Minimize surface roughness'

                        hold(app.UIAxes7, 'on');
                        slider_value_ra = app.MaximumRamSlider_2.Value;
                        app.EditField_2.Value = slider_value_ra;  
                        app.EditField_2.Enable = 'off';
                    
                        n_Ra_points = [];
                        ap_Ra_points = [];
                        Ra_values = [];
                        
                        for in = 1:nx
                            for iap = 1:ny
                                Ra_value = app.Ra(in, iap);                         
                                if Ra_value <= slider_value_ra
                                    n_Ra_points = [n_Ra_points; n_s(in)];
                                    ap_Ra_points = [ap_Ra_points; ap_s(iap)];
                                    Ra_values = [Ra_values; Ra_value];
                                end
                            end
                        end
                    
                        if ~isempty(n_Ra_points) && ~isempty(ap_Ra_points)
                            app.RaMarkers = plot(app.UIAxes7, n_Ra_points, ap_Ra_points, 'mo', 'MarkerSize', 3, 'MarkerFaceColor', 'm');
                            
                            tolerance_ra = 0.5;
                            valid_indices = find(abs(Ra_values - slider_value_ra) <= tolerance_ra);
                    
                            if ~isempty(valid_indices)
                                valid_n_Ra_points = n_Ra_points(valid_indices);
                                valid_ap_Ra_points = ap_Ra_points(valid_indices);
                                valid_Ra_values = Ra_values(valid_indices);
                    
                                [~, max_idx] = max(valid_n_Ra_points .* valid_ap_Ra_points);
                                optimal_in = valid_n_Ra_points(max_idx);
                                optimal_iap = valid_ap_Ra_points(max_idx);
                                min_Ra = valid_Ra_values(max_idx);
                    
                                if isprop(app, 'FeedpertoothfzmmZrevEditField_5')
                                    app.FeedpertoothfzmmZrevEditField_5.Value = app.FeedpertoothfzmmZrevEditField_4.Value;
                                end
                                if isprop(app, 'RadialdepthofcutaemmEditField_5')
                                    app.RadialdepthofcutaemmEditField_5.Value = app.RadialdepthofcutaemmEditField_4.Value;
                                end
                                if isprop(app, 'AxialdepthofcutapmmEditField_7')
                                    app.AxialdepthofcutapmmEditField_7.Value = optimal_iap;
                                    app.AxialdepthofcutapmmEditField_7.Visible = 'on';
                                end
                                if isprop(app, 'SpindlespeednrpmEditField_7')
                                    app.SpindlespeednrpmEditField_7.Value = optimal_in;
                                    app.SpindlespeednrpmEditField_7.Visible = 'on';
                                end

                                s_start = (Smin + Smax) / 2;
                                ap_start = (apmin + apmax) / 2;                        
                                [~, idx_n] = min(abs(n_s - s_start));
                                [~, idx_ap] = min(abs(ap_s - ap_start)); 
                                Ra_start = app.Ra(idx_n, idx_ap);
    
                                if isprop(app, 'SurfroughnessimprovedbyEditField')
                                    app.SurfroughnessimprovedbyEditField.Value = ((Ra_start - min_Ra) / Ra_start) * 100;
                                    app.SurfroughnessimprovedbyEditField.Visible = 'on';
                                end
                        
                                if isprop(app, 'ProductivityimprovedbyEditField_4')
                                    app.ProductivityimprovedbyEditField_4.Value = Ra_start;
                                    app.ProductivityimprovedbyEditField_4.Visible = 'on';
                                end
                                if isprop(app, 'ProductivityimprovedbyEditField_5')
                                    app.ProductivityimprovedbyEditField_5.Value = min_Ra;
                                    app.ProductivityimprovedbyEditField_5.Visible = 'on';
                                end
                        
                                if ~isnan(optimal_in) && ~isnan(optimal_iap)
                                    app.SurfroughnessimprovedbyEditField.Visible = true;
                                    app.ProductivityimprovedbyEditField_4.Visible = true;
                                    app.ProductivityimprovedbyEditField_5.Visible = true;
                                end
    
                                app.StartingPointMarker = plot(app.UIAxes7, app.SpindlespeednrpmEditField_6.Value, app.AxialdepthofcutapmmEditField_6.Value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                                app.StartingPointText = text(app.UIAxes7, app.SpindlespeednrpmEditField_6.Value, app.AxialdepthofcutapmmEditField_6.Value, 'Starting Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');
                                app.OptimalPointMarker = plot(app.UIAxes7, optimal_in, optimal_iap, 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
                                app.OptimalPointText = text(app.UIAxes7, optimal_in, optimal_iap, 'Optimal Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');
                                app.Panel_5.Visible = 'off';
                                app.Panel_6.Visible = 'on';
                                app.Panel_7.Visible = 'off';

                            else
                                uialert(app.UIFigure, 'No points found with Ra within the tolerance range. Refine the mesh', 'Plot Alert');
                            end
                    
                        else
                            uialert(app.UIFigure, 'No points found with Ra lower or equal to the slider value.', 'Plot Alert');
                        end

                       case 'Productivity and surface accuracy'

                            hold(app.UIAxes7, 'on');

                            [n_s_grid, ap_s_grid] = meshgrid(n_s, ap_s);
                            Qc2 = app.FeedpertoothfzmmZrevEditField_4.Value .* zn .* ap_s_grid .* app.RadialdepthofcutaemmEditField_4.Value .* n_s_grid;
                            MRR_threshold = app.MinimumMRRmmminSlider.Value;
                            idx_MRR = Qc2 < MRR_threshold;
                            if any(idx_MRR, 'all')
                                [row_MRR, col_MRR] = find(idx_MRR);
                                app.ForceRatioMarkers = plot(app.UIAxes7, n_s(col_MRR), ap_s(row_MRR), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
                            else
                                uialert(app.UIFigure, 'No points found with MRR below the selected value', 'Plot Alert');
                            end
                        
                            Ra_threshold = app.MaximumRamSlider.Value;
                            Ra_tr = app.Ra';
                        
                            idx_Ra = Ra_tr <= Ra_threshold;
                            if any(idx_Ra, 'all')
                                [row_Ra, col_Ra] = find(idx_Ra);
                                app.SurfaceRoughnessMarkers = plot(app.UIAxes7, n_s(col_Ra), ap_s(row_Ra), 'mo', 'MarkerSize', 3, 'MarkerFaceColor', 'm');
                            else
                                uialert(app.UIFigure, 'No points found with Ra below the selected value', 'Plot Alert');
                            end
                        
                            common_idx = idx_MRR & idx_Ra;  
                            if any(common_idx, 'all')
                                [row_common, col_common] = find(common_idx);
                                n_ap_product = n_s(col_common) .* ap_s(row_common);
                                [~, max_idx] = max(n_ap_product);
                            
                                optimal_n_s = n_s(col_common(max_idx));
                                optimal_ap_s = ap_s(row_common(max_idx));
                                optimal_MRR = Qc2(row_common(max_idx), col_common(max_idx));
                                optimal_Ra = Ra_tr(row_common(max_idx), col_common(max_idx));
                            
                                app.OptimalPointMarker = plot(app.UIAxes7, optimal_n_s, optimal_ap_s, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
                                app.OptimalPointText = text(app.UIAxes7, optimal_n_s, optimal_ap_s, 'Optimal Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');
                            
                                app.ProductivityimprovedbyEditField_9.Value = 1e-3 * double(optimal_MRR);
                                app.ProductivityimprovedbyEditField_7.Value = double(optimal_Ra);
                                
                                app.Panel_5.Visible = 'off';
                                app.Panel_6.Visible = 'off';
                                app.Panel_7.Visible = 'on';

                                if isprop(app, 'FeedpertoothfzmmZrevEditField_5')
                                    app.FeedpertoothfzmmZrevEditField_5.Value = app.FeedpertoothfzmmZrevEditField_4.Value;
                                end
                                if isprop(app, 'RadialdepthofcutaemmEditField_5')
                                    app.RadialdepthofcutaemmEditField_5.Value = app.RadialdepthofcutaemmEditField_4.Value;
                                end
                                if isprop(app, 'AxialdepthofcutapmmEditField_7')
                                    app.AxialdepthofcutapmmEditField_7.Value = optimal_ap_s;
                                    app.AxialdepthofcutapmmEditField_7.Visible = 'on';
                                end
                                if isprop(app, 'SpindlespeednrpmEditField_7')
                                    app.SpindlespeednrpmEditField_7.Value = optimal_n_s;
                                    app.SpindlespeednrpmEditField_7.Visible = 'on';
                                end
                            else
                                uialert(app.UIFigure, 'No points found that meet both MRR and Ra conditions', 'Optimization Alert');
                            end                        
                            
                            s_start = (Smin + Smax) / 2;
                            ap_start = (apmin + apmax) / 2;
                        
                            [~, idx_s_start] = min(abs(n_s - s_start));
                            [~, idx_ap_start] = min(abs(ap_s - ap_start));
                        
                            Qc_start = Qc2(idx_s_start, idx_ap_start);
                            Ra_start = Ra_tr(idx_s_start, idx_ap_start);
                        
                            if isscalar(Qc_start) && isnumeric(Qc_start)
                                app.ProductivityimprovedbyEditField_8.Value = 1e-3 * double(Qc_start);
                            else
                                uialert(app.UIFigure, 'Starting MRR is not a valid scalar value.', 'Error');
                            end
                            if isscalar(Ra_start) && isnumeric(Ra_start)
                                app.ProductivityimprovedbyEditField_6.Value = double(Ra_start);
                            else
                                uialert(app.UIFigure, 'Starting Ra is not a valid scalar value.', 'Error');
                            end
                                                   
                            app.StartingPointMarker = plot(app.UIAxes7, s_start, ap_start, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                            app.StartingPointText = text(app.UIAxes7, s_start, ap_start, 'Starting Point', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold');                
                        end
        end

        % Value changed function: FxymeanFxySlider
        function FxymeanFxySliderValueChanged(app, event)
            value = app.FxymeanFxySlider.Value;
            app.EditField.Value = app.FxymeanFxySlider.Value; 
        end

        % Value changed function: RadialdepthofcutaemmEditField_4
        function RadialdepthofcutaemmEditField_4ValueChanged(app, event)
            value = app.RadialdepthofcutaemmEditField_4.Value;
        end

        % Value changed function: AxialdepthofcutapmmEditField_6
        function AxialdepthofcutapmmEditField_6ValueChanged(app, event)
            value = app.AxialdepthofcutapmmEditField_6.Value;
        end

        % Value changed function: SpindlespeednrpmEditField_6
        function SpindlespeednrpmEditField_6ValueChanged(app, event)
            value = app.SpindlespeednrpmEditField_6.Value;
        end

        % Value changed function: FeedpertoothfzmmZrevEditField_5
        function FeedpertoothfzmmZrevEditField_5ValueChanged(app, event)
            value = app.FeedpertoothfzmmZrevEditField_5.Value;
        end

        % Value changed function: RadialdepthofcutaemmEditField_5
        function RadialdepthofcutaemmEditField_5ValueChanged(app, event)
            value = app.RadialdepthofcutaemmEditField_5.Value;
        end

        % Value changed function: AxialdepthofcutapmmEditField_7
        function AxialdepthofcutapmmEditField_7ValueChanged(app, event)
            value = app.AxialdepthofcutapmmEditField_7.Value;
        end

        % Value changed function: SpindlespeednrpmEditField_7
        function SpindlespeednrpmEditField_7ValueChanged(app, event)
            value = app.SpindlespeednrpmEditField_7.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField
        function ProductivityimprovedbyEditFieldValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField.Value;
        end

        % Value changed function: FeedpertoothfzmmZrevEditField_4
        function FeedpertoothfzmmZrevEditField_4ValueChanged(app, event)
            value = app.FeedpertoothfzmmZrevEditField_4.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField_2
        function ProductivityimprovedbyEditField_2ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_2.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField_3
        function ProductivityimprovedbyEditField_3ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_3.Value;
        end

        % Size changed function: Panel_5
        function Panel_5SizeChanged(app, event)
            position = app.Panel_5.Position;
        end

        % Value changed function: MaximumRamSlider_2
        function MaximumRamSlider_2ValueChanged(app, event)
            value = app.MaximumRamSlider_2.Value;
            app.EditField_2.Value = app.MaximumRamSlider_2.Value;           
        end

        % Size changed function: Panel_6
        function Panel_6SizeChanged(app, event)
            position = app.Panel_6.Position;         
        end

        % Value changed function: SurfroughnessimprovedbyEditField
        function SurfroughnessimprovedbyEditFieldValueChanged(app, event)
            value = app.SurfroughnessimprovedbyEditField.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField_4
        function ProductivityimprovedbyEditField_4ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_4.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField_5
        function ProductivityimprovedbyEditField_5ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_5.Value;
        end

        % Value changed function: EditField
        function EditFieldValueChanged(app, event)
            value = app.EditField.Value;
        end

        % Value changed function: EditField_2
        function EditField_2ValueChanged(app, event)
            value = app.EditField_2.Value;       
        end

        % Button down function: UIAxes1
        function UIAxes1ButtonDown(app, event)
            
        end

        % Button down function: UIAxes2
        function UIAxes2ButtonDown(app, event)
            
        end

        % Value changed function: EditField_3
        function EditField_3ValueChanged(app, event)
            value = app.EditField_3.Value;
            
        end

        % Value changed function: EditField_4
        function EditField_4ValueChanged(app, event)
            value = app.EditField_4.Value;
            
        end

        % Value changed function: MinimumMRRmmminSlider
        function MinimumMRRmmminSliderValueChanged(app, event)
            value = app.MinimumMRRmmminSlider.Value;
            app.EditField_3.Value = app.MinimumMRRmmminSlider.Value;
            
        end

        % Value changed function: MaximumRamSlider
        function MaximumRamSliderValueChanged(app, event)
            value = app.MaximumRamSlider.Value;
            app.EditField_4.Value = app.MaximumRamSlider.Value;
        end

        % Value changed function: ProductivityimprovedbyEditField_9
        function ProductivityimprovedbyEditField_9ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_9.Value;
            
        end

        % Value changed function: ProductivityimprovedbyEditField_7
        function ProductivityimprovedbyEditField_7ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_7.Value;
            
        end

        % Value changed function: ProductivityimprovedbyEditField_6
        function ProductivityimprovedbyEditField_6ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_6.Value;
            
        end

        % Value changed function: ProductivityimprovedbyEditField_8
        function ProductivityimprovedbyEditField_8ValueChanged(app, event)
            value = app.ProductivityimprovedbyEditField_8.Value;
            
        end

        % Size changed function: Panel_7
        function Panel_7SizeChanged(app, event)
            position = app.Panel_7.Position;
        end

        % Button pushed function: ClearfigureButton
        function ClearfigureButtonPushed(app, event)
            % Limpiar los elementos dibujados en la figura
            cla(app.UIAxes7);
            if isprop(app, 'AxialdepthofcutapmmEditField_7')
                app.AxialdepthofcutapmmEditField_7.Value = 0;
            end
            if isprop(app, 'SpindlespeednrpmEditField_7')
                app.SpindlespeednrpmEditField_7.Value = 0;
            end
        end

        % Callback function
        function DropDownValueChanged(app, event)

        end

        % Size changed function: MILLTab
        function MILLTabSizeChanged(app, event)
            position = app.MILLTab.Position;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [0 40 1400 700];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.Layout.Row = [1 9];
            app.TabGroup.Layout.Column = [1 4];

            % Create MILLTab
            app.MILLTab = uitab(app.TabGroup);
            app.MILLTab.SizeChangedFcn = createCallbackFcn(app, @MILLTabSizeChanged, true);
            app.MILLTab.Title = 'MILL+';
            app.MILLTab.BackgroundColor = [0.1804 0.1647 0.1647];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.MILLTab);
            app.GridLayout2.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x'};
            app.GridLayout2.BackgroundColor = [0.0627 0.0706 0.0706];

            % Create Image5
            app.Image5 = uiimage(app.GridLayout2);
            app.Image5.Layout.Row = [1 3];
            app.Image5.Layout.Column = [1 2];
            app.Image5.ImageSource = fullfile(pathToMLAPP, 'software_logo5.jpg');

            % Create Label
            app.Label = uilabel(app.GridLayout2);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontName = 'Century';
            app.Label.FontSize = 48;
            app.Label.FontWeight = 'bold';
            app.Label.FontColor = [0.0588 1 1];
            app.Label.Layout.Row = 1;
            app.Label.Layout.Column = 3;
            app.Label.Text = 'MILL+ v1.0';

            % Create Label_4
            app.Label_4 = uilabel(app.GridLayout2);
            app.Label_4.HorizontalAlignment = 'center';
            app.Label_4.FontName = 'Century';
            app.Label_4.FontSize = 30;
            app.Label_4.FontWeight = 'bold';
            app.Label_4.FontColor = [0.0588 1 1];
            app.Label_4.Layout.Row = 2;
            app.Label_4.Layout.Column = [3 4];
            app.Label_4.Text = 'CHATTER AND SURFACE ROUGHNESS';

            % Create Label_3
            app.Label_3 = uilabel(app.GridLayout2);
            app.Label_3.HorizontalAlignment = 'center';
            app.Label_3.FontName = 'Century';
            app.Label_3.FontSize = 30;
            app.Label_3.FontWeight = 'bold';
            app.Label_3.FontColor = [0.0588 1 1];
            app.Label_3.Layout.Row = 3;
            app.Label_3.Layout.Column = [3 4];
            app.Label_3.Text = 'FOR  ADVANCED  MILLING                       ';

            % Create Label_2
            app.Label_2 = uilabel(app.GridLayout2);
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.FontName = 'Century';
            app.Label_2.FontSize = 18;
            app.Label_2.FontWeight = 'bold';
            app.Label_2.FontColor = [0.0588 1 1];
            app.Label_2.Layout.Row = 3;
            app.Label_2.Layout.Column = 3;
            app.Label_2.Text = '';

            % Create CraftedbyDptofMechEngineeringUPVEHULabel
            app.CraftedbyDptofMechEngineeringUPVEHULabel = uilabel(app.GridLayout2);
            app.CraftedbyDptofMechEngineeringUPVEHULabel.HorizontalAlignment = 'center';
            app.CraftedbyDptofMechEngineeringUPVEHULabel.FontName = 'Century';
            app.CraftedbyDptofMechEngineeringUPVEHULabel.FontSize = 24;
            app.CraftedbyDptofMechEngineeringUPVEHULabel.FontWeight = 'bold';
            app.CraftedbyDptofMechEngineeringUPVEHULabel.FontColor = [0.0588 1 1];
            app.CraftedbyDptofMechEngineeringUPVEHULabel.Layout.Row = 4;
            app.CraftedbyDptofMechEngineeringUPVEHULabel.Layout.Column = [1 4];
            app.CraftedbyDptofMechEngineeringUPVEHULabel.Text = 'Crafted by Dpt. of Mech. Engineering, UPV/EHU and Mechanical Eng. and Adv. Mat. Dpt. of Tecnológico de Monterrey';

            % Create ToolmaterialandmodalparametersTab
            app.ToolmaterialandmodalparametersTab = uitab(app.TabGroup);
            app.ToolmaterialandmodalparametersTab.Title = 'Tool-material and modal parameters';

            % Create ToolparametersPanel
            app.ToolparametersPanel = uipanel(app.ToolmaterialandmodalparametersTab);
            app.ToolparametersPanel.ForegroundColor = [0.0588 1 1];
            app.ToolparametersPanel.Title = 'Tool parameters';
            app.ToolparametersPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.ToolparametersPanel.FontName = 'Franklin Gothic Medium';
            app.ToolparametersPanel.FontWeight = 'bold';
            app.ToolparametersPanel.Position = [12 11 572 632];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.ToolparametersPanel);
            title(app.UIAxes4, 'Title')
            xlabel(app.UIAxes4, 'X')
            ylabel(app.UIAxes4, 'Y')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.ButtonDownFcn = createCallbackFcn(app, @UIAxes4ButtonDown, true);
            app.UIAxes4.Position = [63 33 452 414];

            % Create TooldiameterDmmEditFieldLabel
            app.TooldiameterDmmEditFieldLabel = uilabel(app.ToolparametersPanel);
            app.TooldiameterDmmEditFieldLabel.HorizontalAlignment = 'right';
            app.TooldiameterDmmEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.TooldiameterDmmEditFieldLabel.FontColor = [0.0588 1 1];
            app.TooldiameterDmmEditFieldLabel.Position = [11 567 120 22];
            app.TooldiameterDmmEditFieldLabel.Text = 'Tool diameter D [mm]';

            % Create TooldiameterDmmEditField
            app.TooldiameterDmmEditField = uieditfield(app.ToolparametersPanel, 'numeric');
            app.TooldiameterDmmEditField.ValueChangedFcn = createCallbackFcn(app, @TooldiameterDmmEditFieldValueChanged, true);
            app.TooldiameterDmmEditField.FontName = 'Franklin Gothic Medium';
            app.TooldiameterDmmEditField.Position = [150 560 45 35];
            app.TooldiameterDmmEditField.Value = 12;

            % Create NofflutesZ10EditFieldLabel
            app.NofflutesZ10EditFieldLabel = uilabel(app.ToolparametersPanel);
            app.NofflutesZ10EditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.NofflutesZ10EditFieldLabel.FontColor = [0.0588 1 1];
            app.NofflutesZ10EditFieldLabel.Position = [13 525 121 22];
            app.NofflutesZ10EditFieldLabel.Text = '  Nº of flutes Z [-] (<10)';

            % Create NofflutesZ10EditField
            app.NofflutesZ10EditField = uieditfield(app.ToolparametersPanel, 'numeric');
            app.NofflutesZ10EditField.ValueChangedFcn = createCallbackFcn(app, @NofflutesZ10EditFieldValueChanged, true);
            app.NofflutesZ10EditField.FontName = 'Franklin Gothic Medium';
            app.NofflutesZ10EditField.Position = [150 518 45 35];
            app.NofflutesZ10EditField.Value = 4;

            % Create PLOTMILLButton
            app.PLOTMILLButton = uibutton(app.ToolparametersPanel, 'push');
            app.PLOTMILLButton.ButtonPushedFcn = createCallbackFcn(app, @PLOTMILLButtonPushed, true);
            app.PLOTMILLButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.PLOTMILLButton.FontName = 'Franklin Gothic Medium';
            app.PLOTMILLButton.FontSize = 18;
            app.PLOTMILLButton.FontColor = [0.0588 1 1];
            app.PLOTMILLButton.Position = [364 509 107 54];
            app.PLOTMILLButton.Text = 'PLOT MILL';

            % Create Helixangle60EditFieldLabel
            app.Helixangle60EditFieldLabel = uilabel(app.ToolparametersPanel);
            app.Helixangle60EditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.Helixangle60EditFieldLabel.FontColor = [0.0588 1 1];
            app.Helixangle60EditFieldLabel.Position = [13 484 127 22];
            app.Helixangle60EditFieldLabel.Text = '  Helix angle β [º] (<60º)';

            % Create Helixangle60EditField
            app.Helixangle60EditField = uieditfield(app.ToolparametersPanel, 'numeric');
            app.Helixangle60EditField.ValueChangedFcn = createCallbackFcn(app, @Helixangle60EditFieldValueChanged, true);
            app.Helixangle60EditField.FontName = 'Franklin Gothic Medium';
            app.Helixangle60EditField.Position = [150 477 45 35];
            app.Helixangle60EditField.Value = 30;

            % Create CuttingcoefficientsPanel
            app.CuttingcoefficientsPanel = uipanel(app.ToolmaterialandmodalparametersTab);
            app.CuttingcoefficientsPanel.ForegroundColor = [0.0588 1 1];
            app.CuttingcoefficientsPanel.Title = 'Cutting coefficients';
            app.CuttingcoefficientsPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.CuttingcoefficientsPanel.FontName = 'Franklin Gothic Medium';
            app.CuttingcoefficientsPanel.FontWeight = 'bold';
            app.CuttingcoefficientsPanel.Position = [600 496 769 147];

            % Create KtMPaEditFieldLabel
            app.KtMPaEditFieldLabel = uilabel(app.CuttingcoefficientsPanel);
            app.KtMPaEditFieldLabel.HorizontalAlignment = 'right';
            app.KtMPaEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.KtMPaEditFieldLabel.FontColor = [0.0588 1 1];
            app.KtMPaEditFieldLabel.Position = [18 79 52 22];
            app.KtMPaEditFieldLabel.Text = 'Kt  [MPa]';

            % Create KtMPaEditField
            app.KtMPaEditField = uieditfield(app.CuttingcoefficientsPanel, 'numeric');
            app.KtMPaEditField.ValueChangedFcn = createCallbackFcn(app, @KtMPaEditFieldValueChanged, true);
            app.KtMPaEditField.FontName = 'Franklin Gothic Medium';
            app.KtMPaEditField.Position = [89 73 58 33];
            app.KtMPaEditField.Value = 700;

            % Create MaterialslibraryLabel
            app.MaterialslibraryLabel = uilabel(app.CuttingcoefficientsPanel);
            app.MaterialslibraryLabel.HorizontalAlignment = 'center';
            app.MaterialslibraryLabel.FontName = 'Franklin Gothic Medium';
            app.MaterialslibraryLabel.FontSize = 24;
            app.MaterialslibraryLabel.Position = [186 30 413 80];
            app.MaterialslibraryLabel.Text = 'Materials library';

            % Create KnMPaEditFieldLabel
            app.KnMPaEditFieldLabel = uilabel(app.CuttingcoefficientsPanel);
            app.KnMPaEditFieldLabel.HorizontalAlignment = 'right';
            app.KnMPaEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.KnMPaEditFieldLabel.FontColor = [0.0588 1 1];
            app.KnMPaEditFieldLabel.Position = [18 32 55 22];
            app.KnMPaEditFieldLabel.Text = 'Kn  [MPa]';

            % Create KnMPaEditField
            app.KnMPaEditField = uieditfield(app.CuttingcoefficientsPanel, 'numeric');
            app.KnMPaEditField.ValueChangedFcn = createCallbackFcn(app, @KnMPaEditFieldValueChanged, true);
            app.KnMPaEditField.FontName = 'Franklin Gothic Medium';
            app.KnMPaEditField.Position = [89 26 58 33];
            app.KnMPaEditField.Value = 200;

            % Create ModalparameteresPanel
            app.ModalparameteresPanel = uipanel(app.ToolmaterialandmodalparametersTab);
            app.ModalparameteresPanel.ForegroundColor = [0.0588 1 1];
            app.ModalparameteresPanel.Title = 'Modal parameteres';
            app.ModalparameteresPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.ModalparameteresPanel.FontName = 'Franklin Gothic Medium';
            app.ModalparameteresPanel.FontWeight = 'bold';
            app.ModalparameteresPanel.Position = [600 11 769 470];

            % Create ModesinXPanel
            app.ModesinXPanel = uipanel(app.ModalparameteresPanel);
            app.ModesinXPanel.ForegroundColor = [0.0588 1 1];
            app.ModesinXPanel.Title = 'Modes in X';
            app.ModesinXPanel.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ModesinXPanel.FontName = 'Franklin Gothic Medium';
            app.ModesinXPanel.FontWeight = 'bold';
            app.ModesinXPanel.Position = [11 14 363 427];

            % Create FrequencyfHzEditField_3Label
            app.FrequencyfHzEditField_3Label = uilabel(app.ModesinXPanel);
            app.FrequencyfHzEditField_3Label.FontName = 'Franklin Gothic Medium';
            app.FrequencyfHzEditField_3Label.FontColor = [0.0588 1 1];
            app.FrequencyfHzEditField_3Label.Position = [35 367 87 22];
            app.FrequencyfHzEditField_3Label.Text = 'Frequency f [Hz]';

            % Create FrequencyfHzEditField_3
            app.FrequencyfHzEditField_3 = uieditfield(app.ModesinXPanel, 'numeric');
            app.FrequencyfHzEditField_3.ValueChangedFcn = createCallbackFcn(app, @FrequencyfHzEditFieldValueChanged, true);
            app.FrequencyfHzEditField_3.FontName = 'Franklin Gothic Medium';
            app.FrequencyfHzEditField_3.Position = [139 367 54 22];
            app.FrequencyfHzEditField_3.Value = 500;

            % Create StiffnesskNmEditField_3Label
            app.StiffnesskNmEditField_3Label = uilabel(app.ModesinXPanel);
            app.StiffnesskNmEditField_3Label.HorizontalAlignment = 'right';
            app.StiffnesskNmEditField_3Label.FontName = 'Franklin Gothic Medium';
            app.StiffnesskNmEditField_3Label.FontColor = [0.0588 1 1];
            app.StiffnesskNmEditField_3Label.Position = [29 336 94 22];
            app.StiffnesskNmEditField_3Label.Text = 'Stiffness k [N/m]';

            % Create StiffnesskNmEditField_3
            app.StiffnesskNmEditField_3 = uieditfield(app.ModesinXPanel, 'numeric');
            app.StiffnesskNmEditField_3.ValueChangedFcn = createCallbackFcn(app, @StiffnesskNmEditFieldValueChanged, true);
            app.StiffnesskNmEditField_3.FontName = 'Franklin Gothic Medium';
            app.StiffnesskNmEditField_3.Position = [139 336 54 22];
            app.StiffnesskNmEditField_3.Value = 10000000;

            % Create DampingEditField_3Label
            app.DampingEditField_3Label = uilabel(app.ModesinXPanel);
            app.DampingEditField_3Label.FontName = 'Franklin Gothic Medium';
            app.DampingEditField_3Label.FontColor = [0.0588 1 1];
            app.DampingEditField_3Label.Position = [32 305 92 22];
            app.DampingEditField_3Label.Text = ' Damping ζ [-]';

            % Create DampingEditField_3
            app.DampingEditField_3 = uieditfield(app.ModesinXPanel, 'numeric');
            app.DampingEditField_3.ValueChangedFcn = createCallbackFcn(app, @DampingEditFieldValueChanged, true);
            app.DampingEditField_3.FontName = 'Franklin Gothic Medium';
            app.DampingEditField_3.Position = [139 305 54 22];
            app.DampingEditField_3.Value = 0.01;

            % Create Image5_3
            app.Image5_3 = uiimage(app.ModesinXPanel);
            app.Image5_3.Position = [-2 69 373 202];
            app.Image5_3.ImageSource = fullfile(pathToMLAPP, 'ModoX.png');

            % Create ModesinYPanel
            app.ModesinYPanel = uipanel(app.ToolmaterialandmodalparametersTab);
            app.ModesinYPanel.ForegroundColor = [0.0588 1 1];
            app.ModesinYPanel.Title = 'Modes in Y';
            app.ModesinYPanel.BackgroundColor = [0.1255 0.3098 0.3098];
            app.ModesinYPanel.FontName = 'Franklin Gothic Medium';
            app.ModesinYPanel.FontWeight = 'bold';
            app.ModesinYPanel.Position = [993 26 363 426];

            % Create FrequencyfHzEditField_4Label
            app.FrequencyfHzEditField_4Label = uilabel(app.ModesinYPanel);
            app.FrequencyfHzEditField_4Label.FontName = 'Franklin Gothic Medium';
            app.FrequencyfHzEditField_4Label.FontColor = [0.0588 1 1];
            app.FrequencyfHzEditField_4Label.Position = [29 366 93 22];
            app.FrequencyfHzEditField_4Label.Text = '  Frequency f [Hz]';

            % Create FrequencyfHzEditField_4
            app.FrequencyfHzEditField_4 = uieditfield(app.ModesinYPanel, 'numeric');
            app.FrequencyfHzEditField_4.ValueChangedFcn = createCallbackFcn(app, @FrequencyfHzEditField_2ValueChanged, true);
            app.FrequencyfHzEditField_4.FontName = 'Franklin Gothic Medium';
            app.FrequencyfHzEditField_4.Position = [139 366 54 22];
            app.FrequencyfHzEditField_4.Value = 500;

            % Create StiffnesskNmEditField_4Label
            app.StiffnesskNmEditField_4Label = uilabel(app.ModesinYPanel);
            app.StiffnesskNmEditField_4Label.HorizontalAlignment = 'right';
            app.StiffnesskNmEditField_4Label.FontName = 'Franklin Gothic Medium';
            app.StiffnesskNmEditField_4Label.FontColor = [0.0588 1 1];
            app.StiffnesskNmEditField_4Label.Position = [29 335 94 22];
            app.StiffnesskNmEditField_4Label.Text = 'Stiffness k [N/m]';

            % Create StiffnesskNmEditField_4
            app.StiffnesskNmEditField_4 = uieditfield(app.ModesinYPanel, 'numeric');
            app.StiffnesskNmEditField_4.ValueChangedFcn = createCallbackFcn(app, @StiffnesskNmEditField_2ValueChanged, true);
            app.StiffnesskNmEditField_4.FontName = 'Franklin Gothic Medium';
            app.StiffnesskNmEditField_4.Position = [139 335 54 22];
            app.StiffnesskNmEditField_4.Value = 10000000;

            % Create DampingEditField_4Label
            app.DampingEditField_4Label = uilabel(app.ModesinYPanel);
            app.DampingEditField_4Label.FontName = 'Franklin Gothic Medium';
            app.DampingEditField_4Label.FontColor = [0.0588 1 1];
            app.DampingEditField_4Label.Position = [32 304 92 22];
            app.DampingEditField_4Label.Text = ' Damping ζ [-]';

            % Create DampingEditField_4
            app.DampingEditField_4 = uieditfield(app.ModesinYPanel, 'numeric');
            app.DampingEditField_4.ValueChangedFcn = createCallbackFcn(app, @DampingEditField_2ValueChanged, true);
            app.DampingEditField_4.FontName = 'Franklin Gothic Medium';
            app.DampingEditField_4.Position = [139 304 54 22];
            app.DampingEditField_4.Value = 0.01;

            % Create Image5_4
            app.Image5_4 = uiimage(app.ModesinYPanel);
            app.Image5_4.Position = [-10 68 376 204];
            app.Image5_4.ImageSource = fullfile(pathToMLAPP, 'ModoY.png');

            % Create MillingoperationTab
            app.MillingoperationTab = uitab(app.TabGroup);
            app.MillingoperationTab.Title = 'Milling operation';
            app.MillingoperationTab.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MillingoperationTab.ForegroundColor = [0.2314 0.2157 0.2157];

            % Create Panel
            app.Panel = uipanel(app.MillingoperationTab);
            app.Panel.BackgroundColor = [0.2314 0.502 0.502];
            app.Panel.Position = [1100 17 260 618];

            % Create Lamp
            app.Lamp = uilamp(app.Panel);
            app.Lamp.Position = [74 422 112 112];

            % Create RUNButton
            app.RUNButton = uibutton(app.Panel, 'push');
            app.RUNButton.ButtonPushedFcn = createCallbackFcn(app, @RUNButtonPushed, true);
            app.RUNButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.RUNButton.FontName = 'Franklin Gothic Medium';
            app.RUNButton.FontSize = 60;
            app.RUNButton.FontWeight = 'bold';
            app.RUNButton.FontColor = [0.0588 1 1];
            app.RUNButton.Position = [49 262 163 86];
            app.RUNButton.Text = 'RUN';

            % Create TypeofmillingButtonGroup
            app.TypeofmillingButtonGroup = uibuttongroup(app.MillingoperationTab);
            app.TypeofmillingButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @TypeofmillingButtonGroupSelectionChanged, true);
            app.TypeofmillingButtonGroup.ForegroundColor = [0.0588 1 1];
            app.TypeofmillingButtonGroup.Title = 'Type of milling';
            app.TypeofmillingButtonGroup.BackgroundColor = [0.2314 0.502 0.502];
            app.TypeofmillingButtonGroup.FontName = 'Franklin Gothic Medium';
            app.TypeofmillingButtonGroup.FontWeight = 'bold';
            app.TypeofmillingButtonGroup.Position = [14 323 567 311];

            % Create UpmillingButton_2
            app.UpmillingButton_2 = uiradiobutton(app.TypeofmillingButtonGroup);
            app.UpmillingButton_2.Text = 'Up-milling';
            app.UpmillingButton_2.FontName = 'Franklin Gothic Medium';
            app.UpmillingButton_2.FontColor = [0.0588 1 1];
            app.UpmillingButton_2.Position = [291 239 73 22];
            app.UpmillingButton_2.Value = true;

            % Create DownmillingButton
            app.DownmillingButton = uiradiobutton(app.TypeofmillingButtonGroup);
            app.DownmillingButton.Text = 'Down-milling';
            app.DownmillingButton.FontName = 'Franklin Gothic Medium';
            app.DownmillingButton.FontColor = [0.0588 1 1];
            app.DownmillingButton.Position = [23 239 89 22];

            % Create Image
            app.Image = uiimage(app.TypeofmillingButtonGroup);
            app.Image.Position = [19 15 255 219];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'Downmilling3.jpg');

            % Create Image_2
            app.Image_2 = uiimage(app.TypeofmillingButtonGroup);
            app.Image_2.Position = [291 15 255 219];
            app.Image_2.ImageSource = fullfile(pathToMLAPP, 'Upmilling3.jpg');

            % Create CuttingparametersButtonGroup
            app.CuttingparametersButtonGroup = uibuttongroup(app.MillingoperationTab);
            app.CuttingparametersButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @updateInmersionLevel, true);
            app.CuttingparametersButtonGroup.ForegroundColor = [0.0588 1 1];
            app.CuttingparametersButtonGroup.Title = 'Cutting parameters';
            app.CuttingparametersButtonGroup.BackgroundColor = [0.2314 0.502 0.502];
            app.CuttingparametersButtonGroup.FontName = 'Franklin Gothic Medium';
            app.CuttingparametersButtonGroup.FontWeight = 'bold';
            app.CuttingparametersButtonGroup.Position = [597 323 486 312];

            % Create FeedpertoothfzmmZrevEditFieldLabel
            app.FeedpertoothfzmmZrevEditFieldLabel = uilabel(app.CuttingparametersButtonGroup);
            app.FeedpertoothfzmmZrevEditFieldLabel.HorizontalAlignment = 'right';
            app.FeedpertoothfzmmZrevEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditFieldLabel.FontColor = [0.0588 1 1];
            app.FeedpertoothfzmmZrevEditFieldLabel.Position = [9 251 155 22];
            app.FeedpertoothfzmmZrevEditFieldLabel.Text = 'Feed per tooth fz [mm/Z/rev]';

            % Create FeedpertoothfzmmZrevEditField
            app.FeedpertoothfzmmZrevEditField = uieditfield(app.CuttingparametersButtonGroup, 'numeric');
            app.FeedpertoothfzmmZrevEditField.ValueChangedFcn = createCallbackFcn(app, @FeedpertoothfzmmZrevEditFieldValueChanged, true);
            app.FeedpertoothfzmmZrevEditField.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField.Position = [180 245 35 35];
            app.FeedpertoothfzmmZrevEditField.Value = 0.1;

            % Create RadialdepthofcutaemmEditFieldLabel
            app.RadialdepthofcutaemmEditFieldLabel = uilabel(app.CuttingparametersButtonGroup);
            app.RadialdepthofcutaemmEditFieldLabel.HorizontalAlignment = 'right';
            app.RadialdepthofcutaemmEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditFieldLabel.FontColor = [0.0588 1 1];
            app.RadialdepthofcutaemmEditFieldLabel.Position = [10 206 149 22];
            app.RadialdepthofcutaemmEditFieldLabel.Text = 'Radial depth of cut ae [mm]';

            % Create RadialdepthofcutaemmEditField
            app.RadialdepthofcutaemmEditField = uieditfield(app.CuttingparametersButtonGroup, 'numeric');
            app.RadialdepthofcutaemmEditField.ValueChangedFcn = createCallbackFcn(app, @RadialdepthofcutaemmEditFieldValueChanged, true);
            app.RadialdepthofcutaemmEditField.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField.Position = [180 200 35 35];
            app.RadialdepthofcutaemmEditField.Value = 6;

            % Create AxialdepthofcutapmmEditFieldLabel
            app.AxialdepthofcutapmmEditFieldLabel = uilabel(app.CuttingparametersButtonGroup);
            app.AxialdepthofcutapmmEditFieldLabel.HorizontalAlignment = 'right';
            app.AxialdepthofcutapmmEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditFieldLabel.FontColor = [0.0588 1 1];
            app.AxialdepthofcutapmmEditFieldLabel.Position = [11 159 141 22];
            app.AxialdepthofcutapmmEditFieldLabel.Text = 'Axial depth of cut ap [mm]';

            % Create AxialdepthofcutapmmEditField
            app.AxialdepthofcutapmmEditField = uieditfield(app.CuttingparametersButtonGroup, 'text');
            app.AxialdepthofcutapmmEditField.Editable = 'off';
            app.AxialdepthofcutapmmEditField.HorizontalAlignment = 'right';
            app.AxialdepthofcutapmmEditField.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField.Enable = 'off';
            app.AxialdepthofcutapmmEditField.Position = [180 149 35 41];
            app.AxialdepthofcutapmmEditField.Value = '[-]';

            % Create SpindlespeednrpmEditField
            app.SpindlespeednrpmEditField = uieditfield(app.CuttingparametersButtonGroup, 'text');
            app.SpindlespeednrpmEditField.Editable = 'off';
            app.SpindlespeednrpmEditField.HorizontalAlignment = 'right';
            app.SpindlespeednrpmEditField.Enable = 'off';
            app.SpindlespeednrpmEditField.Position = [180 98 35 41];
            app.SpindlespeednrpmEditField.Value = '[-]';

            % Create vcDn1000Label
            app.vcDn1000Label = uilabel(app.CuttingparametersButtonGroup);
            app.vcDn1000Label.FontColor = [0.0588 1 1];
            app.vcDn1000Label.Position = [234 89 124 58];
            app.vcDn1000Label.Text = '( vc = (π∙D∙n) / 1000 )';

            % Create SpindlespeednrpmEditFieldLabel
            app.SpindlespeednrpmEditFieldLabel = uilabel(app.CuttingparametersButtonGroup);
            app.SpindlespeednrpmEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditFieldLabel.FontColor = [0.0588 1 1];
            app.SpindlespeednrpmEditFieldLabel.Position = [15 108 118 22];
            app.SpindlespeednrpmEditFieldLabel.Text = 'Spindle speed n [rpm]';

            % Create SaveoptionsinAStabilitychartsPanel
            app.SaveoptionsinAStabilitychartsPanel = uipanel(app.MillingoperationTab);
            app.SaveoptionsinAStabilitychartsPanel.ForegroundColor = [0.0588 1 1];
            app.SaveoptionsinAStabilitychartsPanel.Title = 'Save options in A) Stability charts';
            app.SaveoptionsinAStabilitychartsPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.SaveoptionsinAStabilitychartsPanel.FontName = 'Franklin Gothic Medium';
            app.SaveoptionsinAStabilitychartsPanel.FontWeight = 'bold';
            app.SaveoptionsinAStabilitychartsPanel.Position = [16 17 566 287];

            % Create SavesimulationdataButton
            app.SavesimulationdataButton = uibutton(app.SaveoptionsinAStabilitychartsPanel, 'push');
            app.SavesimulationdataButton.ButtonPushedFcn = createCallbackFcn(app, @SavesimulationdataButtonPushed, true);
            app.SavesimulationdataButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.SavesimulationdataButton.FontName = 'Franklin Gothic Medium';
            app.SavesimulationdataButton.FontSize = 14;
            app.SavesimulationdataButton.FontWeight = 'bold';
            app.SavesimulationdataButton.FontColor = [0.0588 1 1];
            app.SavesimulationdataButton.Position = [21 209 180 40];
            app.SavesimulationdataButton.Text = 'Save simulation data';

            % Create ClearandreleaseButton_3
            app.ClearandreleaseButton_3 = uibutton(app.SaveoptionsinAStabilitychartsPanel, 'push');
            app.ClearandreleaseButton_3.ButtonPushedFcn = createCallbackFcn(app, @ClearandreleaseButton_3Pushed, true);
            app.ClearandreleaseButton_3.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ClearandreleaseButton_3.FontName = 'Franklin Gothic Medium';
            app.ClearandreleaseButton_3.FontSize = 14;
            app.ClearandreleaseButton_3.FontWeight = 'bold';
            app.ClearandreleaseButton_3.FontColor = [0.0588 1 1];
            app.ClearandreleaseButton_3.Position = [21 153 180 40];
            app.ClearandreleaseButton_3.Text = 'Clear and release';

            % Create GeneratereportButton
            app.GeneratereportButton = uibutton(app.SaveoptionsinAStabilitychartsPanel, 'push');
            app.GeneratereportButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.GeneratereportButton.FontName = 'Franklin Gothic Medium';
            app.GeneratereportButton.FontSize = 14;
            app.GeneratereportButton.FontWeight = 'bold';
            app.GeneratereportButton.FontColor = [0.0588 1 1];
            app.GeneratereportButton.Position = [21 97 180 40];
            app.GeneratereportButton.Text = 'Generate report';

            % Create SimulationparametersPanel
            app.SimulationparametersPanel = uipanel(app.MillingoperationTab);
            app.SimulationparametersPanel.ForegroundColor = [0.0588 1 1];
            app.SimulationparametersPanel.Title = 'Simulation parameters';
            app.SimulationparametersPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.SimulationparametersPanel.FontName = 'Franklin Gothic Medium';
            app.SimulationparametersPanel.FontWeight = 'bold';
            app.SimulationparametersPanel.Position = [597 17 486 287];

            % Create SpindlespeedrangerpmLabel
            app.SpindlespeedrangerpmLabel = uilabel(app.SimulationparametersPanel);
            app.SpindlespeedrangerpmLabel.FontName = 'Franklin Gothic Medium';
            app.SpindlespeedrangerpmLabel.FontColor = [0.0588 1 1];
            app.SpindlespeedrangerpmLabel.Position = [14 237 141 22];
            app.SpindlespeedrangerpmLabel.Text = 'Spindle speed range [rpm]';

            % Create S_minEditFieldLabel
            app.S_minEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.S_minEditFieldLabel.HorizontalAlignment = 'right';
            app.S_minEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.S_minEditFieldLabel.FontColor = [0.0588 1 1];
            app.S_minEditFieldLabel.Position = [11 206 38 22];
            app.S_minEditFieldLabel.Text = 'S_min';

            % Create S_minEditField
            app.S_minEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.S_minEditField.ValueChangedFcn = createCallbackFcn(app, @S_minEditFieldValueChanged, true);
            app.S_minEditField.FontName = 'Franklin Gothic Medium';
            app.S_minEditField.Position = [65 202 49 29];
            app.S_minEditField.Value = 3000;

            % Create S_maxEditFieldLabel
            app.S_maxEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.S_maxEditFieldLabel.HorizontalAlignment = 'right';
            app.S_maxEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.S_maxEditFieldLabel.FontColor = [0.0588 1 1];
            app.S_maxEditFieldLabel.Position = [131 206 40 22];
            app.S_maxEditFieldLabel.Text = 'S_max';

            % Create S_maxEditField
            app.S_maxEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.S_maxEditField.ValueChangedFcn = createCallbackFcn(app, @S_maxEditFieldValueChanged, true);
            app.S_maxEditField.FontName = 'Franklin Gothic Medium';
            app.S_maxEditField.Position = [187 202 49 29];
            app.S_maxEditField.Value = 10000;

            % Create MESHLabel
            app.MESHLabel = uilabel(app.SimulationparametersPanel);
            app.MESHLabel.HorizontalAlignment = 'center';
            app.MESHLabel.FontName = 'Franklin Gothic Medium';
            app.MESHLabel.FontColor = [0.0588 1 1];
            app.MESHLabel.Position = [389 211 36 22];
            app.MESHLabel.Text = 'MESH';

            % Create AxialdepthofcutmmLabel
            app.AxialdepthofcutmmLabel = uilabel(app.SimulationparametersPanel);
            app.AxialdepthofcutmmLabel.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutmmLabel.FontColor = [0.0588 1 1];
            app.AxialdepthofcutmmLabel.Position = [14 169 125 22];
            app.AxialdepthofcutmmLabel.Text = 'Axial depth of cut [mm]';

            % Create nxEditField
            app.nxEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.nxEditField.ValueChangedFcn = createCallbackFcn(app, @nxEditFieldValueChanged, true);
            app.nxEditField.FontName = 'Franklin Gothic Medium';
            app.nxEditField.Position = [358 176 33 29];
            app.nxEditField.Value = 21;

            % Create XLabel
            app.XLabel = uilabel(app.SimulationparametersPanel);
            app.XLabel.HorizontalAlignment = 'center';
            app.XLabel.FontName = 'Franklin Gothic Medium';
            app.XLabel.FontWeight = 'bold';
            app.XLabel.FontColor = [0.0588 1 1];
            app.XLabel.Position = [394 179 25 22];
            app.XLabel.Text = 'X';

            % Create nyEditField
            app.nyEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.nyEditField.ValueChangedFcn = createCallbackFcn(app, @nyEditFieldValueChanged, true);
            app.nyEditField.FontName = 'Franklin Gothic Medium';
            app.nyEditField.Position = [423 176 33 29];
            app.nyEditField.Value = 21;

            % Create ap_minEditFieldLabel
            app.ap_minEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.ap_minEditFieldLabel.HorizontalAlignment = 'right';
            app.ap_minEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.ap_minEditFieldLabel.FontColor = [0.0588 1 1];
            app.ap_minEditFieldLabel.Position = [10 138 44 22];
            app.ap_minEditFieldLabel.Text = 'ap_min';

            % Create ap_minEditField
            app.ap_minEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.ap_minEditField.ValueChangedFcn = createCallbackFcn(app, @ap_minEditFieldValueChanged, true);
            app.ap_minEditField.Editable = 'off';
            app.ap_minEditField.FontName = 'Franklin Gothic Medium';
            app.ap_minEditField.Position = [67 135 47 29];
            app.ap_minEditField.Value = 0.001;

            % Create ap_maxEditFieldLabel
            app.ap_maxEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.ap_maxEditFieldLabel.HorizontalAlignment = 'right';
            app.ap_maxEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.ap_maxEditFieldLabel.FontColor = [0.0588 1 1];
            app.ap_maxEditFieldLabel.Position = [130 138 46 22];
            app.ap_maxEditFieldLabel.Text = 'ap_max';

            % Create ap_maxEditField
            app.ap_maxEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.ap_maxEditField.ValueChangedFcn = createCallbackFcn(app, @ap_maxEditFieldValueChanged, true);
            app.ap_maxEditField.FontName = 'Franklin Gothic Medium';
            app.ap_maxEditField.Position = [189 135 47 29];
            app.ap_maxEditField.Value = 6;

            % Create nxEditFieldLabel
            app.nxEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.nxEditFieldLabel.HorizontalAlignment = 'center';
            app.nxEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.nxEditFieldLabel.FontColor = [0.0588 1 1];
            app.nxEditFieldLabel.Position = [362 152 25 22];
            app.nxEditFieldLabel.Text = '(nx)';

            % Create nyEditFieldLabel
            app.nyEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.nyEditFieldLabel.HorizontalAlignment = 'center';
            app.nyEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.nyEditFieldLabel.FontColor = [0.0588 1 1];
            app.nyEditFieldLabel.Position = [428 152 25 22];
            app.nyEditFieldLabel.Text = '(ny)';

            % Create OthersimulationparametersLabel
            app.OthersimulationparametersLabel = uilabel(app.SimulationparametersPanel);
            app.OthersimulationparametersLabel.FontName = 'Franklin Gothic Medium';
            app.OthersimulationparametersLabel.FontColor = [0.0588 1 1];
            app.OthersimulationparametersLabel.Position = [13 75 153 22];
            app.OthersimulationparametersLabel.Text = 'Other simulation parameters';

            % Create IsolinesEditFieldLabel
            app.IsolinesEditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.IsolinesEditFieldLabel.HorizontalAlignment = 'right';
            app.IsolinesEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.IsolinesEditFieldLabel.FontColor = [0.0588 1 1];
            app.IsolinesEditFieldLabel.Position = [13 40 44 22];
            app.IsolinesEditFieldLabel.Text = 'Isolines';

            % Create IsolinesEditField
            app.IsolinesEditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.IsolinesEditField.ValueChangedFcn = createCallbackFcn(app, @IsolinesEditFieldValueChanged, true);
            app.IsolinesEditField.FontName = 'Franklin Gothic Medium';
            app.IsolinesEditField.Position = [64 33 39 36];
            app.IsolinesEditField.Value = 100;

            % Create Periodstosimulate10EditFieldLabel
            app.Periodstosimulate10EditFieldLabel = uilabel(app.SimulationparametersPanel);
            app.Periodstosimulate10EditFieldLabel.HorizontalAlignment = 'right';
            app.Periodstosimulate10EditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.Periodstosimulate10EditFieldLabel.FontColor = [0.0588 1 1];
            app.Periodstosimulate10EditFieldLabel.Position = [121 40 135 22];
            app.Periodstosimulate10EditFieldLabel.Text = 'Periods to simulate (>10)';

            % Create Periodstosimulate10EditField
            app.Periodstosimulate10EditField = uieditfield(app.SimulationparametersPanel, 'numeric');
            app.Periodstosimulate10EditField.ValueChangedFcn = createCallbackFcn(app, @Periodstosimulate10EditFieldValueChanged, true);
            app.Periodstosimulate10EditField.FontName = 'Franklin Gothic Medium';
            app.Periodstosimulate10EditField.Position = [267 33 39 36];
            app.Periodstosimulate10EditField.Value = 30;

            % Create SimulationprogressGauge
            app.SimulationprogressGauge = uigauge(app.MillingoperationTab, 'circular');
            app.SimulationprogressGauge.FontName = 'Franklin Gothic Medium';
            app.SimulationprogressGauge.Position = [1145 73 172 172];

            % Create SimulationLabel
            app.SimulationLabel = uilabel(app.MillingoperationTab);
            app.SimulationLabel.HorizontalAlignment = 'center';
            app.SimulationLabel.FontName = 'Franklin Gothic Medium';
            app.SimulationLabel.FontSize = 18;
            app.SimulationLabel.FontWeight = 'bold';
            app.SimulationLabel.FontColor = [0.0588 1 1];
            app.SimulationLabel.Position = [1142 43 178 24];
            app.SimulationLabel.Text = 'Simulation progress %';

            % Create AStabilitychartsTab
            app.AStabilitychartsTab = uitab(app.TabGroup);
            app.AStabilitychartsTab.Title = 'A) Stability charts';
            app.AStabilitychartsTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create Panel_2
            app.Panel_2 = uipanel(app.AStabilitychartsTab);
            app.Panel_2.Position = [16 11 697 631];

            % Create UIAxes1
            app.UIAxes1 = uiaxes(app.Panel_2);
            title(app.UIAxes1, 'Title')
            xlabel(app.UIAxes1, 'X')
            ylabel(app.UIAxes1, 'Y')
            zlabel(app.UIAxes1, 'Z')
            app.UIAxes1.ButtonDownFcn = createCallbackFcn(app, @UIAxes1ButtonDown, true);
            app.UIAxes1.Position = [17 24 587 575];

            % Create p2pXButton
            app.p2pXButton = uibutton(app.AStabilitychartsTab, 'push');
            app.p2pXButton.ButtonPushedFcn = createCallbackFcn(app, @p2pXButtonPushed, true);
            app.p2pXButton.BackgroundColor = [0.2314 0.502 0.502];
            app.p2pXButton.FontName = 'Franklin Gothic Medium';
            app.p2pXButton.FontColor = [0.0588 1 1];
            app.p2pXButton.Position = [642 202 58 30];
            app.p2pXButton.Text = 'p2pX';

            % Create p2pYButton
            app.p2pYButton = uibutton(app.AStabilitychartsTab, 'push');
            app.p2pYButton.ButtonPushedFcn = createCallbackFcn(app, @p2pYButtonPushed, true);
            app.p2pYButton.BackgroundColor = [0.2314 0.502 0.502];
            app.p2pYButton.FontName = 'Franklin Gothic Medium';
            app.p2pYButton.FontColor = [0.0588 1 1];
            app.p2pYButton.Position = [642 154 58 30];
            app.p2pYButton.Text = 'p2pY';

            % Create p2pXYButton
            app.p2pXYButton = uibutton(app.AStabilitychartsTab, 'push');
            app.p2pXYButton.ButtonPushedFcn = createCallbackFcn(app, @p2pXYButtonPushed, true);
            app.p2pXYButton.BackgroundColor = [0.2314 0.502 0.502];
            app.p2pXYButton.FontName = 'Franklin Gothic Medium';
            app.p2pXYButton.FontColor = [0 1 1];
            app.p2pXYButton.Position = [642 106 58 30];
            app.p2pXYButton.Text = 'p2pXY';

            % Create RaButton
            app.RaButton = uibutton(app.AStabilitychartsTab, 'push');
            app.RaButton.ButtonPushedFcn = createCallbackFcn(app, @RaButtonPushed, true);
            app.RaButton.BackgroundColor = [0.2314 0.502 0.502];
            app.RaButton.FontName = 'Franklin Gothic Medium';
            app.RaButton.FontColor = [0.0588 1 1];
            app.RaButton.Position = [642 54 58 30];
            app.RaButton.Text = 'Ra';

            % Create Panel_3
            app.Panel_3 = uipanel(app.AStabilitychartsTab);
            app.Panel_3.Position = [726 11 633 631];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.Panel_3);
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.ButtonDownFcn = createCallbackFcn(app, @UIAxes2ButtonDown, true);
            app.UIAxes2.Position = [15 24 585 574];

            % Create BDynamicforcesTab
            app.BDynamicforcesTab = uitab(app.TabGroup);
            app.BDynamicforcesTab.Title = 'B) Dynamic forces';
            app.BDynamicforcesTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create UIAxes3_2
            app.UIAxes3_2 = uiaxes(app.BDynamicforcesTab);
            xlabel(app.UIAxes3_2, 'X')
            ylabel(app.UIAxes3_2, 'Y')
            zlabel(app.UIAxes3_2, 'Z')
            app.UIAxes3_2.Position = [42 177 969 170];

            % Create UIAxes3_3
            app.UIAxes3_3 = uiaxes(app.BDynamicforcesTab);
            xlabel(app.UIAxes3_3, 'X')
            ylabel(app.UIAxes3_3, 'Y')
            zlabel(app.UIAxes3_3, 'Z')
            app.UIAxes3_3.Position = [41 14 969 166];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.BDynamicforcesTab);
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Position = [42 336 969 309];

            % Create CuttingandsimulationparametersPanel
            app.CuttingandsimulationparametersPanel = uipanel(app.BDynamicforcesTab);
            app.CuttingandsimulationparametersPanel.ForegroundColor = [0.0588 1 1];
            app.CuttingandsimulationparametersPanel.Title = 'Cutting and simulation parameters';
            app.CuttingandsimulationparametersPanel.BackgroundColor = [0.2314 0.502 0.502];
            app.CuttingandsimulationparametersPanel.FontName = 'Franklin Gothic Medium';
            app.CuttingandsimulationparametersPanel.FontWeight = 'bold';
            app.CuttingandsimulationparametersPanel.Position = [1053 24 310 603];

            % Create FeedpertoothfzmmZrevEditField_2Label
            app.FeedpertoothfzmmZrevEditField_2Label = uilabel(app.CuttingandsimulationparametersPanel);
            app.FeedpertoothfzmmZrevEditField_2Label.HorizontalAlignment = 'right';
            app.FeedpertoothfzmmZrevEditField_2Label.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_2Label.FontColor = [0.0588 1 1];
            app.FeedpertoothfzmmZrevEditField_2Label.Position = [40 538 155 22];
            app.FeedpertoothfzmmZrevEditField_2Label.Text = 'Feed per tooth fz [mm/Z/rev]';

            % Create FeedpertoothfzmmZrevEditField_2
            app.FeedpertoothfzmmZrevEditField_2 = uieditfield(app.CuttingandsimulationparametersPanel, 'numeric');
            app.FeedpertoothfzmmZrevEditField_2.ValueChangedFcn = createCallbackFcn(app, @FeedpertoothfzmmZrevEditField_2ValueChanged, true);
            app.FeedpertoothfzmmZrevEditField_2.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_2.Position = [221 532 39 35];
            app.FeedpertoothfzmmZrevEditField_2.Value = 0.1;

            % Create RadialdepthofcutaemmEditField_2Label
            app.RadialdepthofcutaemmEditField_2Label = uilabel(app.CuttingandsimulationparametersPanel);
            app.RadialdepthofcutaemmEditField_2Label.HorizontalAlignment = 'right';
            app.RadialdepthofcutaemmEditField_2Label.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_2Label.FontColor = [0.0588 1 1];
            app.RadialdepthofcutaemmEditField_2Label.Position = [40 481 149 22];
            app.RadialdepthofcutaemmEditField_2Label.Text = 'Radial depth of cut ae [mm]';

            % Create RadialdepthofcutaemmEditField_2
            app.RadialdepthofcutaemmEditField_2 = uieditfield(app.CuttingandsimulationparametersPanel, 'numeric');
            app.RadialdepthofcutaemmEditField_2.ValueChangedFcn = createCallbackFcn(app, @RadialdepthofcutaemmEditField_2ValueChanged, true);
            app.RadialdepthofcutaemmEditField_2.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_2.Position = [221 475 39 35];
            app.RadialdepthofcutaemmEditField_2.Value = 6;

            % Create AxialdepthofcutapmmEditField_2Label
            app.AxialdepthofcutapmmEditField_2Label = uilabel(app.CuttingandsimulationparametersPanel);
            app.AxialdepthofcutapmmEditField_2Label.HorizontalAlignment = 'right';
            app.AxialdepthofcutapmmEditField_2Label.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_2Label.FontColor = [0.0588 1 1];
            app.AxialdepthofcutapmmEditField_2Label.Position = [40 424 141 22];
            app.AxialdepthofcutapmmEditField_2Label.Text = 'Axial depth of cut ap [mm]';

            % Create AxialdepthofcutapmmEditField_2
            app.AxialdepthofcutapmmEditField_2 = uieditfield(app.CuttingandsimulationparametersPanel, 'numeric');
            app.AxialdepthofcutapmmEditField_2.ValueChangedFcn = createCallbackFcn(app, @AxialdepthofcutapmmEditField_2ValueChanged, true);
            app.AxialdepthofcutapmmEditField_2.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_2.Position = [221 418 39 35];
            app.AxialdepthofcutapmmEditField_2.Value = 6;

            % Create SpindlespeednrpmEditField_2Label
            app.SpindlespeednrpmEditField_2Label = uilabel(app.CuttingandsimulationparametersPanel);
            app.SpindlespeednrpmEditField_2Label.HorizontalAlignment = 'right';
            app.SpindlespeednrpmEditField_2Label.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_2Label.FontColor = [0.0588 1 1];
            app.SpindlespeednrpmEditField_2Label.Position = [40 368 118 22];
            app.SpindlespeednrpmEditField_2Label.Text = 'Spindle speed n [rpm]';

            % Create SpindlespeednrpmEditField_2
            app.SpindlespeednrpmEditField_2 = uieditfield(app.CuttingandsimulationparametersPanel, 'numeric');
            app.SpindlespeednrpmEditField_2.ValueChangedFcn = createCallbackFcn(app, @SpindlespeednrpmEditField_2ValueChanged, true);
            app.SpindlespeednrpmEditField_2.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_2.Position = [221 361 39 35];
            app.SpindlespeednrpmEditField_2.Value = 8000;

            % Create PeriodstosimulateT10EditFieldLabel
            app.PeriodstosimulateT10EditFieldLabel = uilabel(app.CuttingandsimulationparametersPanel);
            app.PeriodstosimulateT10EditFieldLabel.HorizontalAlignment = 'right';
            app.PeriodstosimulateT10EditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.PeriodstosimulateT10EditFieldLabel.FontColor = [0.0588 1 1];
            app.PeriodstosimulateT10EditFieldLabel.Position = [40 306 144 22];
            app.PeriodstosimulateT10EditFieldLabel.Text = 'Periods to simulate T (>10)';

            % Create PeriodstosimulateT10EditField
            app.PeriodstosimulateT10EditField = uieditfield(app.CuttingandsimulationparametersPanel, 'numeric');
            app.PeriodstosimulateT10EditField.ValueChangedFcn = createCallbackFcn(app, @PeriodstosimulateT10EditFieldValueChanged, true);
            app.PeriodstosimulateT10EditField.FontName = 'Franklin Gothic Medium';
            app.PeriodstosimulateT10EditField.Position = [221 299 39 36];
            app.PeriodstosimulateT10EditField.Value = 30;

            % Create RUNButton_2
            app.RUNButton_2 = uibutton(app.CuttingandsimulationparametersPanel, 'push');
            app.RUNButton_2.ButtonPushedFcn = createCallbackFcn(app, @RUNButton_2Pushed, true);
            app.RUNButton_2.BackgroundColor = [0.1294 0.3098 0.3098];
            app.RUNButton_2.FontName = 'Franklin Gothic Medium';
            app.RUNButton_2.FontSize = 48;
            app.RUNButton_2.FontWeight = 'bold';
            app.RUNButton_2.FontColor = [0.0588 1 1];
            app.RUNButton_2.Position = [87 174 141 79];
            app.RUNButton_2.Text = 'RUN';

            % Create ChatterButton
            app.ChatterButton = uibutton(app.CuttingandsimulationparametersPanel, 'push');
            app.ChatterButton.ButtonPushedFcn = createCallbackFcn(app, @ChatterButtonPushed, true);
            app.ChatterButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ChatterButton.FontName = 'Franklin Gothic Medium';
            app.ChatterButton.FontSize = 24;
            app.ChatterButton.FontWeight = 'bold';
            app.ChatterButton.FontColor = [0.0588 1 1];
            app.ChatterButton.Position = [98 90 121 50];
            app.ChatterButton.Text = 'Chatter ?';

            % Create ClearandreleaseButton
            app.ClearandreleaseButton = uibutton(app.CuttingandsimulationparametersPanel, 'push');
            app.ClearandreleaseButton.ButtonPushedFcn = createCallbackFcn(app, @ClearandreleaseButtonPushed, true);
            app.ClearandreleaseButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ClearandreleaseButton.FontName = 'Franklin Gothic Medium';
            app.ClearandreleaseButton.FontSize = 24;
            app.ClearandreleaseButton.FontWeight = 'bold';
            app.ClearandreleaseButton.FontColor = [0.0588 1 1];
            app.ClearandreleaseButton.Position = [17 13 276 41];
            app.ClearandreleaseButton.Text = 'Clear and release';

            % Create COptimizationtoolTab
            app.COptimizationtoolTab = uitab(app.TabGroup);
            app.COptimizationtoolTab.SizeChangedFcn = createCallbackFcn(app, @COptimizationtoolTabSizeChanged, true);
            app.COptimizationtoolTab.Title = 'C) Optimization tool';

            % Create OptimizationobjectivePanel
            app.OptimizationobjectivePanel = uipanel(app.COptimizationtoolTab);
            app.OptimizationobjectivePanel.ForegroundColor = [0.0588 1 1];
            app.OptimizationobjectivePanel.Title = 'Optimization objective';
            app.OptimizationobjectivePanel.BackgroundColor = [0.2314 0.502 0.502];
            app.OptimizationobjectivePanel.FontWeight = 'bold';
            app.OptimizationobjectivePanel.Position = [15 10 545 632];

            % Create StartingpointDropDown
            app.StartingpointDropDown = uidropdown(app.OptimizationobjectivePanel);
            app.StartingpointDropDown.Items = {'Maximize productivity', 'Minimize surface roughness', 'Productivity and surface accuracy'};
            app.StartingpointDropDown.ValueChangedFcn = createCallbackFcn(app, @StartingpointDropDownValueChanged, true);
            app.StartingpointDropDown.FontName = 'Franklin Gothic Medium';
            app.StartingpointDropDown.FontSize = 14;
            app.StartingpointDropDown.Position = [8 577 206 25];
            app.StartingpointDropDown.Value = 'Maximize productivity';

            % Create Label_7
            app.Label_7 = uilabel(app.OptimizationobjectivePanel);
            app.Label_7.FontName = 'Franklin Gothic Medium';
            app.Label_7.FontColor = [0.0588 1 1];
            app.Label_7.Position = [13 587 75 25];
            app.Label_7.Text = '';

            % Create OPTIMUMButton
            app.OPTIMUMButton = uibutton(app.OptimizationobjectivePanel, 'push');
            app.OPTIMUMButton.ButtonPushedFcn = createCallbackFcn(app, @OPTIMUMButtonPushed, true);
            app.OPTIMUMButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.OPTIMUMButton.FontSize = 36;
            app.OPTIMUMButton.FontWeight = 'bold';
            app.OPTIMUMButton.FontColor = [0.0588 1 1];
            app.OPTIMUMButton.Position = [42 233 199 66];
            app.OPTIMUMButton.Text = 'OPTIMUM';

            % Create OptimumpointPanel
            app.OptimumpointPanel = uipanel(app.OptimizationobjectivePanel);
            app.OptimumpointPanel.AutoResizeChildren = 'off';
            app.OptimumpointPanel.ForegroundColor = [0 1 0];
            app.OptimumpointPanel.Title = '2. Optimum point';
            app.OptimumpointPanel.BackgroundColor = [0.1294 0.3098 0.3098];
            app.OptimumpointPanel.FontName = 'Franklin Gothic Medium';
            app.OptimumpointPanel.FontWeight = 'bold';
            app.OptimumpointPanel.Position = [288 7 233 211];

            % Create FeedpertoothfzmmZrevEditField_5Label
            app.FeedpertoothfzmmZrevEditField_5Label = uilabel(app.OptimumpointPanel);
            app.FeedpertoothfzmmZrevEditField_5Label.HorizontalAlignment = 'right';
            app.FeedpertoothfzmmZrevEditField_5Label.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_5Label.FontColor = [0 1 0];
            app.FeedpertoothfzmmZrevEditField_5Label.Position = [5 148 155 22];
            app.FeedpertoothfzmmZrevEditField_5Label.Text = 'Feed per tooth fz [mm/Z/rev]';

            % Create FeedpertoothfzmmZrevEditField_5
            app.FeedpertoothfzmmZrevEditField_5 = uieditfield(app.OptimumpointPanel, 'numeric');
            app.FeedpertoothfzmmZrevEditField_5.ValueChangedFcn = createCallbackFcn(app, @FeedpertoothfzmmZrevEditField_5ValueChanged, true);
            app.FeedpertoothfzmmZrevEditField_5.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_5.FontSize = 18;
            app.FeedpertoothfzmmZrevEditField_5.FontColor = [0 1 0];
            app.FeedpertoothfzmmZrevEditField_5.BackgroundColor = [0.1294 0.3098 0.3098];
            app.FeedpertoothfzmmZrevEditField_5.Enable = 'off';
            app.FeedpertoothfzmmZrevEditField_5.Position = [171 142 53 35];

            % Create RadialdepthofcutaemmEditField_5Label
            app.RadialdepthofcutaemmEditField_5Label = uilabel(app.OptimumpointPanel);
            app.RadialdepthofcutaemmEditField_5Label.HorizontalAlignment = 'right';
            app.RadialdepthofcutaemmEditField_5Label.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_5Label.FontColor = [0 1 0];
            app.RadialdepthofcutaemmEditField_5Label.Position = [6 103 149 22];
            app.RadialdepthofcutaemmEditField_5Label.Text = 'Radial depth of cut ae [mm]';

            % Create RadialdepthofcutaemmEditField_5
            app.RadialdepthofcutaemmEditField_5 = uieditfield(app.OptimumpointPanel, 'numeric');
            app.RadialdepthofcutaemmEditField_5.ValueChangedFcn = createCallbackFcn(app, @RadialdepthofcutaemmEditField_5ValueChanged, true);
            app.RadialdepthofcutaemmEditField_5.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_5.FontSize = 18;
            app.RadialdepthofcutaemmEditField_5.FontColor = [0 1 0];
            app.RadialdepthofcutaemmEditField_5.BackgroundColor = [0.1294 0.3098 0.3098];
            app.RadialdepthofcutaemmEditField_5.Enable = 'off';
            app.RadialdepthofcutaemmEditField_5.Position = [171 97 53 35];

            % Create AxialdepthofcutapmmEditField_7Label
            app.AxialdepthofcutapmmEditField_7Label = uilabel(app.OptimumpointPanel);
            app.AxialdepthofcutapmmEditField_7Label.HorizontalAlignment = 'right';
            app.AxialdepthofcutapmmEditField_7Label.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_7Label.FontColor = [0 1 0];
            app.AxialdepthofcutapmmEditField_7Label.Position = [7 58 141 22];
            app.AxialdepthofcutapmmEditField_7Label.Text = 'Axial depth of cut ap [mm]';

            % Create AxialdepthofcutapmmEditField_7
            app.AxialdepthofcutapmmEditField_7 = uieditfield(app.OptimumpointPanel, 'numeric');
            app.AxialdepthofcutapmmEditField_7.ValueChangedFcn = createCallbackFcn(app, @AxialdepthofcutapmmEditField_7ValueChanged, true);
            app.AxialdepthofcutapmmEditField_7.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_7.FontSize = 18;
            app.AxialdepthofcutapmmEditField_7.FontColor = [0 1 0];
            app.AxialdepthofcutapmmEditField_7.BackgroundColor = [0.1294 0.3098 0.3098];
            app.AxialdepthofcutapmmEditField_7.Enable = 'off';
            app.AxialdepthofcutapmmEditField_7.Position = [171 52 53 35];

            % Create SpindlespeednrpmEditField_7Label
            app.SpindlespeednrpmEditField_7Label = uilabel(app.OptimumpointPanel);
            app.SpindlespeednrpmEditField_7Label.HorizontalAlignment = 'right';
            app.SpindlespeednrpmEditField_7Label.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_7Label.FontColor = [0 1 0];
            app.SpindlespeednrpmEditField_7Label.Position = [6 13 118 22];
            app.SpindlespeednrpmEditField_7Label.Text = 'Spindle speed n [rpm]';

            % Create SpindlespeednrpmEditField_7
            app.SpindlespeednrpmEditField_7 = uieditfield(app.OptimumpointPanel, 'numeric');
            app.SpindlespeednrpmEditField_7.ValueChangedFcn = createCallbackFcn(app, @SpindlespeednrpmEditField_7ValueChanged, true);
            app.SpindlespeednrpmEditField_7.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_7.FontSize = 18;
            app.SpindlespeednrpmEditField_7.FontColor = [0 1 0];
            app.SpindlespeednrpmEditField_7.BackgroundColor = [0.1294 0.3098 0.3098];
            app.SpindlespeednrpmEditField_7.Enable = 'off';
            app.SpindlespeednrpmEditField_7.Position = [171 7 53 35];

            % Create StartingpointPanel
            app.StartingpointPanel = uipanel(app.OptimizationobjectivePanel);
            app.StartingpointPanel.AutoResizeChildren = 'off';
            app.StartingpointPanel.ForegroundColor = [1 0.4118 0.1608];
            app.StartingpointPanel.Title = '1. Starting point';
            app.StartingpointPanel.BackgroundColor = [0.1294 0.3098 0.3098];
            app.StartingpointPanel.FontName = 'Franklin Gothic Medium';
            app.StartingpointPanel.FontWeight = 'bold';
            app.StartingpointPanel.Position = [25 7 233 211];

            % Create FeedpertoothfzmmZrevEditField_4Label
            app.FeedpertoothfzmmZrevEditField_4Label = uilabel(app.StartingpointPanel);
            app.FeedpertoothfzmmZrevEditField_4Label.HorizontalAlignment = 'right';
            app.FeedpertoothfzmmZrevEditField_4Label.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_4Label.FontColor = [1 0.4118 0.1608];
            app.FeedpertoothfzmmZrevEditField_4Label.Position = [5 147 155 22];
            app.FeedpertoothfzmmZrevEditField_4Label.Text = 'Feed per tooth fz [mm/Z/rev]';

            % Create FeedpertoothfzmmZrevEditField_4
            app.FeedpertoothfzmmZrevEditField_4 = uieditfield(app.StartingpointPanel, 'numeric');
            app.FeedpertoothfzmmZrevEditField_4.ValueChangedFcn = createCallbackFcn(app, @FeedpertoothfzmmZrevEditField_4ValueChanged, true);
            app.FeedpertoothfzmmZrevEditField_4.FontName = 'Franklin Gothic Medium';
            app.FeedpertoothfzmmZrevEditField_4.FontSize = 18;
            app.FeedpertoothfzmmZrevEditField_4.FontColor = [1 0.4118 0.1608];
            app.FeedpertoothfzmmZrevEditField_4.BackgroundColor = [0.1294 0.3098 0.3098];
            app.FeedpertoothfzmmZrevEditField_4.Enable = 'off';
            app.FeedpertoothfzmmZrevEditField_4.Position = [168 141 56 35];
            app.FeedpertoothfzmmZrevEditField_4.Value = 0.1;

            % Create RadialdepthofcutaemmEditField_4Label
            app.RadialdepthofcutaemmEditField_4Label = uilabel(app.StartingpointPanel);
            app.RadialdepthofcutaemmEditField_4Label.HorizontalAlignment = 'right';
            app.RadialdepthofcutaemmEditField_4Label.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_4Label.FontColor = [1 0.4118 0.1608];
            app.RadialdepthofcutaemmEditField_4Label.Position = [6 102 149 22];
            app.RadialdepthofcutaemmEditField_4Label.Text = 'Radial depth of cut ae [mm]';

            % Create RadialdepthofcutaemmEditField_4
            app.RadialdepthofcutaemmEditField_4 = uieditfield(app.StartingpointPanel, 'numeric');
            app.RadialdepthofcutaemmEditField_4.ValueChangedFcn = createCallbackFcn(app, @RadialdepthofcutaemmEditField_4ValueChanged, true);
            app.RadialdepthofcutaemmEditField_4.FontName = 'Franklin Gothic Medium';
            app.RadialdepthofcutaemmEditField_4.FontSize = 18;
            app.RadialdepthofcutaemmEditField_4.FontColor = [1 0.4118 0.1608];
            app.RadialdepthofcutaemmEditField_4.BackgroundColor = [0.1294 0.3098 0.3098];
            app.RadialdepthofcutaemmEditField_4.Enable = 'off';
            app.RadialdepthofcutaemmEditField_4.Position = [168 96 56 35];
            app.RadialdepthofcutaemmEditField_4.Value = 6;

            % Create AxialdepthofcutapmmEditField_6Label
            app.AxialdepthofcutapmmEditField_6Label = uilabel(app.StartingpointPanel);
            app.AxialdepthofcutapmmEditField_6Label.HorizontalAlignment = 'right';
            app.AxialdepthofcutapmmEditField_6Label.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_6Label.FontColor = [1 0.4118 0.1608];
            app.AxialdepthofcutapmmEditField_6Label.Position = [7 57 141 22];
            app.AxialdepthofcutapmmEditField_6Label.Text = 'Axial depth of cut ap [mm]';

            % Create AxialdepthofcutapmmEditField_6
            app.AxialdepthofcutapmmEditField_6 = uieditfield(app.StartingpointPanel, 'numeric');
            app.AxialdepthofcutapmmEditField_6.ValueChangedFcn = createCallbackFcn(app, @AxialdepthofcutapmmEditField_6ValueChanged, true);
            app.AxialdepthofcutapmmEditField_6.FontName = 'Franklin Gothic Medium';
            app.AxialdepthofcutapmmEditField_6.FontSize = 18;
            app.AxialdepthofcutapmmEditField_6.FontColor = [1 0.4118 0.1608];
            app.AxialdepthofcutapmmEditField_6.BackgroundColor = [0.1294 0.3098 0.3098];
            app.AxialdepthofcutapmmEditField_6.Enable = 'off';
            app.AxialdepthofcutapmmEditField_6.Position = [167 51 56 35];

            % Create SpindlespeednrpmEditField_6Label
            app.SpindlespeednrpmEditField_6Label = uilabel(app.StartingpointPanel);
            app.SpindlespeednrpmEditField_6Label.HorizontalAlignment = 'right';
            app.SpindlespeednrpmEditField_6Label.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_6Label.FontColor = [1 0.4118 0.1608];
            app.SpindlespeednrpmEditField_6Label.Position = [8 12 118 22];
            app.SpindlespeednrpmEditField_6Label.Text = 'Spindle speed n [rpm]';

            % Create SpindlespeednrpmEditField_6
            app.SpindlespeednrpmEditField_6 = uieditfield(app.StartingpointPanel, 'numeric');
            app.SpindlespeednrpmEditField_6.ValueChangedFcn = createCallbackFcn(app, @SpindlespeednrpmEditField_6ValueChanged, true);
            app.SpindlespeednrpmEditField_6.FontName = 'Franklin Gothic Medium';
            app.SpindlespeednrpmEditField_6.FontSize = 18;
            app.SpindlespeednrpmEditField_6.FontColor = [1 0.4118 0.1608];
            app.SpindlespeednrpmEditField_6.BackgroundColor = [0.1294 0.3098 0.3098];
            app.SpindlespeednrpmEditField_6.Enable = 'off';
            app.SpindlespeednrpmEditField_6.Position = [167 6 56 35];

            % Create FxymeanFxySliderLabel
            app.FxymeanFxySliderLabel = uilabel(app.OptimizationobjectivePanel);
            app.FxymeanFxySliderLabel.HorizontalAlignment = 'right';
            app.FxymeanFxySliderLabel.FontName = 'Franklin Gothic Medium';
            app.FxymeanFxySliderLabel.FontSize = 10;
            app.FxymeanFxySliderLabel.FontColor = [0.302 0.749 0.9294];
            app.FxymeanFxySliderLabel.Position = [8 552 82 22];
            app.FxymeanFxySliderLabel.Text = '%∆Fxy/mean(Fxy)';

            % Create FxymeanFxySlider
            app.FxymeanFxySlider = uislider(app.OptimizationobjectivePanel);
            app.FxymeanFxySlider.Limits = [0 1000];
            app.FxymeanFxySlider.ValueChangedFcn = createCallbackFcn(app, @FxymeanFxySliderValueChanged, true);
            app.FxymeanFxySlider.FontName = 'Franklin Gothic Medium';
            app.FxymeanFxySlider.FontSize = 9;
            app.FxymeanFxySlider.FontColor = [0.302 0.749 0.9294];
            app.FxymeanFxySlider.Position = [18 542 443 7];

            % Create MaximumalowedRamumLabel
            app.MaximumalowedRamumLabel = uilabel(app.OptimizationobjectivePanel);
            app.MaximumalowedRamumLabel.HorizontalAlignment = 'right';
            app.MaximumalowedRamumLabel.FontName = 'Franklin Gothic Medium';
            app.MaximumalowedRamumLabel.FontSize = 10;
            app.MaximumalowedRamumLabel.FontColor = [0.8 0.6 0.902];
            app.MaximumalowedRamumLabel.Position = [8 492 85 22];
            app.MaximumalowedRamumLabel.Text = 'Maximum Ra [μm]';

            % Create MaximumRamSlider_2
            app.MaximumRamSlider_2 = uislider(app.OptimizationobjectivePanel);
            app.MaximumRamSlider_2.Limits = [0 50];
            app.MaximumRamSlider_2.ValueChangedFcn = createCallbackFcn(app, @MaximumRamSlider_2ValueChanged, true);
            app.MaximumRamSlider_2.FontName = 'Franklin Gothic Medium';
            app.MaximumRamSlider_2.FontSize = 10;
            app.MaximumRamSlider_2.FontColor = [0.8 0.6 0.902];
            app.MaximumRamSlider_2.Position = [18 482 442 7];

            % Create EditField
            app.EditField = uieditfield(app.OptimizationobjectivePanel, 'numeric');
            app.EditField.ValueChangedFcn = createCallbackFcn(app, @EditFieldValueChanged, true);
            app.EditField.Editable = 'off';
            app.EditField.FontName = 'Franklin Gothic Medium';
            app.EditField.FontSize = 18;
            app.EditField.FontColor = [0.302 0.7451 0.9333];
            app.EditField.BackgroundColor = [0.2314 0.502 0.502];
            app.EditField.Position = [471 531 69 35];

            % Create EditField_2
            app.EditField_2 = uieditfield(app.OptimizationobjectivePanel, 'numeric');
            app.EditField_2.ValueChangedFcn = createCallbackFcn(app, @EditField_2ValueChanged, true);
            app.EditField_2.Editable = 'off';
            app.EditField_2.FontName = 'Franklin Gothic Medium';
            app.EditField_2.FontSize = 18;
            app.EditField_2.FontColor = [0.8 0.6 0.902];
            app.EditField_2.BackgroundColor = [0.2314 0.502 0.502];
            app.EditField_2.Position = [471 473 69 35];

            % Create EditField_3
            app.EditField_3 = uieditfield(app.OptimizationobjectivePanel, 'numeric');
            app.EditField_3.ValueChangedFcn = createCallbackFcn(app, @EditField_3ValueChanged, true);
            app.EditField_3.Editable = 'off';
            app.EditField_3.FontName = 'Franklin Gothic Medium';
            app.EditField_3.FontSize = 18;
            app.EditField_3.FontColor = [0.302 0.749 0.9294];
            app.EditField_3.BackgroundColor = [0.2314 0.502 0.502];
            app.EditField_3.Position = [471 409 69 35];

            % Create EditField_4
            app.EditField_4 = uieditfield(app.OptimizationobjectivePanel, 'numeric');
            app.EditField_4.ValueChangedFcn = createCallbackFcn(app, @EditField_4ValueChanged, true);
            app.EditField_4.Editable = 'off';
            app.EditField_4.FontName = 'Franklin Gothic Medium';
            app.EditField_4.FontSize = 18;
            app.EditField_4.FontColor = [0.8 0.6 0.902];
            app.EditField_4.BackgroundColor = [0.2314 0.502 0.502];
            app.EditField_4.Position = [471 343 69 35];

            % Create MRRmm3minLabel
            app.MRRmm3minLabel = uilabel(app.OptimizationobjectivePanel);
            app.MRRmm3minLabel.HorizontalAlignment = 'right';
            app.MRRmm3minLabel.FontName = 'Franklin Gothic Medium';
            app.MRRmm3minLabel.FontSize = 10;
            app.MRRmm3minLabel.FontColor = [0.302 0.749 0.9294];
            app.MRRmm3minLabel.Position = [9 430 120 22];
            app.MRRmm3minLabel.Text = 'Minimum MRR [mm³/min]';

            % Create MinimumMRRmmminSlider
            app.MinimumMRRmmminSlider = uislider(app.OptimizationobjectivePanel);
            app.MinimumMRRmmminSlider.Limits = [0 100000];
            app.MinimumMRRmmminSlider.ValueChangedFcn = createCallbackFcn(app, @MinimumMRRmmminSliderValueChanged, true);
            app.MinimumMRRmmminSlider.FontName = 'Franklin Gothic Medium';
            app.MinimumMRRmmminSlider.FontSize = 10;
            app.MinimumMRRmmminSlider.FontColor = [0.302 0.749 0.9294];
            app.MinimumMRRmmminSlider.Position = [18 420 443 7];

            % Create MaximumRamSliderLabel
            app.MaximumRamSliderLabel = uilabel(app.OptimizationobjectivePanel);
            app.MaximumRamSliderLabel.HorizontalAlignment = 'right';
            app.MaximumRamSliderLabel.FontName = 'Franklin Gothic Medium';
            app.MaximumRamSliderLabel.FontSize = 10;
            app.MaximumRamSliderLabel.FontColor = [0.8 0.6 0.902];
            app.MaximumRamSliderLabel.Position = [7 365 85 22];
            app.MaximumRamSliderLabel.Text = 'Maximum Ra [μm]';

            % Create MaximumRamSlider
            app.MaximumRamSlider = uislider(app.OptimizationobjectivePanel);
            app.MaximumRamSlider.Limits = [0 50];
            app.MaximumRamSlider.ValueChangedFcn = createCallbackFcn(app, @MaximumRamSliderValueChanged, true);
            app.MaximumRamSlider.FontName = 'Franklin Gothic Medium';
            app.MaximumRamSlider.FontSize = 10;
            app.MaximumRamSlider.FontColor = [0.8 0.6 0.902];
            app.MaximumRamSlider.Position = [17 355 443 7];

            % Create ClearfigureButton
            app.ClearfigureButton = uibutton(app.OptimizationobjectivePanel, 'push');
            app.ClearfigureButton.ButtonPushedFcn = createCallbackFcn(app, @ClearfigureButtonPushed, true);
            app.ClearfigureButton.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ClearfigureButton.FontName = 'Franklin Gothic Medium';
            app.ClearfigureButton.FontSize = 14;
            app.ClearfigureButton.FontWeight = 'bold';
            app.ClearfigureButton.FontColor = [0.0588 1 1];
            app.ClearfigureButton.Position = [314 247 180 40];
            app.ClearfigureButton.Text = 'Clear figure';

            % Create ChoosecriterionLabel
            app.ChoosecriterionLabel = uilabel(app.OptimizationobjectivePanel);
            app.ChoosecriterionLabel.FontAngle = 'italic';
            app.ChoosecriterionLabel.FontColor = [0.0588 1 1];
            app.ChoosecriterionLabel.Position = [219 579 92 22];
            app.ChoosecriterionLabel.Text = 'Choose criterion';

            % Create Panel_4
            app.Panel_4 = uipanel(app.COptimizationtoolTab);
            app.Panel_4.Position = [575 11 788 631];

            % Create UIAxes7
            app.UIAxes7 = uiaxes(app.Panel_4);
            xlabel(app.UIAxes7, 'X')
            ylabel(app.UIAxes7, 'Y')
            zlabel(app.UIAxes7, 'Z')
            app.UIAxes7.ButtonDownFcn = createCallbackFcn(app, @UIAxes7ButtonDown, true);
            app.UIAxes7.Position = [25 151 727 450];

            % Create Panel_5
            app.Panel_5 = uipanel(app.Panel_4);
            app.Panel_5.BackgroundColor = [0.1294 0.3098 0.3098];
            app.Panel_5.SizeChangedFcn = createCallbackFcn(app, @Panel_5SizeChanged, true);
            app.Panel_5.Position = [48 9 692 136];

            % Create Label_5
            app.Label_5 = uilabel(app.Panel_5);
            app.Label_5.FontSize = 48;
            app.Label_5.FontWeight = 'bold';
            app.Label_5.FontColor = [0.0588 1 1];
            app.Label_5.Position = [640 64 48 63];
            app.Label_5.Text = '%';

            % Create ProductivityimprovedbyEditFieldLabel_2
            app.ProductivityimprovedbyEditFieldLabel_2 = uilabel(app.Panel_5);
            app.ProductivityimprovedbyEditFieldLabel_2.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_2.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_2.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_2.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_2.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_2.Position = [13 14 190 47];
            app.ProductivityimprovedbyEditFieldLabel_2.Text = 'MRR   from:';

            % Create ProductivityimprovedbyEditField_2
            app.ProductivityimprovedbyEditField_2 = uieditfield(app.Panel_5, 'numeric');
            app.ProductivityimprovedbyEditField_2.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_2ValueChanged, true);
            app.ProductivityimprovedbyEditField_2.FontSize = 36;
            app.ProductivityimprovedbyEditField_2.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_2.FontColor = [1 0.4118 0.1608];
            app.ProductivityimprovedbyEditField_2.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_2.Position = [208 8 167 49];

            % Create ProductivityimprovedbyEditFieldLabel_3
            app.ProductivityimprovedbyEditFieldLabel_3 = uilabel(app.Panel_5);
            app.ProductivityimprovedbyEditFieldLabel_3.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_3.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_3.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_3.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_3.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_3.Position = [372 13 44 47];
            app.ProductivityimprovedbyEditFieldLabel_3.Text = 'to:';

            % Create ProductivityimprovedbyEditField_3
            app.ProductivityimprovedbyEditField_3 = uieditfield(app.Panel_5, 'numeric');
            app.ProductivityimprovedbyEditField_3.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_3ValueChanged, true);
            app.ProductivityimprovedbyEditField_3.FontSize = 36;
            app.ProductivityimprovedbyEditField_3.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_3.FontColor = [0 1 0];
            app.ProductivityimprovedbyEditField_3.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_3.Position = [418 8 167 49];

            % Create cmminLabel
            app.cmminLabel = uilabel(app.Panel_5);
            app.cmminLabel.FontName = 'Franklin Gothic Medium';
            app.cmminLabel.FontSize = 20;
            app.cmminLabel.FontAngle = 'italic';
            app.cmminLabel.FontColor = [0.0588 1 1];
            app.cmminLabel.Position = [592 19 81 26];
            app.cmminLabel.Text = 'cm³/min';

            % Create ProductivityimprovedbyEditFieldLabel
            app.ProductivityimprovedbyEditFieldLabel = uilabel(app.Panel_5);
            app.ProductivityimprovedbyEditFieldLabel.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditFieldLabel.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel.Position = [40 69 398 47];
            app.ProductivityimprovedbyEditFieldLabel.Text = 'Productivity improved by:';

            % Create ProductivityimprovedbyEditField
            app.ProductivityimprovedbyEditField = uieditfield(app.Panel_5, 'numeric');
            app.ProductivityimprovedbyEditField.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditFieldValueChanged, true);
            app.ProductivityimprovedbyEditField.FontSize = 50;
            app.ProductivityimprovedbyEditField.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditField.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField.Position = [447 60 186 64];

            % Create Panel_6
            app.Panel_6 = uipanel(app.Panel_4);
            app.Panel_6.BackgroundColor = [0.1294 0.3098 0.3098];
            app.Panel_6.SizeChangedFcn = createCallbackFcn(app, @Panel_6SizeChanged, true);
            app.Panel_6.Position = [48 9 692 136];

            % Create Label_6
            app.Label_6 = uilabel(app.Panel_6);
            app.Label_6.FontSize = 48;
            app.Label_6.FontWeight = 'bold';
            app.Label_6.FontColor = [0.0588 1 1];
            app.Label_6.Position = [640 65 48 63];
            app.Label_6.Text = '%';

            % Create ProductivityimprovedbyEditFieldLabel_4
            app.ProductivityimprovedbyEditFieldLabel_4 = uilabel(app.Panel_6);
            app.ProductivityimprovedbyEditFieldLabel_4.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_4.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_4.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_4.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_4.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_4.Position = [47 14 156 47];
            app.ProductivityimprovedbyEditFieldLabel_4.Text = 'Ra   from:';

            % Create ProductivityimprovedbyEditField_4
            app.ProductivityimprovedbyEditField_4 = uieditfield(app.Panel_6, 'numeric');
            app.ProductivityimprovedbyEditField_4.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_4ValueChanged, true);
            app.ProductivityimprovedbyEditField_4.FontSize = 36;
            app.ProductivityimprovedbyEditField_4.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_4.FontColor = [1 0.4118 0.1608];
            app.ProductivityimprovedbyEditField_4.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_4.Position = [208 8 167 49];

            % Create ProductivityimprovedbyEditFieldLabel_5
            app.ProductivityimprovedbyEditFieldLabel_5 = uilabel(app.Panel_6);
            app.ProductivityimprovedbyEditFieldLabel_5.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_5.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_5.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_5.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_5.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_5.Position = [372 13 44 47];
            app.ProductivityimprovedbyEditFieldLabel_5.Text = 'to:';

            % Create ProductivityimprovedbyEditField_5
            app.ProductivityimprovedbyEditField_5 = uieditfield(app.Panel_6, 'numeric');
            app.ProductivityimprovedbyEditField_5.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_5ValueChanged, true);
            app.ProductivityimprovedbyEditField_5.FontSize = 36;
            app.ProductivityimprovedbyEditField_5.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_5.FontColor = [0 1 0];
            app.ProductivityimprovedbyEditField_5.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_5.Position = [418 8 167 49];

            % Create mLabel
            app.mLabel = uilabel(app.Panel_6);
            app.mLabel.FontName = 'Franklin Gothic Medium';
            app.mLabel.FontSize = 20;
            app.mLabel.FontAngle = 'italic';
            app.mLabel.FontColor = [0.0588 1 1];
            app.mLabel.Position = [592 19 33 26];
            app.mLabel.Text = 'μm';

            % Create SurfroughnessimprovedbyEditFieldLabel
            app.SurfroughnessimprovedbyEditFieldLabel = uilabel(app.Panel_6);
            app.SurfroughnessimprovedbyEditFieldLabel.BackgroundColor = [0.1294 0.3098 0.3098];
            app.SurfroughnessimprovedbyEditFieldLabel.HorizontalAlignment = 'right';
            app.SurfroughnessimprovedbyEditFieldLabel.FontName = 'Franklin Gothic Medium';
            app.SurfroughnessimprovedbyEditFieldLabel.FontSize = 36;
            app.SurfroughnessimprovedbyEditFieldLabel.FontAngle = 'italic';
            app.SurfroughnessimprovedbyEditFieldLabel.FontColor = [0.0588 1 1];
            app.SurfroughnessimprovedbyEditFieldLabel.Position = [7 71 458 47];
            app.SurfroughnessimprovedbyEditFieldLabel.Text = 'Surf. roughness improved by:';

            % Create SurfroughnessimprovedbyEditField
            app.SurfroughnessimprovedbyEditField = uieditfield(app.Panel_6, 'numeric');
            app.SurfroughnessimprovedbyEditField.ValueChangedFcn = createCallbackFcn(app, @SurfroughnessimprovedbyEditFieldValueChanged, true);
            app.SurfroughnessimprovedbyEditField.FontSize = 50;
            app.SurfroughnessimprovedbyEditField.FontAngle = 'italic';
            app.SurfroughnessimprovedbyEditField.FontColor = [0.0588 1 1];
            app.SurfroughnessimprovedbyEditField.BackgroundColor = [0.1294 0.3098 0.3098];
            app.SurfroughnessimprovedbyEditField.Position = [468 61 165 64];

            % Create Panel_7
            app.Panel_7 = uipanel(app.Panel_4);
            app.Panel_7.BackgroundColor = [0.1294 0.3098 0.3098];
            app.Panel_7.SizeChangedFcn = createCallbackFcn(app, @Panel_7SizeChanged, true);
            app.Panel_7.Position = [48 8 692 136];

            % Create ProductivityimprovedbyEditFieldLabel_6
            app.ProductivityimprovedbyEditFieldLabel_6 = uilabel(app.Panel_7);
            app.ProductivityimprovedbyEditFieldLabel_6.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_6.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_6.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_6.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_6.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_6.Position = [47 14 156 47];
            app.ProductivityimprovedbyEditFieldLabel_6.Text = 'Ra   from:';

            % Create ProductivityimprovedbyEditField_6
            app.ProductivityimprovedbyEditField_6 = uieditfield(app.Panel_7, 'numeric');
            app.ProductivityimprovedbyEditField_6.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_6ValueChanged, true);
            app.ProductivityimprovedbyEditField_6.FontSize = 36;
            app.ProductivityimprovedbyEditField_6.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_6.FontColor = [1 0.4118 0.1608];
            app.ProductivityimprovedbyEditField_6.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_6.Position = [208 5 167 58];

            % Create ProductivityimprovedbyEditFieldLabel_7
            app.ProductivityimprovedbyEditFieldLabel_7 = uilabel(app.Panel_7);
            app.ProductivityimprovedbyEditFieldLabel_7.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_7.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_7.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_7.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_7.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_7.Position = [378 13 44 47];
            app.ProductivityimprovedbyEditFieldLabel_7.Text = 'to:';

            % Create ProductivityimprovedbyEditField_7
            app.ProductivityimprovedbyEditField_7 = uieditfield(app.Panel_7, 'numeric');
            app.ProductivityimprovedbyEditField_7.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_7ValueChanged, true);
            app.ProductivityimprovedbyEditField_7.FontSize = 36;
            app.ProductivityimprovedbyEditField_7.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_7.FontColor = [0 1 0];
            app.ProductivityimprovedbyEditField_7.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_7.Position = [424 5 167 58];

            % Create mLabel_2
            app.mLabel_2 = uilabel(app.Panel_7);
            app.mLabel_2.FontName = 'Franklin Gothic Medium';
            app.mLabel_2.FontSize = 20;
            app.mLabel_2.FontAngle = 'italic';
            app.mLabel_2.FontColor = [0.0588 1 1];
            app.mLabel_2.Position = [601 21 33 26];
            app.mLabel_2.Text = 'μm';

            % Create ProductivityimprovedbyEditFieldLabel_8
            app.ProductivityimprovedbyEditFieldLabel_8 = uilabel(app.Panel_7);
            app.ProductivityimprovedbyEditFieldLabel_8.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_8.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_8.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_8.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_8.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_8.Position = [13 73 190 47];
            app.ProductivityimprovedbyEditFieldLabel_8.Text = 'MRR   from:';

            % Create ProductivityimprovedbyEditField_8
            app.ProductivityimprovedbyEditField_8 = uieditfield(app.Panel_7, 'numeric');
            app.ProductivityimprovedbyEditField_8.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_8ValueChanged, true);
            app.ProductivityimprovedbyEditField_8.FontSize = 36;
            app.ProductivityimprovedbyEditField_8.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_8.FontColor = [1 0.4118 0.1608];
            app.ProductivityimprovedbyEditField_8.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_8.Position = [208 64 167 58];

            % Create ProductivityimprovedbyEditFieldLabel_9
            app.ProductivityimprovedbyEditFieldLabel_9 = uilabel(app.Panel_7);
            app.ProductivityimprovedbyEditFieldLabel_9.HorizontalAlignment = 'right';
            app.ProductivityimprovedbyEditFieldLabel_9.FontName = 'Franklin Gothic Medium';
            app.ProductivityimprovedbyEditFieldLabel_9.FontSize = 36;
            app.ProductivityimprovedbyEditFieldLabel_9.FontAngle = 'italic';
            app.ProductivityimprovedbyEditFieldLabel_9.FontColor = [0.0588 1 1];
            app.ProductivityimprovedbyEditFieldLabel_9.Position = [378 72 44 47];
            app.ProductivityimprovedbyEditFieldLabel_9.Text = 'to:';

            % Create ProductivityimprovedbyEditField_9
            app.ProductivityimprovedbyEditField_9 = uieditfield(app.Panel_7, 'numeric');
            app.ProductivityimprovedbyEditField_9.ValueChangedFcn = createCallbackFcn(app, @ProductivityimprovedbyEditField_9ValueChanged, true);
            app.ProductivityimprovedbyEditField_9.FontSize = 36;
            app.ProductivityimprovedbyEditField_9.FontAngle = 'italic';
            app.ProductivityimprovedbyEditField_9.FontColor = [0 1 0];
            app.ProductivityimprovedbyEditField_9.BackgroundColor = [0.1294 0.3098 0.3098];
            app.ProductivityimprovedbyEditField_9.Position = [424 64 167 58];

            % Create cmminLabel_2
            app.cmminLabel_2 = uilabel(app.Panel_7);
            app.cmminLabel_2.FontName = 'Franklin Gothic Medium';
            app.cmminLabel_2.FontSize = 20;
            app.cmminLabel_2.FontAngle = 'italic';
            app.cmminLabel_2.FontColor = [0.0588 1 1];
            app.cmminLabel_2.Position = [601 80 81 26];
            app.cmminLabel_2.Text = 'cm³/min';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Millplus

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
