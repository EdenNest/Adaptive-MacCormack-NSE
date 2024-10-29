function UserParameters()
%% USER INPUT FOR VISUALIZATION



    % Plot every # snapshot :
         % Insert your desired gap between animation snapshots
         % 0 = plotting only the last snapshot
      Animation_plot_gap = 15;


    % Plot convergence animation (1 = active, 0 = inactive)
        Plot_convergence = 1; 




    % Activate CFL condition  (1 = active, 0 = in-active (constant dt))  
        CFLdt = 1;
    
 
    
    % Mach Number
        Mach = 4;
     
    % 

    M = Mach;
    nx = 75;
    ny = 80;
    L = 1e-5;
    H = 8e-6;
    dt_cnst = 2.35e-11/3; 
  

   % standard sea level
    p0 = 101300;
    rho0 = 1.225; 
    R = 287;
    cp = 1005;
    cv = 718;
    gamma = 1.4;
    Pr = 0.71;
    mu0 = 1.735e-5; 
    S1 = 110.4;
    T0 = 288.15;
    u_inf = M*sqrt(gamma*R*T0);
  

% Coarse Domain

  
        npts = nx*ny; % Total number of initial grid points
        x_pos = linspace(0,L,nx); % Array of x positions
        y_pos = linspace(0,H,ny); % Array of y positions
        dx = x_pos(2)-x_pos(1); % Calculating horizontal step size
        dy = y_pos(2)-y_pos(1); % Calculating vertical step sizes
        % 
        save('parameters')%,'p0','R','rho0','cp','cv','gamma','Pr','mu0','S1','T0'...
        %     ,'u_inf','M','nx','ny','L','H','dt_cnst','dx','dy','npts',...
        %     'x_pos','y_pos','iteration','CFLdt','Adiabatic','SchlierenDensity'...
        %     ,'Plot_convergence','Animation_plot_gap')

end
