% Eden - Anthony 
% 290C - Spring 2024 Final
% Adaptive Mesh On MacCormack

% RUN THIS CODE
% this code works Standalone.

% Use UserParameters.m if you want to suppress animations of Solution or
% change the iteration-step for which you would like to see it. Or to
% supress the convregence plots or turning off the CFL




clc;clear;close all

global T0 dx dy
UserParameters()
load('parameters.mat')
addpath('functions')


iteration = 1500;
refining_interval = 50;


%% set up grid Book Keeping

        ptspos = zeros(npts,2); % Initializing array of node coordinates
        point = 1; % Index for number of nodes being created
        for yindex = 1:ny
            for xindex = 1:nx
                ptspos(point,:) = [x_pos(xindex),y_pos(yindex)]; % Create node coords.
                point = point + 1;
            end
        end
    
        ptsnei = zeros(npts,5); % Initializing array of node neighbors
        point = 1; % Index for finding the neighbors of this node index
        for xindex = 1:nx
            for yindex = 1:ny
                if mod(point,nx) == 0 % If node is on the right boundary
                    ptsnei(point,2) = 0;
                    ptsnei(point,4) = point - 1;
                elseif mod(point,nx) == 1 % If node is on the left boundary
                    ptsnei(point,2) = point + 1;
                    ptsnei(point,4) = 0;
                else % All other nodes
                    ptsnei(point,2) = point + 1;
                    ptsnei(point,4) = point - 1;
                end
                if point<=nx % If node is on the bottom boundary
                    ptsnei(point,1) = point + nx;
                    ptsnei(point,3) = 0;
                elseif point>nx*(ny-1) % If node is on the top boundary
                    ptsnei(point,1) = 0;
                    ptsnei(point,3) = point - nx;
                else % All other nodes
                    ptsnei(point,1) = point + nx;
                    ptsnei(point,3) = point - nx;
                end
                point = point + 1;
            end
            
        end
        clear point xindex yindex
        % Fill column 5 of these nodes with level value of 0 (coarse grid)
        ptsnei(:,5) = 0;


%% INITIAL CONDITION
    
    % 3D matrices for Primitive u v p T With 3rd dimention as time
 
    
    u = ones(npts,1) .* u_inf; v = zeros(npts,1);
    T = ones(npts,1) .* T0;    p = ones(npts,1) .* p0;
    rho = ones(npts,1).*rho0;  e = cv*T;
    Et = rho.*(e+1/2*(u.^2+v.^2));
    U = prim2cons(u,v,T,rho,cv);
    num_of_refined = [];

   
    MainSol(1,:) = { [u,v,T,p,rho,e,Et],ptsnei,ptspos};%[u,v,T,p,e,rho],
 
    
%% SOLVER STARTER

    % plotting related parameters 
        time = 0;
        sum_u = sum(u,'all'); % convergence variable 1
        absdudt = sum(abs(u)/dt_cnst,'all')/npts;  % convergence variable 2
        dtarray = dt_cnst; % array of every dt calculated in each iteration
      
        plotRefinement(2 , 0 ,npts,ptsnei, ptspos);drawnow

 %% SOLVER LOOP %%%%%%%

    for n=1:iteration

      fprintf('computing iteration %d / %d \n' , n, iteration)



        % calculating time step     - with CFL or not 
      
            if CFLdt == 1; dt = CFL(u,v,T,rho);
            else;          dt = dt_cnst;
            end
            
     

      [u,v,T,p,rho,e,Et,U] = ...  
      MacCormack(u,v,T,p,rho,U,dt ,ptsnei,ptspos);
      
      
      if n > 20 && mod(n,refining_interval)==21
          disp('REFINING')
          percent = 90;
                      lightup = DetectGradient(ptsnei,ptspos,v,u,T,percent,dx,dy);
                      if ~isempty(lightup) 
                           fprintf('Refined %d nodes \n' , length(lightup))
                           
                      
                      [newptspos,newptsnei] = ...
                          refinement_v3(lightup,ptspos,ptsnei,dx,dy,npts);
                   
                        
                      u = intpRefined(newptspos,ptspos,u);
                      v = intpRefined(newptspos,ptspos,v);
                      T = intpRefined(newptspos,ptspos,T);
                      p = intpRefined(newptspos,ptspos,p);
                      rho = intpRefined(newptspos,ptspos,rho);
                      e = intpRefined(newptspos,ptspos,e);
                      Et = intpRefined(newptspos,ptspos,Et);
                      for i=1:4
                      U_new(i,:)=intpRefined(newptspos,ptspos,U(i,:)');
                      end
                      U = U_new;
                      clear U_new;


                      dummy = -length(ptsnei)+length(newptsnei);
                      num_of_refined=cat(2,num_of_refined,dummy);
                    
                      plotRefinement(2 , n ,dummy,newptsnei, newptspos)
                      % exportgraphics(figure(2),'refinement.gif','Append',true);
                      ptsnei = newptsnei ; ptspos = newptspos;


                      else
                       disp('Could not refine')
                      end
      end

          

          MainSol(n+1,:) = { [u,v,T,p,rho,e,Et],ptsnei,ptspos};
            

                
        %%%% Visualization %%%%%%



            % saving the refinement animation        
                if Plot_convergence ~= 0 && mod(n-1,10)==0 
                    figure(3)
                    ConvergencePlot(sum_u,absdudt,dtarray,n)
                end
                if Animation_plot_gap ~= 0 && mod(n-1,Animation_plot_gap)==0 && n > 0
                    
                    fig = figure(1);
                    if n == 1;set(fig ,'WindowState','maximized');end
                    visualize(n,ptspos,u,v,T,p,e,rho,time(end),iteration)
                end
               % Updating convergence variables and time arrays for plotting
                 sum_u = cat(2,sum_u,sum(u(1:6e3)));%/length(ptsnei));
                 AA = MainSol{n};
                 absdudt = cat(2,absdudt,     mean(abs(u(1:6e3))-mean(abs(AA(1:6e3,1)))/dt,'all'));%/length(ptsnei)     );
                 time = cat( 2 , time , time(end)+dt) ;
                 dtarray = cat( 2 , dtarray , dt) ;
      
    end
       %%%%%%% End of calculation loop %%%%%%





%%

% Plotting only the last snapshot if desired
fig = figure();
if Animation_plot_gap == 0
    sol = MainSol{length(MainSol),1};
    ptsposs = MainSol{length(MainSol),3};
 
    set(fig ,'WindowState','maximized')
    visualize(length(MainSol),ptsposs,sol(:,1),...
        sol(:,2),sol(:,3),sol(:,4),sol(:,6),sol(:,5),time(i),1500)

end



% save('convergencestuff','absdudt','sum_u','dtarray','time')
% save('refinedstuff','num_of_refined' )
% save('mainsol','MailSol')
   

%% %%%% FUNCTIONS %%%%%%%%%



 function dt = CFL(u,v,T,rho)
   % global gamma Pr R T0 S1 mu0 dx dy
    dy = 8e-6/79;
    dx=1e-5/74;
    R = 287;
   
    gamma = 1.4;
    Pr = 0.71;
    mu0 = 1.735e-5; 
    S1 = 110.4;
    T0 = 288.15;

   mu = mu0*((T/T0).^(3/2)).*((T0+S1)./(T+S1));
   a = @(T) sqrt(gamma*R*T);

   v_prime = max(4/3.*mu.*(gamma*mu/Pr)./rho);
   dtCFL = [abs(u)./dx+abs(v)./dy+a(T).*sqrt(1/dx^2+1/dy^2)...
       +2*ones(size(u)).*v_prime.*(1/dx^2+1/dy^2)].^(-1);

   dt = min(0.2*dtCFL,[],"all");

 end


