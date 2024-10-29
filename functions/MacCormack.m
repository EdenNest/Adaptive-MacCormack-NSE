function [u,v,T,p,rho,e,Et,U]=...
    MacCormack(u,v,T,p,rho,U,dt ,ptsnei,ptspos)

    
    addpath('functions')

%%
    global dx dy 
    global p0 T0 u_inf 
    global R Pr cp cv S1 mu0 gamma


   

    % Adiabatic Choice   ( 0 = constant wall temp | 1 = adiabatic wall)
        Adiabatic = 0;  

    % Activate CFL condition  (1 = active, 0 = in-active (constant dt))  
        CFLdt = 0;
    
    % Iteration 
        % iteration = 50 ; 
    
    % Mach Number
        Mach = 4;
     
    % 

    M = Mach;
    nx = 75;
    ny = 80;
    L = 1e-5;
    H = 8e-6;
    dt_cnst = 2.35e-11; 
    
    % Adiabatic doesn't converge with CFL. Enforce constant dt
    if Adiabatic==1
          CFLdt = 0;
    end
  

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

%%

 % global dx dy mu0 cp Pr T0 S1  Adiabatic ptsnei


    % U_n = prim2cons(u,v,T,rho,cv);


       %%%%%% PREDICTOR %%%%%%%%

            % disp('FD biasing ')
            % 

            % disp('Imposing BC')
            [u,v,T,p,rho] = imposeBC(u,v,T,p,rho);

            [E,F] = GenerateBiasedEF( 'fwd' , u, v, T ,p, rho);

            % disp(' main predictor time march')
            U_bar = predictor(U,E,F,dt) ;

            % disp('calculating primitives from U bar')
            [u_bar,v_bar,T_bar,p_bar,rho_bar,~,~] = ...
                cons2prim(U_bar,R,cv);
       
           
            % disp('Imposing Boundary Conditions on Primitives (also updates rho)')
            [u_bar,v_bar,T_bar,p_bar,rho_bar] = ...
                imposeBC(u_bar,v_bar,T_bar,p_bar,rho_bar);

            % disp('Update U bar')
            U_bar = prim2cons(u_bar,v_bar,T_bar,rho_bar,cv);



        % disp('%%%%% CORRECTOR %%%%%%')
% 
            % disp('FD biasing')
                [E,F] = GenerateBiasedEF('bwd',u_bar,v_bar,T_bar,p_bar,rho_bar);

            % disp('main corrector time march')
                U = corrector(U_bar,U,E,F,dt);

            % disp('Calculating primitives from U')
                [u,v,T,p,rho,~,~] = cons2prim(U,R,cv);

            % disp('Imposing BC')
                [u,v,T,p,rho] = imposeBC(u,v,T,p,rho);

            % Update U
                U = prim2cons(u,v,T,rho,cv);

            % disp('getting the final e')
                Et=squeeze(U(4,:))';
                e=(Et./rho)-1/2*(u.^2+v.^2);




  %%% inner functions
  function [txx,tyy,txy,tyx,qx_dot,qy_dot]=InnerDerivatives...
    (xdir,ydir,u,v,T)

    % Calculates stresses and fluxes using the given xdir and ydir as the
    % desired direction for x derivative and y derivative. (Doesn't do
    % biasing. it just simply recieves the directions)

    % global dx dy mu0 cp Pr T0 S1  Adiabatic ptsnei

    % precalculates derivatives according to the input direction

    switch xdir
        case 'bwd'
            ux = ddx_bwd(ptsnei,ptspos,u,dx);
            vx = ddx_bwd(ptsnei,ptspos,v,dx);
            Tx = ddx_bwd(ptsnei,ptspos,T,dx);
        case 'fwd'
            ux = ddx_fwd(ptsnei,ptspos,u,dx);
            vx = ddx_fwd(ptsnei,ptspos,v,dx);
            Tx = ddx_fwd(ptsnei,ptspos,T,dx);
        case 'ctr'
            ux = ddx_central(ptsnei,ptspos,u,dx);
            vx = ddx_central(ptsnei,ptspos,v,dx);
            Tx = ddx_central(ptsnei,ptspos,T,dx);

    end
    % 
    

    % 
    switch ydir
        case 'bwd'
            uy = ddy_bwd(ptsnei,ptspos,u,dy);
            vy = ddy_bwd(ptsnei,ptspos,v,dy);
            Ty = ddy_bwd(ptsnei,ptspos,T,dy);
        case 'fwd'
            uy = ddy_fwd(ptsnei,ptspos,u,dy);
            vy = ddy_fwd(ptsnei,ptspos,v,dy);
            Ty = ddy_fwd(ptsnei,ptspos,T,dy);
        case 'ctr'
            uy = ddy_central(ptsnei,ptspos,u,dy);
            vy = ddy_central(ptsnei,ptspos,v,dy);
            Ty = ddy_central(ptsnei,ptspos,T,dy);
    end
 % disp('checkpoint')
     % Bind = find(ptsnei(:,3)==0);  % Bottom border index
     if Adiabatic == 1
        Ty(ptsnei(:,3)==0) = 0;
     end


    % Get Temperature dependent paramteres

    mu_now = mu0*((T/T0).^(3/2)).*((T0+S1)./(T+S1));
    k_now = cp/Pr.*mu_now;

    % normal and shear stresses 
    txx = 2*mu_now.*(ux - 1/3*(ux + vy));
    tyy = 2*mu_now.*(vy - 1/3*(ux + vy));
    txy = mu_now.*(uy + vx);
    tyx = txy ;

    % heat fluxes 
    qx_dot = - k_now .* Tx;
    qy_dot = - k_now .* Ty;


end






function [E,F] = GenerateBiasedEF(direction, u , v , T , p , rho)

    % Computes E and F using the InnerDerivatives function. It
    % automatically finds the correct biasing to insert into
    % InnerDerivatives function based on "direction" input. "direction"
    % should be the direction of the outer derivative. for example it
    % should be forward for predictor level and backward for the corrector
    % level.
    % example: ydirE = direction for y derivative inside of E.



    % Biasing Logic Block
        ydirE = 'ctr';    
        xdirF = 'ctr';
        switch direction
          case 'fwd'                 % means predictor
              xdirE = 'bwd';
              ydirF = 'bwd';
          case 'bwd'                 % means corrector
              xdirE = 'fwd';   
              ydirF = 'fwd';
        end


    % preallocate E and F
    E = zeros([4,length(u)]);  % this will be 4 x npts
    F = zeros([4,length(u)]);
   

    % Get inner derivatives for F
        [~,tyy,txy,~,~,qy_dot]=...
            InnerDerivatives(xdirF,ydirF,u,v,T);


    % global cv 
    e = cv*T;
    Et = rho.*(e+1/2*(u.^2+v.^2));

    % Compute F
        F(1,:) = rho.*v; 
        F(2,:) = rho.*u.*v - txy; 
        F(3,:) = rho.*v.^2 + p - tyy; 
        F(4,:) = (Et+p).*v - u.*txy -v.*tyy + qy_dot;

    % Get inner derivatives for E
        [txx,~,txy,~,qx_dot,~]=...
            InnerDerivatives(xdirE,ydirE,u,v,T);

    % Compute E
        E(1,:) = rho.*u; 
        E(2,:) = rho.*u.^2 + p - txx;
        E(3,:) = rho.*u.*v - txy; 
        E(4,:) = (Et+p).*u - u.*txx -v.*txy + qx_dot;

end


function [u,v,T,p,rho]=imposeBC(u,v,T,p,rho)

    % Imposes boundary Conditions and updates rho

    % global T0 p0 u_inf R Adiabatic ptsnei
    err = 5e-8;
    Rind = find(ptspos(:,1)>L-err);
    Lind = find(ptspos(:,1)< err);
    Tind = find(ptspos(:,2)>H-err);
    Bind = find(ptspos(:,2)<err);

    % check
%     J = [Rind;Bind;Lind; Tind];
% for j = J;for i = 1: length(J)
% plusplot(J(i),ptsnei,ptspos);drawnow;hold on ;end;end

    % Rind = find(ptsnei(:,2)==0);  % right border index
    % Lind = find(ptsnei(:,4)==0);  % Left border index
    % Tind = find(ptsnei(:,1)==0);  % Top border index
    % Bind = find(ptsnei(:,3)==0);  % Bottom border index
    % 


    % In-flow
    u(Lind) = u_inf;
    v(Lind) = 0;
    T(Lind) = T0;
    p(Lind) = p0;

    % Farfield
    u(Tind) = u_inf;
    v(Tind) = 0;
    T(Tind) = T0;
    p(Tind) = p0;

    % Out-flow (right) 

    left_2nei_idx = slnn(Rind,4,2,ptsnei,ptspos,dx);

    u(Rind) = (4*u(ptsnei(Rind,4)) - u(left_2nei_idx))/3;
    p(Rind) = (4*p(ptsnei(Rind,4)) - p(left_2nei_idx))/3;
    T(Rind) = (4*T(ptsnei(Rind,4)) - T(left_2nei_idx))/3;
    v(Rind) = (4*v(ptsnei(Rind,4)) - v(left_2nei_idx))/3;


    % v(Rind) = (4*v(end-1,:) - v(end-2,:))/3;
    % p(Rind) = (4*p(end-1,:) - p(end-2,:))/3;
    % T(Rind) = (4*T(end-1,:) - T(end-2,:))/3;


    % Wall


    top_2nei_idx = slnn(Bind,1,2,ptsnei,ptspos,dy);
    u(Bind) = 0;
    v(Bind) = 0;
    p(Bind) = (4*p(ptsnei(Bind,1)) - p(top_2nei_idx))/3;
    T(Bind) = T0; 

    % enforce ADIABATIC 
            if Adiabatic == 1
                 T(Bind) = (T(top_2nei_idx)-4*(T(ptsnei(Bind,1)))) / (-3);
                % T(:,1) = (T(:,3)-4*(T(:,2))) / (-3);

            end

    % Leading Edge
    [LEindx ,~] = find(ptspos(:,2)<err & ptspos(:,1)<err);
    u(LEindx) = 0;
    v(LEindx) = 0;
    p(LEindx) = p0;
    T(LEindx) = T0; 


    % Update rho (efficiently)
    rho(Lind) = p(Lind) ./ T(Lind) ./ R;
    rho(Bind) = p(Bind) ./ T(Bind) ./ R;
    rho(Tind) = p(Tind) ./ T(Tind) ./ R;
    rho(Rind) = p(Rind) ./ T(Rind) ./ R;
    % rho(1,:) = p(1,:) ./ T(1,:) ./ R;
    % rho(:,1) = p(:,1) ./ T(:,1) ./ R;
    % rho(:,end) = p(:,end) ./ T(:,end) ./ R;
    % rho(end,:) = p(end,:) ./ T(end,:) ./ R;

   

end




function U_bar = predictor(U,E,F,dt) 
    U_bar=zeros(size(U));
    % global dx dy
    for i=1:4
        Ei=squeeze(E(i,:));
        Fi=squeeze(F(i,:)); 
        Ui = squeeze(U(i,:)) ;
        U_bar(i,:) = Ui - dt  .* ( ddx_fwd(ptsnei,ptspos,Ei,dx) +...
            ddy_fwd(ptsnei,ptspos,Fi,dy ));
    end 

end

function newU = corrector(U_bar,U,E_bar,F_bar,dt)
    % global dx dy
    newU = zeros(size(U));
    for i=1:4
        Ei = squeeze(E_bar(i,:));
        Fi = squeeze(F_bar(i,:));
        Ui = squeeze(U(i,:));
      
        Ubari = squeeze(U_bar(i,:));
      

        newU(i,:) = 1/2 .* (Ui + Ubari ...
        - dt.*( ddx_bwd(ptsnei,ptspos,Ei,dx) + ddy_bwd(ptsnei,ptspos,Fi,dy)));
    end
end

end