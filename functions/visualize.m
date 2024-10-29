function visualize(n,ptspos,u,v,T,p,e,rho,timenow,iteration)
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');


    tcl = tiledlayout(2,3);

    titlepos = [0.5e-5,8.7e-6,0];
    nexttile
      scatteredpcolor(ptspos(:,1),ptspos(:,2),u);
    axis tight equal
    % pcolor(X,Y,u(:,:,n));shading interp;
    axis equal tight;
    xlabel('x');ylabel('y');title('x-velocity ($u$)','Position',titlepos);
    colorbarEden('u $\left[\frac{m}{s}\right]$');

    nexttile
      scatteredpcolor(ptspos(:,1),ptspos(:,2),v);
    % pcolor(X,Y,v(:,:,n));shading interp;
    axis equal tight;
    xlabel('x');ylabel('y');title('y-velocity ($v$)','Position',titlepos);
    colorbarEden('v $\left[\frac{m}{s}\right]$');

    nexttile
      scatteredpcolor(ptspos(:,1),ptspos(:,2),p);
    % pcolor(X,Y,p(:,:,n));shading interp;
    axis equal tight;
    xlabel('x');ylabel('y');title('Pressure','Position',titlepos);
    colorbarEden('p $\left[Pa\right]$');

    nexttile
      scatteredpcolor(ptspos(:,1),ptspos(:,2),T);
    % pcolor(X,Y,T(:,:,n));shading interp;
    axis equal tight;
    xlabel('x');ylabel('y');title('Temperature','Position',titlepos);
    colorbarEden('T $\left[K\right]$');

    nexttile();
      scatteredpcolor(ptspos(:,1),ptspos(:,2),rho);
    % pcolor(X,Y,rho(:,:,n));shading interp;
    axis equal tight;
    colorbarEden('$\rho\ \left[\frac{kg}{m^3}\right]$');


    shading interp;axis equal tight;
    xlabel('x');ylabel('y');
    title('Density\ ($\rho$) ','Position',titlepos);

    nexttile
    scatteredpcolor(ptspos(:,1),ptspos(:,2),e);axis equal tight;
 
    
    xlabel('x');ylabel('y');title('\ \ \ \  Internal\ Energy\ ($e$)','Position',titlepos);
    colorbarEden('$e\ \left[ \frac{J}{kg} \right]$');
    title(tcl,{sprintf(' t = %0.3g  ( %d / %d )', timenow , n , iteration ), '  '},...
        'FontSize',19,'Interpreter','latex');
    drawnow


    % 
    % 
    % tcl = tiledlayout(2,3);
    % 
    % titlepos = [0.5e-5,8.7e-6,0];
    % nexttile
    % scatter(ptspos(:,1),ptspos(:,2),[],u,"filled");
    % axis tight equal
    % % pcolor(X,Y,u(:,:,n));shading interp;
    % axis equal tight;
    % xlabel('x');ylabel('y');title('x-velocity ($u$)','Position',titlepos);
    % colorbarEden('u $\left[\frac{m}{s}\right]$');
    % 
    % nexttile
    % scatter(ptspos(:,1),ptspos(:,2),[],v,"filled");
    % % pcolor(X,Y,v(:,:,n));shading interp;
    % axis equal tight;
    % xlabel('x');ylabel('y');title('y-velocity ($v$)','Position',titlepos);
    % colorbarEden('v $\left[\frac{m}{s}\right]$');
    % 
    % nexttile
    % scatter(ptspos(:,1),ptspos(:,2),[],p,"filled");
    % % pcolor(X,Y,p(:,:,n));shading interp;
    % axis equal tight;
    % xlabel('x');ylabel('y');title('Pressure','Position',titlepos);
    % colorbarEden('p $\left[Pa\right]$');
    % 
    % nexttile
    % scatter(ptspos(:,1),ptspos(:,2),[],T,"filled");
    % % pcolor(X,Y,T(:,:,n));shading interp;
    % axis equal tight;
    % xlabel('x');ylabel('y');title('Temperature','Position',titlepos);
    % colorbarEden('T $\left[K\right]$');
    % 
    % nexttile();
    % scatter(ptspos(:,1),ptspos(:,2),[],rho,"filled");
    % % pcolor(X,Y,rho(:,:,n));shading interp;
    % axis equal tight;
    % colorbarEden('$\rho\ \left[\frac{kg}{m^3}\right]$');
    % 
    % 
    % shading interp;axis equal tight;
    % xlabel('x');ylabel('y');
    % title('Density\ ($\rho$) ','Position',titlepos);
    % 
    % nexttile
    % scatter(ptspos(:,1),ptspos(:,2),[],e,"filled");axis equal tight;
    % % pcolor(X,Y,e(:,:,n));shading interp;axis equal tight;
    % xlabel('x');ylabel('y');title('\ \ \ \  Internal\ Energy\ ($e$)','Position',titlepos);
    % colorbarEden('$e\ \left[ \frac{J}{kg} \right]$');
    % title(tcl,{sprintf(' t = %0.3g  ( %d / %d )', timenow , n , iteration ), '  '},...
    %     'FontSize',19,'Interpreter','latex');
    % drawnow
function scatteredpcolor(x, y, c)
    triT = delaunay(x, y);
    trisurf(triT, x, y, zeros(size(x)), c, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(0, 90);
end
end