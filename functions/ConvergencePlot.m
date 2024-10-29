

function ConvergencePlot(sum_u,absdudt,time,nn)
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    % nn=length(sum_u);

    % convergence block 
    converge=gcf;
    subplot(3,1,1)
    plot(3:nn,sum_u(3:end),'LineWidth',4);title('$ \sum{u}_{ij}$');
    ylabel('$\frac{m}{s}$')

    subplot(3,1,2)
    plot(3:nn,absdudt(3:end),'LineWidth',4);
    title('$ \sum{ \frac{|\Delta u|}{\Delta t}  }$');
    ylabel('$\frac{m}{s}$')
  
    subplot(3,1,3)
    plot(3:nn,time(3:end),'LineWidth',4);title('$ \Delta$t CFL');
    xlabel('iteration');ylabel('s')

    sgtitle('converging variables and $\Delta$t plots ')
    set(findall(converge,'-property','FontSize'),'FontSize',15)
    drawnow
     
end
