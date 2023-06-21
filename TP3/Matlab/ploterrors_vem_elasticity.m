%%%%%%%%%%%%% L^2 errors %%%%%%%%%%%%%

PLOT_PEPSILON_L2_ERRORS = true;
% PLOT_PEPSILON_L2_ERRORS = false;


PLOT_SIGMA_L2_ERRORS = true;
% PLOT_SIGMA_L2_ERRORS = false;


% PLOT_l2_ON_EDGES_ERRORS = true;
PLOT_l2_ON_EDGES_ERRORS = false;
% graphs

DRAW_MARKER = true;
%DRAW_MARKER = false;

LINE_WIDTH = 2;

% choose HMAX o HMEAN
%
hv = [errors.hmean];
%hv = [errors.hmax];

slope = [
    'for i=1:length(xv)-1;'...
    ' x = sqrt(xv(i)*xv(i+1));'...
    ' y = sqrt(yv(i)*yv(i+1));'...
    ' m = (log(yv(i+1))-log(yv(i)))/(log(xv(i+1))-log(xv(i)));'...
    ' text(x,y,sprintf(''%.2f'',(m)),''HorizontalAlignment'',''center'',''BackgroundColor'',''y'',''FontSize'',6);'...
    'end'
    ];

%
%
%%%%%%%%%%%%%%%%%%% errori L2 %%%%%%%%%%%%%%%%%%%%%%%
%
figure(10)
%
clear mylegendL2;
ilegend = 0;
%
if PLOT_PEPSILON_L2_ERRORS
    %
    % P0 L^2
    %
    xv = hv;
    yv = [errors.errL2PiEpsilonsq]./[errors.normL2Epsilonsq];
    %
    if DRAW_MARKER
        loglog(xv,yv,'r-d','LineWidth',LINE_WIDTH)
    else
        loglog(xv,yv,'r-','LineWidth',LINE_WIDTH)
    end
    %
    eval(slope);
    %
    ilegend = ilegend+1;
    mylegendL2{ilegend} = '||\Pi\epsilon(v_h) - \epsilon(v)||_{0,\Omega}';
    %
    if not(ishold)
        hold on
    end
end
%
grid on
%
title(['k=' num2str(errors(1).k) ' - L^2 errors']);
%
legend(mylegendL2,'Location','NorthWest')
%
drawnow

%%%%%%%%%%%%%%%%%%% errori L2 %%%%%%%%%%%%%%%%%%%%%%%
%
figure(11)
%
clear mylegendL2;
ilegend = 0;
%
if PLOT_SIGMA_L2_ERRORS
    %
    % P0 L^2
    %
    xv = hv;
    yv = [errors.errL2Sigma_hsq]./[errors.normL2Sigmasq];
    %
    if DRAW_MARKER
        loglog(xv,yv,'r-d','LineWidth',LINE_WIDTH)
    else
        loglog(xv,yv,'r-','LineWidth',LINE_WIDTH)
    end
    %
    eval(slope);
    %
    ilegend = ilegend+1;
    mylegendL2{ilegend} = '||C\Pi\epsilon(v_h) - \sigma||_{0,\Omega}';
    %
    if not(ishold)
        hold on
    end
end
%
grid on
%
title(['k=' num2str(errors(1).k) ' - L^2 errors']);
%
legend(mylegendL2,'Location','NorthWest')
%
drawnow

%%%%%%%%%%%%%%%%%%% errori l2 %%%%%%%%%%%%%%%%%%%%%%%
%
if PLOT_l2_ON_EDGES_ERRORS
%
figure(12)
clear mylegendl2;
ilegend = 0;
    %
    % P0 L^2
    %
    xv = hv;
    yv = [errors.errl2]./[errors.norml2ue];
    %
    if DRAW_MARKER
        loglog(xv,yv,'r-d','LineWidth',LINE_WIDTH)
    else
        loglog(xv,yv,'r-','LineWidth',LINE_WIDTH)
    end
    %
    eval(slope);
    %
    ilegend = ilegend+1;
    mylegendl2{ilegend} = '||u_h - u_e||_{l^2}';
    %
    if not(ishold)
        hold on
    end
end


% %
% grid on
% %
% title(['k=' num2str(errors(1).k) ' - l^2 errors']);
% %
% legend(mylegendl2,'Location','NorthWest')
% %
% drawnow


