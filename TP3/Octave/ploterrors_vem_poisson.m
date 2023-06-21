%%%%%%%%%%%%% L^2 errors %%%%%%%%%%%%%

PLOT_P0k_L2_ERRORS = true;
%PLOT_P0k_L2_ERRORS = false;

PLOT_PNk_L2_ERRORS = true;
%PLOT_PNk_L2_ERRORS = false;

%%%%%%%%%%%% H^1 errors %%%%%%%%%%%%%

PLOT_P0k_H1_ERRORS = true;
%PLOT_P0k_H1_ERRORS = false;

PLOT_PNk_H1_ERRORS = true;
%PLOT_PN_H1_ERRORS = false;

%%%%%%%%%%%% l^2 errors %%%%%%%%%%%%%%

%PLOT_l2_ON_EDGES_ERRORS = true;
PLOT_l2_ON_EDGES_ERRORS = false;

%%%%%%%%%%%% max errors %%%%%%%%%%%%%%

PLOT_MAX_ON_EDGES_ERRORS = true;
%PLOT_MAX_ON_EDGES_ERRORS = false;

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
%%%%%%%%%%%%%%%%%%% errori H1 %%%%%%%%%%%%%%%%%%%%%%%
%
figure(10)
%
clear mylegendH1;
ilegend = 0;
%
if PLOT_P0k_H1_ERRORS
    %
    % P0k H^1
    %
    xv = hv;
    yv = [errors.errH1P0k]./[errors.normH1ue];
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
    %
    mylegendH1{ilegend} = '|\Pi^0_k u_h - u_e|_{1,\Omega}';
    %
    if not(ishold)
        hold on
    end
    %
end
%
if PLOT_PNk_H1_ERRORS
    %
    % PN H^1
    %
    xv = hv;
    yv = [errors.errH1PNk]./[errors.normH1ue];
    %
    if DRAW_MARKER
        loglog(xv,yv,'b-d','LineWidth',LINE_WIDTH)
    else
        loglog(xv,yv,'b-','LineWidth',LINE_WIDTH)
    end
    %
    eval(slope);
    %
    ilegend = ilegend+1;
    mylegendH1{ilegend} = '|\Pi^\nabla_k u_h - u_e|_{1,\Omega}';
    %
    if not(ishold)
        hold on
    end
    %
end
%
grid on
%
title(['k=' num2str(errors(1).k) ' - H^1 errors']);
%
legend(mylegendH1,'Location','NorthWest')
%
drawnow
%
%%%%%%%%%%%%%%%%%%% errori L2 %%%%%%%%%%%%%%%%%%%%%%%
%
figure(11)
%
clear mylegendL2;
ilegend = 0;
%
if PLOT_P0k_L2_ERRORS
    %
    % P0 L^2
    %
    xv = hv;
    yv = [errors.errL2P0k]./[errors.normL2ue];
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
    mylegendL2{ilegend} = '||\Pi^0_k u_h - u_e||_{0,\Omega}';
    %
    if not(ishold)
        hold on
    end
end
%
if PLOT_PNk_L2_ERRORS
    %
    % PN L^2
    %
    xv = hv;
    yv = [errors.errL2PNk]./[errors.normmaxue];
    %
    if DRAW_MARKER
        loglog(xv,yv,'b-d','LineWidth',LINE_WIDTH)
    else
        loglog(xv,yv,'b-','LineWidth',LINE_WIDTH)
    end
    %
    eval(slope);
    %
    ilegend = ilegend+1;
    mylegendL2{ilegend} = '||\Pi^\nabla_k u_h - u_e||_{0,\Omega}';
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

% if PLOT_l2_ON_EDGES_ERRORS
% %
% figure(14)
% clear mylegendl2;
% ilegend = 0;
%     %
%     % PN L^2
%     %
%     xv = hv;
%     yv = [errors.errl2]./[errors.norml2ue];
%     %
%     if DRAW_MARKER
%         loglog(xv,yv,'b-d','LineWidth',LINE_WIDTH)
%     else
%         loglog(xv,yv,'b-','LineWidth',LINE_WIDTH)
%     end
%     %
%     eval(slope);
%     %
%     ilegend = ilegend+1;
%     mylegendl2{ilegend} = '||u_h - u_e||_{l^2}';
%     %
%     if not(ishold)
%         hold on
%     end
% end
% %
% grid on
% %
% title(['k=' num2str(errors(1).k) ' - l^2 errors']);
% %
% legend(mylegendl2,'Location','NorthWest')
% %
% drawnow