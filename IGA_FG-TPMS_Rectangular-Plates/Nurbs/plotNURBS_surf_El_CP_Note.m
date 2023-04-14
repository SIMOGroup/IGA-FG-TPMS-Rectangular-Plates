function plotNURBS_surf_El_CP_Note(p,q,U,V,CP)
% plots the surface, elements and control points

[X,Y,Z] = create_surf(p,q,U,V,CP);

% geometry
figure('Color','w','Units','normalized','Outerposition',[0 0 0.33 0.91],'visible','on'); axis equal
% figure('Color','w','Units','normalized','Outerposition',[0 0 0.51 0.95],'visible','on'); axis equal
surf(X,Y,Z,'FaceColor','none','EdgeColor','none','LineStyle','--', 'LineWidth',1.2);
% xlabel('x'); ylabel('y'); zlabel('z');
hold on;

%% Rectangular Plate
view(2)
xlim([0, 1]); ylim([0, 1]);
set(gca, 'YTick', [1.48], 'XTick', [0.98], 'TickLength', [0.001 0.001], 'TickLabelInterpreter', 'latex')
pgca = gca;
pgca.YTickLabel = ['$b$'];
pgca.XTickLabel = ['$a$'];
pgca.XRuler.TickLabelGapOffset = 0;
pgca.YRuler.TickLabelGapOffset = 10;
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
set(gca, 'box', 'off', 'GridLineStyle', 'none')
axp = get(gca,'Position');
annotation('arrow', [axp(1) axp(1)+axp(3)+0.07], [axp(2) axp(2)]+0.026, 'linewidth', 1);
annotation('arrow', [axp(1) axp(1)]+0.03, [axp(2) axp(2)+axp(4)+0.07], 'linewidth', 1);
annotation('textbox', [0.92 0.02 0.1 0.1],'EdgeColor','none','String','\textbf{$x$}','FitBoxToText','on','fontsize',18,'Interpreter','latex');
annotation('textbox', [0.07 0.88 0.1 0.1],'EdgeColor','none','String','\textbf{$y$}','FitBoxToText','on','fontsize',18,'Interpreter','latex');
annotation('textbox', [0.07 0.02 0.1 0.1],'EdgeColor','none','String','\textbf{$O$}','FitBoxToText','on','fontsize',18,'Interpreter','latex');
set(gca,'fontsize', 18)

%% Circular Plate
% view(2)
% xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);
% set(gca, 'YTick', [-1.1, 1.1], 'XTick', [-1.13, 1.1], 'TickLength', [0.001 0.001], 'TickLabelInterpreter', 'latex')
% pgca = gca;
% pgca.YTickLabel = ["$-R$", "$R$"];
% pgca.XTickLabel = ["$-R$", "$R$"];
% pgca.XRuler.TickLabelGapOffset = 6;
% pgca.YRuler.TickLabelGapOffset = 10;
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% set(gca, 'box', 'off', 'GridLineStyle', 'none')
% axp = get(gca,'Position');
% annotation('arrow', [axp(1) axp(1)+axp(3)+0.03], [axp(2)+axp(4)/2 axp(2)+axp(4)/2], 'linewidth', 1);
% annotation('arrow', [axp(1)+axp(3)/2 axp(1)+axp(3)/2], [axp(2) axp(2)+axp(4)+0.03], 'linewidth', 1);
% annotation('textbox', [0.88 0.41 0.1 0.1],'EdgeColor','none','String','\textbf{$x$}','FitBoxToText','on','fontsize',16,'Interpreter','latex');
% annotation('textbox', [0.46 0.85 0.1 0.1],'EdgeColor','none','String','\textbf{$y$}','FitBoxToText','on','fontsize',16,'Interpreter','latex');
% annotation('textbox', [0.46 0.41 0.1 0.1],'EdgeColor','none','String','\textbf{$O$}','FitBoxToText','on','fontsize',16,'Interpreter','latex');
% set(gca,'fontsize', 14)

% element edges
create_el_edges(p,q,U,V,CP)

axis equal;

hold off;

end