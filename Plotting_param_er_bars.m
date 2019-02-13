cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files');
% load('Comp_matrix.mat');

ColorSet = parula(7);
ColorSet = ColorSet(1:6,:);

for k = 1:6
    for p = 1:4
        a(k,p) = Comp_matrix{k}.a(p);
        erbar(k,p) = Comp_matrix{k}.se(p)*1.96;
        if k==2 && p==1
            a(k,p) = a(k,p)/100
        end
    end
end
figure(1);clf(1); hold on
for k = 1:4
    subplot(2,2,k);hold on;
    x = [1:6];
    if k == 1
        erbar(2,1) = 15;
    end
    for p = 1:6
        h(p) = errorbar(x(p),a(p,k),erbar(p,k),'o','Color',ColorSet(p,:),'linewidth',2);
    end
    xlim([0 7]);
    y_str = sprintf('$a_%g$',k);
    ylabel(y_str,'Interpreter','latex');
    set(gca,'xtick',[1:6],...
        'xticklabel',{'Torque','Force','Active State', 'Neural Drive', 'Energetic','M Power'});
    xtickangle(35);
    pbaspect([1 1 1]);
%     legend(h, {'Torque','Force','Active State', 'Neural Drive', 'Energetic','M Power'});
end
beautifyfig;