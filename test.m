clear
figure(1);clf(1);hold on
for drive = [1,.6,.3]
    for k = 1:200
        u(k) = drive;
        if k > 160
            u(k) = 0;
        end    
    end
    a(1) = 0;
    for k = 2:200
        a(k) = rk4(a(k-1),u(k));
    end
    plot(a,'Linewidth',3)
end
legend({'u = 1','u = .6','u = .3'})
set(gca,'Xtick',[0,40,80,120,160,200],'XtickLabels',[0,40,80,120,160,200]*0.005)
set(gca,'Ytick',[0,.2,.4,.6,.8,1])
xlabel('Time (s)');
ylabel('Active State');