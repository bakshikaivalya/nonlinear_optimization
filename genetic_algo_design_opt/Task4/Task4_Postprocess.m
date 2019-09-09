figure(1)
generations = 50;
subplot(2,1,1);
xlim = [1,generations+1];
plot(1:1:generations+1,avgdist,'*');
xlabel('Generations','Fontsize',12)
ylabel('Average Distance','Fontsize',12)
grid on
hold on
title('Average distance averaged over 10 optimized samples','Fontsize',12)
subplot(2,1,2);
xlim = [1,generations+1];
% plot(1:1:generations+1,avgf);
% errorbar(1:1:generations+1,avgf,minf,maxf,'-');
shadedErrorBar(1:1:generations+1, avgf, [maxf; minf],'-',1);
xlabel('Generations','Fontsize',12)
ylabel('f','Fontsize',12)
grid on
hold on
title('Average f averaged over 10 optimized samples','Fontsize',12)