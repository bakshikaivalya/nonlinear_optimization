clc

figure(1)
generations = 50;
subplot(2,1,1);
xlim = [1,generations+1];
plot(1:1:generations+1,avgd,'*');
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

% subplot(2,1,1)
% legend('Elite count = 2','Elite count = 10')
% set(legend,'FontSize',12);
% subplot(2,1,2)
% legend('Elite count = 2', '', '', 'Elite count 2 Average', 'Elite count = 10', '', '', 'Elite count 10 Average Value')
% set(legend,'FontSize',12);
% saveas(gca,'elite10avg.fig')


% subplot(2,1,1)
% legend('Crossover fraction = 0.75','Crossover fraction = 0.99')
% set(legend,'FontSize',12);
% subplot(2,1,2)
% legend('Crossover fraction = 0.75', '', '', 'Crossover fraction 0.75 Average', 'Crossover fraction = 0.99', '',...
% '', 'Crossover fraction 0.99 Average Value')
% set(legend,'FontSize',12);
% saveas(gca,'cross0_99avg.fig')


% subplot(2,1,1)
% legend('Mutation rate = 0.01','Mutation rate = 0.1')
% set(legend,'FontSize',12);
% subplot(2,1,2)
% legend('Mutation rate = 0.01', '', '', 'Mutation rate 0.01 Average', 'Mutation rate = 0.1', '',...
% '', 'Mutation rate 0.1 Average Value')
% set(legend,'FontSize',12);
% saveas(gca,'mut0_1avg.fig')


% subplot(2,1,1)
% legend('Population size = 50','Population size 100')
% set(legend,'FontSize',12);
% subplot(2,1,2)
% legend('Population size = 50', '', '', 'Population size 50 Average', 'Population size = 100', '',...
% '', 'Population size 100 Average Value')
% set(legend,'FontSize',12);
% saveas(gca,'popsize100avg.fig')