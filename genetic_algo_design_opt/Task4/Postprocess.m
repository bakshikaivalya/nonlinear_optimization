figure(1)
subplot(2,1,1);
xlim = [1,generations+1];
plot(1:1:generations+1,avgd,'*');
xlabel('Generations')
ylabel('Average Distance')
grid on
hold on
title('Average distance averaged over 10 optimized samples')
subplot(2,1,2);
xlim = [1,generations+1];
% plot(1:1:generations+1,avgf);
% errorbar(1:1:generations+1,avgf,minf,maxf,'-');
shadedErrorBar(1:1:generations+1, avgf, [maxf; minf],'-',1);
xlabel('Generations')
ylabel('f')
grid on
hold on
title('Average f averaged over 10 optimized samples')