%% Sort Data from less to more nodes
[n_sorted, I] = sort(n);
RAM_sorted = vertcat(RAM(I).bytes);
sRAM_sorted = vertcat(sRAM(I).bytes);
t_sorted = t(I);
st_sorted = st(I);

%% Plots
figure (2)
plot(n_sorted, RAM_sorted)
hold on
plot(n_sorted, sRAM_sorted)
xlabel('n'); ylabel('storage [MB]');
legend('full', 'sparce', 'Location','northwest');
set(gcf, 'Position',[10,150,620,550]);

figure (3)
plot(n_sorted, t_sorted)
hold on
plot(n_sorted, st_sorted)
xlabel('x'); ylabel('time [sec]');
legend('full', 'sparce', 'Location','northwest');
set(gcf, 'Position',[640,150,620,550]);