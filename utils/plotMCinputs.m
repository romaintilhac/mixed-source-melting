% Plot randomized inputs for Monte Carlo simulations

if pass==0;fig_distrib=figure('Name','Monte Carlo randomized input ranges'); end

figure(fig_distrib)
ncols=round(MonteCarlo_save.Nrand_inputs/2);

for i=1:MonteCarlo_save.Nrand_inputs
    subplot(2,ncols,i)
    histogram(MonteCarlo_save.(str_pass)(:,i));
    xmin=min(MonteCarlo_save.(strcat('pass0'))(:,i)); %force to pass0 xlim
    xmax=max(MonteCarlo_save.(strcat('pass0'))(:,i));  %force to pass0 xlim
    if xmin~=xmax; xlim([xmin, xmax]); end
    title(MonteCarlo_save.rand_inputs(i),'Interpreter','none')
end