function [best_p,fit_em,fit_vm,exit_flag] = fit_Nyquist_plot_bcwt_2020_02_07(varargin)

% Handle inputs
params = inputParser;
addRequired(params,'freq');
addRequired(params,'elastic_mod');
addRequired(params,'viscous_mod');
addOptional(params,'scaling_vector',[1e6 0.1 1e6 1 1e6 10]);
%addOptional(params,'scaling_vector',[1e5 0.1 1e6 1 1e6 10]);
%addOptional(params,'initial_p',[1 1 1 1 1 1]);
addOptional(params,'initial_p',[0.1 1.5 0.2 3 0.5 1.2]);
addOptional(params,'lower_bounds',zeros(1,6));
%addOptional(params,'lower_bounds',[10 0 100 0.1 100 1]);
addOptional(params,'upper_bounds',inf*ones(1,6));
%addOptional(params,'upper_bounds',[inf 1 inf 10 inf 100]);
addOptional(params,'figure_number',0);
parse(params,varargin{:})
params = params.Results;


% Code
expt_data.freq = params.freq;
expt_data.elastic_mod = params.elastic_mod;
expt_data.viscous_mod = params.viscous_mod;
    
[best_p, feval, exit_flag] = fminsearchbnd( ...
    @return_Nyquist_fit, ...
    params.initial_p, ...
    params.lower_bounds,params.upper_bounds, ...
    optimset('MaxFunEvals',15000), ...
    expt_data, ...
    params.scaling_vector, ...
    0);

for in_exit=1:10
    if  exit_flag == 0

    [best_p, feval, exit_flag] = fminsearchbnd( ...
    @return_Nyquist_fit, ...
    best_p, ...
    params.lower_bounds,params.upper_bounds, ...
    optimset('MaxFunEvals',15000), ...
    expt_data, ...
    params.scaling_vector, ...
    0);
    end
end

disp(['feval = ' num2str(feval)])
disp(['exit flag = ' num2str(exit_flag)])

Nyquist_values = return_Nyquist_values(best_p, expt_data.freq, params.scaling_vector);
fit_em = real(Nyquist_values);
fit_vm = imag(Nyquist_values);

plot_fit(params.figure_number, expt_data, [], best_p, params.scaling_vector);

% Store for output
best_p = best_p .* params.scaling_vector;

end

function Nyquist_values = return_Nyquist_values(p, freq, scaling_vector)

    x_data = 2*pi*1i*freq;
    
    p = p.*scaling_vector;

    Nyquist_values = ...
        p(1) * (x_data .^ p(2)) - ...
        p(3) * (x_data ./ (2*pi*p(4) + x_data)) + ...
        p(5) * (x_data ./ (2*pi*p(6) + x_data));
end

function error_value = return_Nyquist_fit(p, expt_data, scaling_vector, figure_display)
    
    % Calculate Nyquist_values
    Nyquist_values = return_Nyquist_values(p,expt_data.freq, scaling_vector);
    
    % Split into elastic and viscous components
    fit_data.elastic_mod = real(Nyquist_values);
    fit_data.viscous_mod = imag(Nyquist_values);
    
    % Calculate sum of squares
    error_value = ...
        sum((fit_data.elastic_mod - expt_data.elastic_mod).^2) + ...
        sum((fit_data.viscous_mod - expt_data.viscous_mod).^2);

    % Display if required
    if (figure_display)
        plot_fit(figure_display, expt_data, fit_data);
    end
       
end

function plot_fit(figure_display, expt_data, fit_data, p, scaling_vector)

    figure(figure_display);
    cla;
    hold on;
    plot(expt_data.elastic_mod, expt_data.viscous_mod,'k^');
    
    if (isempty(fit_data))
         % Calculate Nyquist_values
        Nyquist_values = return_Nyquist_values(p,expt_data.freq, scaling_vector);
        % Split into elastic and viscous components
        fit_data.elastic_mod = real(Nyquist_values);
        fit_data.viscous_mod = imag(Nyquist_values);
    end
        
    plot(fit_data.elastic_mod, fit_data.viscous_mod,'bo-');
%     for i=1:numel(fit_data.elastic_mod)
%         plot([fit_data.elastic_mod(i) expt_data.elastic_mod(i)], ...
%             [fit_data.viscous_mod(i) expt_data.viscous_mod(i)],'r-');
%     end
end
    
    


