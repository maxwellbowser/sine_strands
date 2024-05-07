function analyze_Nyquist_plots

close all
clc

% Change this to the directory containing the modulus file
input_dir = 'C:\Users\panch\Desktop\Important\Gregorio Lab\Sinosoidal Perturbration\MATLAB modulus';

% Insert the modulus file name here!!!
filename = 'Full_Data_Summary-Gerrie-Test.xlsx';
input_excel_file_string = [input_dir filesep filename];
    sl_conditions = [200 225 250 0];
    treatment_conditions = {'Treated', 'Untreated', 'Unknown'};
    pCa_conditions = [9 5.87 5.82];

    % This needs to be changed depending on the file
    % For some reason this cannot be multipl options
    heart_sample = {'MouseWT'};
    

    min_freq = 0;
    max_freq = 150;

    % This is where the data will be written for the fits
    output_excel_file_string = append('output-', filename);
    % This is where the figures for each fit will be written
    output_fit_image = 'output\';
    
    % Read table
    d = readtable(input_excel_file_string);
    vn = d.Properties.VariableNames;
    
    unique_fibers = unique(d.Filename);
    
    % Code
    data_counter = 0;
    
    for fiber_counter = 1:numel(unique_fibers)
    for sl_counter = 1:numel(sl_conditions)
        for treatment_counter = 1:numel(treatment_conditions)
            for pCa_counter = 1:numel(pCa_conditions)
                vi = find( ...
                    (strcmp(d.Filename,unique_fibers{fiber_counter})) & ...
                    (strcmp(d.HeartSample,heart_sample)) & ...
                    (d.SL == sl_conditions(sl_counter)) & ...
                    (strcmp(d.Treatment,treatment_conditions{treatment_counter})) & ...
                    (d.pCa == pCa_conditions(pCa_counter)));
                
                if (~isempty(vi))
                    
                    % Pull off data
                    % Sort on min freq cutoff, if not want lowest
                    %vi2 = vi(find(d.Freq_Hz_(vi) >= min_freq));  
                    % Sort on max freq cutoff, if not want highest
                    vi2 = vi(d.Freq_Hz_(vi) <= max_freq);
                    f = d.Freq_Hz_(vi2);
                    em = 1000 * d.Em_kPa_(vi2);
                    vm = 1000 * d.Vm_kPa_(vi2);

                    [best_p,fit_em,fit_vm, exit_flag] = fit_Nyquist_plot_bcwt_2020_02_07(f,em,vm, ...
                        'figure_number',1);                    

                    
                    % Store data
                    data_counter = data_counter + 1;
                    od.Hashcode{data_counter} = d.HashCode{vi2(1)};
                    od.Filename{data_counter} = d.Filename{vi2(1)};
                    od.Treatment{data_counter} = ...
                        strrep(d.Treatment{vi2(1)},' ','_');
                    od.SL{data_counter} = ...
                        strrep(sprintf('%.1f',d.SL(vi2(1))),'.','_');
                    od.pCa(data_counter) = d.pCa(vi2(1));
                    od.A(data_counter) = best_p(1);
                    od.k(data_counter) = best_p(2);
                    od.B(data_counter) = best_p(3);
                    od.rate_b(data_counter) = best_p(4);
                    od.C(data_counter) = best_p(5);
                    od.rate_c(data_counter) = best_p(6);
                    od.r_squared(data_counter) = ...
                        calculate_r_squared([em ; vm],[fit_em ; fit_vm]);
                    
                    
                    % Make figure
                    figure(2);
                    clf
                    subplot(2,2,1)
                    hold on;
                    plot(em,vm,'ks');
                    plot(fit_em,fit_vm,'bo-');
                    plot([0 1e6],[0 0],'k:');
                    xlim([0 1e6]);
                    ylim([-3e5 1e6]);
                    title(['R^2 = ' num2str(od.r_squared(data_counter))])
                    
                    % Make figure with Em and Vm split vs. frequency
                    %figure(3);
                    %clf
                    subplot(2,2,2)
                    semilogx(f, em,'ks'); hold on;
                    semilogx(f, fit_em, 'bo-');
                    xlim([0.1 250]);
                    title(['Fit Exit Flag = ' num2str(exit_flag)])
                    
                    subplot(2,1,2)
                    semilogx(f, vm,'ks'); hold on;
                    semilogx(f, fit_vm,'bo-');
                    xlim([0.1 250]);
                    title(['Fiber = ' d.Filename{vi2(1)}])
                    
                    
                    pause(1)
                    
                    text(5e5,1e6, ...
                        sprintf('Filename: %s\nSL: %s\nTreatment: %s \nRsq: %s', ...
                        od.Filename{data_counter}, ...
                        od.SL{data_counter}, ...
                        od.Treatment{data_counter}, ...
                        od.r_squared(data_counter) ), ...
                        'Interpreter','none');


                    saveas(gcf, ...
                        sprintf('%s%s_%s_%s.png', ...
                        output_fit_image, ...
                        od.Filename{data_counter}, ...
                        od.SL{data_counter}, ...
                        od.Treatment{data_counter}),...
                        'png' );
                   
                end
            end
        end
    end
    end
    
    % Output
    fn = fieldnames(od);
    for i=1:numel(fn)
        od.(fn{i}) = od.(fn{i})';
    end
    writetable(struct2table(od),output_excel_file_string);
    
