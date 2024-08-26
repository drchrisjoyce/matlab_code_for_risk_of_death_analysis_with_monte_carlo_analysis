function MonteCarloAnalysisFunction

% This script uses individual risk predictions to generate 100 estimates of the probablity distribution of the predicted number of
% deaths. Each of these estimates is 'Run'. A cumulative probabilty distribution is calculated for each of estimates of the probablity distribution.

% The number of simulations within each 'Run' is initially set at 1000. The 100 estimates of the probablity distribution are then generated.
% Two tests are applied to these 100 estimates.
%      a) Test 1 ... The highest number of deaths below the lower 95% confidence limit (0.025) must be the same for all the runs.
%      b) Test 2 ... The lowest number of deaths above the upper 95% confidence limit (0.975) must be the same for all the runs.
% If either of the tests fails, the number of simulations is deemed to be too small and the process is repeated,
% with the number of simulations within each Run increased by a factor of 10. If both tests pass, the program outputs the results and terminates.

% Input:
%        A csv file containing individual risk of death predications for each patient, as a fraction of 1.
% Output:
%        A bar graph of the probability distribution, with 95% confidence intervals, is generated and output as a .fig file.
%        Bars representing a number of deaths outside the 95% confidence intervals are red, while those within the 95% confidence intervals are green.

%********************************************************************************************************

% The file with the input data is selected by the user. This data contains the predicted mortality for
% each patient as a fraction of 1. The data must be a single row and with the data elements separated by commas.
% This data is read into the array of double "Predicted_Mortalities".

[infilename, outfilename, Figoutfilename, Predicted_Mortalities] = UserInput;
tic

% For the input predicted mortalities, the probabilty distribution of the number of deaths is calculated using a Monte Carlo technique.
% Multiple simulations (typically 10 000 to 1 000 000) are required to produce a probablity distribution with stable 95% confidence intervals.
% Each set of simulations generating a probablity distibution was termed a 'Run'.
% This process was replicated 100 times, with each replication being a separate 'Run'.

if infilename ~= 0
    tic
    simulations = 1000;
    Runs = 100;
    too_few_simulations = true;
    while too_few_simulations
        fprintf("%i simulations running\n", simulations);
        %Set up waitbar
        D = parallel.pool.DataQueue;
        h = waitbar(0, '0.0 % complete');
        afterEach(D, @UpdateWaitbar)
        p = 1;

        n = size(Predicted_Mortalities,2);
    
        % Preallocate space for the data matrix
        Death_distribution = zeros(Runs,n+1);
        % This for loop is the key section of the program. It is the only part that takes any significant time.
        for kk = 1:Runs
            %Find the number of deaths in each simulation of the run and store them in a vector containing the frequency distribution of the number of death in the run.
            DeathsForRun = zeros(1,n+1);
            for jj = 1:simulations
                DeathsForSimulation = sum(rand(1,n) < Predicted_Mortalities);
                DeathsForRun(DeathsForSimulation+1) =DeathsForRun(DeathsForSimulation+1)+1;
            end
    
            %Store the frequency distribution of the number of deaths in the run in a matrix containing this data for all the runs.
            Death_distribution(kk,:) = DeathsForRun;
    
            %Update the wait bar
            send(D, kk)
        end
        % Find the probability distributions for each of the runs
        Fractional_death_distribution = Death_distribution ./ sum(Death_distribution,2);
    
        %Apply the first test ...  The highest number of deaths below the lower 95% confidence limit (0.025), must be the same for all the runs.
        cpdf = cumsum(Fractional_death_distribution,2);
        cpdLowerLimitTEST = cpdf < 0.025;
        Test1failed = logical(nnz(sum(cpdLowerLimitTEST) ~= 0 & sum(cpdLowerLimitTEST) ~= Runs));
    
        %Apply the second test ...  The lowest number of deaths above the upper 95% confidence limit (0.975), must be the same for all the runs.
        Reverse_cpdf = cumsum(Fractional_death_distribution,2,'reverse');
        cpdUpperLimitTEST = Reverse_cpdf < 0.025;
        Test2failed = logical(nnz(sum(cpdUpperLimitTEST) ~= 0 & sum(cpdUpperLimitTEST) ~= Runs));
    
        %Calculate the final fractional probablity distribution, combining all the runs into a single run for which the number of simulations = (Runs x simulations)
        RunsCombined_death_distribution =sum(Death_distribution);
        RunsCombined_Fractional_death_distribution = RunsCombined_death_distribution ./ sum(RunsCombined_death_distribution,2);
    
        %Output results if tests have been passed. If not, increase the number of simulations in  aru.
        if (Test1failed||Test2failed)
            fprintf('As at least one of the tests failed, the number of simulations needs to be increased.\n');
        else
            GraphData(outfilename, RunsCombined_Fractional_death_distribution,n);
            savefig(Figoutfilename);
        end
        close(h);
        if ~(Test1failed || Test2failed)
            too_few_simulations = false;
        else
            simulations = simulations * 10;
        end
    end
end
    function UpdateWaitbar(~)
        waitbar(p/Runs, h, [num2str(100*p/Runs), ' % complete'])
        p = p +1;
    end
toc
end
%************************************************************************************************



function [infilename, outfilename, Figoutfilename, Predicted_Mortalities] = UserInput
% The file with the input data is selected by the user. This data contains the predicted mortality for
% each patient as a fraction of 1. The data is in on a single row and the data elements are comma separated.
% This data is read into the array of double "Predicted_Mortalities".


input('\nPress "enter" to select the .csv file containing the predicted mortalities\n')
infilename = uigetfile('*.csv', 'Select the .csv file containing the predicted mortalities\n');
if infilename == 0
    fprintf('You cancelled the program\n\n')
    Predicted_Mortalities = [];
    outfilename = 0;
    Figoutfilename = 0;
else
    Predicted_Mortalities = readmatrix(infilename);
    outfilename = [extractBefore(infilename,'.csv'),'_RESULTS.csv'];
    Figoutfilename = [extractBefore(outfilename,'.csv'),'.fig'];
end
end


%************************************************************************************************

function GraphData(outfilename,RunsCombined_Fractional_death_distribution,n)
%The distribution of the number of deaths is output to a green bar graph in a figure
figure('name',outfilename);
DeathRange = 0:n;
bar(DeathRange,RunsCombined_Fractional_death_distribution(1:n+1),'g');

%The figure is modified so the bars representing numbers of deaths that lie outside the 95 confidence intervals are changed to a red colour.
hold on;

BelowLCI = (cumsum(RunsCombined_Fractional_death_distribution) < 0.025);
AboveUCI = (cumsum(RunsCombined_Fractional_death_distribution,'reverse') < 0.025);
OutsideCI = BelowLCI | AboveUCI;
bar(DeathRange(OutsideCI),RunsCombined_Fractional_death_distribution(OutsideCI),'r');

%The figure is formatted
Fig = gcf;
Fig.Position = [0 98 1193 774];
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
ax.XLim =[0,n+0.5];
ax.XTickLabelRotation = 0;
ax.Box = 'off';
ax.TickDir = 'out';
ax.XLabel.String = 'Number of deaths';
ax.YLabel.String = 'Probability';
ax.Title.String = 'Probability distribution of number of deaths, from individual risk of death predictions';
t = annotation('textbox');
t.String = {'If the observed number of deaths', string(['is in the range ',num2str(nnz(BelowLCI)),' to ',num2str(n - nnz(AboveUCI)), ' (green bars), it lies within']), 'the 95% confidence intervals.'};
t.FontSize = 18;
t.LineStyle = 'none';
t.Position(1:2) = [0.6,0.8];

hold off
end
