function TenureAreaPlots()
%
% This function loads sandfish count data from 'Sandfish_count_data.xlsx'
% and produces violin plots comparing sandfish abundances in two management 
% areas: LMMA Areas (TA) and Open Access Areas (OA). The analysis focuses 
% on four groups:
%   1. Timoenai  - Split from the Timoenai station where the site does NOT
%                  contain the string "Tawi" in column 4.
%   2. Tawi      - Split from the Timoenai station where the site DOES contain
%                  the string "Tawi" (e.g., "Tawi 1").
%   3. Pere      - Where the station (column 3) equals "Pere".
%   4. Mbunai    - Where the station (column 3) equals "Mbunai".
%
% For each group, the function compares sandfish abundances over three time 
% periods ('Before', 'After1', 'After2') for each management area (OA and TA).
% Bootstrap resampling (with NBS resamples) is used to estimate the mean 
% density. Violin plots are generated for each combination.
%
% In the final plot, the fixed x-axis positions are:
%   1 for Timoenai, 2 for Tawi, 3 for Pere, and 4 for Mbunai.
%
% For LMMA Areas (TA; plotted in green) the top of each violin plot is 
% annotated vertically with a year corresponding to the time period:
%   'Before'  -> Feb 2017
%   'After1'  -> Oct 2017
%   'After2'  -> Feb 2018
%
% The analysis parameter "Adult" controls whether only adult sandfish are 
% analyzed (Adult = 1) or all sandfish (Adult = 0).
%
% Author: <Your Name>
% Date: <Today's Date>

% Parameters and Data Loading

% Set analysis flag: 1 = adult sandfish only, 0 = all sandfish
Adult = 1;

% Load data from Excel file
[D, T] = xlsread('Sandfish_count_data.xlsx');

% Extract headers and metadata
Labels = T(1, :);   % Column headers (for reference)
Info = T(2:end, :); % Metadata (columns assumed as below)
% Assumptions for column indices:
%   Column 3: Station name (e.g., Timoenai, Pere, Mbunai)
%   Column 4: Site information (e.g., "Tawi 1", which distinguishes Tawi subset)
%   Column 5: Area type ('OA' or 'TA')
%   Column 6: Adult sandfish counts
%   Column 8: Overall sandfish counts and/or sampling period
%   Column 9: Sampling period (e.g., 'Before', 'After1', 'After2')

% Choose count data:
if Adult == 1
    CountData = 100 * D(:, 6);   % Use adult counts scaled by 100 (density per hectare)
else
    CountData = 100 * D(:, 8);     % Use overall counts scaled by 100
end

% Plotting Parameters and Group Definitions

% Plot formatting settings
FS = 14;           % Font size for axis labels and title
FS2 = 10;          % Font size for text annotations
NBS = 5000;        % Number of bootstrap resamples
% CLR = [0.4 0.8 0.3];   % Color for LMMA Areas (TA, greenish)
CLR(1,:) = [0.64 0.88 0.58];
CLR(2,:) = [0.32 0.64 0.24];
CLR(3,:) = [0.24 0.48 0.18];
COL = ones(1, 3) * 0.6; % Color for Open Access Areas (OA, grey)
width = 0.25;           % Width factor for violin plots
BN = 15;               % Number of bins for histogram used in violin plots

% Define group names and fixed x-axis positions:
%   Group names: 'Timoenai', 'Tawi', 'Pere', 'Mbunai'
groupNames = {'Timoenai', 'Tawi', 'Pere', 'Mbunai'};
xPositions = [1, 2, 3, 4]*1.2;

% Define the time periods and management area types:
timePeriods = {'Before', 'After1', 'After2'};
areaTypes   = {'OA', 'TA'};

% Prepare figure:
cla;
hold on;
box on;
pp = patch([-10 10 10 -10],[-10 -10 50 50],[0.6 0 0]);
set(pp,'edgecolor','none','facealpha',0.1)

% Analysis Loop over Groups, Areas, and Time Periods

for g = 1:length(groupNames)
    group = groupNames{g};
    base_x = xPositions(g);  % Fixed x-axis position for this group
    
    % Select rows based on group.
    idxGroup = contains(Info(:,3), group);

    % Loop over management areas (OA and TA)
    for a = 1:length(areaTypes)
        area = areaTypes{a};
        
        % Loop over time periods: 'Before', 'After1', 'After2'
        for t = 1:length(timePeriods)
            tp = timePeriods{t};
            
            % Logical index for rows matching current group, area, and time period.
            % Column 5 is assumed to denote area type and column 9 the sampling period.
            idx = idxGroup & strcmp(Info(:,4), area) & strcmp(Info(:,8), tp);
            Data = CountData(idx);

            % Bootstrap resampling to estimate the mean density.
            bootMeans = zeros(1, NBS);
            for b = 1:NBS
                bootMeans(b) = mean(randsample(Data, length(Data), true));
            end
            
            % Select plot color based on area type.
            if strcmp(area, 'OA')
                plotColor = COL;
            else  % area is 'TA'
                plotColor = CLR(t,:);
            end
            
            % For visual separation, add a small horizontal offset:
            % Offsets: -width, 0, +width for the three time periods respectively.
            offset = (-width) + (t - 1) * width;
            current_x = base_x + offset;
            
            % Plot the violin plot for the current combination.
            if strcmp(area,'OA') == 1
                YP = Plot_violin(current_x, bootMeans, 2*BN/3, width, plotColor, 5);
            else
                YP = Plot_violin(current_x, bootMeans, BN, width, plotColor, 5);
            end
            
        end  % end time period loop
    end  % end area type loop
end  % end group loop

% Final Plot Settings
set(gca, 'XTick', xPositions, 'XTickLabel', groupNames,'fontsize',FS-2);
xlabel('Group','fontsize',FS);
ylabel('Adults per hectare','fontsize',FS+2)
xlim([0.5 5.5])
ylim([0 2000])
hold off;

% Create offscreen patch handles for custom legend entries
h1 = patch('XData', nan, 'YData', nan, 'FaceColor', CLR(1,:), 'EdgeColor', 'none');
h2 = patch('XData', nan, 'YData', nan, 'FaceColor', CLR(2,:), 'EdgeColor', 'none');
h3 = patch('XData', nan, 'YData', nan, 'FaceColor', CLR(3,:), 'EdgeColor', 'none');
h4 = patch('XData', nan, 'YData', nan, 'FaceColor', COL, 'EdgeColor', 'none');

% Create the legend using the patch handles
legend([h1, h2, h3, h4], {'LMMA Mar 2017', 'LMMA Oct 2017', 'LMMA Sep 2019', 'Open access'}, ...
       'Location', 'northeast', 'Box', 'off','fontsize',FS);

end

% Helper Function: Plot_violin
function YP = Plot_violin(x, data, bins, width, CL, THR)
% PLOT_VIOLIN Creates a violin plot at the specified x-coordinate.
%
%   YP = Plot_violin(x, data, bins, width, CL, THR) generates a violin 
%   plot representing the distribution of 'data' using a normalized histogram.
%
%   Inputs:
%     x     - x-axis coordinate for the center of the violin plot.
%     data  - Vector of bootstrap mean estimates.
%     bins  - Number of bins for constructing the histogram.
%     width - Scaling factor for horizontal spread.
%     CL    - Color for the plot (patch fill and marker).
%     THR   - (Unused, kept for compatibility).
%
%   Output:
%     YP    - Representative y-coordinate (mean of histogram edges) used for
%             positioning text labels.
%
% Compute histogram edges and normalized counts.
Edges = linspace(min(data), max(data), bins);
H = histc(data, Edges);
H = H / sum(H);  % Normalize counts

% Create a patch for each histogram bin.
for e = 1:length(Edges)-1
    patch(x + 1.2.*width .* [-H(e) H(e) H(e) -H(e)], Edges([e e e+1 e+1]), CL, 'EdgeColor', 'none');
end

% Plot the mean as a point.
plot(x, mean(data), '.', 'MarkerSize', 12, 'Color', CL);

% Return a representative y-position.
YP = mean(Edges);
end