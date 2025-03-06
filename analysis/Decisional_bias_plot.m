%% Decisional bais to restructure or repeat negative self-cognitions

% -----------------------------------------------------------------------
% This script reproduces findings on participant negative self-cognition
% endorsment and choice biases reported in the Supplementary Results and
% Supplementary Figure 6.
% 
% Participants' pre-task endorsement of negative self-cognition statements 
% and their choice to restructure or repeate each statement during the 
% cognitive restructuring task were used to calculate bias scores and to
% perform k-mean clustering. These allowed the examination of participants' 
% tendency to repeat or restructure statements based on their overall 
% level of negative self-cognition endorsement.
% -----------------------------------------------------------------------

clear

%% 1: Load Data & create variables

data1 = readtable('../data/CNBTQ_discovery.csv'); % Discovery sample (n = 45)
data2 = readtable('../data/CNBTQ_replication.csv'); % Replication sample (n = 54)

responses1 = data1{:, 1:24};
ratings1 = data1{:, 25:48};
responses2 = data2{:, 1:16};
ratings2 = data2{:, 17:32};

chal_ratings = [];
rep_ratings = [];
bias_scores = [];
avg_ratings = [];

%% 2: Extract Participant Data

% Discovery
for i = 1:size(responses1, 1)
    if any(isnan(responses1(i, :))) || any(isnan(ratings1(i, :)))
        continue;
    end
    participant_responses = responses1(i, :);
    participant_ratings = ratings1(i, :);
    
    chal_idx = participant_responses == 1;
    rep_idx = participant_responses == 2 | participant_responses == 3;
    
    chal_ratings = [chal_ratings, participant_ratings(chal_idx)];
    rep_ratings = [rep_ratings, participant_ratings(rep_idx)];
    
    avg_chal = mean(participant_ratings(chal_idx));
    avg_rep = mean(participant_ratings(rep_idx));
    bias_score = avg_rep - avg_chal; % higher values, greater bias to repeat
    bias_scores = [bias_scores, bias_score];
    
    avg_rating = mean(participant_ratings);
    avg_ratings = [avg_ratings, avg_rating];
end

% Replication
for i = 1:size(responses2, 1)
    if any(isnan(responses2(i, :))) || any(isnan(ratings2(i, :)))
        continue;
    end
    participant_responses = responses2(i, :);
    participant_ratings = ratings2(i, :);
    
    chal_idx = participant_responses == 1;
    rep_idx = participant_responses == 2 | participant_responses == 3;
    
    chal_ratings = [chal_ratings, participant_ratings(chal_idx)];
    rep_ratings = [rep_ratings, participant_ratings(rep_idx)];
    
    avg_chal = mean(participant_ratings(chal_idx));
    avg_rep = mean(participant_ratings(rep_idx));
    bias_score = avg_rep - avg_chal; % higher values, greater bias to repeat
    bias_scores = [bias_scores, bias_score];
    
    avg_rating = mean(participant_ratings);
    avg_ratings = [avg_ratings, avg_rating];
end

%% 4: Endorsement ratings Avg and SD (pre-task)
avg_chal_rating = mean(chal_ratings);
avg_rep_rating = mean(rep_ratings);

std_chal_rating = std(chal_ratings);
std_rep_rating = std(rep_ratings);

disp(['Average CHAL rating: ', num2str(avg_chal_rating)]);
disp(['Average REP rating: ', num2str(avg_rep_rating)]);

%% 5: Save Results in .mat
save('ratings_merged.mat', 'avg_chal_rating', 'avg_rep_rating', 'std_chal_rating', 'std_rep_rating', 'bias_scores', 'avg_ratings');

%% 6: Combined Figure with Three Panels
figure('Units', 'inches', 'Position', [0, 0, 28, 6]);
set(groot, 'DefaultAxesFontName', 'Arial');

subplot(1, 3, 1);
bar_values = [avg_chal_rating, avg_rep_rating];
bar_plot = bar(1:2, bar_values, 'FaceColor', 'flat', 'FaceAlpha', 0.8);
bar_plot.CData(1, :) = [191 044 035]/256; 
bar_plot.CData(2, :) = [058 103 180]/256; 
hold on;
errorbar(1, avg_chal_rating, std_chal_rating, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(2, avg_rep_rating, std_rep_rating, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Challenge', 'Repeat'},'FontSize', 14);
xlabel('Choice Type', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Pre-task endorsement ratings', 'FontSize', 16, 'FontWeight', 'bold');
title('Negative self-cognition endorsement', 'FontSize', 18, 'FontWeight', 'bold');
ylim([0, 7]); 
yticks(0:0.5:7);
hold off;

subplot(1, 3, 2);
[kde_density, kde_x] = ksdensity(bias_scores, 'Bandwidth', 0.2);
fill([kde_x, fliplr(kde_x)], [kde_density, zeros(size(kde_density))], [191 44 35]/256, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
set(gca, 'XDir')%, 'reverse');
set(gca, 'FontSize', 14);
xlabel('Bias score','FontSize', 16, 'FontWeight', 'bold');
ylabel('Density','FontSize', 16, 'FontWeight', 'bold');
title('Decisional bias score','FontSize', 18, 'FontWeight', 'bold');
hold off;

data_for_clustering = [avg_ratings(:), bias_scores(:)];
num_clusters = 2;
[idx, C] = kmeans(data_for_clustering, num_clusters);

subplot(1, 3, 3);
scatter(avg_ratings(idx == 1), bias_scores(idx == 1), 150, [064 176 166]/256, 'filled');
hold on;
scatter(avg_ratings(idx == 2), bias_scores(idx == 2), 150, [255 190 106]/256, 'filled');
set(gca, 'FontSize', 14);
xlabel('Pre-task endorsement ratings', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Bias score', 'FontSize', 16, 'FontWeight', 'bold');
title('Decisional bias clustering','FontSize', 18, 'FontWeight', 'bold');
coeffs = polyfit(avg_ratings, bias_scores, 1);
x_fit = linspace(min(avg_ratings), max(avg_ratings), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k--', 'LineWidth', 1.5);
hold off;

%% 7: Report Statistical Differences
[~, p_avg_rating] = ttest2(avg_ratings(idx == 1), avg_ratings(idx == 2));
[~, p_bias_score] = ttest2(bias_scores(idx == 1), bias_scores(idx == 2));

[corr_coeff, p_corr] = corr(avg_ratings(:), bias_scores(:), 'Type', 'Pearson');

disp(['p-value for Average Rating difference between clusters: ', num2str(p_avg_rating)]);
disp(['p-value for Bias Score difference between clusters: ', num2str(p_bias_score)]);
disp(['Pearson correlation coefficient (Bias Scores vs. Avg Ratings): ', num2str(corr_coeff)]);
disp(['p-value for Pearson correlation: ', num2str(p_corr)]);

%% 8: save figure
exportgraphics(gcf, 'FigRatings.jpeg', 'Resolution', 300);
exportgraphics(gcf, 'FigRatings.tiff', 'Resolution', 300);
