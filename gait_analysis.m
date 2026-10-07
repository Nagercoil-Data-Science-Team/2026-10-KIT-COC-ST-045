clc;
clear;
close all;
%% ============================================================
%  STEP 1: DATASET FOLDER
% =============================================================
dataFolder = 'data';
%% ============================================================
%  STEP 2: FIND ALL FILES
% =============================================================
markerFiles = dir(fullfile(dataFolder, '*_marker.csv'));
forceFiles = dir(fullfile(dataFolder, '*_force.csv'));
overstepFiles = dir(fullfile(dataFolder, '*_oversteps.csv'));
%% ============================================================
%  STEP 3: DISPLAY TOTAL NUMBER OF FILES
% =============================================================
fprintf('========================================\n');
fprintf('          DATASET INFORMATION\n');
fprintf('========================================\n');
fprintf('Marker files   : %d\n', length(markerFiles));
fprintf('Force files    : %d\n', length(forceFiles));
fprintf('Overstep files : %d\n', length(overstepFiles));
fprintf('Total files    : %d\n', length(markerFiles) + length(forceFiles) + length(overstepFiles));
%% ============================================================
%  STEP 4: LOAD ALL MARKER FILES
% =============================================================
MarkerData = cell(length(markerFiles), 1);
for i = 1:length(markerFiles)
    filePath = fullfile(dataFolder, markerFiles(i).name);
    MarkerData{i} = readtable(filePath);
end
%% ============================================================
%  STEP 5: LOAD ALL FORCE FILES
% =============================================================
ForceData = cell(length(forceFiles), 1);
for i = 1:length(forceFiles)
    filePath = fullfile(dataFolder, forceFiles(i).name);
    ForceData{i} = readtable(filePath);
end
%% ============================================================
%  STEP 6: LOAD ALL OVERSTEP FILES
% =============================================================
OverstepData = cell(length(overstepFiles), 1);
for i = 1:length(overstepFiles)
    filePath = fullfile(dataFolder, overstepFiles(i).name);
    OverstepData{i} = readtable(filePath);
end
fprintf('\n========================================\n');
fprintf('          ORGANIZING TRIALS\n');
fprintf('========================================\n');
%% ============================================================
%  STEP 7: CREATE TRIAL NAMES
% =============================================================
trialNames = cell(length(markerFiles), 1);
for i = 1:length(markerFiles)
    name = markerFiles(i).name;
    trialNames{i} = erase(name, '_marker.csv');
end
fprintf('\nFirst 10 trials:\n');
for i = 1:min(10, length(trialNames))
    fprintf('%d. %s\n', i, trialNames{i});
end
%% ============================================================
%  STEP 8: CREATE ORGANIZED TRIAL STRUCTURE
% ============================================================
Trials = struct();
for i = 1:length(trialNames)
    trialName = trialNames{i};
    Trials(i).Name = trialName;
    Trials(i).MarkerFile = [trialName '_marker.csv'];
    Trials(i).ForceFile = [trialName '_force.csv'];
    Trials(i).OverstepFile = [trialName '_oversteps.csv'];
    % Marker data
    Trials(i).Marker = MarkerData{i};
    % Find corresponding force file
    forceIndex = find(strcmp({forceFiles.name}, [trialName '_force.csv']), 1);
    if ~isempty(forceIndex)
        Trials(i).Force = ForceData{forceIndex};
    else
        Trials(i).Force = [];
    end
    % Find corresponding overstep file
    overstepIndex = find(strcmp({overstepFiles.name}, [trialName '_oversteps.csv']), 1);
    if ~isempty(overstepIndex)
        Trials(i).Overstep = OverstepData{overstepIndex};
    else
        Trials(i).Overstep = [];
    end
end
%% ============================================================
%  STEP 9: SAMPLING FREQUENCIES
% ============================================================
Fs_marker = 200;       % Marker data
Fs_force  = 1000;      % Force data
fprintf('\n========================================\n');
fprintf('       SAMPLING FREQUENCIES\n');
fprintf('========================================\n');
fprintf('Marker sampling frequency : %d Hz\n', Fs_marker);
fprintf('Force sampling frequency  : %d Hz\n', Fs_force);
%% ============================================================
%  STEP 10: CREATE TIME VECTOR
% ============================================================
for i = 1:length(Trials)
    % Marker time
    nMarker = height(Trials(i).Marker);
    Trials(i).TimeMarker = (0:nMarker-1)' / Fs_marker;
    % Force time
    if ~isempty(Trials(i).Force)
        nForce = height(Trials(i).Force);
        Trials(i).TimeForce = (0:nForce-1)' / Fs_force;
    else
        Trials(i).TimeForce = [];
    end
    % Overstep event time
    if ~isempty(Trials(i).Overstep)
        Trials(i).TimeOverstep = Trials(i).Overstep.Var1;
    else
        Trials(i).TimeOverstep = [];
    end
end
%% ============================================================
%  STEP 11: DISPLAY ONE ORGANIZED TRIAL
% ============================================================
fprintf('\n========================================\n');
fprintf('          ORGANIZED TRIAL\n');
fprintf('========================================\n');
fprintf('Trial name: %s\n', Trials(1).Name);
fprintf('\nMarker file:\n');
fprintf('%s\n', Trials(1).MarkerFile);
fprintf('\nForce file:\n');
fprintf('%s\n', Trials(1).ForceFile);
fprintf('\nOverstep file:\n');
fprintf('%s\n', Trials(1).OverstepFile);
%% ============================================================
%  STEP 12: DISPLAY DATA SIZE
% ============================================================
fprintf('\n========================================\n');
fprintf('          DATA SIZE\n');
fprintf('========================================\n');
fprintf('Marker samples   : %d\n', height(Trials(1).Marker));
fprintf('Marker columns   : %d\n', width(Trials(1).Marker));
fprintf('Force samples    : %d\n', height(Trials(1).Force));
fprintf('Force columns    : %d\n', width(Trials(1).Force));
fprintf('Overstep events  : %d\n', height(Trials(1).Overstep));
%% ============================================================
%  STEP 13: DISPLAY TIME INFORMATION
% ============================================================
fprintf('\n========================================\n');
fprintf('          MARKER TIME VECTOR\n');
fprintf('========================================\n');
disp(Trials(1).TimeMarker(1:min(10,end)));
fprintf('\n========================================\n');
fprintf('          FORCE TIME VECTOR\n');
fprintf('========================================\n');
disp(Trials(1).TimeForce(1:min(10,end)));
fprintf('\n========================================\n');
fprintf('          OVERSTEP TIMES\n');
fprintf('========================================\n');
disp(Trials(1).TimeOverstep);
%% ============================================================
%  STEP 14: DISPLAY MARKER DATA
% ============================================================
fprintf('\n========================================\n');
fprintf('          FIRST MARKER DATA\n');
fprintf('========================================\n');
disp(Trials(1).Marker(1:min(10, height(Trials(1).Marker)), :));
%% ============================================================
%  STEP 15: DISPLAY FORCE DATA
% ============================================================
fprintf('\n========================================\n');
fprintf('          FIRST FORCE DATA\n');
fprintf('========================================\n');
disp(Trials(1).Force(1:min(10, height(Trials(1).Force)), :));
%% ============================================================
%  STEP 16: DISPLAY OVERSTEP DATA
% ============================================================
fprintf('\n========================================\n');
fprintf('          OVERSTEP DATA\n');
fprintf('========================================\n');
disp(Trials(1).Overstep);
%% ============================================================
%  STEP 17: DISPLAY COLUMN NAMES
% ============================================================
fprintf('\n========================================\n');
fprintf('          COLUMN INFORMATION\n');
fprintf('========================================\n');
fprintf('\nMARKER COLUMNS:\n');
disp(Trials(1).Marker.Properties.VariableNames');
fprintf('\nFORCE COLUMNS:\n');
disp(Trials(1).Force.Properties.VariableNames');
fprintf('\nOVERSTEP COLUMNS:\n');
disp(Trials(1).Overstep.Properties.VariableNames');
fprintf('\n========================================\n');
fprintf('       LOADING COMPLETED\n');
fprintf('========================================\n');
%% ============================================================
%  STEP 18: SELECT TRIAL FOR GAIT ANALYSIS
% ============================================================
trialID = 1;
marker = Trials(trialID).Marker;
t = Trials(trialID).TimeMarker;
markerRaw = table2array(marker);
Fs = Fs_marker;
% Plots that check the processing (set to false to hide them)
showValidationPlots = true;
%% ============================================================
%  STEP 18.1: DETECT MARKER UNITS AND CONVERT TO METRES
% ============================================================
maxRange = max(max(markerRaw, [], 1) - min(markerRaw, [], 1));
if maxRange > 20
    unitScale = 0.001;
    unitName = 'millimetres';
else
    unitScale = 1;
    unitName = 'metres';
end
markerRaw = markerRaw * unitScale;
fprintf('\n========================================\n');
fprintf('       MARKER DATA PREPROCESSING\n');
fprintf('========================================\n');
fprintf('Trial                : %s\n', Trials(trialID).Name);
fprintf('Number of samples    : %d\n', size(markerRaw,1));
fprintf('Number of columns    : %d\n', size(markerRaw,2));
fprintf('Detected marker units: %s (all data converted to metres)\n', unitName);
%% ============================================================
%  STEP 19.1: CHECK MISSING VALUES
% ============================================================
missingBefore = sum(isnan(markerRaw), 'all');
fprintf('\nMissing values before preprocessing : %d\n', missingBefore);
%% ============================================================
%  STEP 19.2: FILL MISSING VALUES
% ============================================================
markerClean = markerRaw;
for c = 1:size(markerClean,2)
    markerClean(:,c) = fillmissing(markerClean(:,c), 'linear');
    % Handle missing values at beginning/end
    markerClean(:,c) = fillmissing(markerClean(:,c), 'nearest');
end
missingAfter = sum(isnan(markerClean), 'all');
fprintf('Missing values after interpolation  : %d\n', missingAfter);
fprintf('\nDetecting marker spikes...\n');
windowSamples = 11;      % moving-median window (samples)
nSigma = 6;              % robust sigma multiplier
minSpikeSize = 0.010;    % minimum spike size (m) = 10 mm
maxSpikeWidth = 3;       % spikes are at most this many samples wide
markerOutlierCorrected = markerClean;
outlierMask = false(size(markerClean));
sampleIndex = (1:size(markerClean,1))';
for c = 1:size(markerClean,2)
    x = markerClean(:,c);
    localMedian = movmedian(x, windowSamples);
    residual = x - localMedian;
    robustSigma = 1.4826 * median(abs(residual - median(residual)));
    threshold = max(nSigma * robustSigma, minSpikeSize);
    spikeIndex = abs(residual) > threshold;
    spikeIndex = removeLongRuns(spikeIndex, maxSpikeWidth);
    outlierMask(:,c) = spikeIndex;
    if any(spikeIndex)
        goodIndex = ~spikeIndex;
        x(spikeIndex) = interp1(sampleIndex(goodIndex), x(goodIndex), sampleIndex(spikeIndex), 'linear', 'extrap');
        markerOutlierCorrected(:,c) = x;
    end
end
outlierCount = sum(outlierMask, 'all');
totalMarkerValues = numel(markerRaw);
outlierPercentage = 100 * outlierCount / totalMarkerValues;
fprintf('Outlier values detected            : %d\n', outlierCount);
fprintf('Total marker values                : %d\n', totalMarkerValues);
fprintf('Outlier percentage                 : %.3f %%\n', outlierPercentage);
fprintf('\nSpikes per marker column:\n');
markerNames = marker.Properties.VariableNames;
for c = 1:size(markerClean,2)
    fprintf('%-10s : %5d (%.3f %%)\n', markerNames{c}, sum(outlierMask(:,c)), 100 * sum(outlierMask(:,c)) / size(markerClean,1));
end
if outlierPercentage > 1
    fprintf('WARNING: more than 1 %% of the values were corrected. Inspect the validation plot.\n');
end
%% ============================================================
%  STEP 19.4: LOW-PASS FILTERING
% ============================================================
Fc = 6;
fprintf('\nApplying low-pass Butterworth filter...\n');
% Second-order Butterworth applied forwards and backwards (zero phase)
[b,a] = butter(2, Fc/(Fs/2), 'low');
markerFiltered = zeros(size(markerOutlierCorrected));
for c = 1:size(markerOutlierCorrected,2)
    markerFiltered(:,c) = filtfilt(b, a, markerOutlierCorrected(:,c));
end
fprintf('Low-pass filtering completed.\n');
fprintf('Sampling frequency : %d Hz\n', Fs);
fprintf('Cut-off frequency  : %d Hz\n', Fc);
fprintf('Filter order       : 2 (zero phase)\n');
%% ============================================================
%  STEP 19.5: FINAL VALIDATION
% ============================================================
remainingNaN = sum(isnan(markerFiltered), 'all');
fprintf('\nFinal preprocessing validation:\n');
fprintf('NaN values after preprocessing : %d\n', remainingNaN);
%% ============================================================
%  STEP 19.6: PREPROCESSING VALIDATION PLOT
%  Raw vs corrected vs filtered
% ============================================================
if showValidationPlots
    [~, worstColumn] = max(sum(outlierMask, 1));
    spikeSamples = find(outlierMask(:,worstColumn));
    if isempty(spikeSamples)
        zoomCentre = round(size(markerRaw,1) / 2);
    else
        zoomCentre = spikeSamples(1);
    end
    zoomIndex = max(1, zoomCentre - 3*Fs) : min(size(markerRaw,1), zoomCentre + 3*Fs);
    figure('Name', 'Preprocessing Validation');
    subplot(2,1,1);
    plot(t, markerRaw(:,worstColumn), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    hold on;
    plot(t, markerOutlierCorrected(:,worstColumn), 'b', 'LineWidth', 1.0);
    plot(t, markerFiltered(:,worstColumn), 'r', 'LineWidth', 1.5);
    plot(t(spikeSamples), markerRaw(spikeSamples,worstColumn), 'ko', 'MarkerSize', 4);
    xlabel('Time (s)');
    ylabel('Position (m)');
    title(['Marker ' markerNames{worstColumn} ' (column with most corrected values)'], 'Interpreter', 'none');
    legend('Raw', 'Outlier corrected', 'Filtered', 'Corrected samples', 'Location', 'best');
    grid on;
    subplot(2,1,2);
    plot(t(zoomIndex), markerRaw(zoomIndex,worstColumn), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    hold on;
    plot(t(zoomIndex), markerOutlierCorrected(zoomIndex,worstColumn), 'b', 'LineWidth', 1.0);
    plot(t(zoomIndex), markerFiltered(zoomIndex,worstColumn), 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Position (m)');
    title('Zoomed view');
    grid on;
end
%% ============================================================
%  STEP 20: CREATE FILTERED MARKER TABLE
% ============================================================
markerFilteredTable = array2table(markerFiltered, 'VariableNames', marker.Properties.VariableNames);
%% ============================================================
%  STEP 21: SEPARATE LEFT AND RIGHT FOOT MARKERS
% ============================================================
fprintf('\n========================================\n');
fprintf('       SEPARATING FOOT MARKERS\n');
fprintf('========================================\n');
L_FCC = [markerFilteredTable.L_FCC_x, markerFilteredTable.L_FCC_y, markerFilteredTable.L_FCC_z];
L_FM1 = [markerFilteredTable.L_FM1_x, markerFilteredTable.L_FM1_y, markerFilteredTable.L_FM1_z];
L_FM2 = [markerFilteredTable.L_FM2_x, markerFilteredTable.L_FM2_y, markerFilteredTable.L_FM2_z];
L_FM5 = [markerFilteredTable.L_FM5_x, markerFilteredTable.L_FM5_y, markerFilteredTable.L_FM5_z];
R_FCC = [markerFilteredTable.R_FCC_x, markerFilteredTable.R_FCC_y, markerFilteredTable.R_FCC_z];
R_FM1 = [markerFilteredTable.R_FM1_x, markerFilteredTable.R_FM1_y, markerFilteredTable.R_FM1_z];
R_FM2 = [markerFilteredTable.R_FM2_x, markerFilteredTable.R_FM2_y, markerFilteredTable.R_FM2_z];
R_FM5 = [markerFilteredTable.R_FM5_x, markerFilteredTable.R_FM5_y, markerFilteredTable.R_FM5_z];
fprintf('Left foot markers separated  : FCC, FM1, FM2, FM5\n');
fprintf('Right foot markers separated : FCC, FM1, FM2, FM5\n');
% Representative foot position (centroid of the four foot markers)
L_Foot = (L_FCC + L_FM1 + L_FM2 + L_FM5) / 4;
R_Foot = (R_FCC + R_FM1 + R_FM2 + R_FM5) / 4;
% Heel marker and forefoot centroid (used for gait events and foot angle)
L_Heel = L_FCC;
R_Heel = R_FCC;
L_Fore = (L_FM1 + L_FM2 + L_FM5) / 3;
R_Fore = (R_FM1 + R_FM2 + R_FM5) / 3;
%% ============================================================
%  STEP 22: FOOT TRAJECTORY RANGE (METRES)
% ============================================================
fprintf('\nLeft Foot Trajectory (m):\n');
fprintf('X range : %.4f to %.4f\n', min(L_Foot(:,1)), max(L_Foot(:,1)));
fprintf('Y range : %.4f to %.4f\n', min(L_Foot(:,2)), max(L_Foot(:,2)));
fprintf('Z range : %.4f to %.4f\n', min(L_Foot(:,3)), max(L_Foot(:,3)));
fprintf('\nRight Foot Trajectory (m):\n');
fprintf('X range : %.4f to %.4f\n', min(R_Foot(:,1)), max(R_Foot(:,1)));
fprintf('Y range : %.4f to %.4f\n', min(R_Foot(:,2)), max(R_Foot(:,2)));
fprintf('Z range : %.4f to %.4f\n', min(R_Foot(:,3)), max(R_Foot(:,3)));
%% ============================================================
%  STEP 23: FOOT VELOCITY AND SPEED (m/s)
%  Each column is differentiated separately over the time samples.
% ============================================================
L_Velocity = zeros(size(L_Foot));
R_Velocity = zeros(size(R_Foot));
for c = 1:3
    L_Velocity(:,c) = gradient(L_Foot(:,c), 1/Fs);
    R_Velocity(:,c) = gradient(R_Foot(:,c), 1/Fs);
end
L_Speed = sqrt(sum(L_Velocity.^2, 2));
R_Speed = sqrt(sum(R_Velocity.^2, 2));
L_MaxSpeed = max(L_Speed);
R_MaxSpeed = max(R_Speed);
fprintf('\nLeft foot  speed (m/s) : mean %.3f | 99th percentile %.3f | max %.3f\n', mean(L_Speed), pctl(L_Speed, 99), L_MaxSpeed);
fprintf('Right foot speed (m/s) : mean %.3f | 99th percentile %.3f | max %.3f\n', mean(R_Speed), pctl(R_Speed, 99), R_MaxSpeed);
if max(L_MaxSpeed, R_MaxSpeed) > 12
    fprintf('WARNING: foot speed above 12 m/s is not physical. Check units and spikes.\n');
end
%% ============================================================
%  STEP 24: OVERSTEP EVENTS
%  NOT used as gait-cycle boundaries. They are only used to FLAG
%  cycles that contain an event. Confirm the meaning of Var1 in
%  the dataset documentation.
% ============================================================
fprintf('\n========================================\n');
fprintf('       OVERSTEP EVENTS\n');
fprintf('========================================\n');
overstepTimes = Trials(trialID).TimeOverstep;
overstepTimes = sort(overstepTimes(~isnan(overstepTimes)));
overstepTimes = overstepTimes(:);
fprintf('Overstep events in trial : %d\n', length(overstepTimes));
for k = 1:length(overstepTimes)
    fprintf('Event %d : %.3f s\n', k, overstepTimes(k));
end
%% ============================================================
%  STEP 25: WALKING DIRECTION AND SAGITTAL COORDINATES
%  The forward axis is taken from the direction the foot points
%  (heel to forefoot). Works for treadmill and overground walking.
% ============================================================
fprintf('\n========================================\n');
fprintf('       WALKING DIRECTION\n');
fprintf('========================================\n');
footVector = [L_Fore - L_Heel; R_Fore - R_Heel];
[~, fwdAxis] = max(mean(abs(footVector(:,1:2)), 1));
latAxis = setdiff([1 2], fwdAxis);
fwdSign = sign(median(footVector(:,fwdAxis)));
if fwdSign == 0
    fwdSign = 1;
end
axisLabels = {'X', 'Y', 'Z'};
fprintf('Forward axis  : %s (sign %+d)\n', axisLabels{fwdAxis}, fwdSign);
fprintf('Lateral axis  : %s\n', axisLabels{latAxis});
fprintf('Vertical axis : Z\n');
% Forward (f) and vertical (z) coordinates, forward positive
L_f_heel = fwdSign * L_Heel(:,fwdAxis);
R_f_heel = fwdSign * R_Heel(:,fwdAxis);
L_f_fore = fwdSign * L_Fore(:,fwdAxis);
R_f_fore = fwdSign * R_Fore(:,fwdAxis);
L_z_heel = L_Heel(:,3);
R_z_heel = R_Heel(:,3);
L_z_fore = L_Fore(:,3);
R_z_fore = R_Fore(:,3);
%% ============================================================
%  STEP 26: HEEL-STRIKE DETECTION (KINEMATIC)
%  Heel strike = most forward position of the heel relative to
%  the body (hip proxy). The hip proxy is a moving average of the
%  feet over one stride, so it works on a treadmill and overground.
% ============================================================
fprintf('\n========================================\n');
fprintf('       HEEL-STRIKE DETECTION\n');
fprintf('========================================\n');
strideTime = 1.1;
for pass = 1:2
    strideWindow = max(round(strideTime * Fs), 3);
    hipF = movmean(0.5 * (L_f_heel + R_f_heel), strideWindow);
    hsL = detectHeelStrikes(L_f_heel - hipF, Fs);
    hsR = detectHeelStrikes(R_f_heel - hipF, Fs);
    if length(hsL) >= 3 && length(hsR) >= 3
        strideTime = median([diff(hsL); diff(hsR)]) / Fs;
    end
end
if length(hsL) < 3 || length(hsR) < 3
    error('Heel-strike detection failed. Check the marker data and the forward axis.');
end
fprintf('Left heel strikes   : %d\n', length(hsL));
fprintf('Right heel strikes  : %d\n', length(hsR));
fprintf('Median stride time  : %.3f s\n', strideTime);
fprintf('Left stride time    : %.3f +/- %.3f s\n', mean(diff(hsL))/Fs, std(diff(hsL))/Fs);
fprintf('Right stride time   : %.3f +/- %.3f s\n', mean(diff(hsR))/Fs, std(diff(hsR))/Fs);
if showValidationPlots
    figure('Name', 'Heel-Strike Detection Check');
    showIndex = 1:min(length(t), 10*Fs);
    plot(t(showIndex), L_f_heel(showIndex) - hipF(showIndex), 'b', 'LineWidth', 1.2);
    hold on;
    plot(t(showIndex), R_f_heel(showIndex) - hipF(showIndex), 'r', 'LineWidth', 1.2);
    hsLshow = hsL(hsL <= showIndex(end));
    hsRshow = hsR(hsR <= showIndex(end));
    plot(t(hsLshow), L_f_heel(hsLshow) - hipF(hsLshow), 'bv', 'MarkerFaceColor', 'b');
    plot(t(hsRshow), R_f_heel(hsRshow) - hipF(hsRshow), 'rv', 'MarkerFaceColor', 'r');
    for k = 1:length(overstepTimes)
        if overstepTimes(k) <= t(showIndex(end))
            xline(overstepTimes(k), '--k');
        end
    end
    xlabel('Time (s)');
    ylabel('Heel position relative to body (m)');
    title('Detected heel strikes (triangles) - first 10 s');
    legend('Left heel', 'Right heel', 'Left HS', 'Right HS', 'Location', 'best');
    grid on;
end
%% ============================================================
%  STEP 27: BIOMECHANICAL MODEL PARAMETERS
%  Joint angles are MODEL-BASED ESTIMATES (inverse kinematics),
%  NOT measured joint angles. There are no hip or knee markers.
% ============================================================
fprintf('\n========================================\n');
fprintf('       BIOMECHANICAL MODEL\n');
fprintf('========================================\n');
% --- Edit these values for your subject ----------------------
subjectHeight = 1.75;                   % subject height (m)
thighLength   = 0.245 * subjectHeight;  % hip to knee (m)
shankLength   = 0.246 * subjectHeight;  % knee to ankle (m)
ankleHeight   = 0.039 * subjectHeight;  % ankle joint above the heel-forefoot line (m)
ankleFwdRatio = 0.30;                   % ankle position along heel to forefoot line
% --------------------------------------------------------------
fprintf('Thigh length  : %.3f m\n', thighLength);
fprintf('Shank length  : %.3f m\n', shankLength);
fprintf('Ankle height  : %.3f m\n', ankleHeight);
%% ============================================================
%  STEP 28: INVERSE KINEMATICS (CONTINUOUS, WHOLE TRIAL)
%  Hip: moves forward with the body (moving average), constant
%  height, chosen automatically so both legs stay within reach.
% ============================================================
fprintf('\n========================================\n');
fprintf('       INVERSE KINEMATICS\n');
fprintf('========================================\n');
% Ankle joint centre from the foot markers
[L_fa, L_za] = ankleFromFoot(L_f_heel, L_z_heel, L_f_fore, L_z_fore, ankleFwdRatio, ankleHeight);
[R_fa, R_za] = ankleFromFoot(R_f_heel, R_z_heel, R_f_fore, R_z_fore, ankleFwdRatio, ankleHeight);
% Hip height: highest value that keeps both legs reachable (2nd percentile of the reach limit)
reach = 0.998 * (thighLength + shankLength);
zCapL = L_za + sqrt(max(reach^2 - (L_fa - hipF).^2, 0));
zCapR = R_za + sqrt(max(reach^2 - (R_fa - hipF).^2, 0));
hipZ = pctl([zCapL; zCapR], 2);
fprintf('Hip height (z)          : %.3f m (%.1f %% of subject height)\n', hipZ, 100 * hipZ / subjectHeight);
fprintf('Typical hip height is about 50-56 %% of body height (if the floor is z = 0).\n');
if hipZ / subjectHeight < 0.45 || hipZ / subjectHeight > 0.62
    fprintf('WARNING: hip height is unusual. Check subjectHeight and the segment lengths.\n');
end
% Two-link inverse kinematics for each leg
[L_hipFlex, L_kneeFlex, L_shank, L_unreach] = legIK(L_fa, L_za, hipF, hipZ, thighLength, shankLength);
[R_hipFlex, R_kneeFlex, R_shank, R_unreach] = legIK(R_fa, R_za, hipF, hipZ, thighLength, shankLength);
fprintf('Unreachable samples - left  : %d\n', sum(L_unreach));
fprintf('Unreachable samples - right : %d\n', sum(R_unreach));
% Foot pitch (+ toes up), zero when the foot is flat on the ground
L_pitch = atan2(L_z_fore - L_z_heel, L_f_fore - L_f_heel);
R_pitch = atan2(R_z_fore - R_z_heel, R_f_fore - R_f_heel);
L_pitchFlat = flatFootPitch(L_pitch, L_z_heel, L_z_fore);
R_pitchFlat = flatFootPitch(R_pitch, R_z_heel, R_z_fore);
fprintf('Flat-foot pitch offset  : left %.2f deg | right %.2f deg (removed)\n', rad2deg(L_pitchFlat), rad2deg(R_pitchFlat));
% Ankle dorsiflexion = foot pitch - shank angle (0 = neutral, + = dorsiflexion)
L_ankleDorsi = (L_pitch - L_pitchFlat) - L_shank;
R_ankleDorsi = (R_pitch - R_pitchFlat) - R_shank;
% Continuous signals used for the gait cycles:
% [hip, knee, ankle (deg), foot forward relative to hip (m), foot lateral (m), foot vertical (m)]
L_streams = [rad2deg(L_hipFlex), rad2deg(L_kneeFlex), rad2deg(L_ankleDorsi), fwdSign * L_Foot(:,fwdAxis) - hipF, L_Foot(:,latAxis), L_Foot(:,3)];
R_streams = [rad2deg(R_hipFlex), rad2deg(R_kneeFlex), rad2deg(R_ankleDorsi), fwdSign * R_Foot(:,fwdAxis) - hipF, R_Foot(:,latAxis), R_Foot(:,3)];
fprintf('NaN in joint signals    : left %d | right %d\n', sum(isnan(L_streams(:,1:3)), 'all'), sum(isnan(R_streams(:,1:3)), 'all'));
if any(isnan(L_streams(:,1:3)), 'all') || any(isnan(R_streams(:,1:3)), 'all')
    error('NaN values in the joint-angle signals. Check the flat-foot pitch offset and the marker data.');
end
%% ============================================================
%  STEP 29: EXTRACT AND NORMALIZE GAIT CYCLES (0-100 %)
%  Cycle = heel strike to the next heel strike of the SAME foot.
% ============================================================
nPoints = 101;
x_new = linspace(0, 100, nPoints);
[L_Cyc, L_dur, L_hasUnreach] = buildCycles(hsL, L_streams, L_unreach, Fs, nPoints);
[R_Cyc, R_dur, R_hasUnreach] = buildCycles(hsR, R_streams, R_unreach, Fs, nPoints);
L_tStart = t(hsL(1:end-1));
L_tEnd   = t(hsL(2:end));
R_tStart = t(hsR(1:end-1));
R_tEnd   = t(hsR(2:end));
fprintf('\n========================================\n');
fprintf('       GAIT CYCLES\n');
fprintf('========================================\n');
fprintf('Left cycles  : %d (duration %.3f to %.3f s)\n', length(L_dur), min(L_dur), max(L_dur));
fprintf('Right cycles : %d (duration %.3f to %.3f s)\n', length(R_dur), min(R_dur), max(R_dur));
%% ============================================================
%  STEP 30: SELECT VALID CYCLES (with automatic relaxation)
% ============================================================
% Overstep times outside the trial cannot be seconds: try frames / ms
if ~isempty(overstepTimes) && max(overstepTimes) > t(end)
    if max(overstepTimes) / Fs <= t(end) * 1.01
        fprintf('NOTE: overstep times look like frames, converting to seconds.\n');
        overstepTimes = overstepTimes / Fs;
    elseif max(overstepTimes) / 1000 <= t(end) * 1.01
        fprintf('NOTE: overstep times look like milliseconds, converting to seconds.\n');
        overstepTimes = overstepTimes / 1000;
    else
        fprintf('WARNING: overstep times exceed the trial length, ignoring them.\n');
        overstepTimes = [];
    end
end
P = struct();
P.minDur = 0.7;                  % shortest accepted cycle (s)
P.maxDur = 1.8;                  % longest accepted cycle (s)
P.excludeOversteps = true;       % skip cycles that contain an overstep event
P.margin = 0.25;                 % time margin around each overstep event (s)
P.maxJump = 20;                  % largest joint step between 1 % points (deg)
P.maxUnreachFrac = 0.05;         % max fraction of unreachable samples per cycle
P.shapeSigma = 3.5;              % shape outlier limit (robust sigma)
P.minCyclesForShape = 5;
P.minValid = 3;                  % minimum cycles needed per side
levelNames = {'strict', 'overstep exclusion off', 'overstep off + loose jump/reach limits', 'duration check only'};
for level = 1:4
    Q = P;
    if level >= 2
        Q.excludeOversteps = false;
    end
    if level >= 3
        Q.maxJump = 45;
        Q.maxUnreachFrac = 0.25;
    end
    if level >= 4
        Q.maxJump = Inf;
        Q.maxUnreachFrac = 1;
        Q.shapeSigma = Inf;
    end
    [validL, infoL] = selectValidCycles(L_Cyc, L_dur, L_hasUnreach, L_tStart, L_tEnd, overstepTimes, Q);
    [validR, infoR] = selectValidCycles(R_Cyc, R_dur, R_hasUnreach, R_tStart, R_tEnd, overstepTimes, Q);
    if length(validL) >= P.minValid && length(validR) >= P.minValid
        break;
    end
end
fprintf('\n========================================\n');
fprintf('       CYCLE SELECTION\n');
fprintf('========================================\n');
fprintf('Selection level used    : %d (%s)\n', level, levelNames{level});
if level > 1
    fprintf('WARNING: strict selection left too few cycles, limits were relaxed.\n');
end
fprintf('                          Left   Right\n');
fprintf('Total cycles            : %4d   %4d\n', length(L_dur), length(R_dur));
fprintf('Bad duration            : %4d   %4d\n', infoL.badDuration, infoR.badDuration);
fprintf('Contain overstep event  : %4d   %4d\n', infoL.overstep, infoR.overstep);
fprintf('Too many unreachable    : %4d   %4d\n', infoL.unreachable, infoR.unreachable);
fprintf('Contains NaN            : %4d   %4d\n', infoL.nan, infoR.nan);
fprintf('Discontinuous (jump)    : %4d   %4d\n', infoL.jump, infoR.jump);
fprintf('Shape outlier           : %4d   %4d\n', infoL.shape, infoR.shape);
fprintf('VALID cycles            : %4d   %4d\n', length(validL), length(validR));
if isempty(validL) || isempty(validR)
    error('No valid gait cycles left. Cycle durations are outside %.1f-%.1f s: check heel-strike detection and the forward axis.', P.minDur, P.maxDur);
end
%% ============================================================
%  STEP 31: AVERAGE GAIT CYCLE (REFERENCE TRAJECTORY)
%  q_ref = (1/N) * sum(q_i)
% ============================================================
L_Sel = L_Cyc(validL,:,:);
R_Sel = R_Cyc(validR,:,:);
% Lateral foot position relative to the cycle mean
L_Sel(:,:,5) = L_Sel(:,:,5) - mean(L_Sel(:,:,5), 2);
R_Sel(:,:,5) = R_Sel(:,:,5) - mean(R_Sel(:,:,5), 2);
L_Ref = reshape(mean(L_Sel, 1), nPoints, 6);
R_Ref = reshape(mean(R_Sel, 1), nPoints, 6);
L_Std = reshape(std(L_Sel, 0, 1), nPoints, 6);
R_Std = reshape(std(R_Sel, 0, 1), nPoints, 6);
L_Hip_ref = L_Ref(:,1);
L_Knee_ref = L_Ref(:,2);
L_Ankle_ref = L_Ref(:,3);
R_Hip_ref = R_Ref(:,1);
R_Knee_ref = R_Ref(:,2);
R_Ankle_ref = R_Ref(:,3);
fprintf('\nAverage cycle computed from %d left and %d right valid cycles (%d points each)\n', length(validL), length(validR), nPoints);
%% ============================================================
%  STEP 32: LEFT-RIGHT PHASE CHECK
%  The opposite heel strike should occur near 50 % of the cycle.
% ============================================================
phaseL = contralateralPhase(hsL, hsR, validL);
phaseR = contralateralPhase(hsR, hsL, validR);
fprintf('\n========================================\n');
fprintf('       LEFT-RIGHT PHASE CHECK\n');
fprintf('========================================\n');
fprintf('Right heel strike inside left cycle  : %.1f %% of cycle (expected about 50)\n', median(phaseL, 'omitnan'));
fprintf('Left heel strike inside right cycle  : %.1f %% of cycle (expected about 50)\n', median(phaseR, 'omitnan'));
if abs(median(phaseL, 'omitnan') - 50) > 10 || abs(median(phaseR, 'omitnan') - 50) > 10
    fprintf('WARNING: left/right phase is not near 50 %%. Check the heel-strike detection.\n');
end
%% ============================================================
%  STEP 33: RANGE-OF-MOTION CHECK (PLAUSIBILITY)
%  Approximate limits for normal walking. These are checks only.
% ============================================================
fprintf('\n========================================\n');
fprintf('       JOINT RANGE CHECK (deg)\n');
fprintf('========================================\n');
romLimits = [-30 50; -5 80; -35 30];   % hip, knee, ankle
jointLabels = {'Hip', 'Knee', 'Ankle'};
for j = 1:3
    lo = min(L_Ref(:,j));
    hi = max(L_Ref(:,j));
    if lo >= romLimits(j,1) && hi <= romLimits(j,2)
        resultL = 'OK';
    else
        resultL = 'CHECK';
    end
    lo2 = min(R_Ref(:,j));
    hi2 = max(R_Ref(:,j));
    if lo2 >= romLimits(j,1) && hi2 <= romLimits(j,2)
        resultR = 'OK';
    else
        resultR = 'CHECK';
    end
    fprintf('%-6s left  : %7.2f to %7.2f  [%s]   right : %7.2f to %7.2f  [%s]   limit %d to %d\n', jointLabels{j}, lo, hi, resultL, lo2, hi2, resultR, romLimits(j,1), romLimits(j,2));
end
%% ============================================================
%  STEP 34: PLOT AVERAGE FOOT TRAJECTORY
% ============================================================
figure('Name', 'Average Foot Trajectory');
footNames = {'Forward position relative to hip (m)', 'Lateral position (m)', 'Vertical position (m)'};
for p = 1:3
    subplot(3,1,p);
    hold on;
    plotMeanStd(x_new, L_Ref(:,3+p), L_Std(:,3+p), [0 0.447 0.741]);
    plotMeanStd(x_new, R_Ref(:,3+p), R_Std(:,3+p), [0.850 0.325 0.098]);
    xlabel('Gait Cycle (%)');
    ylabel(footNames{p});
    title(['Average foot trajectory - ' footNames{p}]);
    if p == 1
        legend('Left \pm SD', 'Left', 'Right \pm SD', 'Right', 'Location', 'best');
    end
    grid on;
end
%% ============================================================
%  STEP 35: PLOT AVERAGE HIP, KNEE AND ANKLE TRAJECTORIES
% ============================================================
jointNames = {'Hip flexion (+)', 'Knee flexion (+)', 'Ankle dorsiflexion (+)'};
jointFiles = {'hip', 'knee', 'ankle'};
colL = [0 0.447 0.741];
colR = [0.850 0.325 0.098];
fprintf('\n========================================\n');
fprintf('       JOINT TRAJECTORY STATISTICS (deg)\n');
fprintf('========================================\n');
fprintf('%-7s %-6s %8s %8s %8s %12s %12s\n', 'Joint', 'Side', 'Min', 'Max', 'ROM', 'Min at (%)', 'Max at (%)');
for j = 1:3
    Lc = reshape(L_Sel(:,:,j), size(L_Sel,1), nPoints);
    Rc = reshape(R_Sel(:,:,j), size(R_Sel,1), nPoints);
    fig = figure('Name', ['Gait Trajectory - ' jointFiles{j}], 'Position', [100 100 1300 420]);
    % Left leg: individual cycles and mean +/- SD
    subplot(1,3,1);
    hold on;
    plot(x_new, Lc', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    plotMeanStd(x_new, L_Ref(:,j), L_Std(:,j), colL);
    xlim([0 100]);
    xlabel('Gait Cycle (%)');
    ylabel('Angle (deg)');
    title(['Left ' jointNames{j} ' (' num2str(size(Lc,1)) ' cycles)']);
    grid on;
    % Right leg: individual cycles and mean +/- SD
    subplot(1,3,2);
    hold on;
    plot(x_new, Rc', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    plotMeanStd(x_new, R_Ref(:,j), R_Std(:,j), colR);
    xlim([0 100]);
    xlabel('Gait Cycle (%)');
    ylabel('Angle (deg)');
    title(['Right ' jointNames{j} ' (' num2str(size(Rc,1)) ' cycles)']);
    grid on;
    % Left vs right comparison
    subplot(1,3,3);
    hold on;
    plotMeanStd(x_new, L_Ref(:,j), L_Std(:,j), colL);
    plotMeanStd(x_new, R_Ref(:,j), R_Std(:,j), colR);
    xlim([0 100]);
    xlabel('Gait Cycle (%)');
    ylabel('Angle (deg)');
    title('Left vs Right (mean \pm SD)');
    legend('Left \pm SD', 'Left', 'Right \pm SD', 'Right', 'Location', 'best');
    grid on;
    sgtitle([jointNames{j} ' - average gait cycle (model-based estimate)']);
    saveas(fig, ['gait_' jointFiles{j} '_trajectory.png']);
    [lmin, lminIdx] = min(L_Ref(:,j));
    [lmax, lmaxIdx] = max(L_Ref(:,j));
    [rmin, rminIdx] = min(R_Ref(:,j));
    [rmax, rmaxIdx] = max(R_Ref(:,j));
    fprintf('%-7s %-6s %8.2f %8.2f %8.2f %12.0f %12.0f\n', jointFiles{j}, 'Left', lmin, lmax, lmax - lmin, x_new(lminIdx), x_new(lmaxIdx));
    fprintf('%-7s %-6s %8.2f %8.2f %8.2f %12.0f %12.0f\n', jointFiles{j}, 'Right', rmin, rmax, rmax - rmin, x_new(rminIdx), x_new(rmaxIdx));
end
%% ============================================================
%  STEP 36: SAVE REFERENCE GAIT TRAJECTORY
% ============================================================
GaitReference = struct();
GaitReference.GaitCyclePercent = x_new';
GaitReference.NumValidCycles = [length(validL), length(validR)];
GaitReference.SelectionLevel = level;
GaitReference.Note = 'Model-based joint angles (inverse kinematics), not measured joint angles';
GaitReference.Left.Hip = L_Hip_ref;
GaitReference.Left.Knee = L_Knee_ref;
GaitReference.Left.Ankle = L_Ankle_ref;
GaitReference.Left.Std = L_Std;
GaitReference.Right.Hip = R_Hip_ref;
GaitReference.Right.Knee = R_Knee_ref;
GaitReference.Right.Ankle = R_Ankle_ref;
GaitReference.Right.Std = R_Std;
GaitReference.Model = struct('SubjectHeight', subjectHeight, 'Thigh', thighLength, 'Shank', shankLength, 'AnkleHeight', ankleHeight, 'AnkleFwdRatio', ankleFwdRatio, 'HipHeight', hipZ);
referenceTable = table(x_new', L_Hip_ref, L_Knee_ref, L_Ankle_ref, R_Hip_ref, R_Knee_ref, R_Ankle_ref, 'VariableNames', {'GaitCycle_percent', 'L_Hip_deg', 'L_Knee_deg', 'L_Ankle_deg', 'R_Hip_deg', 'R_Knee_deg', 'R_Ankle_deg'});
save('gait_reference.mat', 'GaitReference', 'referenceTable');
writetable(referenceTable, 'gait_reference_joint_angles.csv');
fprintf('\nSaved: gait_reference.mat\n');
fprintf('Saved: gait_reference_joint_angles.csv\n');
%% ============================================================
%  STEP 36.1: ANGULAR VELOCITY AND ACCELERATION
%  Differentiation amplifies noise, so each joint angle is
%  low-pass filtered again (zero phase) before differentiating.
%  Velocity and acceleration are computed against REAL time
%  (deg/s and deg/s^2) on the continuous signals, then cut into
%  the same valid heel-strike cycles and normalised to 0-100 %.
% ============================================================
fprintf('\n========================================\n');
fprintf('       ANGULAR VELOCITY AND ACCELERATION\n');
fprintf('========================================\n');
[bJ, aJ] = butter(2, Fc/(Fs/2), 'low');
L_dq = zeros(size(L_streams,1), 3);
R_dq = zeros(size(R_streams,1), 3);
L_ddq = zeros(size(L_streams,1), 3);
R_ddq = zeros(size(R_streams,1), 3);
for c = 1:3
    L_qs = filtfilt(bJ, aJ, L_streams(:,c));
    R_qs = filtfilt(bJ, aJ, R_streams(:,c));
    L_dq(:,c) = gradient(L_qs, 1/Fs);
    R_dq(:,c) = gradient(R_qs, 1/Fs);
    L_ddq(:,c) = gradient(L_dq(:,c), 1/Fs);
    R_ddq(:,c) = gradient(R_dq(:,c), 1/Fs);
end
% Cut into cycles and keep only the valid ones
[L_dqCyc, ~, ~] = buildCycles(hsL, L_dq, L_unreach, Fs, nPoints);
[R_dqCyc, ~, ~] = buildCycles(hsR, R_dq, R_unreach, Fs, nPoints);
[L_ddqCyc, ~, ~] = buildCycles(hsL, L_ddq, L_unreach, Fs, nPoints);
[R_ddqCyc, ~, ~] = buildCycles(hsR, R_ddq, R_unreach, Fs, nPoints);
L_qSel = L_Sel(:,:,1:3);
R_qSel = R_Sel(:,:,1:3);
L_dqSel = L_dqCyc(validL,:,:);
R_dqSel = R_dqCyc(validR,:,:);
L_ddqSel = L_ddqCyc(validL,:,:);
R_ddqSel = R_ddqCyc(validR,:,:);
% Average (reference) curves and standard deviations, 101 x 3
L_qRef = L_Ref(:,1:3);
R_qRef = R_Ref(:,1:3);
L_dqRef = reshape(mean(L_dqSel, 1), nPoints, 3);
R_dqRef = reshape(mean(R_dqSel, 1), nPoints, 3);
L_ddqRef = reshape(mean(L_ddqSel, 1), nPoints, 3);
R_ddqRef = reshape(mean(R_ddqSel, 1), nPoints, 3);
L_dqStd = reshape(std(L_dqSel, 0, 1), nPoints, 3);
R_dqStd = reshape(std(R_dqSel, 0, 1), nPoints, 3);
L_ddqStd = reshape(std(L_ddqSel, 0, 1), nPoints, 3);
R_ddqStd = reshape(std(R_ddqSel, 0, 1), nPoints, 3);
% Real-time axis of the average cycle
L_meanDur = mean(L_dur(validL));
R_meanDur = mean(R_dur(validR));
tCycL = linspace(0, L_meanDur, nPoints)';
tCycR = linspace(0, R_meanDur, nPoints)';
gaitFreqL = 1 / L_meanDur;
gaitFreqR = 1 / R_meanDur;
fprintf('Mean cycle duration : left %.3f s | right %.3f s\n', L_meanDur, R_meanDur);
fprintf('Gait frequency      : left %.3f Hz | right %.3f Hz\n', gaitFreqL, gaitFreqR);
%% ============================================================
%  STEP 36.2: MAX AND RMS OF ANGLE, VELOCITY AND ACCELERATION
%  "Ref"      : computed from the average gait cycle
%  "AllCycles": largest value found in any valid cycle
%  "CyclePeak": mean over the valid cycles of each cycle's peak
% ============================================================
jointNames3 = {'Hip', 'Knee', 'Ankle'};
sideNames = {'Left', 'Right'};
statRows = {};
for sd = 1:2
    if sd == 1
        dqR = L_dqRef; ddqR = L_ddqRef; dqS = L_dqSel; ddqS = L_ddqSel; nCyc = length(validL);
    else
        dqR = R_dqRef; ddqR = R_ddqRef; dqS = R_dqSel; ddqS = R_ddqSel; nCyc = length(validR);
    end
    for j = 1:3
        v = dqR(:,j);
        ac = ddqR(:,j);
        vAll = reshape(dqS(:,:,j), size(dqS,1), nPoints);
        aAll = reshape(ddqS(:,:,j), size(ddqS,1), nPoints);
        statRows(end+1,:) = {jointNames3{j}, sideNames{sd}, nCyc, max(v), min(v), max(abs(v)), sqrt(mean(v.^2)), max(abs(vAll), [], 'all'), mean(max(abs(vAll), [], 2)), max(ac), min(ac), max(abs(ac)), sqrt(mean(ac.^2)), max(abs(aAll), [], 'all'), mean(max(abs(aAll), [], 2))};
    end
end
statVars = {'Joint', 'Side', 'ValidCycles', 'Vel_Max_deg_s', 'Vel_Min_deg_s', 'Vel_MaxAbs_Ref_deg_s', 'Vel_RMS_Ref_deg_s', 'Vel_MaxAbs_AllCycles_deg_s', 'Vel_MeanCyclePeak_deg_s', 'Acc_Max_deg_s2', 'Acc_Min_deg_s2', 'Acc_MaxAbs_Ref_deg_s2', 'Acc_RMS_Ref_deg_s2', 'Acc_MaxAbs_AllCycles_deg_s2', 'Acc_MeanCyclePeak_deg_s2'};
motionStats = cell2table(statRows, 'VariableNames', statVars);
% Same table in radians (for the dynamics / actuator calculation)
motionStatsRad = motionStats;
motionStatsRad{:,4:end} = motionStats{:,4:end} * pi / 180;
radVars = strrep(strrep(statVars, 'deg_s2', 'rad_s2'), 'deg_s', 'rad_s');
motionStatsRad.Properties.VariableNames = radVars;
fprintf('\nVelocity and acceleration statistics (deg/s, deg/s^2):\n');
disp(motionStats);
if max(motionStats.Vel_MaxAbs_AllCycles_deg_s) > 800
    fprintf('WARNING: joint speed above 800 deg/s is unusual for walking. Check noise and spikes.\n');
end
%% ============================================================
%  STEP 36.3: COMPLETE 101-POINT GAIT-REFERENCE TABLES
%  Degrees / deg/s / deg/s^2 and radians / rad/s / rad/s^2
% ============================================================
refFull = table(x_new', tCycL, tCycR, 'VariableNames', {'GaitCycle_percent', 'Time_Left_s', 'Time_Right_s'});
refRad = refFull;
for sd = 1:2
    if sd == 1
        pre = 'L_'; qR = L_qRef; dqR = L_dqRef; ddqR = L_ddqRef;
    else
        pre = 'R_'; qR = R_qRef; dqR = R_dqRef; ddqR = R_ddqRef;
    end
    for j = 1:3
        jn = jointNames3{j};
        refFull.([pre jn '_q_deg']) = qR(:,j);
        refFull.([pre jn '_dq_deg_s']) = dqR(:,j);
        refFull.([pre jn '_ddq_deg_s2']) = ddqR(:,j);
        refRad.([pre jn '_q_rad']) = deg2rad(qR(:,j));
        refRad.([pre jn '_dq_rad_s']) = deg2rad(dqR(:,j));
        refRad.([pre jn '_ddq_rad_s2']) = deg2rad(ddqR(:,j));
    end
end
% Standard deviation across the valid cycles (degrees)
stdFull = table(x_new', 'VariableNames', {'GaitCycle_percent'});
for sd = 1:2
    if sd == 1
        pre = 'L_'; qS = L_Std(:,1:3); dqS2 = L_dqStd; ddqS2 = L_ddqStd;
    else
        pre = 'R_'; qS = R_Std(:,1:3); dqS2 = R_dqStd; ddqS2 = R_ddqStd;
    end
    for j = 1:3
        jn = jointNames3{j};
        stdFull.([pre jn '_q_SD_deg']) = qS(:,j);
        stdFull.([pre jn '_dq_SD_deg_s']) = dqS2(:,j);
        stdFull.([pre jn '_ddq_SD_deg_s2']) = ddqS2(:,j);
    end
end
% Joint range of motion and location of min / max
romRows = {};
for sd = 1:2
    if sd == 1
        qR = L_qRef;
    else
        qR = R_qRef;
    end
    for j = 1:3
        [qMin, iMin] = min(qR(:,j));
        [qMax, iMax] = max(qR(:,j));
        romRows(end+1,:) = {jointNames3{j}, sideNames{sd}, qMin, qMax, qMax - qMin, x_new(iMin), x_new(iMax)};
    end
end
romTable = cell2table(romRows, 'VariableNames', {'Joint', 'Side', 'Min_deg', 'Max_deg', 'ROM_deg', 'MinAt_percent', 'MaxAt_percent'});
% Model and gait parameters
paramRows = {'SubjectHeight', subjectHeight, 'm'; 'ThighLength', thighLength, 'm'; 'ShankLength', shankLength, 'm'; 'AnkleHeight', ankleHeight, 'm'; 'HipHeight', hipZ, 'm'; 'AnkleFwdRatio', ankleFwdRatio, '-'; 'SamplingFrequency', Fs, 'Hz'; 'PointsPerCycle', nPoints, '-'; 'ValidCycles_Left', length(validL), '-'; 'ValidCycles_Right', length(validR), '-'; 'MeanCycleDuration_Left', L_meanDur, 's'; 'MeanCycleDuration_Right', R_meanDur, 's'; 'MedianStrideTime', strideTime, 's'; 'GaitFrequency_Left', gaitFreqL, 'Hz'; 'GaitFrequency_Right', gaitFreqR, 'Hz'; 'ContralateralPhase_Left', median(phaseL, 'omitnan'), '% of cycle'; 'ContralateralPhase_Right', median(phaseR, 'omitnan'), '% of cycle'; 'SelectionLevel', level, '-'};
paramTable = cell2table(paramRows, 'VariableNames', {'Parameter', 'Value', 'Unit'});
%% ============================================================
%  STEP 36.4: PLOT ANGULAR VELOCITY AND ACCELERATION
% ============================================================
figure('Name', 'Angular Velocity and Acceleration', 'Position', [100 100 1200 800]);
for j = 1:3
    subplot(3,2,2*j-1);
    hold on;
    plotMeanStd(x_new, L_dqRef(:,j), L_dqStd(:,j), [0 0.447 0.741]);
    plotMeanStd(x_new, R_dqRef(:,j), R_dqStd(:,j), [0.850 0.325 0.098]);
    xlim([0 100]);
    xlabel('Gait Cycle (%)');
    ylabel('deg/s');
    title([jointNames3{j} ' angular velocity']);
    grid on;
    if j == 1
        legend('Left \pm SD', 'Left', 'Right \pm SD', 'Right', 'Location', 'best');
    end
    subplot(3,2,2*j);
    hold on;
    plotMeanStd(x_new, L_ddqRef(:,j), L_ddqStd(:,j), [0 0.447 0.741]);
    plotMeanStd(x_new, R_ddqRef(:,j), R_ddqStd(:,j), [0.850 0.325 0.098]);
    xlim([0 100]);
    xlabel('Gait Cycle (%)');
    ylabel('deg/s^2');
    title([jointNames3{j} ' angular acceleration']);
    grid on;
end
sgtitle('Average angular velocity and acceleration (model-based estimate)');
saveas(gcf, 'gait_velocity_acceleration.png');
%% ============================================================
%  STEP 36.5: EXPORT TO EXCEL, CSV AND MAT
% ============================================================
xlsxFile = 'gait_reference_full.xlsx';
if isfile(xlsxFile)
    delete(xlsxFile);
end
try
    writetable(paramTable, xlsxFile, 'Sheet', 'Model_Parameters');
    writetable(romTable, xlsxFile, 'Sheet', 'Joint_ROM');
    writetable(motionStats, xlsxFile, 'Sheet', 'Vel_Acc_Stats_deg');
    writetable(motionStatsRad, xlsxFile, 'Sheet', 'Vel_Acc_Stats_rad');
    writetable(refFull, xlsxFile, 'Sheet', 'Reference_101pts_deg');
    writetable(refRad, xlsxFile, 'Sheet', 'Reference_101pts_rad');
    writetable(stdFull, xlsxFile, 'Sheet', 'SD_101pts_deg');
    fprintf('\nSaved: %s (7 sheets)\n', xlsxFile);
catch excelError
    fprintf('WARNING: Excel export failed (%s). Close the file if it is open in Excel. CSV files are still saved.\n', excelError.message);
end
writetable(refFull, 'gait_reference_101pts_deg.csv');
writetable(refRad, 'gait_reference_101pts_rad.csv');
writetable(stdFull, 'gait_reference_101pts_SD.csv');
writetable(motionStats, 'gait_motion_statistics_deg.csv');
writetable(motionStatsRad, 'gait_motion_statistics_rad.csv');
writetable(romTable, 'gait_joint_ROM.csv');
writetable(paramTable, 'gait_model_parameters.csv');
GaitFull = struct();
GaitFull.Note = 'Model-based joint angles (inverse kinematics), not measured joint angles';
GaitFull.GaitCyclePercent = x_new';
GaitFull.TimeLeft = tCycL;
GaitFull.TimeRight = tCycR;
GaitFull.Left = struct('q', L_qRef, 'dq', L_dqRef, 'ddq', L_ddqRef, 'qStd', L_Std(:,1:3), 'dqStd', L_dqStd, 'ddqStd', L_ddqStd);
GaitFull.Right = struct('q', R_qRef, 'dq', R_dqRef, 'ddq', R_ddqRef, 'qStd', R_Std(:,1:3), 'dqStd', R_dqStd, 'ddqStd', R_ddqStd);
GaitFull.Units = 'q [deg], dq [deg/s], ddq [deg/s^2]; columns = hip, knee, ankle';
GaitFull.MotionStatistics = motionStats;
GaitFull.MotionStatisticsRad = motionStatsRad;
GaitFull.ROM = romTable;
GaitFull.Parameters = paramTable;
GaitFull.ReferenceDeg = refFull;
GaitFull.ReferenceRad = refRad;
save('gait_reference_full.mat', 'GaitFull');
fprintf('Saved: gait_reference_full.mat\n');
fprintf('Saved: gait_reference_101pts_deg.csv, gait_reference_101pts_rad.csv, gait_reference_101pts_SD.csv\n');
fprintf('Saved: gait_motion_statistics_deg.csv, gait_motion_statistics_rad.csv\n');
fprintf('Saved: gait_joint_ROM.csv, gait_model_parameters.csv, gait_velocity_acceleration.png\n');
%% ============================================================
%  STEP 37: FINAL SUMMARY
% ============================================================
fprintf('\n========================================\n');
fprintf('       HUMAN GAIT ANALYSIS SUMMARY\n');
fprintf('========================================\n');
fprintf('Trial                  : %s\n', Trials(trialID).Name);
fprintf('Samples                : %d\n', length(t));
fprintf('Duration               : %.3f seconds\n', t(end));
fprintf('Sampling frequency     : %d Hz\n', Fs);
fprintf('Missing values before  : %d\n', missingBefore);
fprintf('Missing values after   : %d\n', missingAfter);
fprintf('Spikes corrected       : %d (%.3f %%)\n', outlierCount, outlierPercentage);
fprintf('Valid cycles (L / R)   : %d / %d\n', length(validL), length(validR));
fprintf('Selection level        : %d (%s)\n', level, levelNames{level});
fprintf('Points per cycle       : %d\n', nPoints);
fprintf('\nLeft  hip / knee / ankle range (deg)  : %.1f to %.1f / %.1f to %.1f / %.1f to %.1f\n', min(L_Hip_ref), max(L_Hip_ref), min(L_Knee_ref), max(L_Knee_ref), min(L_Ankle_ref), max(L_Ankle_ref));
fprintf('Right hip / knee / ankle range (deg)  : %.1f to %.1f / %.1f to %.1f / %.1f to %.1f\n', min(R_Hip_ref), max(R_Hip_ref), min(R_Knee_ref), max(R_Knee_ref), min(R_Ankle_ref), max(R_Ankle_ref));
fprintf('\nNOTE: joint angles are model-based estimates (inverse kinematics).\n');
fprintf('Verify the results before sending them to a controller.\n');
fprintf('\n========================================\n');
fprintf('       ANALYSIS COMPLETED\n');
fprintf('========================================\n');
%% ============================================================
%  LOCAL FUNCTIONS (must stay at the end of the script)
% ============================================================
function mask = removeLongRuns(mask, maxLen)
% Keep only short runs of true values (isolated spikes)
    mask = mask(:);
    d = diff([0; mask; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    for k = 1:length(starts)
        if (ends(k) - starts(k) + 1) > maxLen
            mask(starts(k):ends(k)) = false;
        end
    end
end
function v = pctl(x, p)
% Percentile without the Statistics Toolbox
    x = sort(x(~isnan(x)));
    n = length(x);
    idx = min(max(round(p / 100 * n), 1), n);
    v = x(idx);
end
function hs = detectHeelStrikes(relHeel, Fs)
% Heel strike = peak of the heel position relative to the body
    span = pctl(relHeel, 95) - pctl(relHeel, 5);
    [~, hs] = findpeaks(relHeel, 'MinPeakDistance', round(0.6 * Fs), 'MinPeakProminence', 0.4 * span);
    hs = hs(:);
end
function [fAnkle, zAnkle] = ankleFromFoot(fHeel, zHeel, fFore, zFore, fwdRatio, ankleH)
% Ankle joint centre in the foot frame: along the heel-forefoot
% line plus a perpendicular offset (rotates with the foot)
    df = fFore - fHeel;
    dz = zFore - zHeel;
    len = sqrt(df.^2 + dz.^2);
    len(len < eps) = eps;
    uF = df ./ len;
    uZ = dz ./ len;
    % Perpendicular pointing up
    nF = -uZ;
    nZ = uF;
    fAnkle = fHeel + fwdRatio * len .* uF + ankleH * nF;
    zAnkle = zHeel + fwdRatio * len .* uZ + ankleH * nZ;
end
function [hipFlex, kneeFlex, shankAngle, unreachable] = legIK(fAnkle, zAnkle, hipF, hipZ, L1, L2)
% Planar two-link inverse kinematics (radians)
%   hipFlex    : thigh angle from vertical, + forward
%   kneeFlex   : 0 = straight leg, + = flexed
%   shankAngle : knee-to-ankle direction from vertical, + forward
    dx = fAnkle - hipF;
    dz = zAnkle - hipZ;
    D = sqrt(dx.^2 + dz.^2);
    % 3 % tolerance: small overshoots are clamped below, not rejected
    unreachable = D > 1.03 * (L1 + L2) | D < abs(L1 - L2);
    D = min(max(D, abs(L1 - L2) + 1e-6), L1 + L2 - 1e-6);
    kneeInterior = acos(max(min((L1^2 + L2^2 - D.^2) / (2 * L1 * L2), 1), -1));
    kneeFlex = pi - kneeInterior;
    gam = acos(max(min((L1^2 + D.^2 - L2^2) ./ (2 * L1 * D), 1), -1));
    alpha = atan2(dx, -dz);
    hipFlex = alpha + gam;
    shankAngle = hipFlex - kneeFlex;
end
function [cyc, dur, hasFlag] = buildCycles(hs, streams, flagStream, Fs, nPoints)
% Cut heel-strike to heel-strike cycles and resample to nPoints.
% hasFlag = fraction of flagged (unreachable) samples in each cycle
    nC = length(hs) - 1;
    nS = size(streams, 2);
    cyc = nan(nC, nPoints, nS);
    dur = zeros(nC, 1);
    hasFlag = zeros(nC, 1);
    xNew = linspace(0, 100, nPoints);
    for k = 1:nC
        s = hs(k);
        e = hs(k+1);
        seg = streams(s:e, :);
        xOld = linspace(0, 100, size(seg,1));
        cyc(k,:,:) = reshape(interp1(xOld, seg, xNew), 1, nPoints, nS);
        dur(k) = (e - s) / Fs;
        hasFlag(k) = mean(double(flagStream(s:e)));
    end
end
function [validIdx, info] = selectValidCycles(cyc, dur, hasUnreach, tStart, tEnd, overstepTimes, P)
% Reject cycles by duration, overstep events, unreachable samples,
% discontinuities and abnormal shape
    nC = size(cyc, 1);
    okDur = dur >= P.minDur & dur <= P.maxDur;
    okOver = true(nC, 1);
    if P.excludeOversteps && ~isempty(overstepTimes)
        for k = 1:nC
            okOver(k) = ~any(overstepTimes >= tStart(k) - P.margin & overstepTimes <= tEnd(k) + P.margin);
        end
    end
    okReach = hasUnreach <= P.maxUnreachFrac;
    joints = cyc(:,:,1:3);
    okNaN = true(nC, 1);
    okJump = true(nC, 1);
    for k = 1:nC
        jk = reshape(joints(k,:,:), size(joints,2), 3);
        okNaN(k) = ~any(isnan(jk), 'all');
        d = abs(diff(jk, 1, 1));
        okJump(k) = ~any(d(:) > P.maxJump);
    end
    base = okDur & okOver & okReach & okNaN & okJump;
    okShape = true(nC, 1);
    idx = find(base);
    if ~isinf(P.shapeSigma) && length(idx) >= P.minCyclesForShape
        J = reshape(joints(idx,:,:), length(idx), []);
        medianCurve = median(J, 1);
        dist = sqrt(mean((J - medianCurve).^2, 2));
        madValue = median(abs(dist - median(dist)));
        limit = median(dist) + P.shapeSigma * 1.4826 * madValue;
        okShape(idx) = dist <= max(limit, eps);
    end
    validIdx = find(base & okShape);
    info.badDuration = sum(~okDur);
    info.overstep = sum(~okOver);
    info.unreachable = sum(~okReach);
    info.nan = sum(~okNaN);
    info.jump = sum(~okJump);
    info.shape = sum(~okShape);
end
function pct = contralateralPhase(hsOwn, hsOther, validIdx)
% Position (% of cycle) of the opposite heel strike inside each valid cycle
    pct = nan(length(validIdx), 1);
    for k = 1:length(validIdx)
        s = hsOwn(validIdx(k));
        e = hsOwn(validIdx(k) + 1);
        h = hsOther(hsOther > s & hsOther < e);
        if ~isempty(h)
            pct(k) = 100 * (h(1) - s) / (e - s);
        end
    end
end
function plotMeanStd(x, meanValue, stdValue, colorValue)
% Mean curve with a shaded +/- 1 SD band
    x = x(:)';
    meanValue = meanValue(:)';
    stdValue = stdValue(:)';
    fill([x fliplr(x)], [meanValue - stdValue fliplr(meanValue + stdValue)], colorValue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(x, meanValue, 'Color', colorValue, 'LineWidth', 2);
end
function offset = flatFootPitch(pitch, zHeel, zFore)
% Median foot pitch while the foot is flat on the ground.
% The foot is taken as flat when its mean height is within 2 cm of its
% lowest height. If that gives too few samples, the lowest 25 % is used.
% The result is never NaN (falls back to 0 with a warning).
    zFoot = 0.5 * (zHeel + zFore);
    flat = zFoot <= pctl(zFoot, 5) + 0.02;
    if sum(flat) < 10
        flat = zFoot <= pctl(zFoot, 25);
    end
    offset = median(pitch(flat), 'omitnan');
    if isnan(offset)
        fprintf('WARNING: flat-foot pitch could not be estimated, offset set to 0.\n');
        offset = 0;
    end
end