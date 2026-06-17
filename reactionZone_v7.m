%% reactionZoneV7.m
% Wavefront detection, alignment, and reaction zone structure analysis
% Author: Raj Punna
% Date: June 2026
% Requires: Image Processing Toolbox (medfilt2, imgaussfilt, rgb2gray, graythresh, otsuthresh, bwareaopen)
%
% WORKFLOW:
%   1. Select folder containing BMP sequence (frame 001 = background)
%   2. Background subtract and normalize each frame
%   3. Segment structural pixels via folded Otsu (abs intensity → Otsu); optional speckle removal
%   4. Detect wavefront as rightmost segmented pixel per row + smoothing
%   5. Align images to detected front (front near right edge to maximize burnt gas coverage)
%   6. Track per-pixel coverage to correct for alignment padding
%   7. Detect test boundaries (front jumps leftward between sequential frames)
%   8. Filter frames by front position (exclude FOV boundary frames)
%   9. Build coverage-corrected probability heatmap using randomized frame accumulation
%  10. Track convergence of P(x = 0.05λ) — probability at fixed cell-width distance
%  11. Build the per-row reaction-zone width distribution (front -> last structural pixel)
%
% V6 CHANGES (from V5) — back-ported from the Python port (detrz-extraction):
%   - BUG FIX (structure threshold space): V5 computed the folded-Otsu threshold on the
%     ENHANCED image but applied it to the NORMALIZED (aligned) image when building the
%     structure mask — a mismatch, since the unsharp mask shifts the intensity scale.
%     The threshold is now computed and applied in the SAME space, controlled by
%     params.segmentationSource ('normalized' default, matching the port). Front detection
%     keeps using the enhanced image (sharper edges). This changes structural probabilities.
%   - BUG FIX (convergence zero-fill bias): the per-step accumulation heatmap in
%     computeConvergence was initialised to 0, so 'omitnan' had nothing to omit and every
%     below-coverage-floor pixel pulled the row-mean toward 0. Now initialised to NaN.
%     The FINAL heatmap/probabilities are unchanged (the final pass already used NaN); only
%     the convergence-vs-frames curve is corrected (lower-biased early in accumulation).
%   - BUG FIX (frontDataRaw): V5 stored the SMOOTHED/outlier-rejected front in frontDataRaw.
%     Front detection is now split into extractFrontRaw (rightmost pixel, untouched) and
%     refineFront (outlier rejection + smoothing); frontDataRaw is now genuinely raw.
%     The refined frontData is byte-identical to V5.
%   - NEW ANALYSIS (reaction-zone width distribution): per-row, per-frame distance from the
%     front to the LAST (leftmost) structural pixel, with FOV censoring. Raw per-row columns
%     are stored so the censor margin can be changed at plot time without re-extraction.
%
%   - NEW (pooled threshold mode): params.thresholdMode = 'pooled' computes ONE folded-Otsu
%     threshold from the whole stack (matching the Python port) instead of one per frame,
%     removing frame-to-frame threshold drift. Implemented as a pre-pass that caches the
%     normalized stack (single precision) and pools its abs-intensity histogram, so the main
%     loop still reads each BMP only once — the cost is holding the normalized stack in
%     memory. Default remains 'per_frame' to preserve existing results; the port defaults
%     to 'pooled'. The threshold is computed in params.segmentationSource space either way.
%
% V7 CHANGES (from V6):
%   - NEW (optional structure-mask denoising): params.structure.denoise applies bwareaopen
%     to the binary structure mask (after coverage masking, before any geometry is read) to
%     strip isolated-pixel speckle that the folded-Otsu threshold admits when the schlieren
%     signal is faint. Because the reaction-zone width is the distance to the LAST (leftmost)
%     structural pixel — an extreme-value statistic — a single detached background pixel sets
%     a row's width, so speckle inflates the trailing edge. Off by default (existing results
%     unchanged). When on, it cleans the single mask used for the width distribution, the
%     probability heatmap, and the Structure overlays, and stores per-frame removed/kept
%     pixel counts for QC (the removed fraction quantifies contamination per dataset).
%     Orthogonal to thresholdMode: the threshold is computed before denoising is applied.
%
% OUTPUT:
%   [analysisName]_Analysis/
%   ├── Fronts/          <- overlay images with detected front
%   ├── Aligned/         <- shifted/aligned images
%   ├── Structure/       <- aligned images with structural mask overlay
%   ├── Stats/
%   │   ├── Heatmap.png
%   │   ├── ProbabilityCurve.png
%   │   ├── Convergence.png
%   │   ├── OtsuStability.png
%   │   ├── CoverageMap.png
%   │   └── ReactionZoneWidth.png
%   └── FrontData.mat    <- all detection and structure data

%% FUNCTION REGISTRY
%
% Internal (defined at end of this script):
%   loadOrCreateTestEntry   - Load test params from database or prompt user for new entry
%   loadAndNormalize        - Read BMP, subtract background, zero-median normalize
%   preprocessForDetection  - Median filter + unsharp mask to improve edge contrast
%   threshFoldedOtsu        - Otsu threshold on abs(intensity) distribution
%   threshFoldedOtsuPooledStack - Pooled folded-Otsu over a whole normalized stack
%   segSourceImage          - Return the structure-segmentation image in the requested space
%   extractFrontRaw         - Rightmost segmented pixel per row (raw, no smoothing)
%   refineFront             - Outlier rejection + median smoothing of a raw front (NaN-aware)
%   alignToFront            - Shift image rows so detected front aligns to a target column
%   buildCoverageMask       - Binary mask of pixels with real (non-padded) data after alignment
%   detectTestBoundaries    - Flag frames where the front jumps backward (new test starts)
%   exportFrameOverlay      - Save front overlay, aligned image, and structure overlay to disk
%   computeConvergence      - Accumulate masks with coverage correction, track P at fixed distances
%   interpProbAtDistance    - Interpolate probability at a specific distance from the curve
%   computeReactionZoneWidths - Per-row front-to-last-structural-pixel widths, with FOV censoring
%   prctileLocal            - Percentile via (i-0.5)/n convention (no Statistics Toolbox)
%   drawFrontLine           - Overlay detected front as colored line on RGB image
%   applyMaskOverlay        - Paint binary mask onto RGB image in a given color
%   toDisplayUint8          - Convert zero-centered double image to uint8 for export

clear; clc; close all;

%% ========== USER PARAMETERS ==========

% --- Alignment ---
% Align front near the right edge to maximize burnt gas (left side) coverage.
% Buffer is the number of pixels between the front and the right image edge.
params.alignBuffer = 10;

% --- Preprocessing ---
params.preprocess.medfiltSize = 3;      % Median filter kernel (odd integer)
params.preprocess.unsharpRadius = 2.5;  % Gaussian sigma for unsharp mask
params.preprocess.unsharpAmount = 0.8;  % Unsharp mask gain

% --- Front Estimation (from segmentation mask) ---
params.front.smoothWindow = 7;          % Median filter width for front smoothing
params.front.outlierThreshPx = 20;      % Deviation from local median to flag as outlier

% --- Structure Detection ---
params.structure.overlayColor = [255, 0, 0];

% Space in which the structure mask is segmented. The folded-Otsu threshold is
% computed and applied in the SAME space (V6 fix). 'normalized' matches the
% Python port default; 'enhanced' reuses the (sharpened) front-detection space.
params.segmentationSource = 'normalized';   % 'normalized' | 'enhanced'

% How the structure-mask threshold is computed. 'per_frame' takes a fresh
% folded-Otsu threshold from every frame (original V6 behaviour); 'pooled'
% computes ONE threshold from the whole stack (matching the Python port),
% which removes frame-to-frame threshold drift but holds the normalized stack
% in memory during a short pre-pass. Both compute it in segmentationSource space.
params.thresholdMode = 'pooled';   % 'per_frame' | 'pooled'

% --- Structure mask denoising (optional) ---
% The reaction-zone width is the distance to the LAST (leftmost) structural
% pixel per row -- an extreme-value statistic, so a single isolated background
% pixel sets that row's width. When the schlieren structure is faint and the
% folded-Otsu threshold drops low enough to admit speckle (e.g. weakly-radiating
% regular mixtures), that speckle inflates the trailing edge. bwareaopen removes
% connected components smaller than minPixels before any geometry is measured.
%
% Off by default (preserves existing results). When on, it cleans the structure
% mask everywhere it is used: width distribution, probability heatmap, and the
% exported Structure overlays -- so the overlays show exactly what was removed.
% Tune per dataset by inspecting the Structure overlays: raise minPixels until
% the detached speckle disappears but genuine trailing structure remains.
% CAUTION: in highly irregular cases, small DETACHED unreacted pockets are real
% signal; too large a minPixels deletes them and biases widths short. Keep it
% small (or off) for clean datasets, larger only for the noisy faint ones.
%
% connectivity 8 is conservative (keeps thin diagonal filaments as one
% component) and identical to 4 for truly isolated pixels; 4 is more aggressive
% (fragments and removes thin diagonal features). 8 recommended.
params.structure.denoise.enable       = true;  % true to apply bwareaopen
params.structure.denoise.minPixels    = 2;      % min connected-component size to keep
params.structure.denoise.connectivity = 8;      % 4 | 8 (8 keeps thin diagonal filaments)

% --- Reaction Zone Width ---
% Per-row distance from the front to the last (leftmost) structural pixel.
params.reactionZoneWidth.censorMarginPx = 2;          % rows whose trailing edge sits within
                                                      % this many px of the real-data boundary
                                                      % are censored (may be truncated by FOV)
params.reactionZoneWidth.percentiles    = [50, 90, 98];
params.reactionZoneWidth.bins           = 60;

% --- Test Boundary Detection ---
params.boundary.jumpFraction = 0.3;

% --- Convergence ---
params.convergence.rngSeed = 42;

% --- Probability Probe Distances ---
params.probeDistances.cellFractions = [0.02, 0.05, 0.10];

% --- Frame Filtering for Accumulation ---
params.frameFilter.frontRangePercent = [0.2, 0.95];

% --- Minimum coverage ---
% Pixels with fewer than this many frames of real data are masked out of the
% final heatmap to avoid noisy estimates from low-coverage regions.
params.minCoverage = 5;

% --- Output ---
params.outputEnhanced = false;
params.overlayColor = [255, 0, 0];
params.overlayLinewidth = 1;
params.outputFormat = 'png';

% --- Paths ---
DATABASE_PATH = '/Users/19rwrp/Desktop/Project/Thesis/Reaction Zone Analysis/ReactionZone_Database.mat';
OUTPUT_ROOT   = '/Users/19rwrp/Desktop/Project/Thesis/Reaction Zone Analysis/Structural Analysis/V7_denoise';

%% ========== SETUP ==========

fprintf('Select folder containing BMP sequence...\n');
inputFolder = uigetdir(pwd, 'Select folder with BMP images');
if inputFolder == 0
    error('No folder selected. Exiting.');
end

[~, params.analysisName] = fileparts(inputFolder);
fprintf('Analysis name (from folder): %s\n', params.analysisName);

% --- Load or create test database entry ---
[testEntry, db] = loadOrCreateTestEntry(DATABASE_PATH, params.analysisName);

params.mixture          = string(testEntry.mixture(1));
params.pressure         = testEntry.pressure(1);
params.px2mm            = testEntry.px2mm(1);
params.inductionLengthM = testEntry.inductionLengthM(1);
params.Ucj              = testEntry.Ucj(1);
params.Ea               = testEntry.Ea(1);
params.cellWidth        = testEntry.cellWidth(1);

params.inductionLengthMm = params.inductionLengthM * 1000;
params.px2induction = params.px2mm / params.inductionLengthMm;

% Convert probe distances from cell fractions to induction lengths
probeDistances_mm = params.probeDistances.cellFractions * params.cellWidth;
probeDistances_di = probeDistances_mm / params.inductionLengthMm;
nProbes = length(params.probeDistances.cellFractions);

fprintf('\nTest Parameters:\n');
fprintf('  Mixture:          %s\n', params.mixture);
fprintf('  Pressure:         %.1f kPa\n', params.pressure);
fprintf('  U_CJ:             %.0f m/s\n', params.Ucj);
fprintf('  Induction length: %.3e m (%.4f mm)\n', params.inductionLengthM, params.inductionLengthMm);
fprintf('  Ea/RT:            %.2f\n', params.Ea);
fprintf('  Cell width:       %.2f mm\n', params.cellWidth);
fprintf('  px2mm:            %.4f mm/px\n', params.px2mm);
fprintf('  px2induction:     %.4f induction lengths/px\n', params.px2induction);

fprintf('\nProbe distances:\n');
for p = 1:nProbes
    fprintf('  %.0f%% cell width = %.4f mm = %.2f δ_i\n', ...
        params.probeDistances.cellFractions(p)*100, probeDistances_mm(p), probeDistances_di(p));
end
fprintf('\n');

% --- Output folder structure ---
analysisFolder = fullfile(OUTPUT_ROOT, [params.analysisName '_Analysis']);
folders.fronts    = fullfile(analysisFolder, 'Fronts');
folders.aligned   = fullfile(analysisFolder, 'Aligned');
folders.structure = fullfile(analysisFolder, 'Structure');
folders.stats     = fullfile(analysisFolder, 'Stats');

folderList = {analysisFolder, folders.fronts, folders.aligned, folders.structure, folders.stats};
for i = 1:length(folderList)
    if ~exist(folderList{i}, 'dir'), mkdir(folderList{i}); end
end

fprintf('Output folder: %s\n', analysisFolder);

% --- Archive script ---
scriptFullPath = [mfilename('fullpath') '.m'];
[~, scriptName, ~] = fileparts(scriptFullPath);
dateStr = datestr(now, 'yyyy-mm-dd');
scriptCopy = fullfile(analysisFolder, sprintf('%s_%s.m', scriptName, dateStr));
copyfile(scriptFullPath, scriptCopy);
fprintf('Script archived: %s\n', scriptCopy);

% --- File discovery ---
bmpFiles = dir(fullfile(inputFolder, '*.bmp'));
if isempty(bmpFiles)
    bmpFiles = dir(fullfile(inputFolder, '*.BMP'));
end
if isempty(bmpFiles)
    error('No BMP files found in selected folder.');
end

nFiles = length(bmpFiles);
fprintf('Found %d BMP files.\n', nFiles);

[~, sortIdx] = sort({bmpFiles.name});
bmpFiles = bmpFiles(sortIdx);

frameNums = zeros(nFiles, 1);
for i = 1:nFiles
    [~, fname, ~] = fileparts(bmpFiles(i).name);
    tokens = regexp(fname, '(\d+)$', 'tokens');
    if ~isempty(tokens)
        frameNums(i) = str2double(tokens{1}{1});
    else
        frameNums(i) = i;
    end
end

bgIdx = find(frameNums == 1, 1);
if isempty(bgIdx)
    warning('Frame 001 not found. Using first file as background.');
    bgIdx = 1;
end

% --- Load background ---
bgPath = fullfile(bmpFiles(bgIdx).folder, bmpFiles(bgIdx).name);
bgImage = double(imread(bgPath));
if ndims(bgImage) == 3
    bgImage = double(rgb2gray(uint8(bgImage)));
end

[imgHeight, imgWidth] = size(bgImage);

% Front aligns near right edge to maximize burnt gas coverage
targetCol = imgWidth - params.alignBuffer;

fprintf('Background: %s  |  Size: %d x %d  |  Target col: %d (buffer=%d)\n', ...
    bmpFiles(bgIdx).name, imgHeight, imgWidth, targetCol, params.alignBuffer);

% --- Initialize data storage ---
processIdx = setdiff(1:nFiles, bgIdx);
nFrames = length(processIdx);

frontData      = NaN(imgHeight, nFrames);
frontDataRaw   = NaN(imgHeight, nFrames);
frameNumbers   = zeros(nFrames, 1);
frameNames     = cell(nFrames, 1);
otsuThresholds = zeros(nFrames, 1);   % structure-mask threshold (drives the statistics)
frontThresholds = zeros(nFrames, 1);  % front-detection threshold (enhanced space)
medianFrontCol = zeros(nFrames, 1);

% Optional speckle-filter bookkeeping (zeros if denoise disabled)
denoiseRemoved = zeros(nFrames, 1);  % structural px removed by speckle filter, per frame
denoiseKept    = zeros(nFrames, 1);  % structural px kept after speckle filter, per frame

% Per-row reaction-zone geometry (raw inputs to the width distribution; stored so
% the censor margin is a plot-time choice that needs no re-extraction)
lastStructureCol = NaN(imgHeight, nFrames);  % leftmost structural column (trailing edge)
coverageLeftCol  = NaN(imgHeight, nFrames);  % leftmost real-data column after alignment

% Store aligned structure masks and coverage masks for convergence
alignedMasks  = false(imgHeight, imgWidth, nFrames);
coverageMasks = false(imgHeight, imgWidth, nFrames);

% Distance axis: positive = behind front (burnt gas side), referenced from targetCol
distanceAxisPx = targetCol - (1:imgWidth);
distanceAxis   = distanceAxisPx * params.px2induction;

%% ========== POOLED THRESHOLD PRE-PASS (pooled mode only) ==========

usePooled = strcmpi(params.thresholdMode, 'pooled');
if ~usePooled && ~strcmpi(params.thresholdMode, 'per_frame')
    error('Unknown params.thresholdMode: %s (use ''per_frame'' or ''pooled'')', ...
        params.thresholdMode);
end

if usePooled
    fprintf('\nPooled threshold mode: caching normalized stack (single precision)...\n');
    % Cache the normalized stack so the main loop reads each BMP only once. Single
    % precision halves the footprint; the values are integer pixel differences plus
    % a scalar median shift, so rounding has no effect on fronts or the threshold.
    normalizedStack = zeros(imgHeight, imgWidth, nFrames, 'single');
    for f = 1:nFrames
        idx = processIdx(f);
        framePath = fullfile(bmpFiles(idx).folder, bmpFiles(idx).name);
        normalizedStack(:, :, f) = single(loadAndNormalize(framePath, bgImage));
    end
    pooledThresh = threshFoldedOtsuPooledStack(normalizedStack, ...
        params.segmentationSource, params.preprocess);
    fprintf('Pooled folded-Otsu threshold = %.2f  (segmentationSource=%s)\n', ...
        pooledThresh, params.segmentationSource);
else
    normalizedStack = [];   % unused in per-frame mode
    pooledThresh    = NaN;  % unused in per-frame mode
end

%% ========== MAIN PROCESSING LOOP ==========

fprintf('\nProcessing %d frames...\n', nFrames);
tic;

for f = 1:nFrames
    idx = processIdx(f);
    framePath = fullfile(bmpFiles(idx).folder, bmpFiles(idx).name);
    [~, fname, ~] = fileparts(bmpFiles(idx).name);
    frameNumbers(f) = frameNums(idx);
    frameNames{f} = fname;

    % --- Detect front (on the sharpened image: better edge localization) ---
    if usePooled
        normalized = double(normalizedStack(:, :, f));   % from the pre-pass cache
    else
        normalized = loadAndNormalize(framePath, bgImage);
    end
    enhanced   = preprocessForDetection(normalized, params.preprocess);

    frontThresh        = threshFoldedOtsu(enhanced);
    frontMask          = abs(enhanced) > frontThresh;
    frontRaw           = extractFrontRaw(frontMask);     % rightmost pixel, untouched
    frontDataRaw(:, f) = frontRaw;                       % V6: now genuinely raw
    frontX             = refineFront(frontRaw, params.front);
    frontData(:, f)    = frontX;
    medianFrontCol(f)  = median(frontX, 'omitnan');
    frontThresholds(f) = frontThresh;

    % --- Align and track coverage ---
    aligned = alignToFront(normalized, frontX, 0, targetCol);
    covMask = buildCoverageMask(frontX, imgHeight, imgWidth, targetCol);

    % --- Structure mask (V6 fix: threshold computed AND applied in the same space) ---
    switch params.segmentationSource
        case 'normalized'
            structImg = aligned;                                      % aligned normalized
        case 'enhanced'
            structImg = alignToFront(enhanced, frontX, 0, targetCol); % aligned enhanced
        otherwise
            error('Unknown params.segmentationSource: %s (use ''normalized'' or ''enhanced'')', ...
                params.segmentationSource);
    end

    if usePooled
        structThresh = pooledThresh;                                 % one threshold for all frames
    elseif strcmp(params.segmentationSource, 'enhanced')
        structThresh = frontThresh;                                  % enhanced-space, per frame
    else
        structThresh = threshFoldedOtsu(normalized);                 % normalized-space, per frame
    end
    otsuThresholds(f) = structThresh;
    alignedMask = abs(structImg) > structThresh;

    % Drop structure in padded regions (belt-and-suspenders with coverage)
    alignedMask = alignedMask & covMask;

    % Optional speckle removal (isolated-pixel denoising). Applied here so the
    % cleaned mask is the single source of truth for width, heatmap, and overlays.
    % Orthogonal to thresholdMode (threshold is already fixed before this point).
    if params.structure.denoise.enable
        nBefore     = nnz(alignedMask);
        alignedMask = bwareaopen(alignedMask, params.structure.denoise.minPixels, ...
            params.structure.denoise.connectivity);
        denoiseKept(f)    = nnz(alignedMask);
        denoiseRemoved(f) = nBefore - denoiseKept(f);
    end

    alignedMasks(:, :, f)  = alignedMask;
    coverageMasks(:, :, f) = covMask;

    % --- Per-row reaction-zone geometry (max finds the first/leftmost true column) ---
    [anyStruct, firstStructCol] = max(alignedMask, [], 2);
    rowsStruct = anyStruct > 0;
    lastStructureCol(rowsStruct, f) = firstStructCol(rowsStruct);

    [anyCov, firstCovCol] = max(covMask, [], 2);
    rowsCov = anyCov > 0;
    coverageLeftCol(rowsCov, f) = firstCovCol(rowsCov);

    % --- Export overlays ---
    dispImage = chooseDisplayImage(normalized, enhanced, params.outputEnhanced);
    exportFrameOverlay(dispImage, aligned, alignedMask, frontX, ...
        frameNumbers(f), folders, params);

    if mod(f, 10) == 0 || f == nFrames
        fprintf('  Processed %d/%d frames\n', f, nFrames);
    end
end

elapsed = toc;
fprintf('Processing complete. %.1f seconds (%.2f fps)\n', elapsed, nFrames / elapsed);

%% ========== TEST BOUNDARY DETECTION ==========

[testID, boundaryFrames] = detectTestBoundaries(medianFrontCol, imgWidth, params.boundary.jumpFraction);
nTests = max(testID);

fprintf('\nTest boundaries detected: %d tests\n', nTests);
for t = 1:nTests
    testFrames = find(testID == t);
    fprintf('  Test %d: frames %d-%d (%d frames)\n', t, ...
        testFrames(1), testFrames(end), length(testFrames));
end
if ~isempty(boundaryFrames)
    fprintf('  Boundary frames (new test starts): %s\n', mat2str(boundaryFrames));
end

%% ========== FRAME FILTERING ==========

frontRangePx = params.frameFilter.frontRangePercent * imgWidth;
validFrames = medianFrontCol >= frontRangePx(1) & medianFrontCol <= frontRangePx(2);

nValid = sum(validFrames);
fprintf('\nFrame filtering (front in %.0f–%.0f%% of FOV):\n', ...
    params.frameFilter.frontRangePercent(1)*100, params.frameFilter.frontRangePercent(2)*100);
fprintf('  Valid: %d / %d frames (%.0f%%)\n', nValid, nFrames, nValid/nFrames*100);

%% ========== RANDOMIZED CONVERGENCE WITH COVERAGE CORRECTION ==========

fprintf('\nComputing convergence with coverage correction...\n');

validFrameIdx = find(validFrames);
nValidFrames = length(validFrameIdx);

rng(params.convergence.rngSeed);
randomOrder = validFrameIdx(randperm(nValidFrames));

[convergence, probabilityHeatmap, coverageMap] = computeConvergence( ...
    alignedMasks, coverageMasks, randomOrder, distanceAxis, probeDistances_di, params.minCoverage);

probVsDistance = mean(probabilityHeatmap, 1, 'omitnan');

fprintf('Convergence computed over %d valid frames (seed=%d).\n', nValidFrames, params.convergence.rngSeed);
for p = 1:nProbes
    fprintf('  P(%.0f%% λ) = %.4f\n', ...
        params.probeDistances.cellFractions(p)*100, convergence.probAtProbe(end, p));
end

%% ========== REACTION ZONE WIDTH DISTRIBUTION ==========

fprintf('\nComputing reaction zone width distribution...\n');

rzWidthStats = computeReactionZoneWidths(lastStructureCol, coverageLeftCol, validFrames, ...
    targetCol, params.px2induction, params.reactionZoneWidth.censorMarginPx);

if isempty(rzWidthStats.widths)
    warning('Reaction zone width: no uncensored measurements; skipping plot.');
else
    fprintf('  N=%d rows | median=%.3g | P90=%.3g | P98=%.3g delta_i\n', ...
        numel(rzWidthStats.widths), prctileLocal(rzWidthStats.widths, 50), ...
        prctileLocal(rzWidthStats.widths, 90), prctileLocal(rzWidthStats.widths, 98));
    if rzWidthStats.censoredFraction > 0.05
        warning('  Censored rows: %d (%.1f%%) — consider raising frameFilter.frontRangePercent lower bound', ...
            rzWidthStats.nCensored, rzWidthStats.censoredFraction * 100);
    else
        fprintf('  Censored rows: %d (%.1f%%)\n', ...
            rzWidthStats.nCensored, rzWidthStats.censoredFraction * 100);
    end
end

%% ========== VISUALIZATION ==========

% --- 1. Probability Heatmap ---
fig1 = figure('Visible', 'off', 'Position', [100 100 800 600]);
imagesc(distanceAxis, 1:imgHeight, probabilityHeatmap);
set(gca, 'XDir', 'reverse');
colormap(hot);
cb = colorbar;
cb.Label.String = 'Structural Probability';
clim([0 1]);
hold on;
xline(0, 'c--', 'LineWidth', 2, 'Label', 'Front');
probeColors = {'g', 'm', 'c'};
for p = 1:nProbes
    xline(probeDistances_di(p), [probeColors{p} '-'], 'LineWidth', 1, ...
        'Label', sprintf('%.0f%% lambda', params.probeDistances.cellFractions(p)*100));
end
hold off;
xlabel('Distance from Front (\Delta_i)');
ylabel('Row (px)');
title(sprintf('Structural Probability Heatmap — Coverage Corrected (N=%d)', nValidFrames));
saveas(fig1, fullfile(folders.stats, 'Heatmap.png'));
close(fig1);

% --- 2. Probability vs Distance ---
fig2 = figure('Visible', 'off', 'Position', [100 100 800 500]);
plot(distanceAxis, probVsDistance, 'b-', 'LineWidth', 1.5);
hold on;
xline(0, 'r--', 'LineWidth', 1.5, 'Label', 'Front');
for p = 1:nProbes
    finalProb = convergence.probAtProbe(end, p);
    xline(probeDistances_di(p), [probeColors{p} '-'], 'LineWidth', 1, ...
        'Label', sprintf('%.0f%% lambda (P=%.3f)', params.probeDistances.cellFractions(p)*100, finalProb));
end
hold off;
set(gca, 'XDir', 'reverse');
xlabel('Distance from Front (\Delta_i)');
ylabel('Structural Probability');
title('Probability vs Distance Behind Front (Coverage Corrected)');
ylim([0 1]);
xlim([min(distanceAxis) max(distanceAxis)]);
grid on;
saveas(fig2, fullfile(folders.stats, 'ProbabilityCurve.png'));
close(fig2);

% --- 3. Convergence ---
fig3 = figure('Visible', 'off', 'Position', [100 100 800 500]);
for p = 1:nProbes
    plot(1:nValidFrames, convergence.probAtProbe(:, p), [probeColors{p} '-'], ...
        'LineWidth', 1.5, 'DisplayName', ...
        sprintf('%.0f%% lambda (P=%.3f)', params.probeDistances.cellFractions(p)*100, ...
        convergence.probAtProbe(end, p)));
    hold on;
end
hold off;
xlabel('Frames Accumulated (randomized order)');
ylabel('Structural Probability');
title('Convergence of P at Fixed Distances');
legend('Location', 'best');
grid on;
saveas(fig3, fullfile(folders.stats, 'Convergence.png'));
close(fig3);

% --- 4. Otsu Threshold Stability ---
fig4 = figure('Visible', 'off', 'Position', [100 100 800 400]);
plot(1:nFrames, otsuThresholds, 'b-', 'LineWidth', 1.25);
hold on;
for b = boundaryFrames(:)'
    xline(b, 'r--', 'LineWidth', 1.0);
end
hold off;
xlabel('Frame Index (sequential)');
ylabel('Otsu Threshold (intensity units)');
title(sprintf('Folded Otsu Threshold Stability — structure mask, %s (CV=%.1f%%)', ...
    params.thresholdMode, std(otsuThresholds)/mean(otsuThresholds)*100));
grid on;
saveas(fig4, fullfile(folders.stats, 'OtsuStability.png'));
close(fig4);

% --- 5. Coverage Map ---
fig5 = figure('Visible', 'off', 'Position', [100 100 800 600]);
imagesc(distanceAxis, 1:imgHeight, coverageMap);
set(gca, 'XDir', 'reverse');
colormap(parula);
cb = colorbar;
cb.Label.String = 'Frame Coverage';
hold on;
xline(0, 'r--', 'LineWidth', 2, 'Label', 'Front');
hold off;
xlabel('Distance from Front (\Delta_i)');
ylabel('Row (px)');
title(sprintf('Coverage Map — Frames with Real Data (max=%d)', nValidFrames));
saveas(fig5, fullfile(folders.stats, 'CoverageMap.png'));
close(fig5);

% --- 6. Reaction Zone Width Distribution ---
if ~isempty(rzWidthStats.widths)
    fig6 = figure('Visible', 'off', 'Position', [100 100 700 500]);
    histogram(rzWidthStats.widths, params.reactionZoneWidth.bins, ...
        'Normalization', 'pdf', 'FaceColor', [0.00 0.45 0.70], ...
        'FaceAlpha', 0.75, 'EdgeColor', 'none');
    hold on;
    pcts = params.reactionZoneWidth.percentiles;
    pctLineColors = {'g', 'm', 'c', 'r', 'b'};
    for q = 1:length(pcts)
        v = prctileLocal(rzWidthStats.widths, pcts(q));
        cIdx = mod(q-1, numel(pctLineColors)) + 1;
        xline(v, [pctLineColors{cIdx} '--'], 'LineWidth', 1.25, ...
            'Label', sprintf('P%d: %.3g', pcts(q), v));
    end
    hold off;
    xlabel('Reaction Zone Width (\delta_i)');
    ylabel('Probability Density');
    title({'Last Structural Pixel Distance from Front', ...
        sprintf('N=%d rows (%d censored = %.1f%%, %d without structure)', ...
        numel(rzWidthStats.widths), rzWidthStats.nCensored, ...
        rzWidthStats.censoredFraction*100, rzWidthStats.nNoStructure)});
    grid on;
    saveas(fig6, fullfile(folders.stats, 'ReactionZoneWidth.png'));
    close(fig6);
end

fprintf('Stats plots saved to: %s\n', folders.stats);

% --- 7. Summary figure (on-screen) ---
figure('Name', 'Analysis Results', 'Position', [50 50 1400 800]);

subplot(2, 2, 1);
imagesc(distanceAxis, 1:imgHeight, probabilityHeatmap);
set(gca, 'XDir', 'reverse');
colormap(gca, hot); colorbar; clim([0 1]);
hold on; xline(0, 'c--', 'LineWidth', 2); hold off;
xlabel('Distance (\Delta_i)'); ylabel('Row (px)'); title('Probability Heatmap');

subplot(2, 2, 2);
plot(distanceAxis, probVsDistance, 'b-', 'LineWidth', 1.5);
hold on; xline(0, 'r--', 'LineWidth', 1.5); hold off;
set(gca, 'XDir', 'reverse');
xlabel('Distance (\Delta_i)'); ylabel('Probability'); title('Probability vs Distance');
ylim([0 1]); grid on;

subplot(2, 2, 3);
for p = 1:nProbes
    plot(1:nValidFrames, convergence.probAtProbe(:, p), [probeColors{p} '-'], 'LineWidth', 1.5);
    hold on;
end
hold off;
xlabel('Frames (randomized)'); ylabel('Probability'); title('Convergence');
legend(arrayfun(@(x) sprintf('%.0f%% lambda', x*100), params.probeDistances.cellFractions, ...
    'UniformOutput', false), 'Location', 'best');
grid on;

subplot(2, 2, 4);
imagesc(distanceAxis, 1:imgHeight, coverageMap);
set(gca, 'XDir', 'reverse');
colormap(gca, parula); colorbar;
hold on; xline(0, 'r--', 'LineWidth', 2); hold off;
xlabel('Distance (\Delta_i)'); ylabel('Row (px)'); title('Coverage Map');

sgtitle(sprintf('%s Analysis Results (V7)', params.analysisName));

%% ========== SAVE & SUMMARY ==========

results.frontData      = frontData;
results.frontDataRaw   = frontDataRaw;
results.frameNumbers   = frameNumbers;
results.frameNames     = frameNames;
results.params         = params;
results.imgSize        = [imgHeight, imgWidth];
results.targetCol      = targetCol;
results.centerCol      = targetCol;  % Backward compatibility with downstream scripts
results.bgFile         = bmpFiles(bgIdx).name;
results.inputFolder    = inputFolder;
results.analysisFolder = analysisFolder;

results.segmentation.otsuThresholds  = otsuThresholds;   % structure-mask threshold
results.segmentation.frontThresholds = frontThresholds;  % front-detection threshold
results.segmentation.method          = 'foldedOtsu';
results.segmentation.segmentationSource = params.segmentationSource;
results.segmentation.thresholdMode   = params.thresholdMode;
if usePooled
    results.segmentation.pooledThreshold = pooledThresh;
end

results.tests.testID         = testID;
results.tests.boundaryFrames = boundaryFrames;
results.tests.nTests         = nTests;

results.structure.probabilityHeatmap = probabilityHeatmap;
results.structure.probVsDistance      = probVsDistance;
results.structure.distanceAxis        = distanceAxis;
results.structure.distanceAxisPx      = distanceAxisPx;
results.structure.convergence         = convergence;
results.structure.coverageMap         = coverageMap;
results.structure.nFrames             = nFrames;
results.structure.nValidFrames        = nValidFrames;
results.structure.validFrames         = validFrames;
results.structure.randomOrder         = randomOrder;

results.structure.probeDistances.cellFractions  = params.probeDistances.cellFractions;
results.structure.probeDistances.distancesDi    = probeDistances_di;
results.structure.probeDistances.distancesMm    = probeDistances_mm;
results.structure.probeDistances.finalProb      = convergence.probAtProbe(end, :);

% Per-row reaction-zone geometry + width distribution (NEW in V6)
results.structure.lastStructureCol = lastStructureCol;
results.structure.coverageLeftCol  = coverageLeftCol;

% Optional speckle-filter settings + per-frame counts (NEW in V7)
results.structure.denoise = params.structure.denoise;
if params.structure.denoise.enable
    results.structure.denoiseRemoved = denoiseRemoved;
    results.structure.denoiseKept    = denoiseKept;
end

results.reactionZoneWidth.stats          = rzWidthStats;
results.reactionZoneWidth.censorMarginPx = params.reactionZoneWidth.censorMarginPx;
results.reactionZoneWidth.percentiles    = params.reactionZoneWidth.percentiles;
results.reactionZoneWidth.unit           = 'inductionLength';

results.scale.px2mm            = params.px2mm;
results.scale.inductionLengthM  = params.inductionLengthM;
results.scale.inductionLengthMm = params.inductionLengthMm;
results.scale.px2induction      = params.px2induction;
results.scale.Ucj               = params.Ucj;
results.scale.Ea                = params.Ea;
results.scale.cellWidth         = params.cellWidth;
results.scale.mixture           = params.mixture;

matPath = fullfile(analysisFolder, 'FrontData.mat');
save(matPath, 'results', '-v7.3');
fprintf('Data saved to: %s\n', matPath);

fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d (%d tests), %d valid for accumulation\n', nFrames, nTests, nValidFrames);
fprintf('Alignment: front at col %d (buffer=%d px from right edge)\n', targetCol, params.alignBuffer);
fprintf('Average detection rate: %.1f%%\n', ...
    mean(sum(~isnan(frontData), 1) / imgHeight * 100));
fprintf('Otsu threshold (structure): mean=%.1f, CV=%.1f%%  |  source=%s, mode=%s\n', ...
    mean(otsuThresholds), std(otsuThresholds)/mean(otsuThresholds)*100, ...
    params.segmentationSource, params.thresholdMode);

if params.structure.denoise.enable
    totRemoved = sum(denoiseRemoved);
    totStruct  = totRemoved + sum(denoiseKept);
    fprintf('Denoise (bwareaopen, minPx=%d, conn=%d): removed %d / %d structural px (%.2f%%)\n', ...
        params.structure.denoise.minPixels, params.structure.denoise.connectivity, ...
        totRemoved, totStruct, 100 * totRemoved / max(1, totStruct));
end

fprintf('\nPhysical Scale:\n');
fprintf('  px2mm:            %.4f mm/px\n', params.px2mm);
fprintf('  Induction length: %.3e m (%.4f mm)\n', params.inductionLengthM, params.inductionLengthMm);
fprintf('  Cell width:       %.2f mm\n', params.cellWidth);
fprintf('  px2induction:     %.4f induction lengths/px\n', params.px2induction);

fprintf('\nCoverage: min=%d, max=%d frames per pixel (behind front)\n', ...
    min(coverageMap(coverageMap > 0)), max(coverageMap(:)));

fprintf('\nStructural Probability at Fixed Distances:\n');
for p = 1:nProbes
    fprintf('  P(%.0f%% lambda) = %.4f  (at %.2f di = %.4f mm)\n', ...
        params.probeDistances.cellFractions(p)*100, ...
        convergence.probAtProbe(end, p), probeDistances_di(p), probeDistances_mm(p));
end

if ~isempty(rzWidthStats.widths)
    fprintf('\nReaction Zone Width (delta_i), source=%s:\n', params.segmentationSource);
    fprintf('  N=%d rows | median=%.3g | P90=%.3g | P98=%.3g\n', ...
        numel(rzWidthStats.widths), prctileLocal(rzWidthStats.widths, 50), ...
        prctileLocal(rzWidthStats.widths, 90), prctileLocal(rzWidthStats.widths, 98));
    fprintf('  Censored: %d (%.1f%%) | No structure: %d rows\n', ...
        rzWidthStats.nCensored, rzWidthStats.censoredFraction*100, rzWidthStats.nNoStructure);
end

fprintf('\nOutputs:\n');
fprintf('  Front overlays:     %s\n', folders.fronts);
fprintf('  Aligned images:     %s\n', folders.aligned);
fprintf('  Structure overlays: %s\n', folders.structure);
fprintf('  Statistics:         %s\n', folders.stats);
fprintf('  Script archive:     %s\n', scriptCopy);
fprintf('==============================\n');

%% ========== LOCAL FUNCTIONS ==========

function [entry, db] = loadOrCreateTestEntry(dbPath, analysisName)
% Load test parameters from database by analysisName, or prompt user to create a new entry.
%
% Inputs:
%   dbPath       - full path to the database .mat file
%   analysisName - string key (derived from input folder name)
%
% Outputs:
%   entry - single-row table with test parameters
%   db    - full database table

    VARS = {'analysisName', 'mixture', 'pressure', 'inductionLengthM', 'Ucj', 'Ea', 'px2mm', 'cellWidth', 'dateCreated'};

    if isfile(dbPath)
        load(dbPath, 'rzDatabase');
    else
        rzDatabase = table(string.empty, string.empty, [], [], [], [], [], [], string.empty, ...
            'VariableNames', VARS);
        fprintf('Creating new reaction zone database: %s\n', dbPath);
    end

    idx = (rzDatabase.analysisName == string(analysisName));

    if any(idx)
        entry = rzDatabase(idx, :);
        fprintf('Loaded existing database entry for "%s" (created %s)\n', ...
            analysisName, entry.dateCreated);
    else
        fprintf('\n--- New dataset: "%s" ---\n', analysisName);
        fprintf('Enter test parameters:\n');
        mixture          = input('  Mixture (e.g. 2H2_O2_2N2): ', 's');
        pressure         = input('  Initial pressure [kPa]: ');
        inductionLengthM = input('  ZND induction length [m]: ');
        Ucj              = input('  CJ velocity [m/s]: ');
        Ea               = input('  Reduced activation energy (Ea/RT): ');
        px2mm            = input('  Pixel scale [mm/px]: ');
        cellWidth        = input('  Cell width [mm]: ');

        newRow = {string(analysisName), string(mixture), pressure, inductionLengthM, ...
                  Ucj, Ea, px2mm, cellWidth, string(datestr(now, 'yyyy-mm-dd'))};
        entry = cell2table(newRow, 'VariableNames', VARS);
        rzDatabase = [rzDatabase; entry];

        dbDir = fileparts(dbPath);
        if ~exist(dbDir, 'dir'), mkdir(dbDir); end

        save(dbPath, 'rzDatabase');
        fprintf('Saved new entry to database.\n');
    end

    db = rzDatabase;
end


function normalized = loadAndNormalize(framePath, bgImage)
% Load a BMP frame, subtract background, and zero-median normalize.
%
% Inputs:
%   framePath - full path to BMP file
%   bgImage   - [H x W] double background image
%
% Output:
%   normalized - [H x W] double, background-subtracted and zero-median

    frameImage = double(imread(framePath));
    if ndims(frameImage) == 3
        frameImage = double(rgb2gray(uint8(frameImage)));
    end
    subtracted = frameImage - bgImage;
    normalized = subtracted - median(subtracted(:));
end


function enhanced = preprocessForDetection(normalized, prepOpts)
% Apply median filter and unsharp mask to sharpen edges for front detection.
%
% Inputs:
%   normalized - [H x W] double normalized image
%   prepOpts   - struct with fields: medfiltSize, unsharpRadius, unsharpAmount
%
% Output:
%   enhanced   - [H x W] double with improved edge contrast

    filtered = medfilt2(normalized, [prepOpts.medfiltSize, prepOpts.medfiltSize]);
    blurred  = imgaussfilt(filtered, prepOpts.unsharpRadius);
    enhanced = filtered + prepOpts.unsharpAmount * (filtered - blurred);
end


function thresh = threshFoldedOtsu(enhanced)
% Otsu threshold on the absolute-value (folded) intensity distribution.
%
% Input:  enhanced - [H x W] double, zero-centered
% Output: thresh   - scalar threshold on abs(intensity) scale

    absVals = abs(enhanced(:));
    maxVal = max(absVals);
    if maxVal == 0
        thresh = 0;
        return;
    end
    otsuLevel = graythresh(absVals / maxVal);
    thresh = otsuLevel * maxVal;
end


function thresh = threshFoldedOtsuPooledStack(normStack, segmentationSource, prepOpts)
% Folded-Otsu threshold pooled over an entire normalized stack.
%
% Builds the abs-intensity histogram of the segmentation-source image across ALL
% frames, then applies Otsu to that single pooled histogram. Two in-memory passes
% (global max, then counts) — no disk I/O. Returns one threshold shared by every
% frame, matching the Python port's pooled mode. Mirrors threshFoldedOtsu's
% normalize-by-max convention so a single-frame stack reproduces the per-frame value.
%
% Inputs:
%   normStack          - [H x W x N] normalized stack (single or double)
%   segmentationSource - 'normalized' | 'enhanced'
%   prepOpts           - preprocessing struct (used only for 'enhanced')
%
% Output:
%   thresh - scalar pooled threshold on the abs(intensity) scale

    NBINS = 256;
    nFrames = size(normStack, 3);

    % Pass A: global max of the folded segmentation source over the stack
    globalMax = 0;
    for f = 1:nFrames
        seg = segSourceImage(double(normStack(:, :, f)), segmentationSource, prepOpts);
        m = max(abs(seg(:)));
        if m > globalMax, globalMax = m; end
    end
    if globalMax == 0
        thresh = 0;
        return;
    end

    % Pass B: accumulate the pooled histogram on common edges, then Otsu on it
    edges = linspace(0, globalMax, NBINS + 1);
    counts = zeros(1, NBINS);
    for f = 1:nFrames
        seg = segSourceImage(double(normStack(:, :, f)), segmentationSource, prepOpts);
        counts = counts + histcounts(abs(seg(:)), edges);
    end

    otsuLevel = otsuthresh(counts);   % normalized [0,1] threshold on the pooled histogram
    thresh = otsuLevel * globalMax;   % map back to abs(intensity) units
end


function seg = segSourceImage(normalized, segmentationSource, prepOpts)
% Return the image used for structure segmentation in the requested space.
%
% 'normalized' returns the input unchanged; 'enhanced' applies the same
% median + unsharp preprocessing used for front detection. Keeps the pooled
% pre-pass consistent with the per-frame structure mask.
%
% Inputs:  normalized - [H x W] normalized image; segmentationSource - space;
%          prepOpts - preprocessing struct (used only for 'enhanced')
% Output:  seg - [H x W] image in the requested space

    switch segmentationSource
        case 'normalized'
            seg = normalized;
        case 'enhanced'
            seg = preprocessForDetection(normalized, prepOpts);
        otherwise
            error('Unknown segmentationSource: %s (use ''normalized'' or ''enhanced'')', ...
                segmentationSource);
    end
end


function frontX = extractFrontRaw(mask)
% Extract the raw wavefront as the rightmost segmented pixel in each row.
%
% No outlier rejection or smoothing — this is the unprocessed detection
% (V6 split so frontDataRaw can store it before refineFront touches it).
%
% Input:  mask   - [H x W] logical segmentation mask
% Output: frontX - [H x 1] rightmost segmented column per row (NaN if none)

    imgHeight = size(mask, 1);
    frontX = NaN(imgHeight, 1);

    for row = 1:imgHeight
        segCols = find(mask(row, :));
        if ~isempty(segCols)
            % Rightmost (max) crossing — the front is the leading edge of the wave
            frontX(row) = max(segCols);
        end
    end
end


function frontX = refineFront(frontX, frontOpts)
% Reject outliers and median-smooth a raw front detection (NaN-aware).
%
% Outliers (deviation from the local median > outlierThreshPx) are replaced
% by the local median, then the whole trace is median-smoothed. Rows that
% were NaN on input remain NaN. Bit-identical to V5's estimateFrontFromMask
% post-extraction logic.
%
% Inputs:
%   frontX    - [H x 1] raw front column per row (NaN = missing)
%   frontOpts - struct with fields: smoothWindow, outlierThreshPx
%
% Output:
%   frontX    - [H x 1] refined front column per row

    validIdx = find(~isnan(frontX));
    if length(validIdx) < 5
        return;
    end

    localMed = movmedian(frontX, frontOpts.smoothWindow, 'omitnan');
    deviation = abs(frontX - localMed);
    isOutlier = deviation > frontOpts.outlierThreshPx & ~isnan(frontX);
    frontX(isOutlier) = localMed(isOutlier);

    nanMask = isnan(frontX);
    win = frontOpts.smoothWindow;
    if mod(win, 2) == 0
        win = win + 1;
    end
    frontX = movmedian(frontX, win, 'omitnan');
    frontX(nanMask) = NaN;
end


function aligned = alignToFront(image, frontX, padValue, targetCol)
% Shift each row so the detected front lands at targetCol.
%
% Rows with missing front data are interpolated/extrapolated from neighbors.
% Pixels shifted outside the frame are filled with padValue.
%
% Inputs:
%   image     - [H x W] grayscale image (double)
%   frontX    - [H x 1] front column position per row (NaN = missing)
%   padValue  - scalar fill value for pixels shifted out of frame
%   targetCol - column where the front should land after alignment
%
% Output:
%   aligned   - [H x W] shifted image with front at targetCol

    [nRows, nCols] = size(image);
    aligned = ones(nRows, nCols) * padValue;

    validIdx = find(~isnan(frontX));
    if isempty(validIdx)
        aligned = image;
        return;
    end

    allRows = (1:nRows)';
    if length(validIdx) < 2
        frontInterp = ones(nRows, 1) * frontX(validIdx(1));
    else
        frontInterp = interp1(validIdx, frontX(validIdx), allRows, 'linear', 'extrap');
    end

    for row = 1:nRows
        shift = round(targetCol - frontInterp(row));
        if shift == 0
            aligned(row, :) = image(row, :);
        elseif shift > 0
            srcEnd = nCols - shift;
            if srcEnd >= 1
                aligned(row, shift+1:nCols) = image(row, 1:srcEnd);
            end
        else
            srcStart = abs(shift) + 1;
            dstEnd = nCols - abs(shift);
            if srcStart <= nCols && dstEnd >= 1
                aligned(row, 1:dstEnd) = image(row, srcStart:nCols);
            end
        end
    end
end


function covMask = buildCoverageMask(frontX, imgHeight, imgWidth, targetCol)
% Binary mask of pixels with real (non-padded) data after alignment.
%
% Computes which pixels would contain real image data (vs padding) when
% the image is aligned so the front lands at targetCol. Equivalent to
% aligning a ones-matrix, but computed analytically for speed.
%
% Inputs:
%   frontX    - [H x 1] front column per row (NaN = missing)
%   imgHeight - image height
%   imgWidth  - image width
%   targetCol - alignment target column
%
% Output:
%   covMask   - [H x W] logical, true = real data after alignment

    covMask = false(imgHeight, imgWidth);

    validIdx = find(~isnan(frontX));
    if isempty(validIdx)
        covMask = true(imgHeight, imgWidth);
        return;
    end

    allRows = (1:imgHeight)';
    if length(validIdx) < 2
        frontInterp = ones(imgHeight, 1) * frontX(validIdx(1));
    else
        frontInterp = interp1(validIdx, frontX(validIdx), allRows, 'linear', 'extrap');
    end

    for row = 1:imgHeight
        shift = round(targetCol - frontInterp(row));
        if shift == 0
            covMask(row, :) = true;
        elseif shift > 0
            % Data lands in columns shift+1 : imgWidth
            srcEnd = imgWidth - shift;
            if srcEnd >= 1
                covMask(row, shift+1:imgWidth) = true;
            end
        else
            % Data lands in columns 1 : imgWidth-|shift|
            dstEnd = imgWidth - abs(shift);
            srcStart = abs(shift) + 1;
            if srcStart <= imgWidth && dstEnd >= 1
                covMask(row, 1:dstEnd) = true;
            end
        end
    end
end


function [testID, boundaryFrames] = detectTestBoundaries(medianFrontCol, imgWidth, jumpFraction)
% Flag frames where the median front position jumps leftward, indicating a new test.
%
% Inputs:
%   medianFrontCol - [N x 1] median front column per frame
%   imgWidth       - image width in pixels
%   jumpFraction   - fraction of image width that constitutes a jump
%
% Outputs:
%   testID         - [N x 1] integer test number per frame
%   boundaryFrames - indices where a new test begins

    nFrames = length(medianFrontCol);
    jumpThreshPx = jumpFraction * imgWidth;

    testID = ones(nFrames, 1);
    boundaryFrames = [];
    currentTest = 1;

    for f = 2:nFrames
        if isnan(medianFrontCol(f)) || isnan(medianFrontCol(f-1))
            testID(f) = currentTest;
            continue;
        end

        leftwardJump = medianFrontCol(f-1) - medianFrontCol(f);

        if leftwardJump > jumpThreshPx
            currentTest = currentTest + 1;
            boundaryFrames = [boundaryFrames; f]; %#ok<AGROW>
        end

        testID(f) = currentTest;
    end
end


function exportFrameOverlay(dispImage, aligned, structureMask, frontX, frameNum, folders, params)
% Write front overlay, aligned image, and structure overlay to disk.

    fmt = params.outputFormat;

    dispRGB = repmat(toDisplayUint8(dispImage), [1, 1, 3]);
    dispRGB = drawFrontLine(dispRGB, frontX, params.overlayColor, params.overlayLinewidth);
    imwrite(dispRGB, fullfile(folders.fronts, sprintf('%03d_Front.%s', frameNum, fmt)));

    alignedOut = toDisplayUint8(aligned);
    imwrite(alignedOut, fullfile(folders.aligned, sprintf('%03d_Aligned.%s', frameNum, fmt)));

    structRGB = repmat(alignedOut, [1, 1, 3]);
    structRGB = applyMaskOverlay(structRGB, structureMask, params.structure.overlayColor);
    imwrite(structRGB, fullfile(folders.structure, sprintf('%03d_Structure.%s', frameNum, fmt)));
end


function [convergence, probabilityHeatmap, coverageMap] = computeConvergence( ...
    alignedMasks, coverageMasks, randomOrder, distanceAxis, probeDistances_di, minCoverage)
% Accumulate structure masks with coverage correction and track convergence.
%
% Normalizes by per-pixel coverage (number of frames with real data) instead
% of total frame count. This gives the true conditional probability
% P(structure | real data exists) at every pixel.
%
% Inputs:
%   alignedMasks     - [H x W x N] logical, aligned structure masks
%   coverageMasks    - [H x W x N] logical, real-data masks per frame
%   randomOrder      - [1 x M] frame indices defining accumulation order
%   distanceAxis     - [1 x W] distance in induction lengths
%   probeDistances_di - [1 x P] probe distances in induction lengths
%   minCoverage      - minimum frame count for a pixel to be included
%
% Outputs:
%   convergence        - struct with probAtProbe [M x P] and randomOrder
%   probabilityHeatmap - [H x W] coverage-corrected probability (NaN where coverage < min)
%   coverageMap        - [H x W] total frame coverage count

    [imgHeight, imgWidth, ~] = size(alignedMasks);
    nAccum = length(randomOrder);
    nProbes = length(probeDistances_di);

    convergence.probAtProbe = NaN(nAccum, nProbes);
    convergence.randomOrder = randomOrder;
    convergence.probeDistances_di = probeDistances_di;

    structureAccum = zeros(imgHeight, imgWidth);
    coverageAccum  = zeros(imgHeight, imgWidth);

    for k = 1:nAccum
        f = randomOrder(k);
        structureAccum = structureAccum + double(alignedMasks(:, :, f));
        coverageAccum  = coverageAccum  + double(coverageMasks(:, :, f));

        % Coverage-corrected heatmap: P = structure / coverage where coverage > 0.
        % V6 fix: init to NaN (not 0) so the 'omitnan' row-mean below actually
        % excludes below-coverage-floor pixels instead of counting them as P=0.
        currentHeatmap = NaN(imgHeight, imgWidth);
        validPx = coverageAccum >= max(1, min(k, minCoverage));
        currentHeatmap(validPx) = structureAccum(validPx) ./ coverageAccum(validPx);

        probVsDist = mean(currentHeatmap, 1, 'omitnan');

        for p = 1:nProbes
            convergence.probAtProbe(k, p) = interpProbAtDistance( ...
                distanceAxis, probVsDist, probeDistances_di(p));
        end
    end

    % Final heatmap with coverage correction
    coverageMap = coverageAccum;
    probabilityHeatmap = NaN(imgHeight, imgWidth);
    validPx = coverageAccum >= minCoverage;
    probabilityHeatmap(validPx) = structureAccum(validPx) ./ coverageAccum(validPx);
end


function prob = interpProbAtDistance(distanceAxis, probVsDist, targetDistance)
% Interpolate the probability value at a specific distance from the front.

    [dSorted, si] = sort(distanceAxis);
    pSorted = probVsDist(si);

    valid = ~isnan(pSorted);
    dSorted = dSorted(valid);
    pSorted = pSorted(valid);

    if length(dSorted) < 2 || targetDistance < min(dSorted) || targetDistance > max(dSorted)
        prob = NaN;
        return;
    end

    prob = interp1(dSorted, pSorted, targetDistance, 'linear');
end


function stats = computeReactionZoneWidths(lastStructureCol, coverageLeftCol, ...
        validFrames, targetCol, px2unit, censorMarginPx)
% Per-row reaction-zone widths across all valid frames, with FOV censoring.
%
% For each row of each valid frame, width = (targetCol - lastStructureCol) in
% pixels, converted to normalization units. Rows whose trailing (leftmost)
% structural pixel sits within censorMarginPx of the real-data boundary are
% censored: the true trailing edge may lie beyond the field of view, so the
% measured width would understate the real one.
%
% Inputs:
%   lastStructureCol - [H x N] leftmost structural column per row/frame (NaN = none)
%   coverageLeftCol  - [H x N] leftmost real-data column per row/frame (NaN = none)
%   validFrames      - [N x 1] logical, frames to include
%   targetCol        - front column after alignment
%   px2unit          - pixel -> normalization-unit scale factor
%   censorMarginPx   - censoring margin in pixels
%
% Output:
%   stats - struct with fields:
%       widths           - [M x 1] uncensored widths (normalization units)
%       nCensored        - rows censored against the FOV boundary
%       censoredFraction - nCensored / (rows with structure)
%       nNoStructure     - rows with no structural pixels

    lsc = lastStructureCol(:, validFrames);
    clc = coverageLeftCol(:, validFrames);
    lsc = lsc(:);
    clc = clc(:);

    hasStructure = ~isnan(lsc);
    nNoStructure = sum(~hasStructure);

    % Trailing edge against the real-data boundary -> possibly truncated by FOV.
    % A NaN coverage column means no real data in that row, so censor it too.
    isCensored = hasStructure & (isnan(clc) | lsc <= clc + censorMarginPx);
    nCensored  = sum(isCensored);

    isUncensored = hasStructure & ~isCensored;
    widths = (targetCol - lsc(isUncensored)) * px2unit;

    nWithStructure = sum(hasStructure);
    if nWithStructure > 0
        censoredFraction = nCensored / nWithStructure;
    else
        censoredFraction = 0;
    end

    stats.widths           = widths;
    stats.nCensored        = nCensored;
    stats.censoredFraction = censoredFraction;
    stats.nNoStructure     = nNoStructure;
end


function v = prctileLocal(x, p)
% Percentile via MATLAB's (i-0.5)/n plotting-position convention.
%
% Avoids the Statistics Toolbox dependency that prctile/quantile carry, so the
% reaction-zone-width stats run on a base + Image Processing Toolbox install.
% Matches prctile: linear interpolation between plotting positions, clamped to
% the extremes outside the sampled range.
%
% Inputs:  x - data vector (NaNs ignored); p - percentile(s) in [0, 100]
% Output:  v - percentile value(s), same shape as p

    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0
        v = NaN(size(p));
        return;
    elseif n == 1
        v = repmat(x, size(p));
        return;
    end

    pos = 100 * ((1:n) - 0.5) / n;        % plotting positions (row vector)
    v = interp1(pos, x(:)', p, 'linear'); % x as row to match pos orientation
    v(p < pos(1))   = x(1);
    v(p > pos(end)) = x(end);
end


function rgbImage = drawFrontLine(rgbImage, frontX, color, lineWidth)
% Overlay detected front as a colored line on an RGB image.

    [~, imgWidth, ~] = size(rgbImage);

    validRows = find(~isnan(frontX));
    if isempty(validRows), return; end

    allRows = validRows(1):validRows(end);
    interpX = interp1(validRows, frontX(validRows), allRows, 'linear');
    halfWidth = floor(lineWidth / 2);

    for i = 1:length(allRows)
        row = allRows(i);
        col = round(interpX(i));
        for dc = -halfWidth:halfWidth
            c = col + dc;
            if c >= 1 && c <= imgWidth
                rgbImage(row, c, 1) = color(1);
                rgbImage(row, c, 2) = color(2);
                rgbImage(row, c, 3) = color(3);
            end
        end
    end
end


function rgbImage = applyMaskOverlay(rgbImage, mask, color)
% Paint binary mask pixels onto RGB image in the given color.

    maskIdx = find(mask);
    if isempty(maskIdx), return; end

    [imgHeight, imgWidth, ~] = size(rgbImage);
    nPixels = imgHeight * imgWidth;

    rgbImage(maskIdx)              = color(1);
    rgbImage(maskIdx + nPixels)    = color(2);
    rgbImage(maskIdx + 2*nPixels)  = color(3);
end


function dispOut = toDisplayUint8(image)
% Convert zero-centered double image to uint8 for display/export.
    DISPLAY_OFFSET = 128;
    dispOut = uint8(max(0, min(255, image + DISPLAY_OFFSET)));
end


function dispImage = chooseDisplayImage(normalized, enhanced, useEnhanced)
% Select which image to use for front overlay export.
    if useEnhanced
        dispImage = enhanced;
    else
        dispImage = normalized;
    end
end