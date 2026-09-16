% Structured-Bessel-beam tractor-field simulation and 3-D visualization.
% 
% This program implements the optical landscape described by Rahman et al.,
% "Tractor beam for fully immersed multiple objects: Long distance pulling,
% trapping, and rotation with a single optical set-up," Annalen der Physik 527,
% 777-793 (2015), DOI 10.1002/andp.201500266.
% 
% Two higher-order Bessel beams with reverse orders produce ``2*abs(m)``
% azimuthal trapping regions. Unequal longitudinal wave numbers produce axial
% modulation, and a small frequency difference moves the interference landscape
% without an externally ramped phase.
%%
clc
clear
close all
projectFolder = fileparts(mfilename('fullpath'));
addpath(projectFolder);

%% 2. Shared options -- control precision, model compatibility, and exports.
options.geometryMode = "legacy"; % "legacy" reproduces the source geometry.
options.storageClass = "double"; % "single" reduces stored-volume memory.
options.storeGasIntensity = false;
options.saveFigures = false;     % Save both .fig and .png when true.
options.saveResults = false;     % Save each result structure as .mat.
options.outputFolder = fullfile(projectFolder, 'results');

%% 3. Figure 1 parameters
config = BesselBeamSimulation.preset("ADP15kHz");
config(1).beam.electricAmplitude = 1.5e6;
config(1).beam.frequencyDifferenceHz = 15000;
config(1).beam.betaDifferenceFraction = 0.16090036886470404;
config(1).beam.incidenceDeg = [45 45];
config(1).beam.orders = [2 -2];
config(1).beam.wavelengthM = 1e-6;
config(1).beam.coneAngleDeg = 40;
config(1).timeS = 0;
config(1).grid.xHalfWidthM = 1e-6;
config(1).grid.yHalfWidthM = 1e-6;
config(1).grid.xStepM = 20e-9;
config(1).grid.yStepM = 20e-9;
config(1).grid.zLimitsM = [-3e-6 3e-6];
config(1).grid.zStepM = 60e-9;
config(1).plot.isoLevel = 0.99;
config(1).particle.displayRadiusM = 150e-9;

%% 4. Figure 2 parameters
config(2) = BesselBeamSimulation.preset("Helical5kHz");
config(2).beam.electricAmplitude = 2.5e6;
config(2).beam.frequencyDifferenceHz = 5000;
config(2).beam.betaDifferenceFraction = 0.45;
config(2).beam.incidenceDeg = [45 45];
config(2).beam.orders = [2 -2];
config(2).beam.wavelengthM = 1e-6;
config(2).beam.coneAngleDeg = 40;
config(2).timeS = 135e-6;         % Original tt=45, t=tt*3e-6.
config(2).grid.xHalfWidthM = 1.5e-6;
config(2).grid.yHalfWidthM = 1.5e-6;
config(2).grid.xStepM = 15e-9;
config(2).grid.yStepM = 15e-9;
config(2).grid.zLimitsM = [-3e-6 3e-6];
config(2).grid.zStepM = 60e-9;
config(2).plot.isoLevel = 0.98;
config(2).particle.displayRadiusM = 180e-9;

%% 5. Run both cases -- each object validates, caches, simulates, and plots.
simulations = cell(1, numel(config));
results = cell(1, numel(config));
figures = gobjects(1, numel(config));
for caseIndex = 1:numel(config)
    config(caseIndex).numerics.geometryMode = options.geometryMode;
    config(caseIndex).numerics.storageClass = options.storageClass;
    config(caseIndex).numerics.storeGasIntensity = options.storeGasIntensity;
    config(caseIndex).output.saveFigures = options.saveFigures;
    config(caseIndex).output.saveResults = options.saveResults;
    config(caseIndex).output.folder = options.outputFolder;
    simulations{caseIndex} = BesselBeamSimulation(config(caseIndex));
    results{caseIndex} = simulations{caseIndex}.run();
    figures(caseIndex) = simulations{caseIndex}.plot();
    simulations{caseIndex}.saveOutputs(figures(caseIndex));
end

%% 6. Summary -- raw |E|^2 and plane-normalized intensity are in results{k}.
fprintf('\nFinished. Figure 1: ADP 15 kHz; Figure 2: helical 5 kHz.\n');
fprintf('Use results{1} and results{2} to inspect grids, intensities, and metrics.\n');
