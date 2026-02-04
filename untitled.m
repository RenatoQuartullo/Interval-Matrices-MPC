clear
clc
close all

files1 = matlab.codetools.requiredFilesAndProducts('main_population_runtimes.m');
toolsFiles1 = files1(contains(files1, [filesep 'tools' filesep]));
files2 = matlab.codetools.requiredFilesAndProducts('main_population.m');
toolsFiles2 = files2(contains(files2, [filesep 'tools' filesep]));
files3 = matlab.codetools.requiredFilesAndProducts('main_coverage.m');
toolsFiles3 = files3(contains(files3, [filesep 'tools' filesep]));
files4 = matlab.codetools.requiredFilesAndProducts('main_feasible_domain_DI.m');
toolsFiles4 = files4(contains(files4, [filesep 'tools' filesep]));
allFiles = [toolsFiles1, toolsFiles2, toolsFiles3, toolsFiles4]';
allFilesUnique = unique(allFiles, 'stable')