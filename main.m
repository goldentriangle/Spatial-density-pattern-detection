clear all; close all; clc;
warning off all;

filename = 'tile.png';
filepath = './data/';
readfile = strcat(filepath,filename);

% set the fourth parameter seedThresh to [0, *] to detect only over-density
% or [*, 255] to detect only under-density.
% for different images, the same image with different resolution, the
% parameters thresVal, seedThresh, minArea, maxAreaRatio may need to be
% adjusted to get optimal results
anomaly_detection(readfile, true, 40, [20,230], 10); 

