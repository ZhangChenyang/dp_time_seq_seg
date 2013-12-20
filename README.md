dp_time_seq_seg
===============

A Time sequence segmentation implementation using Dynamic Programming in Matlab

Usage:

global C; global V;

%% initialize some A \in R^{m,n}
%% specify some K, number of segments you want

>>[e,s,m] = dp_tseg_c(A,K,0);

%% e: total error
%% s: segmentation, K-1 integers
%% m: K means of the segments
