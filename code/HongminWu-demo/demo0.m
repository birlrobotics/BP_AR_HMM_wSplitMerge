clc;
clear all;
close all;

% for mocap6 dataset
%function [] = RunMocap6Experiment( jobID, taskID, infName, initName, TimeLimit )
%INPUT
%  jobID : integer/name of job
%  taskID : integer/name of task
%  infName  : {'Prior','DD','SM'} indicates type of sampler     
%  initName : {'one', 'unique5'}  indicates initialization
%  TimeLimit : # of seconds to allow for MCMC

RunMocap6Experiment('2', '1', 'SM+DD' , 'seq', 6)