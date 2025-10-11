clear
close all
clc
tic;

%% load data

addpath('utils');
data=load('sardata.mat');data=data.sardata;
sarData=data.range_slice;figure;imagesc(abs(sarData));    %range slice echo
xyz=data.slice_xyz;                                       %dx dy and distance z0
params=data.params;

%% imaging parameters

dx=xyz.dx;
dy=xyz.dy;                        % Sampling interval at x,y axis in mm
z0=xyz.distance;                  %distance of radar and target at mm
f0=params.f0;                     % start frequency
c=params.c;                       % speed of light
fslope=params.fslope;             % Slope const (Hz/sec)
fs=params.fs;                     % Sampling rate (sps)
num_sample=params.num_sample;     % Number of sampling points at z axis
nFFTspace=params.nFFTspace;       % Number of FFT points for wave-domain

k = 2*pi*f0/c;                    %wave number
imSize =400;                      %image size at mm
amplitude=-40;                    % amplitude of [-40,0] dB after normalization

%% mask of different sparse sampling

[M,N]= size(sarData);mask=zeros(M,N);truncated_L=floor(N/2);
% uniform row undersample
% masknum=4;mask([1:masknum:M],:)=1;     

% 30% non-uniform row undersample
mask([1,3,10,12,16,20,25,26,30,32,36,39,44,48,50,52,57,61,66,68,70,72,75,79,82,84,87,88,93,97],:)=1;

%sarData=sarData.*mask;figure;imagesc(abs(sarData));     % mask
fprintf('percent=%.6f\n',nnz(sarData)/numel(sarData));  % percent
maxMod=max(abs(sarData(:)));sarData=sarData/maxMod;     %normalization

%% Truncated-DCT  (TD)

sarData=TD(sarData,truncated_L);
figure;imagesc(abs(sarData));

%% Matrix Completion  (TVecA/TCTA)

% proposed TVecA        
%R=80;   % pencil parameter
%mu=2;           %penalty parameter
%e_rank=R/4;     %e_rank
%K=30;           %iterations
%sarData=TVecA(sarData,R,mu,K,e_rank);

% proposed TCTA
% P=10;  % pencil parameter
% Q=P;            % pencil parameter
% mu=2;           %penalty parameter
% e_rank=P+Q;     %e_rank 
% K=30;           %iterations
% sarData=TCTA(sarData,P,Q,mu,K,e_rank);

figure;imagesc(abs(sarData));

%% zero-padding and IDCT(ITD)

sarData=ITD(sarData,N,truncated_L);
figure;imagesc(abs(sarData));

%% rma 2d imaging

rma_2d(dx,dy,k,z0,sarData,nFFTspace,amplitude,imSize); %imaging of one range slice 

elapsedTime = toc ;


