clear
close all
clc
tic;

%% load data

addpath('utils');
data=load('sardata.mat');
data=data.sardata;
sarData=data.range_slice;
figure;imagesc(abs(sarData));axis square;colormap('jet');             % range slice echo
title 'Fully Sampled SAR Echo';
xyz=data.slice_xyz;                                                   % dx dy and distance z0
params=data.params;

%% chose MC methods(TCTA/TVNNA)

IS_TCTA  = 0;
IS_TVNNA = 1;

%% imaging parameters

dx=xyz.dx;
dy=xyz.dy;                        % sampling interval at x,y axis in mm
z0=xyz.distance;                  % distance of radar and target at mm
f0=params.f0;                     % start frequency
c=params.c;                       % speed of light
fslope=params.fslope;             % slope const (Hz/sec)
fs=params.fs;                     % sampling rate (sps)
num_sample=params.num_sample;     % number of sampling points at z axis
nFFTspace=params.nFFTspace;       % number of FFT points for wave-domain
k = 2*pi*f0/c;                    % wave number
imSize =400;                      % image size at mm
amplitude=-40;                    % amplitude of [-40,0] dB after normalization

%% mask of different sparse sampling pattern

[M,N]= size(sarData);mask=zeros(M,N);

%% random sparse

% samp_rate = 0.2; chosen = randperm(M*N,round(samp_rate*M*N)); mask(chosen)=1;                                 % random: 20%

%% structured sparse

% % uniform sparse

% masknum=2;mask([1:masknum:end],:)=1;                                                                          % row: 50%

% masknum=4;mask([1:masknum:end],:)=1;                                                                          % row: 25%

% masknum=2;mask(:,[1:masknum:end])=1;                                                                          % column: 50%

% masknum=4;mask(:,[1:masknum:end])=1;                                                                          % column: 25%


% % non-uniform

% mask([1,3,4,7,9,10,12,14,15,18,20,21,24,26,29,31,32,34,36,39,41,43,45,47,49, ...
%     50,53,55,58,59,61,62,65,67,69,70,72,75,76,79,80,82,84,87,88,90,93,95,97,100],:)=1;                        % row 50%

mask([3,7,10,12,16,20,22,26,31,36,39,42,44,48,53,57,58,61,64,67,70,72,75,79,82,84,87,90,93,97],:)=1;          % row 30%

% mask([7,11,15,22,27,32,36,39,42,45,48,51,53,57,61,65,68,71,75,78,81,84,88,93,96],:)=1;                        % row 25%

% mask([11,15,22,27,33,39,43,48,52,56,60,63,67,70,73,78,81,84,88,93],:)=1;                                      % row 20%

% mask([7,15,25,36,42,48,54,59,64,70,76,81,87,90,96],:)=1;                                                      % row 15%

% mask(:,[1,3,10,12,16,20,25,26,30,32,36,39,44,48,50,52,57,61,66,68,70,72,75,79, ...
%     82,84,87,88,93,97,101,105,109,111,115,119,124,125,129,131,135,138,143,147, ...
%     149,151,156,160,165,167,169,171,174,178,181,183,186,187,192,196])=1;                                      % column 30%

% mask(:,[3,6,12,14,20,25,28,32,39,45,47,57,51,61,66,70,74,79,82,88,93,96,101, ...
%     107,112,118,127,124,136,133,141,147,154,158,163,169,172,177,184,181,190,197])=1;                          % column 20%

% mask(:,[7,15,22,29,36,42,48,57,64,71,78,84,90,96,100,102,109,115,122,127,136, ...
%     142,148,157,164,171,178,184,190,193,197])=1;                                                              % column 15%


% zigzag sparse sample

% masknum=2;mask([1:masknum*2:M,4:masknum*2:M],1:100)=1;mask([2:masknum*2:M,3:masknum*2:M],101:N)=1;            % 50%

% masknum=3;mask([1:masknum*2:M,6:masknum*2:M],1:67)=1;
% mask([2:masknum*2:M,5:masknum*2:M],68:133)=1;mask([3:masknum*2:M,4:masknum*2:M],134:N)=1;                     % 33.3%

% masknum=4;mask([1:masknum*2:M,8:masknum*2:M],1:50)=1;mask([2:masknum*2:M,7:masknum*2:M],51:100)=1;
% mask([3:masknum*2:M,6:masknum*2:M],101:150)=1;mask([4:masknum*2:M,5:masknum*2:M],151:N)=1;                    % 25%


sarData=sarData.*mask;                                       % mask
figure;imagesc(abs(sarData));axis square;colormap('jet');
title 'Sparse SAR Echo';

percent=nnz(sarData)/numel(sarData);                         % percent
fprintf('percent=%.6f\n',percent);
maxMod=max(abs(sarData(:)));sarData=sarData/maxMod;          % normalization

%% matrix completion(TCTA/TVNNA)

if IS_TCTA
    mu0=2;                                               % ADMM penalty parameter
    e_rank=20;                                           % e_rank
    P=e_rank/2;                                          % pencil parameter
    Q=P;                                                 % pencil parameter
    maxIter_K=30;                                        % iterations
    sarData=TCTA(sarData,P,Q,mu0,maxIter_K,e_rank);      % toeplitz-column-toeplitz ADMM (TCTA)
    figure;imagesc(abs(sarData));axis square;colormap('jet');
    title 'Reconstructed SAR Echo After MC(TCTA)';
end

if IS_TVNNA
    lambda=1/percent;                                    % smoothness parameter
    rho=1;                                               % ADMM penalty parameter
    mu1=1.25;                                            % penalty parameter multiplier
    maxIter=100;                                         % iterations
    tol=1e-6;                                            % tolerance
    sarData=TVNNA(sarData,lambda,rho,mu1,maxIter,tol);   % Total Variation and Nuclear Norm ADMM (TVNNA)
    figure;imagesc(abs(sarData));axis square;colormap('jet');
    title 'Reconstructed SAR Echo After MC(TVNNA)';
end

%% rma 2d imaging

RMA_2D(dx,dy,k,z0,sarData,nFFTspace,amplitude,imSize);   % imaging of one range slice

ReconstrcutAndImagingTime = toc;
fprintf('Reconstruction and imaging time=%.6f\n',ReconstrcutAndImagingTime);

