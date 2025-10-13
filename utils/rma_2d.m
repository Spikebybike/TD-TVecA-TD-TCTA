%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function rma_2d for 2D rma imaging
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2025-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rma_2d(dx,dy,k,z0,sarData,nFFTspace,amplitude,imSize)

%% phase creating
    wx = 2*pi/(dx*1e-3);                                    % Sampling frequency for x
    kX = linspace(-(wx/2),(wx/2),nFFTspace);               % kX
    wy = 2*pi/(dy*1e-3);                                    % Sampling frequency for y
    kY = (linspace(-(wy/2),(wy/2),nFFTspace)).';           % kY
    kz = single(sqrt((2*k).^2 - kX.^2 - kY.^2));
    phaseFactor = exp(-1i*z0*kz);
    phaseFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
    phaseFactor1 = kz.*phaseFactor;
    phaseFactor1 = fftshift(fftshift(phaseFactor1,1),2);
    
%% zero-padding
    [yPointM,xPointM] = size(sarData);
    [yPointF,xPointF] = size(phaseFactor1);
    left=zeros(yPointM,floor((xPointF-xPointM)/2));%padding zeros around sardata 0     0     0 from ny nx to nfftspace nfftspace
    right=zeros(yPointM,ceil((xPointF-xPointM)/2));%                             0  sardata  0   
    top=zeros(floor((yPointF-yPointM)/2),xPointF);%                              0     0     0
    bottom=zeros(ceil((yPointF-yPointM)/2),xPointF);
    sarData=cat(2,left,sarData,right);
    sarData=cat(1,top,sarData,bottom);

%% image to center
    yPointT=nFFTspace;
    xPointT=nFFTspace;
    xRangeT = dx * (-(xPointT-1)/2 : (xPointT-1)/2);
    yRangeT = dy * (-(yPointT-1)/2 : (yPointT-1)/2);
    indXpartT = xRangeT>(-imSize/2) & xRangeT<(imSize/2);
    indYpartT = yRangeT>(-imSize/2) & yRangeT<(imSize/2);
    xRangeT = xRangeT(indXpartT);
    yRangeT = yRangeT(indYpartT);

%% rma-2d
    sarDataFFT = fft2(sarData);
    sarImage_2D_RMA = ifft2(sarDataFFT.*phaseFactor1);

%% imaging
    sarImage_2D = abs(sarImage_2D_RMA(indYpartT,indXpartT));
    finalsar=flipud(2*db(sarImage_2D/max(sarImage_2D(:))));
    finalsar=max(min(finalsar, 0), amplitude);
    figure;imagesc(xRangeT,yRangeT,finalsar);
    axis equal xy on;colormap('jet');
    title 'Sar Image-2D RMA'
    xlabel('Horizontal (mm)')
    ylabel('Vertical (mm)')
    
end






