% set up video reader for sample video
% v = VideoReader('sample_video_144.mp4');
% 
% % create empty tensor of dimensions 
% % video_height*3 x video_width x % video_frames
% % having 1:video_height = 1st pixel color, video_height+1:2*video_height
% % 2nd pixel color, 2*video_height+1:3*video_height for 3rd pixel color
% n = v.Width;
% m = v.Height;
% time = v.Duration;
% nFrames = ceil(v.FrameRate*time);
% X = tenzeros([m*3 n nFrames]);
% 
% % fill in video frames
% i = 1;
% while hasFrame(v) && i <= nFrames
%     vid = readFrame(v);
%     X(1:m,1:n,i) = vid(:,:,1);
%     X(m+1:2*m,1:n,i) = vid(:,:,2);
%     X(2*m+1:3*m,1:n,i) = vid(:,:,3);
%     i = i + 1;
% end

m = size(X,1);
n = size(X,2);
nFrames = size(X,3);
L0 = tenzeros([m n nFrames]);
% initialize L0 to mean for each frame
for i=1:nFrames
   L0(1:m,:,i) = mean(X(1:m,:,i));
   L0(m+1:m*2,:,i) = mean(X(m+1:m*2,:,i));
   L0(m*2+1:m*3,:,i) = mean(X(m*2+1:m*3,:,i));    
end
S0 = X - L0;

% set tolerance for difference in solution
epsilon = 5e-5*tensorFroNorm(X);

% ADMM code for PCP from HW2
% need to modify for tensor
% see page 3 of Tensor-SVD paper for calculating nuclear norm of tensor
% see page 7 of Tensor-SVD paper for tensor RPCA
for k=1:maxiter
    iter = k;
    oldL = solL;
    oldS = solS;
    solL = proxNukeNorm(1, solL - (Y + rho*(solL + solS - M)));
    solS = softThresh(lambdaS/rho, M - solL - (1/rho)*Y);
    Y = Y + rho*(solL + solS - M);

    cost = abs((norm(solL + solS, 'fro')-normM));

    if (abs(norm(solL-oldL, 'fro')) + abs(norm(solS-oldS, 'fro')) < epsilon)
        break;
    end
end

