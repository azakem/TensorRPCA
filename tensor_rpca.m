function [X, L, S, iters] = tensor_rpca() % set up video reader for sample video
    v = VideoReader('sample_video_144.mp4');

    % create empty tensor of dimensions 
    % video_height*3 x video_width x % video_frames
    % having 1:video_height = 1st pixel color, video_height+1:2*video_height
    % 2nd pixel color, 2*video_height+1:3*video_height for 3rd pixel color
    % using matlab matrices; tensor toolkit does not support fft
    n = v.Width;
    m = v.Height;
    time = v.Duration;
    nFrames = ceil(v.FrameRate*time);
    X = zeros(m*3, n, nFrames);

    % fill in video frames
    i = 1;
    while hasFrame(v) && i <= nFrames
        vid = readFrame(v);
        X(1:m,:,i) = vid(:,:,1);
        X(m+1:2*m,:,i) = vid(:,:,2);
        X(2*m+1:3*m,:,i) = vid(:,:,3);
        i = i + 1;
    end
    
    % set ADMM parameters
    froNormX = tensorFroNorm(X);
    epsilon = 5e-5*froNormX;
    lambda = 1 / sqrt(max(m*3, n)); % per Zhang paper
    tau = 1; 
    rho = 1; 
    maxiter = 100;
    L = zeros(m*3, n, nFrames);
    % initialize columns of L to column means for each frame
    for i=1:nFrames
       d1mean = mean(X(1:m,:,i));
       for j=1:m
           L(j,:,i) = d1mean;
       end
       d2mean = mean(X(m+1:m*2,:,i));
       for j=1:m
           L(m+j,:,i) = d2mean; 
       end
       d3mean = mean(X(m*2+1:m*3,:,i));
       for j=1:m
           L(m*2+j,:,i) = d3mean;
       end
    end
    % set S to X - L for initial guess
    S = X - L;
    Y = zeros(m*3,n,nFrames);
    
    % transform M to fourier space along 3rd dimension for L update
    Xhat = fft(X,[],3);

    % ADMM based on Zhang paper
    for k=1:maxiter
        iters = k;
        sprintf('Iteration: %k', k);
        oldL = L;
        oldS = S;
        L = L_Update(tau, Xhat,S,rho*Y);
        S = S_Update(lambda, rho, S);
        Y = rho*Y + L + S - X;

        diff = abs(tensorFroNorm(L - oldL)) + abs(tensorFroNorm(S - oldS));
        sprintf('Change: %f', diff);

        % if (diff < epsilon)
        %    break;
        % end
    end
end

function result = L_Update(tau, Mhat, S, Y)
    % first, transform S and Y to fourier space
    Shat = fft(S, [], 3);
    Yhat = fft(Y, [], 3);
    result = zeros(size(S,1), size(S,2), size(S,3));
    
    % calculate SVD for each ith frontal slice Mhat(i) - Shat(i) - What(i)
    % threshold diagonal of Sigma by tau
    % set ith frontal slice of result 
    for i=1:size(S,3)
        [U, Sigma, V] = svd(Mhat(:,:,i) - Shat(:,:,i) - Yhat(:,:,i));
        for j=1:min(size(Sigma,1), size(Sigma,2))
            Sigma(j,j) = subplus(Sigma(j,j) - tau); 
        end
        result(:,:,i) = U * Sigma * V';
    end
    % transform to original space
    result = ifft(result, [], 3);
end

function result = S_Update(lambda, rho, S)
    result = zeros(size(S,1), size(S,2), size(S,3));
   % determine result for each tube of S
    for i=1:size(S,1)
        for j=1:size(S,2)
            fnorm = 0;
            for k=1:size(S,3)
                fnorm = fnorm + abs(S(i,j,k))^2;
            end
            fnorm = sqrt(fnorm);
            inner = 1 - (lambda / (rho * fnorm));
            result(i,j,:) = subplus(inner) * S(i,j,:);
        end
    end

%     % determine result for each frontal slice
%     for i=1:size(S,3)
%         inner = 1 - (lambda / (rho * norm(S(:,:,i), 'fro')));
%         result(:,:,i) = subplus(inner) * S(:,:,i);
%     end
end
   


