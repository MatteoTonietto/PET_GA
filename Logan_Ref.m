function [DVR,stdDVR,int,stdint] = Logan_Ref(TACs,time_PET,ref,frame_length,weights,nrpoints)
% Patlack graphical method
%
% Expects
%--------------------------------------------------------------------------
% TACs                      % tissue TAC(s)
% time_PET                  % time of TACs
% ref                       % reference region TAC
% frame_length              % length of PET frame
% weights                   % weight defined as 1./SD
% nrpoints                  % number of points to use (the last nrpoints)
%
% Evaluates:
%--------------------------------------------------------------------------
% DVR                       % Distribution volume ratio
% stdDVR                    % standard error DVR
% int                       % Logan intercept
% stdint                    % standard error int
%
%__________________________________________________________________________
% Matteo Tonietto

nTAC = size(TACs,2);

N   = length(time_PET);
idx = N - nrpoints + 1 : N;

refint = cumsum(ref.*frame_length);
TACint = cumsum(TACs.*repmat(frame_length,1,nTAC));

DVR        = zeros(1,nTAC);
stdDVR     = zeros(1,nTAC);
if nargout > 2
    int    = zeros(1,nTAC);
    stdint = zeros(1,nTAC);
end

for i = 1 : nTAC
    Y = TACint(idx,i)./TACs(idx,i);
    X = [refint(idx)./TACs(idx,i), ones(nrpoints,1)];
    W = diag(1./weights(idx,i).^2);

    % estimated_par, estimated_standev, gamma (or sigmaq), covariance matrix
    [P,StDev] = lscov(X,Y,W);

    DVR(i)    = P(1,:);
    stdDVR(i) = StDev(1,:);

    if nargout > 2
        int(i)    = P(2,:);
        stdint(i) = StDev(2,:);
    end
    
    % Graphic
%     figure
%     plot(refint(:)./TAC(:,i),TACint(:,i)./TAC(:,i),'or')
%     hold on
%     plot(X(:,1),(P(1)*X(:,1) + P(2)),'.-b')
%     hold off
%     pause
    
end


   