function [Ki,stdKi,int,stdint] = Patlak_modified(TAC,time_PET,Cp,Cpint,weights,nrpoints)
% Patlack graphical method (modified to remove the ratio of Cp) 
%
% Expects
%--------------------------------------------------------------------------
% TAC                       % tissue TAC
% time_PET                  % time of TAC
% Cp                        % input function at times time_TAC
% Cpint                     % cumulative integral of Cp at times time_TAC
% weights                   % weight defined as 1./SD
% nrpoints                  % number of points to use (the last nrpoints)
%
% Evaluates:
%--------------------------------------------------------------------------
% Ki                        % Net influx rate Ki
% stdKi                     % standard error Ki
% int                       % Patlak intercept
% stdint                    % standard error int
%
%__________________________________________________________________________
% Matteo Tonietto

N   = length(time_PET);
idx = N - nrpoints + 1 : N;

X = [Cpint(idx) Cp(idx)];
Y = TAC(idx,:);
W = diag(1./weights(idx).^2);

% estimated_par, estimated_standev, gamma (or sigmaq), covariance matrix
[P,StDev] = lscov(X,Y,W);

Ki    = P(1,:);
stdKi = StDev(1,:);

if nargout > 2
    int    = P(2,:);
    stdint = StDev(2,:);
end

% % Graphic Patlak
% figure(1)
% plot(X(:,1)./Cp(idx),Y./Cp(idx),'or')
% hold on
% plot(X(:,1)./Cp(idx),(P(1)*X(:,1)./Cp(idx) + P(2)),'.-b')
% hold off
   
    
    


