function [Ki,stdKi,int,stdint] = Patlak(TAC,time_TAC,Cp,Cpint,weights,nrpoints)
% Patlack graphical method
%
% Expects
%--------------------------------------------------------------------------
% TAC                       % tissue TAC
% time_TAC                  % time of TAC
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

N = length(time_TAC);

X = [Cpint(N-nrpoints+1:N)./Cp(N-nrpoints+1:N) ones(1,nrpoints)];
Y = TAC(N-nrpoints+1:N)./Cp(N-nrpoints+1:N);
W = diag(1./weights(N-nrpoints+1 : N).^2);

% estimated_par, estimated_standev, gamma (or sigmaq), covariance matrix
[P,StDev] = lscov(X,Y,W);

Ki    = P(1);
stdKi = StDev(1);

if nargout > 2
    int    = P(2);
    stdint = StDev(2);
end

% % Graphic
% figure(3)
% plot(X(:,1),Y,'or')
% hold on
% plot(X(:,1),(P(1)*X(:,1) + P(2)),'.-b')
% hold off
   
    
    


