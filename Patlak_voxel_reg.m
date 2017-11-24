function [mapKi,mapstdKi] = Patlak_voxel_reg(dynPET,mask,time_TAC,infoCp,formatCp,nrpoints)

% TAC processing
TAC     = conv4Dto2D(dynPET,mask);

meanTAC = mean(TAC,2);
frames  = mid2frames(time_TAC);
delta   = frames(:,2) - frames(:,1);
weights = correct_weights(sqrt(delta./abs(meanTAC)));

% Input Function processing
if strcmp(formatCp,'modelled')
    Cp_tTAC    = infoCp.FUN(infoCp.par,infoCp.fixed_par,time_TAC);
    
    Cpint_tTAC = zeros(size(time_TAC));
    for j = 1 : length(time_TAC)
        Cpint_tTAC(j) = integral(@(t)infoCp.FUN(infoCp.par,infoCp.fixed_par,t),0,time_TAC(j),'ArrayValued',true);
    end
    
elseif strcmp(formatCp,'measured')
    indtplasma = find(infoCp.tCp>0);
    tCp        = [0;infoCp.tCp(indtplasma)];
    Cp         = [0;infoCp.Cp(indtplasma)];
    Cp_tTAC    = interp1(tCp,Cp,time_TAC,'linear','extrap');
    
    tv         = [0;logspace(-4,log10(time_TAC(end)),1000)'];
    Cp_tv      = interp1(tCp,Cp,tv);
    Cpint_tv   = cumtrapz(tv,Cp_tv);
    Cpint_tTAC = interp1(tv,Cpint_tv,time_TAC);
    
else
    error('Input function format not recognized')
end

[nT,nTAC] = size(TAC);
TAC_1 = zeros(size(TAC));
for i = 1 : nTAC
    TAC_1(:,i) = (meanTAC*nTAC - TAC(:,i))/(nTAC-1);
end

% Patlak mean
[Ki_all,stdKi_all] = Patlak_modified(meanTAC,time_TAC,Cp_tTAC,Cpint_tTAC,weights,nrpoints);

% Patlak mean -1
[Ki_1,stdKi_1] = Patlak_modified(TAC_1,time_TAC,Cp_tTAC,Cpint_tTAC,weights,nrpoints);

% Patlak
Ki = zeros(size(Ki_1));
for i = 1 : nTAC
    Ki(1,i) = nTAC*Ki_all - (nTAC-1)*Ki_1(i);
end

% Results
idx      = find(mask);
mapKi    = zeros(size(mask));
mapstdKi = mapKi;

mapKi(idx)    = Ki;
mapstdKi(idx) = stdKi_1;
