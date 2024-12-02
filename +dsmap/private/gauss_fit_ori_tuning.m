function [tcFit,stat] = gauss_fit_ori_tuning(rawCurv,oriList)
% fit the tuning curve of neurons orientation selectivity with a wraped gaussian funtion
% G(angle) = C + Rp*exp(wrapA(angle-pA)^2/(2/log(4)*W^2))
% where C is constant
% Rp is response at prefered orientation angle pA
% W is tuning width (half-width at half height)
% wrapA warp angle difference from [-180,180] to [0, 90]
% warped gaussian function then becomes a circular function in [0, 180]
% optimize G with MSE using FMINSEARCHBND with constrains
% 1) C within [-F,F], where F is max reponse obseverd in data
% 2) Rp within [0, 3F]
% 3) W within [S, 45], where S is half step of oriList
% 4) start searching with pA = A where R(A)=F, Rp = F, C = 0
% for details, refer to M. Mazurek, et al. 2014 and N.V.Swindale, 1998

% INPUT
% rawCurv[KxT], response of K neurons to T stimuli with different orientations
% oriList[Tx1], orientations tested, ranged in [0,180]

% OUTPUT
% tcFit [Kx180], fited turinng curves at 0 to 179 degree
% stat [Kx1 struct] contains statistics of fitted turining curve

% wirtten by Ruix.Li in Jul, 2021

K = size(rawCurv,1);
N = 0:1:179;
S = (oriList(2) - oriList(1))/2;
tcFit = zeros(K,180);
mse = zeros(K);
tw = zeros(K);
osi = zeros(K);
pa = zeros(K);
rp = zeros(K);

G = @(N,M) M(2)*exp(-(wrpAngle(N-M(1))).^2/(2*M(3)^2/log(4)))+M(4);

ops = optimset('Display','off','TolX',1e-3); 
%%
for i = 1:K
D = @(X) sum((G(oriList,X) - rawCurv(i,:)').^2);
lowBound = [0,0,S,-max(rawCurv(i,:))];
upBound = [180,3*max(rawCurv(i,:)),45,max(rawCurv(i,:))];
initX = [oriList(atmaxof(rawCurv(i,:))),max(rawCurv(i,:)),S,0];
[x,mse(i)] = fminsearchbnd(D,initX,lowBound,upBound,ops);
tcFit(i,:) = G(N,x);
pa(i) = x(1);
tw(i) = x(3);
osi(i) = (G(x(1),x) - G(x(1)+90,x)) / (G(x(1),x) + G(x(1)+90,x));
rp(i) = x(2);
end

if nargout == 2
    for i = 1:K
        stat(i).mse = mse(i);
        stat(i).osi = osi(i);
        stat(i).turningWith = tw(i);
        stat(i).prefOri = pa(i);
        stat(i).prefOriResp = rp(i);
    end
end
end

function B = wrpAngle(A)
A = A(:);
A = min(abs([A,A-180,A+180]),[],2);
A(A>90) = 180 - A(A>90);
B = A;
end

function b = atmaxof(a)
% return positon of max value
% usage, get value of B where A has max
% v = B(atmaxof(A));
[~,b] = max(a(:));
end
