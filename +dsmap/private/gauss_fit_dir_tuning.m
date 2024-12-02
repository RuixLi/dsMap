function [tcFit,stat] = gauss_fit_dir_tuning(rawCurv,dirList)
% fit the tuning curve of neurons direction selectivity with wraped double gaussian funtion
% G(angle) = C + Rp*exp(wrapA(angle-pA)^2/(2/log(4)*W^2)) + Rn*exp(wrapA(angle+180-pA)^2/(2/log(4)*W^2))
% where C is constant
% Rp is response at prefered orientation angle pA
% Rn is response at prefered orientation angle pA
% W is tuning width (half-width at half height)
% wrapA warp angle difference from [-360,360] to [0, 180]
% warped gaussian function then becomes a circular function in [0, 360]
% optimize G with MSE using FMINSEARCHBND with constrains
% 1) C within [-F,F], where F is max reponse obseverd in data
% 2) Rp within [0, 3F]
% 3) W1 and W2 within [S, 45], where S is half step of oriList
% 4) start searching with pA = A where R(A)=F, Rp = F, C = 0
% for details, refer to M. Mazurek, et al. 2014 and N.V.Swindale, 1998

% INPUT
% rawCurv[KxT], response of K neurons to T stimuli with different orientations
% dirList[Tx1], orientations tested, ranged in [0,180]

% OUTPUT
% tcFit [Kx360], fited turinng curves at 0 to 359 degree
% stat [Kx1 struct] contains statistics of fitted turining curve

% wirtten by Ruix.Li in Jul, 2021

K = size(rawCurv,1);
N = 0:1:359;
S = (dirList(2) - dirList(1))/2;
tcFit = zeros(K,360);
mse = zeros(K);
tw = zeros(K);
dsi = zeros(K);
pa = zeros(K);
rp = zeros(K);
rn = zeros(K);

G = @(N,M) M(2)*exp(-(wrpAngle(N-M(1))).^2/(2*M(3)^2/log(4))) + ...
           M(4)*exp(-(wrpAngle(N+180-M(1))).^2/(2*M(3)^2/log(4))) + M(5);

ops = optimset('Display','off','TolX',1e-3);       
%%
for i = 1:K
D = @(X) sum((G(dirList,X) - rawCurv(i,:)').^2);
lowBound = [0,0,S,0,-max(rawCurv(i,:))];
upBound = [360,3*max(rawCurv(i,:)),45,3*max(rawCurv(i,:)),max(rawCurv(i,:))];
initX = [dirList(atmaxof(rawCurv(i,:))),max(rawCurv(i,:)),S,max(rawCurv(i,:)),0];

[x,mse(i)] = fminsearchbnd(D,initX,lowBound,upBound,ops);
tcFit(i,:) = G(N,x);
pa(i) = N(atmaxof(tcFit(i,:)));
tw(i) = x(3);
dsi(i) = (G(pa(i),x) - G(pa(i)+180,x)) / abs(G(pa(i),x) + G(pa(i)+180,x));
rp(i) = x(2);
rn(i) = x(4);
end

if nargout == 2
    for i = 1:K
        stat(i).mse = mse(i);
        stat(i).dsi = dsi(i);
        stat(i).turningWith = tw(i);
        stat(i).prefDir = pa(i);
        stat(i).prefDirResp = rp(i);
        stat(i).nullDirResp = rn(i);
    end
end
end

function B = wrpAngle(A)
A = A(:);
A = min(abs([A,A-360,A+360]),[],2);
A(A>180) = 360 - A(A>180);
B = A;
end

function b = atmaxof(a)
% return positon of max value
% usage, get value of B where A has max
% v = B(atmaxof(A));
[~,b] = max(a(:));
end
