function [Xlag,Ytrim] = build_lagged_run(X,Y,lags)
% Positive lag means earlier GM predicts later WM: predictor X(t-lag).
max_abs=max(abs(lags));
valid=(1+max_abs):(size(X,1)-max_abs);
Xlag=zeros(numel(valid),size(X,2)*numel(lags),'single');
for j=1:numel(lags)
    source=valid-lags(j);
    cols=(j-1)*size(X,2)+(1:size(X,2));
    Xlag(:,cols)=X(source,:);
end
Ytrim=Y(valid,:);
end
