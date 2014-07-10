% This function uses the parameters for (Fisher's) quadratic discriminant analysis model
% to classify a given dataset of features into the modelled classes.
% Usage: ld = qdalda_classify(featvec, mu0, mu1, C0, C1)
% featvec - input feature matrix of n features (rows) and m observations (columns)
% mu0     - means matrix with n feature rows (class 0)
% mu1     - means matrix with n feature rows (class 1)
% C0      - square n x n covariance matrix for class labels 0
% C1      - square n x n covariance matrix for class labels 1
% l       - log-likelihood ratio values vector between classes 0 and 1 for each observation in featvec
% ld      - class assignments, one per column of featvec
% (c) 2006 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, I.M. Moroz, S.J. Roberts (2006)
% Nonlinear, Biophysically-Informed Speech Pathology Detection
% in Proceedings of IEEE ICASSP 2006, IEEE Publishers: Toulouse, France.
function [l, ld] = qdalda_classify(featvec, mu0, mu1, C0, C1)

Nobs = size(featvec, 2);

C0i = inv(C0);
C1i = inv(C1);

% Calculate conic section boundary from class model parameters
w0 = -0.5*log(det(C0)) + 0.5*log(det(C1)) - 0.5*mu0'*C0i*mu0 + 0.5*mu1'*C1i*mu1;
w1 = mu0'*C0i - mu1'*C1i;
w2 = C0i - C1i;

% Calculate log-likelihood ratio for each observation
for j=1:Nobs
   l(j) = -0.5*featvec(:,j)'*w2*featvec(:,j) + w1*featvec(:,j) + w0;
end

% Apply threshold to find class label for each observation
ld = (l < 0);
