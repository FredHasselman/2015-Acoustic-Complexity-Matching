% This function finds the bootstrap resampled parameters for (Fisher's) quadratic discriminant analysis, finding the
% (hyper) conic sections that minimise the probability of misclassification of the feature vectors into
% the two given classes, given multivariate Gaussian class models.
% Usage: [mu0, mu1, C0, C1] = qdalda_train(featvec, labels, bootN)
% Usage: [mu0, mu1, C0, C1] = qdalda_train(featvec, labels, bootN, lda)
% featvec - input feature matrix of n features (rows) and m observations (columns)
% labels  - row vector of numerical class labels 0 and 1 for each of the m input feature vectors
% bootN   - number of bootstrap resampling trials, set to zero for no bootstrap resampling
% lda     - (optional) set this parameter to 1 in order to perform LDA training rather than QDA (covariance matrices are then identical for each class)
% mu0     - means matrix with n feature rows and bootN columns (class 0)
% mu1     - means matrix with n feature rows and bootN columns (class 1)
% C0      - square n x n x bootN covariance tensor for class labels 0
% C1      - square n x n x bootN covariance tensor for class labels 1
% 
% (c) 2006 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, I.M. Moroz, S.J. Roberts (2006)
% Nonlinear, Biophysically-Informed Speech Pathology Detection
% in Proceedings of IEEE ICASSP 2006, IEEE Publishers: Toulouse, France.
function [mu0, mu1, C0, C1] = qdalda_train(featvec, labels, bootN, varargin)

iclass0 = find(labels == 0);
iclass1 = find(labels == 1);

N = size(featvec,1);

if (length(varargin{:}) > 0)
    lda = 1;
else
    lda = 0;
end

if (bootN > 0)
    C0 = zeros(N,N,bootN);
    C1 = C0;
    mu0 = zeros(N,bootN);
    mu1 = mu0;
    N0 = length(iclass0);
    N1 = length(iclass1);
    for k = 1:bootN
        iselect0  = iclass0(round(rand(N0,1) * (N0-1)) + 1);
        iselect1  = iclass1(round(rand(N1,1) * (N1-1)) + 1);
        if (lda == 1)
            C0(:,:,k) = cov(featvec');
            C1(:,:,k) = C0(:,:,k);
        else
            C0(:,:,k) = cov(featvec(:,iselect0)');
            C1(:,:,k) = cov(featvec(:,iselect1)');
        end
        mu0(:,k)  = mean(featvec(:,iselect0), 2);
        mu1(:,k)  = mean(featvec(:,iselect1), 2);
    end
else
    C0  = cov(featvec(:,iclass0));
    C1  = cov(featvec(:,iclass1));
    mu0 = mean(featvec(:,iclass0), 2);
    mu1 = mean(featvec(:,iclass1), 2);
end
