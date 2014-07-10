% This function finds the bootstrap resampled parameters for (Fisher's) quadratic discriminant analysis,
% of a feature vector, and then tests the classification performance on that feature vector, returning
% the mean confusion matrix and 95% confidence intervals for that matrix.
% Usage: [perform, confidence, mu0E, mu1E, C0E, C1E] = qdalda_traintest(featvec, labels, bootN)
% Usage: [perform, confidence, mu0E, mu1E, C0E, C1E] = qdalda_traintest(featvec, labels, bootN, lda)
% featvec    - input feature matrix of n features (rows) and m observations (columns)
% labels     - row vector of numerical class labels 0 and 1 for each of the m input feature vectors
% bootN      - number of bootstrap resampling trials, set to zero for no bootstrap resampling
% lda        - (optional) set this parameter to 1 to perform LDA rather than QDA (class covariance matrices will be identical)
% perform    - mean performance vector with correct class 0 (row 1), correct class 1 (row 2) and correct overall (row 3) classification fractions
% confidence - 95% confidence intervals for above performance vector, assuming a Gaussian classifier performance model
% mu0E       - average feature vector mean over all bootstrap trials for class 0
% mu01       - average feature vector mean over all bootstrap trials for class 1
% C0E        - average feature vector covariance matrix over all bootstrap trials for class 0
% C1E        - average feature vector covariance matrix over all bootstrap trials for class 1
% (c) 2006 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, I.M. Moroz, S.J. Roberts (2006)
% Nonlinear, Biophysically-Informed Speech Pathology Detection
% in Proceedings of IEEE ICASSP 2006, IEEE Publishers: Toulouse, France.
function [perform, confidence, mu0E, mu1E, C0E, C1E] = qdalda_traintest(featvec, labels, bootN, varargin)

% Return bootstrap trained classifier model parameters
[mu0, mu1, C0, C1] = qdalda_train(featvec, labels, bootN, varargin);

iclass0 = find(labels == 0);
iclass1 = find(labels == 1);
N0   = length(iclass0);
N1   = length(iclass1);
Nobs = length(labels);

% For each bootstrap trial, test the classifier performance against known class labels
if (bootN > 0)
    for k = 1:bootN
        [l, ld] = qdalda_classify(featvec, mu0(:,k), mu1(:,k), C0(:,:,k), C1(:,:,k));
        class_correct = (ld == labels');
        tru0(k)  = sum(class_correct(iclass0))/N0;
        tru1(k)  = sum(class_correct(iclass1))/N1;
        overall(k) = sum(class_correct)/Nobs;
    end

    perform(1) = mean(tru0);
    perform(2) = mean(tru1);
    perform(3) = mean(overall);
    
    confidence(1) = 1.96*std(tru0);
    confidence(2) = 1.96*std(tru1);
    confidence(3) = 1.96*std(overall);
    
    mu0E = mean(mu0,2);
    mu1E = mean(mu1,2);
    C0E  = mean(C0,3);
    C1E  = mean(C1,3);
else
    ld = qda_classify(featvec, mu0, mu1, C0, C1);
    class_correct = (ld == labels);
    perform(1) = sum(class_correct(iclass0))/N0;
    perform(2) = sum(class_correct(iclass1))/N1;
    perform(3) = sum(class_correct)/Nobs;
    
    confidence = zeros(3,1);
    
    mu0E = mu0;
    mu1E = mu1;
    C0E  = C0;
    C1E  = C1;
end
