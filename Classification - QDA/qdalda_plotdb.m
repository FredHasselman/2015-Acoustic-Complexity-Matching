% This function uses the parameters for (Fisher's) quadratic discriminant analysis model
% to plot the decision boundary (for two features) between the modelled classes on the current figure,
% and also to plot the given feature data with different class label markers.
% Usage: qdalda_plotdb(featvec, labels, intvs, mu0, mu1, C0, C1)
% featvec  - input feature matrix of 2 features (rows) and m observations (columns)
% labels   - row vector of numerical class labels 0 and 1 for each of the m input observation feature vectors
% intvs    - number of subdivisions of the axes
% mu0      - means vector with 2 feature rows (class labels 0)
% mu1      - means vector with 2 feature rows (class labels 1)
% C0       - square 2 x 2 covariance matrix for class labels 0
% C1       - square 2 x 2 covariance matrix for class labels 1
% segments - the contour matrix of line segments for further processing
% (c) 2006 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, I.M. Moroz, S.J. Roberts (2006)
% Nonlinear, Biophysically-Informed Speech Pathology Detection
% in Proceedings of IEEE ICASSP 2006, IEEE Publishers: Toulouse, France.
function segments = qdalda_plotdb(featvec, labels, intvs, mu0, mu1, C0, C1)

minx = min(featvec(1,:));
maxx = max(featvec(1,:));
miny = min(featvec(2,:));
maxy = max(featvec(2,:));

x       = linspace(minx, maxx, intvs);
y       = linspace(miny, maxy, intvs);
[X Y]   = meshgrid(x, y);
[z, ld] = qdalda_classify([X(:)'; Y(:)'], mu0, mu1, C0, C1);
z       = reshape(z, intvs, intvs);
segments = contour(x, y, z, [0 0], 'k-');

iclass0 = find(labels == 0);
iclass1 = find(labels == 1);
plot(featvec(1,iclass0), featvec(2,iclass0), 'k.', 'MarkerSize', 12);
plot(featvec(1,iclass1), featvec(2,iclass1), 'bx', 'MarkerSize', 5);

axis([minx maxx miny maxy]);
