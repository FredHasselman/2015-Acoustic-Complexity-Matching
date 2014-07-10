% Demonstration of using QDALDA package to classify features from two different classes

% Setup class 0 data points
% x = randn(2, 200);
% x = [0.5 0.1; 0.6 -0.3]*x + 1;

bak = [round(data(:,2)); round(data(:,5)); round(data(:,8)); round(data(:,11))];

%bakd= [round(data2(:,2)); round(data2(:,5)); round(data2(:,8)); round(data2(:,11))];
% 
% cnt1=0;
% cnt2=0;
% for i=1:length(bak)
%    if bak(i)==0
%     cnt1=cnt1+1;
%     x1(1,cnt1)=Hn(i);
%     x1(2,cnt1)=RPDE(i);
%    else
%     cnt2=cnt2+1;
%     y1(1,cnt2)=Hn(i);
%     y1(2,cnt2)=RPDE(i);
%    end
% end

cnt1=0;
cnt2=0;
for i=1:length(bak)
   if bak(i)==0
    cnt1=cnt1+1;
    x2(1,cnt1)=HNR(i);
    x2(2,cnt1)=NOISY(i);
   else
    cnt2=cnt2+1;
    y2(1,cnt2)=HNR(i);
    y2(2,cnt2)=NOISY(i);
   end
end



% 
% cnt1=0;
% cnt2=0;
% for i=1:length(bakd)
%    if bakd(i)==0
%     cnt1=cnt1+1;
%     x2(1,cnt1)=Hn(i);
%     x2(2,cnt1)=RPDE(i);
%    else
%     cnt2=cnt2+1;
%     y2(1,cnt2)=Hn(i);
%     y2(2,cnt2)=RPDE(i);
%    end
% end

% class = [bak bakd];
% x = [x1 x2];
% y = [y1 y2];

% class = [bak];
% x = [x1];
% y = [y1];


class = [bak];
x = [x2];
y = [y2];

% Setup class 1 data points
% y = randn(2, 200);
% y = [0.85 0.15; -0.15 0.85]*y - 1;

% Augment into feature vector
featvec = [x y];

% Setup class labels
labels = zeros(length(class), 1);
ind = find(class);
labels(ind) = 1;

% Perform bootstrapped resampled training and testing using QDA
[qda_perf, qda_conf, qda_mu0E, qda_mu1E, qda_C0E, qda_C1E] = ...
    qdalda_traintest(featvec, labels, 100);

disp(sprintf('QDA: Class 0 correct: %3.1f%%, class 1 correct: %3.1f%%, total correct: %3.1f%%', ...
    qda_perf(1)*100, qda_perf(2)*100, qda_perf(3)*100));

% % Perform bootstrapped resampled training and testing using LDA
% [lda_perf, lda_conf, lda_mu0E, lda_mu1E, lda_C0E, lda_C1E] = ...
%     qdalda_traintest(featvec, labels, 100, 1);
% 
% disp(sprintf('LDA: Class 0 correct: %3.1f%%, class 1 correct: %3.1f%%, total correct: %3.1f%%', ...
%     lda_perf(1)*100, lda_perf(2)*100, lda_perf(3)*100));

% Plot data and classification boundaries on figures
figure;
%subplot(1,2,1);
hold on;
qdalda_plotdb(featvec, labels, 40, qda_mu0E, qda_mu1E, qda_C0E, qda_C1E);

% subplot(1,2,2);
% hold on;
% qdalda_plotdb(featvec, labels, 40, lda_mu0E, lda_mu1E, lda_C0E, lda_C1E);
