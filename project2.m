N = 10000; % sample
T= 1;
A = [1 T; 0 1];
C = [1 0];
sigma_w = 0.1;
sigma_v = 1;
R1 = [T^4/4 T^3/2; T^3/2 T^2] * sigma_w^2;
R2 = sigma_v^2; 
Pm(:,:,1) = 1e5*eye(2); %Pm(:,:,k) = P(N|N-1)
% kalman gain and estimate uncertainty:
for k = 1: N
    Kf(:, k) = Pm(:, :, k) * C' * (C * Pm(:, :, k) * C' + R2) ^ (-1);
    K(:, k) = A * Kf(:, k);
    P(:, :, k) = Pm(:, :, k) - Pm(:, :, k) * C' * (C * Pm(:, :, k) * C' + R2) ^ (-1) * C * Pm(:, :, k);
    Pm(:, :, k+1) = A * Pm(:, :, k)*A' - K(:, k)*(C*Pm(:, :,k)*C' + R2)*K(:, k)' + R1;
end
w = random('normal', 0, sigma_w, 1, N);
v = random('normal', 0, sigma_v, 1, N);
x(:, 1) = [0; 30];
xhm(:, 1) = [0; 0];
for k = 1: N
    x(:, k+1) = A*x(:, k) + [T^2/2; T] * w(:, k);
    y(k) = C * x(:, k) + v(:, k);
    xh(:, k) = xhm(:, k) + Kf(:, k)*(y(k) - C*xhm(:, k));
    xhm(:, k+1) = A*xhm(:, k) + K(:, k)*(y(:, k) - C*xhm(:, k));
end
% calculate biases:
%x1:
sum_bias_x1 = 0;
for i=1: N
    sum_bias_x1 = sum_bias_x1 + (x(1, i) - xh(1, i));
end
bias_x1 = sum_bias_x1 / N;

%x2:
sum_bias_x2 = 0;
for i=1: N
    sum_bias_x2 = sum_bias_x2 + (x(2, i) - xh(2, i));
end
bias_x2 = sum_bias_x2 / N;

% calculate variance:
sum_variance_x1 = 0;
for i=1: N
    sum_variance_x1 = sum_variance_x1 + (x(1, i) - xh(1, i))^2;
end
variance_x1 = sum_variance_x1 / N;

%x2:
sum_variance_x2 = 0;
for i=1: N
    sum_variance_x2 = sum_variance_x2 + (x(2, i) - xh(2, i))^2;
end
variance_x2 = sum_variance_x2 / N;