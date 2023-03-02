N = 10000; % sample
sigma_w = 0.1;
sigma_v = 1;
% true plant:
A = 1;
B = 0;
C = 1;
R1 = sigma_w^2; % covariance of process noise
R2 = sigma_v^2; % covariance of measurement noise
% Pm(:, k) = P(N|N-1)
Pm(:, 1) = 1e5 * eye(1); % P(0|-1) = 100000 P(1|0) = 100000
% P: covariance of state(uncertainty of an estimate)
% Kalman Gain:
for k = 1: N
    Kf(:, k) = Pm(:, k) * C' * (C * Pm(:, k) * C' + R2)^(-1); % Kf(:, k) = Kf(k)  equation 1.52
    K(:, k) = (A*Pm(:,k)*C')*(C*Pm(:,k)*C'+R2)^(-1); % K(:, k) = K(k)  equation 1.53
    P(:, k) = Pm(:, k) - (Pm(:, k)*C')*(C*Pm(:, k)*C'+R2)^(-1)*C*Pm(:, k); % P(:, k) = P(N|N)  equation 1.56
    Pm(:, k+1) = A*Pm(:, k)*A' + R1 - K(:,k)*(C*Pm(:, k)*C'+R2)*K(:,k)'; % Pm(:, k+1) = P(N+1|N)  equation 1.57
end
% Initialization:
w = random('normal', 0, sigma_w, 1, N);
v = random('normal', 0, sigma_v, 1, N);
xhm = [0]'; % xhm(:, 1) = x(0|-1) or x(1|0)
% x(:, k) = x(N|N)
x(:, 1) = [5]; % x(:, 1) = x(0) or x(1) (matlab's matrix don't have 0 subscript)
for k = 1: N
    x(:, k+1) = A*x(:, k) + w(k); %  x true state
    y(k) = C*x(:,k)+v(k); % measurement
    % Kalman filter:
    xh(:, k) = xhm(:, k) + Kf(:, k) * (y(k) - C * xhm(:, k)); % xh(:, k) = x^hat(k|k)
    xhm(:, k+1) = A*xhm(:, k) + K(:,k)*(y(k) - C*xhm(:, k)); % xhm(:, k+1) = x^hat(k+1|k)
end

% calculate bias:
sum_bias = 0;
for i=1: N
    sum_bias = sum_bias + (x(1, i) - xh(1, i));
end
bias = sum_bias / N;

% calculate variance:
sum_variance = 0;
for i=1: N
    sum_variance = sum_variance + (x(1, i) - xh(1, i))^2;
end
variance = sum_variance / N;
