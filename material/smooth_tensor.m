function tensor_smoothed = smooth_tensor(tensor1, tensor2, rvol1, n12)

% For performance, the error checking code is removed.
% chkarg(istypesizeof(tensor1, 'complex', [Axis.count Axis.count]), ...
% 	'"tensor1" should be %d-by-%d array with complex elements.', Axis.count, Axis.count);
% chkarg(istypesizeof(tensor2, 'complex', [Axis.count Axis.count]), ...
% 	'"tensor2" should be %d-by-%d array with complex elements.', Axis.count, Axis.count);
% chkarg(istypesizeof(rvol1, 'real') && rvol1 >= 0 && rvol1 <= 1, ...
% 	'"rvol1" should be real scalar between 0 and 1.');
% chkarg(istypesizeof(n12, 'real', [1 Axis.count]), ...
% 	'"n12" should be length-%d row vector with real elements.', Axis.count);

n = n12(:) ./ norm(n12);
if any(n == 0)
	h = double(n == 0);
else
	h = [1; 0; 0];
end

h = cross(n, h);
h = h ./ norm(h);

v = cross(n, h);
v = v ./ norm(v);

S = [n, h, v];

tau1 = tau(S.' * tensor1 * S);  % express tensor1 in S coordinates, then transform by tau
tau2 = tau(S.' * tensor2 * S);  % express tensor2 in S coordinates, then transform by tau

tau_mean = tau1 .* rvol1 + tau2 .* (1-rvol1);

tensor_smoothed = S * tau_inv(tau_mean) * S.';

function t = tau(s)
t = zeros(3);
t(2:3, 2:3) = s(2:3, 2:3);
s11 = s(1,1);
s(1,1) = -1;

t = t - (s(:,1) * s(1,:)) ./ s11;

function s = tau_inv(t)
s = zeros(3);
s(2:3, 2:3) = t(2:3, 2:3);
t11 = t(1,1);
t(1,1) = 1;

s = s - (t(:,1) * t(1,:)) ./ t11;
	