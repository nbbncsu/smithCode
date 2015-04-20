function Weights = gauss_weights(N,h)

% determine the Gauss weights for a four point quadrature rule

weights = zeros(4,1);
weights(1) = 49*h/(12*(18 + sqrt(30)));
weights(2) = 49*h/(12*(18 - sqrt(30)));
weights(3) = weights(2);
weights(4) = weights(1);

% copy the weights to form a vector for all N intervals

Weights = weights;
for gct = 1:N-1
    Weights = [Weights; weights];
end;
