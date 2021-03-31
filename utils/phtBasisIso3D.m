function [N,dN]=phtBasisIso3D(u_hat,v_hat,w_hat,Ce,p,q,r,wgts)

%calculate the shape function and first derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_hat - evaluation point in v-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p,q - polynomial degrees
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the shape functions in %       
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[B_u,dB_u] = bernstein_basis(u_hat,p);
[B_v,dB_v] = bernstein_basis(v_hat,q);
[B_w,dB_w] = bernstein_basis(w_hat,r);

B = zeros(1,(p+1)*(q+1)*(r+1));
dBxi = zeros(1,(p+1)*(q+1)*(r+1));
dBeta = zeros(1,(p+1)*(q+1)*(r+1));
dBzeta = zeros(1,(p+1)*(q+1)*(r+1));

basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;            
            %  Bernstein basis functions;
            B(basisCounter) = B_u(i)*B_v(j)*B_w(k);
            dBxi(basisCounter) = dB_u(i)*B_v(j)*B_w(k);
            dBeta(basisCounter) = B_u(i)*dB_v(j)*B_w(k);
            dBzeta(basisCounter) = B_u(i)*B_v(j)*dB_w(k);
        end
    end
end

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';
dBeta = dBeta*Ce';
dBzeta = dBzeta*Ce';

B = B.*wgts';
dBxi = dBxi.*wgts';
dBeta = dBeta.*wgts';
dBzeta = dBzeta.*wgts';

% Compute the sums of B-spline functions
w_sum = sum(B);
dw_xi = sum(dBxi);
dw_eta = sum(dBeta);
dw_zeta = sum(dBzeta);

dBxi = dBxi/w_sum - B*dw_xi/w_sum^2;
dBeta = dBeta/w_sum - B*dw_eta/w_sum^2;
dBzeta = dBzeta/w_sum - B*dw_zeta/w_sum^2;

N = B/w_sum;
dN = [dBxi;dBeta;dBzeta];
end
