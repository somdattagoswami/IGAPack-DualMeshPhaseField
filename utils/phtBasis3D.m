function [N,dN]=phtBasis3D(u_hat,v_hat,w_hat,Ce,p,q,r)

[B_u,dB_u] = bernstein_basis(u_hat,p);
[B_v,dB_v] = bernstein_basis(v_hat,q);
[B_w,dB_w] = bernstein_basis(w_hat,r);

B = zeros(1, (p+1)*(q+1)*(r+1));
dBxi = zeros(1, (p+1)*(q+1)*(r+1));
dBeta = zeros(1, (p+1)*(q+1)*(r+1));
dBzeta = zeros(1, (p+1)*(q+1)*(r+1));

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

N = B;
dN = [dBxi; dBeta; dBzeta];
end
