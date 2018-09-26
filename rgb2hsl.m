function hslImg = rgb2hsl(rgbImg)
%%%% input RGB image, with each component uint8 [0,255]
%%%% output HSL image, with each component double H:[0,360];S:[0,1];L:[0,1]
if((ndims(rgbImg) ~= 3) || (size(rgbImg, 3) ~= 3))
   error('Input image must be RGB.');
end
f = double(rgbImg);
[M,N] = size(f(:,:,3));
R = f(:,:,1)/255;
G = f(:,:,2)/255;
B = f(:,:,3)/255;
H = zeros(M,N);
S = zeros(M,N);
L = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        r = R(ii,jj);
        g = G(ii,jj);
        b = B(ii,jj);
        cmax = max([r,g,b]);
        cmin = min([r,g,b]);
        del = cmax - cmin;
        if(del==0)
            H(ii,jj) = 0;
        elseif(cmax==r)
            H(ii,jj) = 60*(mod((g-b)/del,6));
        elseif(cmax==g)
            H(ii,jj) = 60*((b-r)/del + 2);
        elseif(cmax==b)
            H(ii,jj) = 60*((r-g)/del + 4);
        end
        ll = (cmax+cmin)/2;
        L(ii,jj) = ll;
        if(del==0)
            S(ii,jj) = 0;
        else
            S(ii,jj) = del/(1-abs(2*ll-1));
        end
    end
end
% H(H<0) = 0; H(H>360) = 360;
% S(S<0) = 0; S(S>1) = 1;
% L(L<0) = 0; L(L>1) = 1;
hslImg = cat(3,H,S,L);