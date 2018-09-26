function rgbImg = hsl2rgb(hslImg)
%%%% input HSL image, with each component double H:[0,360];S:[0,1];L:[0,1]
%%%% output RGB image, with each component uint8 [0,255]

if((ndims(hslImg) ~= 3) || (size(hslImg, 3) ~= 3))
   error('Input image must be HSL.');
end
[M,N] = size(hslImg(:,:,1));
H = hslImg(:,:,1);
S = hslImg(:,:,2);
L = hslImg(:,:,3);
R = zeros(M,N);
G = zeros(M,N);
B = zeros(M,N);
% C = (1-abs(2*L-1)) .* S;
% X = C .* (1-abs(mod((H/60),2)-1));
% m = L - C/2;
for ii = 1:M
    for jj = 1:N
        h = H(ii,jj);
        s = S(ii,jj);
        ll = L(ii,jj);
        c = (1-abs(2*ll-1))*s;
        x = c*(1-abs(mod(h/60,2)-1));
        m = ll - c/2;
        if(h<60)
            R(ii,jj) = c;
            G(ii,jj) = x;
            B(ii,jj) = 0;
        elseif(h<120)
            R(ii,jj) = x;
            G(ii,jj) = c;
            B(ii,jj) = 0;
        elseif(h<180)
            R(ii,jj) = 0;
            G(ii,jj) = c;
            B(ii,jj) = x;
        elseif(h<240)
            R(ii,jj) = 0;
            G(ii,jj) = x;
            B(ii,jj) = c;
        elseif(h<300)
            R(ii,jj) = x;
            G(ii,jj) = 0;
            B(ii,jj) = c;
        else%if(h<360)
            R(ii,jj) = c;
            G(ii,jj) = 0;
            B(ii,jj) = x;
        end
        R(ii,jj) = R(ii,jj) + m;
        G(ii,jj) = G(ii,jj) + m;
        B(ii,jj) = B(ii,jj) + m;
        %disp([R(ii,jj),G(ii,jj),B(ii,jj)]);
    end
end
rgbImg = cat(3,R,G,B);
rgbImg = 255*rgbImg;
rgbImg(rgbImg>255) = 255;
rgbImg(rgbImg<0) = 0;
rgbImg = uint8(rgbImg);
            