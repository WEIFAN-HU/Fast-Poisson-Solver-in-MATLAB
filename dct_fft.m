%% DCT-I:IV
% type = 1 : DCT-I
%      = 2 : DCT-II
%      = 3 : DCT-III
%      = 4 : DCT-IV
% dim  = 1 : columnwise transform
% dim  = 2 : rowwise transform
function Y = dct_fft(X,type,dim)

    if nargin < 2, type = 2; end
    if nargin < 3, dim = 1; end

    if dim ~= 1 && dim ~= 2
        error('dim must be 1 or 2.');
    end

    if ~isreal(X)
        Y = dct_fft_real(real(X),type,dim) ...
          + 1i*dct_fft_real(imag(X),type,dim);
    else
        Y = dct_fft_real(X,type,dim);
    end

end

function Y = dct_fft_real(X,type,dim)

    switch type
        case 1  % DCT-I

            if dim == 1
                N = size(X,1);
                if N == 1
                   Y = X;
                   return;
                end

                Z = [X; flipud(X(2:N-1,:))];

                F = fft(Z,[],1);

                Y = F(1:N,:);

                % if ~isreal(Y)
                %    Y = real(Y);
                % end
            else
                N = size(X,2);
                if N == 1
                   Y = X;
                   return;
                end

                Z = [X, fliplr(X(:,2:N-1))];

                F = fft(Z,[],2);

                Y = F(:,1:N);

                % if ~isreal(Y)
                %    Y = real(Y);
                % end
            end

        case 2  % DCT-II

            if dim == 1
                N = size(X,1);
                Z = [X; flipud(X)];
                F = fft(Z,[],1);

                k = (0:N-1)';
                phase = exp(-1i*pi*k/(2*N));

                Y = real(phase .* F(1:N,:));
            else
                N = size(X,2);
                Z = [X, fliplr(X)];
                F = fft(Z,[],2);

                k = 0:N-1;
                phase = exp(-1i*pi*k/(2*N));

                Y = real(F(:,1:N) .* phase);
            end

        case 3  % DCT-III

            if dim == 1
                N = size(X,1);

                C = X;
                if N > 1
                    C(2:N,:) = 2*C(2:N,:);
                end

                n = (0:N-1)';
                phase = exp(1i*pi*n/(2*N));

                Z = zeros(2*N,size(X,2));
                Z(1:N,:) = phase .* C;

                F = ifft(Z,[],1);

                Y = real(2*N*F(1:N,:));
            else
                N = size(X,2);

                C = X;
                if N > 1
                    C(:,2:N) = 2*C(:,2:N);
                end

                n = 0:N-1;
                phase = exp(1i*pi*n/(2*N));

                Z = zeros(size(X,1),2*N);
                Z(:,1:N) = C .* phase;

                F = ifft(Z,[],2);

                Y = real(2*N*F(:,1:N));
            end

        case 4  % DCT-IV

            if dim == 1
                N = size(X,1);
                M = size(X,2);
                Z = zeros(8*N,M);
                idx = 2:2:2*N;
                Z(idx,:) = X;
                Z(8*N - idx + 2,:) = X;

                F = fft(Z,[],1);

                kidx = 2:2:2*N;

                Y = real(F(kidx,:));

             elseif dim == 2
                N = size(X,2);
                M = size(X,1);
                Z = zeros(M,8*N);
                idx = 2:2:2*N;
                Z(:,idx) = X;
                Z(:,8*N - idx + 2) = X;

                F = fft(Z,[],2);

                kidx = 2:2:2*N;
                Y = real(F(:,kidx));
            end

        otherwise
            error('type must be 1, 2, 3, or 4.');
    end

end