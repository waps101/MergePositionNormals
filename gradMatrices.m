function [ Dx,Dy,mask,Scentral ] = gradMatrices( mask,method )
%GRADMATRICES Compute matrices to evaluate gradients of function in mask
%   Method:
%       'SmoothedCentral' - Convolution of central difference kernel with
%       Gaussian smoothing kernel, if not available use central, then
%       forward, then backward differences

rows = size(mask,1);
cols = size(mask,2);

%mask = bwareaopen(mask,2048);

% Pad to avoid boundary problems
mask = pad(mask);

rows = rows+2;
cols = cols+2;

checkermask = false(size(mask));
checkermask(1:2:end,1:2:end)=true;
checkermask(2:2:end,2:2:end)=true;

% The number of usable pixels
npix = sum(mask(:));

invalid = false(rows,cols);

% Build lookup table relating x,y coordinate of valid pixels to index
% position in vectorised representation
count = 0;
indices = zeros(size(mask));
for col=1:cols
    for row=1:rows
        if mask(row,col)
            count=count+1;
            indices(row,col)=count;
        end
    end
end


if strcmp(method,'SmoothedCentral')
    % Preallocate maximum required space
    % This would be if all valid pixels had equations for all 8 neighbours for
    % all possible pairs of images - it will be less in practice
    i = zeros(npix*6,1);
    j = zeros(npix*6,1);
    s = zeros(npix*6,1);
    k=0;
    NumEq=0;
    
    i2 = zeros(npix*3,1);
    j2 = zeros(npix*3,1);
    s2 = zeros(npix*3,1);
    k2=0;
    NumEq2=0;

    % X derivatives
    for col=1:cols
        for row=1:rows
            if mask(row,col)
                if mask(row,col-1) && mask(row,col+1)
                    % Both X neighbours present
                    if mask(row-1,col-1) && mask(row-1,col+1) && mask(row+1,col-1) && mask(row+1,col+1)
                        % All 6 neighbours present - smoothed
                        % central differences
                        NumEq=NumEq+1;
                        % Edge neighbours
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-(4/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col+1); s(k)=(4/12);
                        % Corner neighbours
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col-1); s(k)=-(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col+1); s(k)=(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col-1); s(k)=-(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col+1); s(k)=(1/12);
                        
                        % Add smoothness term because of use of central
                        % differences
                        if checkermask(row,col)
                        NumEq2=NumEq2+1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col-1); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col+1); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col); s2(k2)=2;
                        end
                    else
                        % Central difference
                        NumEq=NumEq+1;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-1/2;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col+1); s(k)=1/2;
                        % Add smoothness term because of use of central
                        % differences
                        if checkermask(row,col)
                        NumEq2=NumEq2+1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col-1); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col+1); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col); s2(k2)=2;
                        end
                    end
                    
                elseif mask(row,col-1)
                    % Only backward in X
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                    
                elseif mask(row,col+1)
                    % Only forward in X
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col+1); s(k)=1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                    
                else
                    invalid(row,col) = true;
                end
            end
        end
        % Finished with a pixel
    end
    i=i(1:k,1);
    j=j(1:k,1);
    s=s(1:k,1);
    Dx = sparse(i,j,s,npix,npix);
    
    % Preallocate maximum required space
    % This would be if all valid pixels had equations for all 8 neighbours for
    % all possible pairs of images - it will be less in practice
    i = zeros(npix*6,1);
    j = zeros(npix*6,1);
    s = zeros(npix*6,1);
    k=0;
    NumEq=0;
    
    
    % Y Derivatives
    for col=1:cols
        for row=1:rows
            if mask(row,col)
                if mask(row-1,col) && mask(row+1,col)
                    % Both Y neighbours present
                    if mask(row-1,col-1) && mask(row-1,col+1) && mask(row+1,col-1) && mask(row+1,col+1)
                        % All 6 neighbours present - smoothed
                        % central differences
                        NumEq=NumEq+1;
                        % Edge neighbours
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col); s(k)=(4/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-(4/12);
                        % Corner neighbours
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col-1); s(k)=(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col+1); s(k)=(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col-1); s(k)=-(1/12);
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col+1); s(k)=-(1/12);
                        % Add smoothness term because of use of central
                        % differences
                        if checkermask(row,col)
                        NumEq2=NumEq2+1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row-1,col); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row+1,col); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col); s2(k2)=2;
                        end
                    else
                        % Central difference
                        NumEq=NumEq+1;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row-1,col); s(k)=1/2;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-1/2;
                        % Add smoothness term because of use of central
                        % differences
                        if checkermask(row,col)
                        NumEq2=NumEq2+1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row-1,col); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row+1,col); s2(k2)=-1;
                        k2=k2+1;
                        i2(k2)=NumEq2; j2(k2)=indices(row,col); s2(k2)=2;
                        end
                    end
                    
                elseif mask(row-1,col)
                    % Only backward in Y
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row-1,col); s(k)=1;
                    
                elseif mask(row+1,col)
                    % Only forward in Y
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-1;
                    
                else
                    invalid(row,col) = true;
                end
            end
        end
        % Finished with a pixel
    end
    i=i(1:k,1);
    j=j(1:k,1);
    s=s(1:k,1);
    Dy = sparse(i,j,s,npix,npix);
    
    i2=i2(1:k2,1);
    j2=j2(1:k2,1);
    s2=s2(1:k2,1);
    Scentral = sparse(i2,j2,s2,NumEq2,npix);

    %elseif strcmp(method,'Central')
elseif strcmp(method,'Backward')
        % Preallocate maximum required space
    % This would be if all valid pixels had equations for all 8 neighbours for
    % all possible pairs of images - it will be less in practice
    i = zeros(npix*2,1);
    j = zeros(npix*2,1);
    s = zeros(npix*2,1);
    k=0;
    NumEq=0;
    
    % X derivatives
    for col=1:cols
        for row=1:rows
            if mask(row,col)
                if mask(row,col-1)
                    % Only backward in X
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                    
                elseif mask(row,col+1)
                    % Only forward in X
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col+1); s(k)=1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                    
                else
                    invalid(row,col) = true;
                end
            end
        end
        % Finished with a pixel
    end
    i=i(1:k,1);
    j=j(1:k,1);
    s=s(1:k,1);
    Dx = sparse(i,j,s,npix,npix);
    
    % Preallocate maximum required space
    % This would be if all valid pixels had equations for all 8 neighbours for
    % all possible pairs of images - it will be less in practice
    i = zeros(npix*2,1);
    j = zeros(npix*2,1);
    s = zeros(npix*2,1);
    k=0;
    NumEq=0;
    
    % Y Derivatives
    for col=1:cols
        for row=1:rows
            if mask(row,col)
                if mask(row-1,col)
                    % Only backward in Y
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row-1,col); s(k)=1;
                    
                elseif mask(row+1,col)
                    % Only forward in Y
                    NumEq=NumEq+1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-1;
                    
                else
                    invalid(row,col) = true;
                end
            end
        end
        % Finished with a pixel
    end
    i=i(1:k,1);
    j=j(1:k,1);
    s=s(1:k,1);
    Dy = sparse(i,j,s,npix,npix);
    
    Scentral = [];
end

% Unpad
mask = mask(2:rows-1,2:cols-1);
invalid = invalid(2:rows-1,2:cols-1);

if sum(invalid(:))>0
    disp(['Removing ' num2str(sum(invalid(:))) ' pixels']);
    [ Dx,Dy,mask,Scentral ] = gradMatrices( mask&~invalid,method );
end

end

function [ X2 ] = pad( X )
%PAD Summary of this function goes here
%   Detailed explanation goes here

[rows,cols]=size(X);
X2 = zeros(rows+2,cols+2);
X2(2:rows+1,2:cols+1)=X(:,:);
if islogical(X)
    X2 = logical(X2);
end

end
