%% 1.read compressed data
clear;
% change the storage path of the compressed data to be read
load('/Users/apple/Desktop/34240/project/good/compressed_image.mat');  
origin_width = double(origin_width);
origin_height = double(origin_height);
resize_width = ceil(origin_width/8) * 8;
resize_height = ceil(origin_height/8) * 8;

%% 2.Huffman decoding to obtain DC coefficients(DPCM encoded) and AC coefficients(RLC encoded)
dcY = myhuffmandeco(dcY_comp, dcY_dict);
dcCb = myhuffmandeco(dcCb_comp, dcCb_dict);
dcCr = myhuffmandeco(dcCr_comp, dcCr_dict);
acY = myhuffmandeco(acY_comp, acY_dict);
acCb = myhuffmandeco(acCb_comp, acCb_dict);
acCr = myhuffmandeco(acCr_comp, acCr_dict);

%% 3.DPCM decoding to get DC coefficients
dcY_rec = dpcm_decode(dcY);
dcCb_rec = dpcm_decode(dcCb);
dcCr_rec = dpcm_decode(dcCr);

%% 4.RLC decoding to get AC coefficients
[~,col_cnt] = size(dcY_rec);
acY_rec = rlc_decode(acY,col_cnt);
acCb_rec = rlc_decode(acCb,col_cnt/4);
acCr_rec = rlc_decode(acCr,col_cnt/4);

%% 5.The AC and DC coefficients are recombined into a complete matrix
zigY = [dcY_rec;acY_rec];
zigCb = [dcCb_rec;acCb_rec];
zigCr = [dcCr_rec;acCr_rec];

%% 6.Inverse Zigzag Transform
qdctY = inv_zigzag(zigY,resize_width);
qdctCb = inv_zigzag(zigCb,resize_width/2);
qdctCr = inv_zigzag(zigCr,resize_width/2);

%% 7.Inverse quantization according to quantization table
[iqdctY,iqdctCb,iqdctCr] = inv_quantify(qdctY,qdctCb,qdctCr);

%% 8.Inverse dct transform
[iy,icb,icr] = inv_blockdct(iqdctY,iqdctCb,iqdctCr);

%% 9.Regenerate the YCbCr image from the result of downsampling
iYCbCr = regenerate_Ycbcr(iy,icb,icr);

%% 10.Convert YCbCr image to RGB image
irgb = ycbcr2rgb(iYCbCr);
irgb_origin = irgb(1:origin_height,1:origin_width,:);
figure(),imshow(irgb_origin),title("ac(i) >= -4 && ac(i) <= 4");
% Enter the path where you want to save the reconstructed image
imwrite(irgb_origin,'/Users/apple/Desktop/34240/project/good/reconstructed_image.bmp'); 
%% function：huffman decoding
function sig = myhuffmandeco(comp,dict)
    % Transform comp to double type again, otherwise the huffmandeco function cannot be called
    comp = double(comp);   
    sig = huffmandeco(comp,dict);
end

%% function：DPCM decoding of DC coefficients
function f_rec = dpcm_decode(en)
    [~,cnt] = size(en);
    f_pre = zeros(1,cnt);
    f_rec = zeros(1,cnt);
    % The first two values of the predicted signal are initialized to the first signal of the original signal
    f_pre(1) = en(1); 
    f_pre(2) = en(1);
    % The first value of the reconstructed signal is also initialized to the first signal of the original signal
    f_rec(1) = en(1); 
    for i=2:cnt
        if(i ~= 2)
            f_pre(i) = (f_rec(i-1) + f_rec(i-2))/2;
        end
        f_rec(i) = f_pre(i) + en(i);        
    end
    % Reconvert to double type after rounding
    f_rec = double(uint8(f_rec));   
end

%% function：Decode AC coefficients(RLC encoded)
function out = rlc_decode(in,col_cnt)
    [~,cnt] = size(in);
    % Transform the input to the form of a two-column rlc table
    rlc_table = reshape(in,cnt/2,2);    
    ac_vec = zeros(1,63 * col_cnt);
    j = 1;
    for i=1:cnt/2
        % If the first number of a row is 0, then directly assign the latter number to ac_vec
        if(rlc_table(i,1) == 0)         
            ac_vec(j) = rlc_table(i,2);
            j = j + 1;
        % If the first number of a row is not 0, then this number is the number of 0s, 
        % and ac_vec is filled with the corresponding number of 0s
        else                            
            for k=1:rlc_table(i,1)
                ac_vec(j) = 0;
                j = j + 1;
            end
            % Then assign the last number of that line to ac_vec
            ac_vec(j) = rlc_table(i,2); 
            j = j + 1;
        end
    end
    % Rescale matrix format into 63 rows of AC coefficients, 
    % which is convenient for subsequent splicing with the DC coefficients
    out = reshape(ac_vec,63,col_cnt);   
end

%% function：Traverse the 8x8 matrix through zigzag to produce a 1x64 vector
function zigx = block_zigzag(x)
    zigzag_table=[
        1 2 9 17 10 3 4 11;
        18 25 33 26 19 12 5 6;
        13 20 27 34 41 49 42 35; 
        28 21 14 7 8 15 22 29;
        36 43 50 57 58 51 44 37;
        30 23 16 24 31 38 45 52;
        59 60 53 46 39 32 40 47;
        54 61 62 55 48 56 63 64];
    % Convert the input 8x8 matrix to a 1x64 vector
    v = reshape(x',1,64); 
    vtable = reshape(zigzag_table',1,64);
    % Simulate zigzag scanning by looking up the table
    zigx = v(vtable); 
end

%% function：Divide the input matrix into 8×8 small blocks and perform zigzag processing respectively
function zigx = zigzag(x)
    fun = @(block_struct) block_zigzag(block_struct.data);
    tzigx = blockproc(x,[8,8],fun);
    [a,b] = size(tzigx);
    % Adjust the matrix format to n×64 after zigzag processing, 
    % which is convenient for the subsequent processing of DC coefficients and AC coefficients.
    zigx = reshape(tzigx',64,a*b/64);   
end

%% function：Transform a 64x1 column vector to an 8x8 matrix by inverse zigzag
function invzigx = block_inv_zigzag(x)
    zigzag_table=[
        1 2 9 17 10 3 4 11;
        18 25 33 26 19 12 5 6;
        13 20 27 34 41 49 42 35; 
        28 21 14 7 8 15 22 29;
        36 43 50 57 58 51 44 37;
        30 23 16 24 31 38 45 52;
        59 60 53 46 39 32 40 47;
        54 61 62 55 48 56 63 64];
    
    vtable = reshape(zigzag_table',1,64);   
    invzigx = zeros(8,8);
    for i=1:64
        % Reverse zigzag scan by look-up table
        invzigx(vtable(i)) = x(i);  
    end
    invzigx = invzigx';
end

%% function：Perform inverse zigzag processing on each column of the input matrix
function invzigx = inv_zigzag(x,resize_width)
    tzigx = reshape(x,resize_width*8,[]);
    tzigx = tzigx';
    fun = @(block_struct) block_inv_zigzag(block_struct.data);
    invzigx = blockproc(tzigx,[1,64],fun);
end

%% fucntion：Inverse quantization of matrix using existing luminance and chroma quantization tables
function [y,cb,cr] = inv_quantify(qy,qcb,qcr)
    % Luminance Quantization Table
    LumiTable=[
        16 11 10 16 24 40 51 61;
        12 12 14 19 26 58 60 55;
        14 13 16 24 40 57 69 56;
        14 17 22 29 51 87 80 62;
        18 22 37 56 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99];
    % Chroma Quantization Table
    ChromiTable=[
        17 18 24 47 99 99 99 99 ;
        18 21 26 66 99 99 99 99 ;
        24 26 56 99 99 99 99 99 ;
        47 66 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99];    
    
    fun1 = @(block_struct) block_struct.data .* LumiTable;
    fun2 = @(block_struct) block_struct.data .* ChromiTable;
    y = blockproc(qy,[8,8],fun1);
    cb = blockproc(qcb,[8,8],fun2);
    cr = blockproc(qcr,[8,8],fun2);
end

%% function：Perform inverse dct transform on the image of 8*8 small blocks and return
function [y,cb,cr] = inv_blockdct(dcty,dctcb,dctcr)
    fun = @(block_struct) idct2(block_struct.data);
    y=blockproc(dcty,[8 8],fun);
    cb=blockproc(dctcb,[8 8],fun);
    cr=blockproc(dctcr,[8 8],fun);
end

%% function：Regenerate the image using the result of downsampling
function down_imgYCbCr = regenerate_Ycbcr(imgY,imgCb,imgCr)
    [YH,YL] = size(imgY);
    [CbH,CbL] = size(imgCb);
    y = uint8(imgY);
    cb = zeros(YH,YL);
    cr = zeros(YH,YL);
    for i=1:CbH
        for j=1:CbL
            % Each value is filled with a small 2x2 block
            [cb(2*i,2*j),cb(2*i-1,2*j),cb(2*i,2*j-1),cb(2*i-1,2*j-1)] = deal(imgCb(i,j));
            [cr(2*i,2*j),cr(2*i-1,2*j),cr(2*i,2*j-1),cr(2*i-1,2*j-1)] = deal(imgCr(i,j));
        end
    end
    cb = uint8(cb);
    cr = uint8(cr);
    down_imgYCbCr = cat(3,y,cb,cr);
end