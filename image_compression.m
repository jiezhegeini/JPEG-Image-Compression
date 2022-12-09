clear;
close all;
clc;

%% 1.Read the bmp bitmap and complete the length and width of the image as a multiple of 8
originBMP = imread('/Users/apple/Desktop/34240/project/good/image.bmp'); 
figure(),imshow(originBMP),title("origin bitmap image");
[origin_height,origin_width,~] = size(originBMP);
resizeBMP = resize(originBMP);

%% 2.Convert RGB to YCbCr and perform color downsampling
imgYCbCr = rgb2ycbcr(resizeBMP);     
[imgY,imgCb,imgCr] = downSample(imgYCbCr);

%% 3.The image is divided into 8*8 small blocks for dct transformation processing
[dctY,dctCb,dctCr] = blockdct(imgY,imgCb,imgCr);

%% 4.Quantize according to the quantization table
[qdctY,qdctCb,qdctCr] = quantify(dctY,dctCb,dctCr);

%% 5.zigzag one-dimensional
zigY = zigzag(qdctY);
zigCb = zigzag(qdctCb);
zigCr = zigzag(qdctCr);

%% 6.DC coefficients are coded using DPCM
dcY = dpcm(zigY);
dcCb = dpcm(zigCb);
dcCr = dpcm(zigCr);

H1 = entropy(dcY);
H2 = entropy(dcCb);
H3 = entropy(dcCr);

%% 7.AC coefficients are coded using RLC
acY = rlc(zigY);
acCb = rlc(zigCb);
acCr = rlc(zigCr);

H4 = entropy(acY);
H5 = entropy(acCb);
H6 = entropy(acCr);

H = H1+H2+H3+H4+H5+H6;

%% 8.Apply Huffman coding to DC and AC coefficients
[dcY_comp, dcY_dict] = huffman(dcY);
[dcCb_comp, dcCb_dict] = huffman(dcCb);
[dcCr_comp, dcCr_dict] = huffman(dcCr);
[acY_comp, acY_dict] = huffman(acY);
[acCb_comp, acCb_dict] = huffman(acCb);
[acCr_comp, acCr_dict] = huffman(acCr);

%% 9.Save the compressed data
sizeofimage = origin_width * origin_height;

origin_width = uint16(origin_width);
origin_height = uint16(origin_height);

save('/Users/apple/Desktop/34240/project/good/compressed_image.mat','dcY_comp','dcCb_comp','dcCr_comp','acY_comp','acCb_comp','acCr_comp','dcY_dict','dcCb_dict','dcCr_dict','acY_dict','acCb_dict','acCr_dict','origin_width','origin_height');

%% 10.code length calculation
[a1,b1] = size(dcY_comp);
[a2,b2] = size(dcCb_comp);
[a3,b3] = size(dcCr_comp);
[a4,b4] = size(acY_comp);
[a5,b5] = size(acCb_comp);
[a6,b6] = size(acCr_comp);

b = b1+b2+b3+b4+b5+b6;

code_length = b/sizeofimage; 
compression = 24/code_length; % compression Ratio = original size / compressed size

%% Function: complete the image to a multiple of 8 in both length and width
function new_img = resize(img)
[m,n,~] = size(img);
new_m = ceil(m/8) * 8; 
new_n = ceil(n/8) * 8;
for i=m+1:new_m
    img(i,:,:)=img(m,:,:);
end
for j=n+1:new_n
    img(:,j,:)=img(:,n,:);
end
new_img = img;
end

%% Function: Perform 4:2:0 color downsampling on the obtained YCbCr
function [y,cb,cr] = downSample(I)
    [imgHeight,imgWidth,~] = size(I);
    y = double(I(:,:,1));                             
    cb = double(I(1:2:imgHeight-1,1:2:imgWidth-1,2));
    cr = double(I(2:2:imgHeight,2:2:imgWidth,3));     
end

%% Function: return the image after DCT transformation in units of 8*8 small blocks
function [y,cb,cr] = blockdct(imgY,imgCb,imgCr)
    fun = @(block_struct) dct2(block_struct.data);
    y=blockproc(imgY,[8 8],fun);
    cb=blockproc(imgCb,[8 8],fun);
    cr=blockproc(imgCr,[8 8],fun);
end

%% Function: Quantize the matrix through the existing luminance and chromaticity quantization tables
function [qy,qcb,qcr] = quantify(y,cb,cr)
    LumiTable=[
        16 11 10 16 24 40 51 61 ;
        12 12 14 19 26 58 60 55 ;
        14 13 16 24 40 57 69 56 ;
        14 17 22 29 51 87 80 62 ;
        18 22 37 56 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99];
    ChromiTable=[
        17 18 24 47 99 99 99 99 ;
        18 21 26 66 99 99 99 99 ;
        24 26 56 99 99 99 99 99 ;
        47 66 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99 ;
        99 99 99 99 99 99 99 99];    
    
    fun1 = @(block_struct) round(block_struct.data ./ LumiTable);
    fun2 = @(block_struct) round(block_struct.data ./ ChromiTable);
    qy=blockproc(y,[8,8],fun1);
    qcb=blockproc(cb,[8,8],fun2);
    qcr=blockproc(cr,[8,8],fun2);
end

%% Function: dequantize the matrix through the existing luminance and chroma quantization tables
function [y,cb,cr] = inv_quantify(qy,qcb,qcr)
    LumiTable=[
        16 11 10 16 24 40 51 61;
        12 12 14 19 26 58 60 55;
        14 13 16 24 40 57 69 56;
        14 17 22 29 51 87 80 62;
        18 22 37 56 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99];
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

%% Function function: Traverse the 8×8 matrix through zigzag to generate a 1×64 vector
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
    
    v = reshape(x',1,64); 
    vtable = reshape(zigzag_table',1,64);
    zigx = v(vtable);
end

%% Function: Divide the input matrix into small blocks of 8×8 and perform zigzag processing respectively
function zigx = zigzag(x)
    fun = @(block_struct) block_zigzag(block_struct.data);
    tzigx = blockproc(x,[8,8],fun);
    [a,b] = size(tzigx);
    zigx = reshape(tzigx',64,a*b/64);  
end

%% Function: DPCM encoding for DC coefficients
function en = dpcm(x)
    f = x(1,:);
    [~,cnt] = size(f);
    f_pre = zeros(1,cnt);
    f_rec = zeros(1,cnt);
    e = zeros(1,cnt);
    en = zeros(1,cnt);
    f_pre(1) = f(1); 
    f_pre(2) = f(1);
    f_rec(1) = f(1);
    en(1) = f(1); 
    for i=2:cnt
        if(i ~= 2)
            f_pre(i) = (f_rec(i-1) + f_rec(i-2))/2;
            %f_pre(i) = (f(i-1) + f(i-2))/2;
        end
        e(i) = f(i) - f_pre(i);
        %en(i) = 16 * floor((255 + e(i))/16) - 256 + 8; 
        %en(i) = 8 * floor((255 + e(i))/8) - 256 + 4; 
        %en(i) = 4 * floor((255 + e(i))/4) - 256 + 2; 
        en(i) = 2 * floor((255 + e(i))/2) - 256 + 1;  %no change
        %en(i) = e(i);  
        f_rec(i) = f_pre(i) + en(i);        
    end
    
end

%% Function: RLC encoding of AC coefficients
function rlc_table = rlc(x)
    ac = x(2:64,:);
    [m,n] = size(ac);
    cnt = m * n;
    zero_cnt = 0;
    rlc_table = [];
    step = [1,2,4];
    for i=1:cnt
        %if ac(i) == 0  
        %if ac(i) >= -1 && ac(i) <= 1
        if ac(i) >= -2 && ac(i) <= 2  %no change
        %if ac(i) >= -4 && ac(i) <= 4
            zero_cnt = zero_cnt + 1;
          else
            rlc_table = [rlc_table;[zero_cnt,ac(i)]];
            zero_cnt = 0;
          end
        end
    rlc_table = [rlc_table;[0,0]];  % 表的末行为[0,0]
    end


%% Function: perform huffman encoding on the incoming parameters and return codewords and dictionaries
function [comp,dict] = huffman(x)
    [m,n] = size(x);
    xx = reshape(x, 1, m*n);     
    table = tabulate(xx(:));     
    symbols = table(:,1);        
    p = table(:,3) ./ 100;       
    dict = huffmandict(symbols,p);  
    comp = huffmanenco(xx,dict);    
    comp = uint8(comp);         
end


%% function: entropy estimation
function H = entropy(I)
  [N,edges] = histcounts(I,100000); 
  edges(N==0) = [];
  N(N==0)=[];
  edges = edges + 1;
  edges = edges(1:length(N));
  H = 0;
  for i = 1:length(N)
      H = H + (N(i)/sum(N)*log2(N(i)/sum(N)));
      H = -H;
  end
end