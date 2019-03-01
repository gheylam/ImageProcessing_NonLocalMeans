%Efficient template matching code 

%Load up a 3by3 test image 
testImg = imread("images/test.jpeg"); 
test2Img = imread("images/alleyReference.png");
test3Img = imread("images/test2.jpeg");
testImgBW = rgb2gray(testImg);
test2ImgBW = rgb2gray(test2Img);
test3ImgBW = rgb2gray(test3Img);
%filtering parameters 
%patch size, search window, sigma, decayH
pSize = 3; 
swSize = pSize * 3;
sigma = 5;
decayH = 0.5; 

%The samplingImg is the image we pick values from, the resultImg is the
%denoised image. 
samplingImg = test2ImgBW;
resultImg = samplingImg;

%Loop through all the pixels and find the intensities and replace the ones
%in the result image 

[imgM, imgN] = size(samplingImg);

% offsets = disFinder(swSize, pSize, 9, 9, samplingImg);
% ssdList = getSSD(offsets, pSize, 9, 9, samplingImg);
for row = 1 : imgM 
    for col = 1 : imgN
        %intImg = integralImageCalc(samplingImg);
        %first we need to check if the patch will exceed the boundaries of
        %the image 
        isOut = isOutBound(imgN, imgM, pSize, col, row);
        if isOut == 0
            
            offsets = disFinder(swSize, pSize, col, row, samplingImg);
            ssdList = getSSD(offsets, pSize, col, row, samplingImg);
            pI = calcIntensity(col, row, offsets, ssdList, samplingImg, sigma, decayH);
            resultImg(row, col) = pI;
        end 
    end
end

figure; 
subplot(2,1,1);
imshow(uint8(samplingImg));
subplot(2,1,2); 
imshow(uint8(resultImg));


%the following function finds the list of displacements for a given point
%and search window 
function displacements = disFinder(swSize, pSize, ptX, ptY, img)
   
    %calculate the correct search window by truncating the pixels that are
    %out of bound from the original image 
    swOffset = int16(swSize - 1) / 2;
    %compute the limits of the search window 
    swMinX = ptX - swOffset; 
    swMinY = ptY - swOffset; 
    swMaxX = ptX + swOffset; 
    swMaxY = ptY + swOffset; 
    %get the limits of the image
    [m, n] = size(img); 
    %resize the limits of the search window 
    if swMinX < 1 
        swMinX = 1; 
    end 
    if swMinY < 1 
        swMinY = 1; 
    end 
    if swMaxX > n
        swMaxX = n; 
    end
    if swMaxY > m
        swMaxY = m;
    end
    %with the search window limits corrected, we can identify all the
    %patches to compare 
    
    %we also initalize the displacement array to store all the
    %displacements 
    disArray = zeros(2, swMaxX * swMaxY);
    
    %loop through all the possible pixels in the corect search window 
    %we will ignore the pixels that create patches that exceed outside the
    %boundaries of the search window 
    
    %set up a displacement array counter
    dCounter = 1; 
    for row = swMinY : swMaxY
        for col = swMinX : swMaxX
            %Check if the patch goes out of bounds of the search window 
            isOut = isOutBound(swMaxX, swMaxY, pSize, col, row);
            if isOut == 0 %if the patch does not exceed the boundaries of the search window
                %we will add its offset into our displacement array 
                disArray(1, dCounter) = col - ptX; %the offset x amount 
                disArray(2, dCounter) = row - ptY; %the offset y amount
                %increment the dCounter 
                dCounter = dCounter + 1; 
            end 
            %if the patch is out of bounds, we do nothing and just loop to
            %the next set of coordinates 
        end 
    end 
    %get rid of the exta zero elements 
    disArray = disArray(:, 1 : dCounter -1);
    
    %output the array containing all the offsets we need to compute for the
    %ssd 
    displacements = disArray; 
end 

%The following function takes all the offsets for the single point and
%computes the array of ssd for each patch 
function ssdArray = getSSD(offsetArray, pSize, ptX, ptY, img) 
    %initlize the ssd array to contain all the ssd's
    [oRow, oCol] = size(offsetArray); 
    ssdA = zeros(1, oCol); 
    %loop through all the offsets and calucalate the SSD 
    for ssdC = 1 : oCol
        %get the offset values 
        offsetX = offsetArray(1, ssdC); 
        offsetY = offsetArray(2, ssdC); 
        %Find the distance ^ 2 image map 
        diffImage = calcDifferenceImage(offsetX, offsetY, img);
        %Next using the difference Image, calculate the integral image
        intImage = integralImageCalc(diffImage);
        %After obtaining the integral image we can find the ssd for the
        %exact patch we have 
        pOffset = (pSize - 1) / 2;
        pMinX = ptX - pOffset; 
        pMaxX = ptX + pOffset;
        pMinY = ptY - pOffset; 
        pMaxY = ptY + pOffset;
        
        %correct the limits to the integral image constraints 
        [intM, intN] = size(intImage);
        if pMinX > intN 
            pMinX = intN;
        end
        if pMinY > intM
            pMinY = intM;
        end
        if pMaxX > intN
            pMaxX = intN;
        end
        if pMaxY > intM 
            pMaxY = intM;
        end
        
        %using the formula we will obtain the sum of this patch 
        l1 = intImage(pMinY, pMinX); 
        l2 = intImage(pMinY, pMaxX); 
        l3 = intImage(pMaxY, pMaxX);
        l4 = intImage(pMaxY, pMinX);
        sumOfPatch = l3 - l2  - l4 + l1; 
        %insert the sum into the array 
        ssdA(1, ssdC) = sumOfPatch;
    end
    %output the ssd of the patches 
    ssdArray = ssdA; 
        
end 

%The following function caluclates the difference image of an image and its offset 
function dImg = calcDifferenceImage(offsetX, offsetY, img) 
    %get the size of the image 
    [mImg, nImg] = size(img); 
    %first we need to identity how the Images will be tructated 
    dMinX = 1; 
    dMinY = 1;
    dMaxX = nImg;
    dMaxY = mImg;
    
    ogMinX = 1; 
    ogMinY = 1;
    ogMaxX = nImg; 
    ogMaxY = mImg; 
    %please note, what ever we do to the displaced image, to apply the
    %opposite operator to the original image
    
    if offsetX > 0 %if the x offset is positive 
        %we will shrink x from the top left
        dMinX = 1 + offsetX;
        %since x offset is positive, the shift for og image will be
        %negative 
        ogMaxX = ogMaxX - offsetX; 
    else
        %if the x offset is negative, we shrink from the bottom right 
        dMaxX = dMaxX + offsetX;
        %since x offset is negative, the shift for og image will be
        %positive 
        ogMinX = 1 - offsetX;
    end 
    if offsetY > 0 %if the y offset is positive
        %we will shrink y from the top left 
        dMinY = 1 + offsetY; 
        ogMaxY = ogMaxY - offsetY;
    else 
        %if the y offset is negative, we will shrink from the bottom right 
        dMaxY = dMaxY + offsetY; 
        ogMinY = 1 - offsetY; 
    end 
    
    %now we construct the two images, dImg and ogImage 
    dImage = double(img(dMinY : dMaxY, dMinX : dMaxX)); 
    ogImage = double(img(ogMinY : ogMaxY, ogMinX : ogMaxX)); 
    
    %the two images should now have the same dimensions for us to caluclate
    %the difference  
    %It is in this moment we also normalise the values 
    diffImage = ((ogImage - dImage) .^ 2) / 255;
    
    %now we return the image 
    dImg = diffImage;
end
   
%the following function calculates the intensity of the pixel we are
%modifying
function pixIntensity = calcIntensity(ptX, ptY, offsetArray, ssdArray, img, sigma, decayH)
    %First we must calc the weights array from the ssd, sigma and decay H 
    [rowC, colC] = size(ssdArray); 
    
    wArray = zeros(1, colC);
    %for every ssd we get add a new weight to the weight array 
    for wC = 1 : colC
        weight = calcWeight(ssdArray(1, wC), sigma, decayH);
        wArray(1, wC) = weight; 
    end 
    
    %Given the weight array we now loop through all the weights and
    %multiply it with the pixel 
    %Then accumulate the values in accumulator 
    accumulator = 0;
    for wC = 1 : colC
        %first we need to discover the pixel intensity 
        pixelX = ptX + offsetArray(1, wC);
        pixelY = ptY + offsetArray(2, wC); 
        pixel = img(pixelY, pixelX);
        pIntensity = pixel * wArray(1, wC); 
        accumulator = uint16(accumulator) + uint16(pIntensity);
    end
    
    pixIntensity = uint8(accumulator / colC);
end 

function value = calcWeight(ssd, sigma, decayH)
    %normalised sigma 
    sigma = sigma / 255;
    temp1 = max(ssd - (2 * sigma), 0); 
    temp2 = temp1 / decayH;
    temp3 = -temp2; 
    value = exp(temp3);
    
end


function out = isOutBound(limitX, limitY, windowSize, ptX, ptY)
    offset = (windowSize - 1) / 2; 
    minX = ptX - offset;
    minY = ptY - offset; 
    maxX = ptX + offset; 
    maxY = ptY + offset; 
    
    if minX < 1 || minY < 1
        out = 1; 
    elseif maxX > limitX || maxY > limitY
        out = 1; 
    else 
        out = 0; 
    end 
end

    
%Integral Image calculator
function image = integralImageCalc(img)
    %create the integral image values holder 
    [m,n] = size(img); 
    integralImg = zeros(m+1,n+1); 
    horizontalImg = zeros(m+1, n+1);
    
    %caluclate the horizontal integral image
    for row = 1 : m+1
        for col = 1: n+1 
            if col == 1 || row == 1
                horizontalImg(row, col) = 0;
            else
                horizontalImg(row, col) = horizontalImg(row, col - 1) + img(row-1, col-1); 
            end
        end
    end
    
    %the sum of the value is the sum of the value above + the rest
    %on its row
    for row = 1 : m +1 
        for col = 1 : n +1 
            %if it is the first element, let the first element by itself 
            if row == 1 || col == 1 
                integralImg(row, col) = 0; 
            else
                %computing the integral value 
                integralImg(row, col) = integralImg(row - 1, col) + horizontalImg(row, col); 
            end
        end
    end
    
    %remove the edge zeros 
    integralImg = integralImg(2:m+1, 2:n+1);
    image = integralImg;
end