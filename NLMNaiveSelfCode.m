%load the image 
ogIMG = imread("images/alleyNoisy_sigma20.png");
test3by3 = imread("images/test.jpeg");
testCrop = imread("images/test2.jpeg");
testCrop2 = imread("images/test3.jpeg");
testCrop3 = imread("images/test4.jpeg");
testFace = imread("images/noise_example.jpg");
wTest = imread("images/weightTest.jpeg");
townTest = imread("images/townNoisy_sigma5.png"); 
wCheck = imread("images/weightCheck.jpeg"); 

%Convert image to greyscale 
bwIMG = rgb2gray(ogIMG);
test3by3 = rgb2gray(test3by3);
testCrop = rgb2gray(testCrop);
testCrop2 = rgb2gray(testCrop2);
testCrop3 = rgb2gray(testCrop3);
testWeight = rgb2gray(wTest);
testFace = rgb2gray(testFace);
townTest = rgb2gray(townTest);
wCheckbw = rgb2gray(wCheck); 


img = ogIMG;
resultIMG = img;
%making the mega for loop to do the non-local means naive filtering 
%for every pixel in the image

%get the size of the image in terms of rows and cols 
[m, n] = size(img);
%get the desired patch size 
pSize = 3;
winSize = 5; 
pOffset = (pSize - 1) / 2; 
for row = 1 : m 
    for col = 1 : n 
        pMinX = col - pOffset; 
        pMinY = row - pOffset; 
        pMaxX = col + pOffset;
        pMaxY = row + pOffset; 
        isOut = outOfBounds(m, n, pMinX, pMinY, pMaxX, pMaxY); 
        if isOut == 0
            referencePatch = makePatch(pSize, col, row, img); 
            foundPatches = findPatches(col, row, pSize, img, winSize); 
            avgPixelFound = computeWeighting(referencePatch, foundPatches, pSize, 20, 0.25);
%             if avgPixelFound == 0
%                figure; 
%                imshow(uint8(resultIMG));
%                row
%                col
%             end
            resultIMG(row, col) = avgPixelFound; 
        end
    end
end
figure; 
subplot(2,1,1);
imshow(uint8(img));
subplot(2,1,2); 
imshow(uint8(resultIMG));


%Function to perform the patch matching in the assigned search window given
%the coordinate of the point we are looking at
function patches = findPatches(ptX, ptY, patchSize, oImg, wSizeMultiplier)


    %get original image size 
    [m, n] = size(oImg);
    %Initate the size of the patch to be a multiplier of the chosen multiplier 
    searchWindowSize = patchSize * wSizeMultiplier;
    %calculate the actual size of the search window that resides inside the
    %image 
    swOffset = uint16((searchWindowSize - 1)/2); 
    tempSWminX = ptX - swOffset;
    tempSWminY = ptY - swOffset;
    tempSWmaxX = ptX + swOffset; 
    tempSWmaxY = ptY + swOffset;
    if tempSWminX < 1 
        tempSWminX = 1; 
    end 
    if tempSWminY < 1 
        tempSWminY = 1; 
    end 
    if tempSWmaxX > n 
        tempSWmaxX = n; 
    end 
    if tempSWmaxY > m
        tempSWmaxY = m; 
    end 
    %Initilize the search window canvas 
    searchWindowRows = uint16(tempSWmaxY - tempSWminY + 1); 
    searchWindowCols = uint16(tempSWmaxX - tempSWminX + 1);
    searchWindow = zeros(searchWindowRows, searchWindowCols); 
    
    %Copy the original image onto the search window 
    for row = 1 : searchWindowRows
        for col = 1 : searchWindowCols
            searchWindow(row, col) = oImg(uint16(tempSWminY) + uint16(row) -1, uint16(tempSWminX) + uint16(col) -1);
        end
    end
   
    
    %initialize the patch array to be the maxiumum number of possible patches in the
    %search window 
    numPatches = searchWindowRows * searchWindowCols; 
    %initalize the array that will hold all the patches 
    tempPatches = zeros(uint16(patchSize), uint16(patchSize), numPatches);
    
    %loop through all the possible patches, each patch isu essentially a
    %shift in a single pixel 
    %we will reject any patch that exceeds outside the search window
    %boundaries via the use of the outOfBound function 
    patchNumber = 1; 
    for row = 1 : searchWindowRows
        for col = 1 : searchWindowCols
            %first check if the patch is out of bounds from the search
            %window 
            pSizeOffset = (patchSize - 1) / 2;
            patchMinX = col - pSizeOffset;
            patchMinY = row - pSizeOffset;
            patchMaxX = col + pSizeOffset;
            patchMaxY = row + pSizeOffset;
            
            outValue = outOfBounds(searchWindowRows, searchWindowCols, patchMinX, patchMinY, patchMaxX, patchMaxY);
            
            if outValue == 0
                tempPatch = makePatch(patchSize, col, row, searchWindow); 
                tempPatches(:, :, patchNumber) = tempPatch; 
                patchNumber = patchNumber + 1;
            end 
        end 
    end
    
%     figure;
%     for p = 1 : patchNumber - 1
%         
%         subplot(patchNumber -1, 1, p);
%         imshow(uint8(tempPatches(:, :, p)));
%     end
%     
    patches = tempPatches(:, :, 1 : patchNumber-1); 
end


function patch = makePatch(patchSize, ptX, ptY, oImg)
    %The patch size offset is the amount of pixels offset from the centre
    pSizeOffset = int8((patchSize - 1) / 2) * -1;
    %initalizing patchfinding coordinate grid 
    grid = zeros(patchSize, patchSize*2);
    for row = 1 : patchSize
        for col = 1 : patchSize 
            grid(row, (col * 2) - 1) = pSizeOffset + col -1;
            grid(row, (col * 2)) = pSizeOffset + row - 1; 
        end
    end
    %initalize the patch 
    tempPatch = zeros(patchSize); 
    for row = 1 : patchSize
        for col = 1 : patchSize
            %find the correct pixel intensity to copy using the grid 
            tempPatch(row, col) = oImg(ptY + grid(row, col * 2), ptX + grid(row, (col * 2) - 1)); 
        end
    end
    patch = tempPatch;
end 


%We need to create a function that detects when the patch is out of bounds
%from the original image 
function out = outOfBounds(canvasR, canvasC, minX, minY, maxX, maxY)
    %first check if the patch is out of bounds from top left corner 
    if minX < 1 || minY < 1
        out = true;
    %then check if the patch is out of bounds from the bottom right corner 
    elseif maxX > canvasC || maxY > canvasR
        out = true; 
    else 
    %If nothing is out of bounds, then return false, it is not out of
    %bounds
        out = false;
    end
end 


%This is the weighting function using the formula, it also contains the
%stack of code that identifies the SSD
function pixelValue = computeWeighting(refPatch, patches, patchSize, sigma, decayH)
    %Calculate the SSD
    %get the number of patches that are not zeros patches 
    [m1, n1, o] = size(patches);
    ssd = zeros(1, o);
    %next we compare all the patches with the reference patch
    for patch = 1 : o
        tempSSD = 0; 
        for row = 1 : m1 
            for col = 1 : n1 
                %normalised values
                rPatchVal = refPatch(row, col) / 255;
                patchVal = patches(row, col, patch) / 255;
                sd = (rPatchVal - patchVal) ^ 2;
                tempSSD = tempSSD + sd; 
            end
        end
        ssd(1, patch) = tempSSD;
    end
    
    %the ssd of all the patches have been calculated, next we produce the
    %weight array 
    weight = zeros(1, o);
    for d = 1 : o 
        weight(1, d) = weightingFormula(ssd(1, d), sigma, decayH);
    end
    
    %multiply the weight of the patch with the centre pixel and get the
    %average
    
    %storing the intensities found for each of the weights and centre pixel
    intensity = zeros(1, o);
    centrePos = (patchSize + 1) / 2;
    total = 0; 
    for w = 1 : o 
        intensity(1, w) = weight(1, w) * patches(centrePos, centrePos, w);
        total = total + intensity(1, w); 
    end
    
    %instead of dividing by the number of patches compared, we should divde
    %by the number of significant patches instead 
    numOfSigInt = 0;
    for w = 1 : o 
        if weight(1, w) > 0.9
            numOfSigInt = numOfSigInt + 1;
        end 
    end 
    
    averageIntensity = uint16(total / numOfSigInt);
    
%     sortedD = sort(ssd);
%     stdDistance = (sqrt(ssd))*255;
%     stdDistance = sort(stdDistance); 
%     sortedWeight = zeros(1, o);
%  
%     
%     figure;
%     for j = 1 : 10
%         for i = 1 : o
%             sortedWeight(1, i) = weightingFormula(sortedD(1, i), sigma - (j*2), decayH/4);
%         end
%         subplot(5,2,j);
%         plot(stdDistance, sortedWeight)
%     end 
    
    pixelValue = averageIntensity; 
    
    
end 

function value = weightingFormula(ssd, sigma, decayH)
    %normalised sigma 
    normalisationDenominator = 2 * 3.14 * sigma; 
    sigma = sigma / normalisationDenominator;
    temp1 = max(ssd - (2 * sigma), 0); 
    temp2 = temp1 / decayH;
    temp3 = -temp2; 
    value = exp(temp3);
    
end


%we need to initially pad out the image 

%after padding out the image we need to gather the template stack 
%how many pads do we get based on a certain search window? 

%The padding should be the searchwindow width + the patch width (I need to
%re implementing 
%for the template patching we shouldn't even consider the pixels outside
%the image. How do i do that?


%The following function returns a new image that is padded by the amount
%required based on the search window size 
function paddedImg = padOut(SearchWindowSize, img)
    %The padding function has been checked and is functional. 
    padAmount = round((SearchWindowSize - 1) / 2); 
    %initalize the new canvas size 
    [m, n] = size(img); 
    pIMG = zeros(m + padAmount*2, n + padAmount*2);
    offSet = padAmount; 
    newM = m + padAmount*2; 
    newN = n + padAmount*2; 
    
    %Copy the image onto the blank canvas centered in the middle of the
    %canvas 
    for row = 1 : m
        for col = 1: n 
            pIMG(row+offSet, col+offSet) = img(row,col);
        end
    end
    
    %Next copy the horizontal components of the image to the edge of the
    %canvas 
    
    %Get the top edge of the image 
    topEdge = img(1, 1 : n); 
    bottomEdge = img(m, 1 : n); 
    
    %copy them to the edge of the canvase 
    for row = 1 : offSet
        pIMG(row, offSet+1 : offSet + n) = topEdge; 
    end
    
    for row = offSet + m : newM 
        pIMG(row, offSet+1 : offSet + n) = bottomEdge; 
    end
    
    %next pad the vertical edges of the padded canvas by first collecting
    %the side edge and copying it 
    
    leftEdge = pIMG( 1 : newM, offSet + 1); 
    rightEdge = pIMG( 1 : newM, offSet + n); %not certin whether to + 1
    
    for col = 1 : offSet
        pIMG(1 : newM, col) = leftEdge; 
    end
    
    for col = offSet + n + 1 : newN
        pIMG(1 : newM, col) = rightEdge; 
    end 
    
    figure; 
    imshow(uint8(pIMG));
    paddedImg = pIMG; 
end 