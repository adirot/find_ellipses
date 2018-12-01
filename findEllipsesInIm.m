function [CC,ellipses,files] = findEllipsesInIm()
prcTh = 90; % prcent threshold for image noise clean (ignore more = higher percent)
minObjectLength = 50; % Delete small objects (they are probably noise)
th = 0.5 ; % threshord for exepting ellipse fit
maxLongAxis = 60; %max ellipse long axis
minShortAxis = 10; %min ellipse long axis
minArea = 1000; %max ellipse area
showProg = false; % show results on image during run

% try
    %load('output.mat');

    %%  get images from folder  
    files = dir( fullfile('*tif*') );
    files = {files.name};
    totNumIm = length(files);
    %%  find ellipes in all pictures  %%%%%
    ellipses = cell(1,totNumIm);
    %parLines = struct;
    for i = 1:totNumIm

        [bw,CC] = getCellsBW(files{1,i},prcTh);
        if showProg
            close all;
            showObjects(CC);
        end
        
        
        %% delete small objects
        objLength = cellfun(@numel,CC.PixelIdxList);
        CC.PixelIdxList(objLength < minObjectLength) = [];
        CC.NumObjects = length(CC.PixelIdxList);
        bw = bwareaopen(bw, 50,4);
        
        if showProg
            close all;
            showObjects(CC);
        end
        
        %% delete child objects 
        bw = delChild(bw);
        bw = bwareaopen(bw, 50,4);
        if showProg
            close all;
            showObjects(CC);
        end
        
        %% sort CC by large to small (more likely to find ellipses in large objs)
        CC = sortCC(CC);
        
        %% find elliptic objects
        [CC,ellipses{1,i}] = findEllipses(CC,100,10,th,showProg,maxLongAxis,minShortAxis,minArea);
        
        %% find cells that don't look like ellipses
%         [lines,point1ind,point2ind] = findParallelLines(bw,CC);
%         parLines(i).lines = lines;
%         parLines(i).point1ind = point1ind;
%         parLines(i).point2ind = point2ind;
        
        
        save('output.mat');
        date_time = datestr(now,'dd-mmm-yyyy HH:MM:SS');
        str = sprintf('finding objects in image number %d of %d images. %s',i,totNumIm,date_time);
        disp(str);
        fileID = fopen('log.txt','w');
        fprintf(fileID,str);
        fclose(fileID);

    end %for i - all files in folder
    
    %% delete small ellipses within bigger ones (slow and stupid code...)
    %newEllipses = delOverlappingEllipses(ellipses,files,0.3); 
    %ellipses = newEllipses;
    %save('output.mat');
    
    date_time = datestr(now,'dd-mmm-yyyy-HH:MM');
    str = sprintf('all %d images done! %s',totNumIm,date_time);
    fileID = fopen('log.txt','w');
    fprintf(fileID,str);
    fclose(fileID);
    disp(str);
 
    
%% Help functions

function [bw,CC] = getCellsBW(fileName,prc)
    %% read image, thrashold it, put found objects in CC

    Im1=imread(fileName);
    [~,~,c] = size(Im1);
    if(c>1)
        Im1=rgb2gray(Im1);
    end

%{
    % delete time labal from image
    xmin = 96; xmax = 493; ymin = 917; ymax = 970;
    Im1(ymin:ymax,xmin:xmax) = 0;
%}

            Im1(:,1:1200) = 0; % make image smaller, for debugging
    %[A ,~ ,~ ,~] = crop4(Im);
    %Im = A;

    h2=[0.1667,     0,      0.6667,       0,      0.1667;...
        0,          0,      0     ,       0,      0     ;...
        0.6667,     0,-3.33333333333333,  0,      0.6667;...
        0,          0,      0,            0,      0     ;...
        0.1667,     0,  0.6667,           0,      0.1667];
    Im_sh2=imfilter(Im1,h2,'replicate'); 

    %figure(3); imshow(Im_sh2,[]);

    % Thershold filered image
    %bw=int16(Im_sh2>20); % mark all values greater than 20 (minima) - for jpg
    AdjustThreshold=prctile(imhist(Im_sh2),prc);
    bw=int16(Im_sh2>AdjustThreshold); %for tif
    %figure(4); imshow(not(bw),[]);

    %Extract from BW image the connected objects. Each object is displayed in  a different
    %color on  result (Fig.  5)
    CC = bwconncomp(bw);
%     L = labelmatrix(CC);
%     RGB2 = label2rgb(L, 'jet', 'k', 'shuffle');
%     figure;imshow(RGB2);

end

% catch err
%    date_time = datestr(now,'dd-mmm-yyyy-HHMM');
%    errFileName = sprintf('errorFile%s',date_time);
%    fid = fopen(errFileName,'a+');
%    fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
%    fclose(fid)
% end

function [] = showObjects(CC)
    L = labelmatrix(CC);
    RGB2 = label2rgb(L, 'jet', 'k', 'shuffle');
    figure;imshow(RGB2);
end

function [CC,final] = findEllipses(CC,runs,stop,th,showProg,maxLongAxis,minShortAxis,minArea)
%th - threshord for exepting ellipse fit
    imageSize = CC.ImageSize;
    ii = 0;
    while(ii < stop)
        %% first try to fit ellipses to objects
        objNum = length(CC.PixelIdxList);
        if (ii == 0)            
            stri = sprintf('blind');
            disp(stri);
            ellipsesA = findEllipBlind(CC,th,maxLongAxis,minShortAxis,minArea);
            final = ellipsesA;
        end
        %% now use algorithem for the rest
        if ( ii > 0 )
            ellipsesA = findSomeEllipses(CC,runs,th,maxLongAxis,minShortAxis,minArea);
            final = [final,ellipsesA];
        end
        %% delete found ellipses from list %%%%
        % also delete unfitted dots from the inside of found ellipses
        % to prevente dubles.
        
        [~,ellipNum] = size(ellipsesA);
        if ( isempty(ellipsesA) && ( ii > 0 ) )
            ii = stop+1;
        end
        if ((ellipNum>0) && (objNum>0) && not(isempty(ellipsesA)))
            jj = 1;
            newobjNum = objNum;
            while( jj <= newobjNum ) 
                for k = 1:ellipNum
                     if ( ~ isempty(ellipsesA{2,k}) )
                         x = ellipsesA{2,k}(:,1);
                         y = ellipsesA{2,k}(:,2);
                         ellipk = ellipsesA(1,k);
                         a = sub2ind(imageSize,x,y);
                         [xcellj,ycellj] = ind2sub(imageSize,CC.PixelIdxList{1,jj});

                         isInsideEllipseAn = isInsideEllipse(xcellj,ycellj,ellipk{1,1});

                         isFoundEllipse = ismember(CC.PixelIdxList{1,jj},a);
                         ind2delete = or(isFoundEllipse,isInsideEllipseAn);

                         CC.PixelIdxList{1,jj}(ind2delete) = []; %delete found ellipse + insides
                         if showProg
                            close all;
                            showObjects(CC);
                         end

                     end
                end
                newobjNum = length(CC.PixelIdxList);
                if(newobjNum == objNum)
                    jj = jj + 1;
                end
            end
        end 
        ii = ii+1;
        if (length(CC.PixelIdxList) == 0)
            ii = stop+1;
        end
    end
    
        function final = findEllipBlind(CC,th,maxLongAxis,minShortAxis,minArea)
            final = {};
            objNum1 = length(CC.PixelIdxList);
            l = 0;
            d = 1.5; % how far should you look for points close to the ellipse?
            for j = 1:objNum1
                [xj,yj] = ind2sub(CC.ImageSize,CC.PixelIdxList{1,j});
                ellipseParam = fit_ellipse(xj,yj);
                isGood = isGoodEllipse(ellipseParam,maxLongAxis,minShortAxis,minArea);
                %% find other dots that are close to the ellipse and set ellipse rank (r)
                if(isGood)
                     l=l+1;
                     final{1,l} = ellipseParam;
                     final{2,l} = ellipseClosePoints(ellipseParam,CC,50,maxLongAxis,minShortAxis,minArea);

                     if (l>0)
                         Sring = 2*d*pi*(ellipseParam.a + ellipseParam.b);
                         r = length(final{2,l})/Sring;
                         if (r>th)
                             final{3,l} = r;
                         end
                         if (r<=th)
                             l = l - 1;
                         end
                     end 
                end %if ellipse is good
            end
            [s,~] = size(final);
            if (s == 2)
                final = [];
            end
        end % findEllipBlind
        
        function final = findSomeEllipses(CC,runs,th,maxLongAxis,minShortAxis,minArea)

            if (nargin < 2)
                runs = 30;
            end 

            finalEllipseInd = 0;
            d = 1.5;


            for i1=1:CC.NumObjects
                l = 0;
                ellipses1 = [];
                ellipsesClosePoints = [];
                celli = CC.PixelIdxList{1,i1};
                imax = length(celli);
                if (imax>50)
                    for k1 = 1:runs
                        str1 = sprintf('run: %d',i1);
                        disp(str1);

                        % choose 5 random points from current object and fit an ellipse to
                        % thease points. the points are with in two ellipses radii
                        % from each other
                        if(imax > 10) %if cell is good
                            bad = true;
                            while bad
                                pointsToFit = UniqRandi(1,imax,5);
                                [x1,y1] = ind2sub(CC.ImageSize,celli(pointsToFit));
                                distMat = ipdm([x1,y1]);

                                % at least 4 of the points are close 
                                bad = distMat > 120;
                                bad = sum(sum(bad)) > 2;
                            end
                            ellipseParam = fit_ellipse(x1,y1);
                            isGood = isGoodEllipse(ellipseParam,maxLongAxis,minShortAxis,minArea);

                            %find other dots that are close to the ellipse
                            dist = 50;
                            if(isGood)
                                l=l+1;
                                ellipses1{1,l} = ellipseParam;
                                ellipsesClosePoints{1,l} = ellipseClosePoints(ellipseParam,CC,dist,maxLongAxis,minShortAxis,minArea);
                            end
                        end %if cell is good 
                    end  %end k runs on cellj
                end
                % find best ellipse
                if(l~=0) % if at least one ellipse was found
                    % d is the width of ellipse ring in which we are looking for supporting
                    % points.
                    % max number of supporting points that the ellipse can have is the area
                    % of the ellipse ring of axis a,b and width d: Sring = 2d pi (a+b)
                    % 
                    % r is a mesure for how good the ellipse fits its surrounding points.
                    % 0<=r<=1. 
                    % r =  number of supporting points / area of a ring of width d around ellipse 

                    n =1;
                    rn = 0;
                    for m=1:length(ellipses1)
                        ellipsem = ellipses1{1,m};
                        isGood = isGoodEllipse(ellipsem,maxLongAxis,minShortAxis,minArea);
                        % Sringm is the area of the m'th ellipse ring
                        Sringm = 2*d*pi*(ellipsem.a + ellipsem.b);
                        rm = length(ellipsesClosePoints{1,m})/Sringm;  
                        isBest = logical(rm>rn);
                        if(isBest && isGood)
                            n=m;
                            rn = rm;
                        end
                    end

                    if(length(ellipsesClosePoints{1,n})>0) %if any supporting points were found
                        finalEllipseInd = finalEllipseInd+1;
                        finalEllipses{1,finalEllipseInd} = ellipses1{1,n};
                        finalEllipses{2,finalEllipseInd} = ellipsesClosePoints{1,n};
                        finalEllipses{3,finalEllipseInd} = rn;
                    else
                        finalEllipses = [];
                    end
                else
                        finalEllipses = [];
                end
            end

            % keep only best ellipses: r>th
            final = finalEllipses;
            if(length(finalEllipses) > 0)
                r=cell2mat(finalEllipses(3,:));
                %delete bad ellipses from final list
                isNotSupported = find(r<th);
                len = length(isNotSupported);
                [~,lenFinal] = size(final);
                if(len >= lenFinal)
                    final = {};
                end

                if ((len>0) && (lenFinal>len))
                    for i3 = 1:len
                        final(:,isNotSupported(len - i3 + 1)) = [];
                    end 
                end
            end
        end
        
        function r=UniqRandi(min,max,len)
            % Produce unique random integers list of length len ranged from a to b
            % linear random unique integers
            if (len > (max - min + 1))
                error('myApp:argChk','len must be lower or equal to (max-min)');
            end 

            clc;
            r=min:max; r=r';
            for i4=1:(max-min)
                rn=round((max-min-1)*rand(1,1))+1;
                tmp=r(i4);
                r(i4)=r(rn);
                r(rn)=tmp;
            end
            r = r(1:len);
        end
    
end % end find ellipses

function CC = sortCC(CC)
    Objs = cellfun(@numel,CC.PixelIdxList);
    [~ , ObjsOrder] = sort(Objs,'descend');
    CC.PixelIdxList = CC.PixelIdxList(ObjsOrder);
end

function newBW = delChild(BW,bg)
% delete "child" objects in a BW image. "child" is an object which is inside
% a hole of another object, and has no holes. 
% bg = 0 or 1: background of the image. default is 0.
    if ( nargin < 2 )
        bg = 0;
    end
    newBW = BW;
    [B,~,~,A] = bwboundaries(BW);
    [a,b] = find(A); 
    % the meaning of a,b: the a'th boundery is within the b'th boundery in B
    
    % if the boundery is enclosed by another object, and there is no other
    % object inclosed in him - the obect is a child. 
    
    isChild = ~ismember(a,b);
    isObject = isObj(B,a,bg,BW);
    isChild = isChild & isObject;
    bounderyNum = a(isChild);
    
    newBW = fillBoundery(B,bounderyNum',newBW);
        function newBW = fillBoundery(B,bounderyNum,bw)
            newBW = bw;
            [m,n] = size(bw);
            for i1 = 1:length(bounderyNum)
                x = B{bounderyNum(i1),1}(:,2);
                y = B{bounderyNum(i1),1}(:,1);
                reg2del = poly2mask(x, y, m, n);
                a1 = newBW(reg2del);
                if( ~isempty(a1) )
                    newBW(reg2del) = deal( ~a1(1,1) );
                end
            end
        end
        
        function isObject = isObj(B,boundaryNum,bg,bw)
            [m,n] = size(bw);
            isObject = false(size(boundaryNum));
            for i2 = 1:length(boundaryNum)
                x = B{boundaryNum(i2),1}(:,2);
                y = B{boundaryNum(i2),1}(:,1);
                obj = poly2mask(x, y, m, n);
                a2 = bw(obj);
                if( ~isempty(a2) && ( a2(1,1) ~= bg ))
                    isObject(i2) = true;
                end
            end
        end
end

function closePoints = ellipseClosePoints(ellipseParam,CC,dist2,maxLongAxis,minShortAxis,minArea)
            d = 1.5;
            isGood = isGoodEllipse(ellipseParam,maxLongAxis,minShortAxis,minArea);

            %find other dots that are close to the ellipse
            if(isGood)
                closePointsX = [];
                closePointsY = [];
                for j = 1:CC.NumObjects
                        cellj = CC.PixelIdxList{1,j};
                        [xcellj,ycellj] = ind2sub(CC.ImageSize,cellj);
                        
                        % is cellj close ("dist" pixells away) from the guessed ellipse
                        dist = dist2+ellipseParam.long_axis;
                        isCellCloseToEllipse=isInsideEllipse(xcellj,ycellj,stretcheEllipse(ellipseParam,dist));

                    if(isCellCloseToEllipse) 
                        % find the dots of cellj that are close (d pixells away)
                        % from the guessed ellipse. add the close points to "closePointsX/Y" 
                        closePoints = isPointCloseToEllipse(ellipseParam,xcellj,ycellj,d);
                        closeInd = find(closePoints);
                        closePointsX = [closePointsX;xcellj(closeInd)];
                        closePointsY = [closePointsY;ycellj(closeInd)];
                    end
                end %end find close points
                closePoints = [closePointsX,closePointsY];

            end %if ellipse is good
            if (~isGood)
                closePoints = [];
            end
            
            % Help functions
            function  isPointCloseToEllipse = isPointCloseToEllipse(ellipseParam,x,y,dist)
                % 
                % isPointCloseToEllipse is true if point/s "x,y" are "dist"
                % away from the ellipse spacified in "ellipseParam"

                % input:    ellipseParam is a structure created with "fit_ellipse.m" function.
                %           x,y are coordinates to check if close to ellipse, can be
                %           scalar or line vectors.
                %           dist is the max dist of points to ellipse.
                bigEllipse =  stretcheEllipse(ellipseParam,dist);
                smallEllipse =  stretcheEllipse(ellipseParam,-dist);
                isInsideBigEllipse = isInsideEllipse(x,y,bigEllipse);
                isOutsideSmallEllipse = not(isInsideEllipse(x,y,smallEllipse)); 
                isPointCloseToEllipse = logical(isInsideBigEllipse & isOutsideSmallEllipse);
            end
            
            
            
            function NewEllipseParam = stretcheEllipse(ellipseParam,stretch)
                NewEllipseParam = ellipseParam;
                NewEllipseParam.a = NewEllipseParam.a + stretch;
                NewEllipseParam.b = NewEllipseParam.b + stretch;
                NewEllipseParam.long_axis = 2*max(NewEllipseParam.a,NewEllipseParam.b);
                NewEllipseParam.short_axis = 2*min(NewEllipseParam.a,NewEllipseParam.b);
            end

end



end