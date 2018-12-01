function plotEllipsesOnIm(ellipses,files,imNum)

    Im1=imread(files{1,imNum});
    imshow(Im1,[]);
    hold on;
    grades = getGrades(ellipses{1,imNum});
    centers = getCenters(ellipses{1,imNum});
    for i = 1:length(ellipses{1,imNum})
        plotEllipse(ellipses{1,imNum}{1,i});
    %         if(~isempty(ellipClosePoints{1,i}))
    %             %plot(ellipClosePoints{1,i}(:,2),ellipClosePoints{1,i}(:,1),'r.');
    %             if ( nargin > 2 )
    %                 str = sprintf('%.1f',grades(i));
    %                 text(centers(i,2),centers(i,1),  str, 'Color', 'g');
    %             end
    %         end
    end

    
    %% help functions
    function grades = getGrades(ellipses)
        for ii = 1:length(ellipses)
            if(~isempty(ellipses{3,ii}))
                grades(ii) = ellipses{3,ii};
            else
                grades(ii) = 0;
            end
        end
    end

    function centers = getCenters(ellipses)
        for ii = 1:length(ellipses)

                centers(ii,1) = ellipses{1,ii}.X0_in;
                centers(ii,2) = ellipses{1,ii}.Y0_in;

        end
    end

    function plotEllipse(ellipseParam)

        if(~isempty(ellipseParam))
            phi = ellipseParam.phi;
            a = ellipseParam.a;
            b = ellipseParam.b;
            x0 = ellipseParam.X0_in;
            y0 = ellipseParam.Y0_in;

            theta = linspace(0,2*pi);

                x = x0+a*cos(theta)*cos(phi)-b*sin(theta)*sin(phi);
                y = y0+a*cos(theta)*sin(phi)+b*sin(theta)*cos(phi);
                plot(y,x);
                hold on;
                plot(y0,x0,'+r');


        end
    end
end
