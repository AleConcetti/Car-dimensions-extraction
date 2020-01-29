function matr = BuildMatrixFrom(nameFigure,figure,lenIm)

    if strcmp(nameFigure, 'conic')
        C1=figure;
        matr=zeros(lenIm,lenIm);
        for i=1:lenIm
            for j=1:lenIm
                matr(i,j)=[j i 1]*C1*[j i 1]';
            end    
        end
        matr=double(matr>-0.00007 & matr<0.00007);
    end

    if strcmp(nameFigure, 'line')
        l1=figure;
        matr=zeros(lenIm,lenIm);
        for i=1:lenIm
            for j=1:lenIm
                matr(i,j)=[j i 1]*l1;
            end
        end
        matr=double(matr<0.002 & matr>-0.002);
    end

    if strcmp(nameFigure, 'point')
        point=figure;
        matr=zeros(lenIm,lenIm);
        matr(int16(point(1))-4:int16(point(1))+4,int16(point(2))-4:int16(point(2))+4)=1;
        matr=matr';
      
    end

end

