%SEGMENTA EN FUNCIÓN DE VORONOI DIRECTO E INVERSO

%% 1: Cargamos mapas de ocupación
load('mapasOccu.mat');
mat=robotnik_small;%willow_garage;

%% 2: PREPROCESAMIENTO
%Zona vacía
tic

J=imbinarize(mat,0.9);
Joriginal=J;
SE = strel('disk',4);
J = imerode(J,SE);
Joriginal=J;
J = medfilt2(J);
J = bwareafilt(J,1);
Joriginal = bwareafilt(Joriginal,1);

%Zona ocupada
Jinv=imbinarize(mat,0.5);
SE = strel('square',5);
Jinv = ~imerode(Jinv,SE);
%Jinv = ~medfilt2(Jinv);

%% 3: ESQUELETIZACIÓN DIRECTA
skel=bwskel(J);%bwmorph(J,'skel',Inf);% %Esqueleto
skel = bwmorph(skel,'spur',10);

localMin=zeros(size(J));
D = bwdist(~J);
branchPt=bwmorph(skel,'branchpoints');
branchPt = imdilate(branchPt,strel('disk',3));
branches=imbinarize(skel-branchPt);
branches=bwareaopen(branches, 10);
[labelBranches, nLabel]=bwlabel(branches);
for i=1:nLabel
    auxImg=labelBranches==i; %Máscara con la label correspondiente
    auxImg=auxImg.*D;
    auxImg(auxImg==0)=Inf;
    minVal=min(min(auxImg)); %Elegir mínimo
    minVal=find((auxImg<(minVal+2)).*(auxImg>=(minVal-2)));%find(auxImg==minVal);
    localMin(minVal)=1;
end
localMin=bwmorph(localMin,'clean');

%% 4: ESQUELETIZACIÓN INVERSA
skelInv=bwmorph(Jinv,'skel',Inf);
skelInv = bwareaopen(skelInv, 20);
invEP=bwmorph(skelInv,'endpoints');
invEP = regionprops(invEP,'centroid');
invEP = round(cat(1,invEP.Centroid));

%% 5: ANCHO DE LAS PAREDES
Dinv = bwdist(~Jinv);
%BdistInv=skelInv.*Dinv;
BdistInv=bwmorph(skelInv,'spur',10).*Dinv;
[row,col,v] = find(BdistInv);
meanWall=max(v); %Distancia media en las ramas
if meanWall>10
    meanWall=round(mean(v));
end
%thresWall=0.5*meanWall; %Distancia máxima entre endpoint y punto ocupado más cercano

%% 6: PUNTOS MÁS CERCANOS A MÍNIMOS LOCALES PERPENDICULARMENTE
localCentroid = regionprops(localMin,'centroid');
localCentroid = round(cat(1,localCentroid.Centroid)); %Coordenadas de los mínimos locales

localOrientation = regionprops(localMin,'orientation');
localOrientation = round(cat(1,localOrientation.Orientation));
localOrientation=localOrientation+90;
tangentes=tand(localOrientation);

perpPoints=[];
% imshow(J+skelInv)
% hold on
% plot(localCentroid(:,1),localCentroid(:,2),'*')
giros=[];
%%
%J = imerode(J,SE);
for i=1:size(tangentes,1)
    localPoints=[];
    if abs(tangentes(i))==Inf %Orientación vertical
        %Hacia arriba
        pointArriba=false;
        deltaY=0;
        while ~pointArriba
            deltaY=deltaY+1;
            if J(localCentroid(i,2)-deltaY,localCentroid(i,1))==0
                pointArriba=true;
            end
        end
        % plot(localCentroid(i,1),localCentroid(i,2)-deltaY,'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
        localPoints=[localPoints,localCentroid(i,2)-deltaY,localCentroid(i,1)];
        %line([localCentroid(i,1) localCentroid(i,1)], [localCentroid(i,2)-deltaY localCentroid(i,2)], 'Color','g','LineWidth',1)

        %Hacia abajo
        pointAbajo=false;
        deltaY=0;
        while ~pointAbajo
            deltaY=deltaY+1;
            if J(localCentroid(i,2)+deltaY,localCentroid(i,1))==0
                pointAbajo=true;
            end
        end
        % plot(localCentroid(i,1),localCentroid(i,2)+deltaY,'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
        localPoints=[localPoints,localCentroid(i,2)+deltaY,localCentroid(i,1)];
        %line([localCentroid(i,1) localCentroid(i,1)], [localCentroid(i,2)+deltaY localCentroid(i,2)], 'Color','g','LineWidth',1)
        perpPoints=[perpPoints;localPoints];
        
    elseif tangentes(i)==0
        %Hacia derecha
        pointDerecha=false;
        deltaX=0;
        while ~pointDerecha
            deltaX=deltaX+1;
            if J(localCentroid(i,2),localCentroid(i,1)+deltaX)==0
                pointDerecha=true;
            end
        end
        % plot(localCentroid(i,1)+deltaX,localCentroid(i,2),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
        localPoints=[localPoints,localCentroid(i,2),localCentroid(i,1)+deltaX];
        %line([localCentroid(i,1)+deltaX localCentroid(i,1)], [localCentroid(i,2) localCentroid(i,2)], 'Color','g','LineWidth',1)

        %Hacia izquierda
        pointIzquierda=false;
        deltaX=0;
        while ~pointIzquierda
            deltaX=deltaX+1;
            if J(localCentroid(i,2),localCentroid(i,1)-deltaX)==0
                pointIzquierda=true;
            end
        end
        % plot(localCentroid(i,1)-deltaX,localCentroid(i,2),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
        localPoints=[localPoints,localCentroid(i,2),localCentroid(i,1)-deltaX];
        %line([localCentroid(i,1)-deltaX localCentroid(i,1)], [localCentroid(i,2) localCentroid(i,2)], 'Color','g','LineWidth',1)
        perpPoints=[perpPoints;localPoints];
        
    else
        tangentes(i)=-tangentes(i);

        %Intersección con los bordes de la imagen:
        upperLim=localCentroid(i,1)+((1-localCentroid(i,2))/tangentes(i));
        lowerLim=localCentroid(i,1)+((size(J,1)-localCentroid(i,2))/tangentes(i));
        leftLim=localCentroid(i,2)+(tangentes(i)*(1-localCentroid(i,1)));
        rightLim=localCentroid(i,2)+(tangentes(i)*(size(J,2)-localCentroid(i,1)));

        giros=[giros;upperLim,lowerLim,leftLim,rightLim];

        if upperLim>=1 && upperLim<=size(J,2) %Corta arriba
            [cx,cy,c]=improfile(~J,[localCentroid(i,1) upperLim],[localCentroid(i,2) 1]);
            desPx = find(c,1);
            % plot(cx(desPx),cy(desPx),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
            localPoints=[localPoints,cy(desPx),cx(desPx)];
            %line([cx(desPx) localCentroid(i,1)], [cy(desPx) localCentroid(i,2)], 'Color','g','LineWidth',1)
        end

        if lowerLim>=1 && lowerLim<=size(J,2) %Corta abajo
            [cx,cy,c]=improfile(~J,[localCentroid(i,1) lowerLim],[localCentroid(i,2) size(J,1)]);
            desPx = find(c,1);
            % plot(cx(desPx),cy(desPx),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
            localPoints=[localPoints,cy(desPx),cx(desPx)];
            %line([cx(desPx) localCentroid(i,1)], [cy(desPx) localCentroid(i,2)], 'Color','g','LineWidth',1)
        end

        if leftLim>=1 && leftLim<=size(J,1) %Corta a la izquierda
            [cx,cy,c]=improfile(~J,[localCentroid(i,1) 1],[localCentroid(i,2) leftLim]);
            desPx = find(c,1);
            % plot(cx(desPx),cy(desPx),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
            localPoints=[localPoints,cy(desPx),cx(desPx)];
            %line([cx(desPx) localCentroid(i,1)], [cy(desPx) localCentroid(i,2)], 'Color','g','LineWidth',1)
        end

        if rightLim>=1 && rightLim<=size(J,1) %Corta a la derecha
            [cx,cy,c]=improfile(~J,[localCentroid(i,1) size(J,2)],[localCentroid(i,2) rightLim]);
            desPx = find(c,1);
            % plot(cx(desPx),cy(desPx),'s','MarkerEdgeColor','b','MarkerFaceColor','c','MarkerSize',5)
            localPoints=[localPoints,cy(desPx),cx(desPx)];
            %line([cx(desPx) localCentroid(i,1)], [cy(desPx) localCentroid(i,2)], 'Color','g','LineWidth',1)
        end

        perpPoints=[perpPoints;localPoints];
        
    end
end

%plot(invEP(:,1),invEP(:,2),'m*','MarkerSize',10)
%plot(localCentroid(:,1),localCentroid(:,2),'o','MarkerEdgeColor','#77AC30','MarkerFaceColor','g','MarkerSize',7)

%% PRUEBA
%Definir punto central para cada par de endpoints
centralEP=[];

[k1,dist1]=dsearchn(invEP,[perpPoints(:,2),perpPoints(:,1)]);
[k2,dist2]=dsearchn(invEP,[perpPoints(:,4),perpPoints(:,3)]);

finalDoors=[];
for i=1:size(k1,1)
    localCentralEP=[(invEP(k1(i),1)+invEP(k2(i),1))/2,(invEP(k1(i),2)+invEP(k2(i),2))/2];

    %line([invEP(k1(i),1) invEP(k2(i),1)],[invEP(k1(i),2) invEP(k2(i),2)]);
    %plot(localCentralEP(1),localCentralEP(2),'x')

    %Comprobar que no corte el Voronoi inverso
    intersectInvSkel1=improfile(skelInv,[invEP(k1(i),1) perpPoints(i,2)],[invEP(k1(i),2) perpPoints(i,1)]);
    intersectInvSkel2=improfile(skelInv,[invEP(k2(i),1) perpPoints(i,4)],[invEP(k2(i),2) perpPoints(i,3)]);

    %Comprobar que la orientación de los dos segmentos es similar
    EPorientation=atan2(abs(invEP(k2(i),2)-invEP(k1(i),2)), abs(invEP(k2(i),1)-invEP(k1(i),1)));
    perpPointsOrientation=atan2(abs(perpPoints(i,1)-perpPoints(i,3)),abs(perpPoints(i,2)-perpPoints(i,4)));
    angleDiff=abs(EPorientation-perpPointsOrientation);

    if pdist2(localCentralEP,localCentroid(i,:))<=meanWall && dist1(i)<=meanWall*2 && dist2(i)<=meanWall*2 && angleDiff<pi*40/180 && nnz(intersectInvSkel1)<3 && nnz(intersectInvSkel2)<3
        finalDoors=[finalDoors;localCentroid(i,:),perpPoints(i,:)];
    end
end

%% 8: IMAGEN FINAL SEGMENTADA
Jseg=Joriginal;
for i=1:size(finalDoors,1)
    Jseg = insertShape(single(Jseg),'Line',[finalDoors(i,1),finalDoors(i,2),finalDoors(i,4),finalDoors(i,3)],'Color','black','Opacity',1,'LineWidth',2);
    Jseg = insertShape(single(Jseg),'Line',[finalDoors(i,1),finalDoors(i,2),finalDoors(i,6),finalDoors(i,5)],'Color','black','Opacity',1,'LineWidth',2);
    Jseg = insertShape(single(Jseg),'FilledCircle',[finalDoors(i,4),finalDoors(i,3),2],'Color','black','Opacity',1);
    Jseg = insertShape(single(Jseg),'FilledCircle',[finalDoors(i,6),finalDoors(i,5),2],'Color','black','Opacity',1);
end
Jseg=im2gray(Jseg);
Jseg=imbinarize(Jseg);
Jseg=bwareaopen(Jseg, 50);
L = bwlabel(Jseg); %Labels de las habitaciones

%Rellena líneas de separación con una etiqueta
for i=1:size(L,1)
    for j=1:size(L,2)
        if L(i,j)==0 && J(i,j)==1
            options=[L(i-1,j-1), L(i-1,j), L(i,j-1), L(i+1,j+1), L(i+1,j), L(i,j+1), L(i+1,j-1), L(i-1,j+1)];
            options=nonzeros(options)';
            newLabel=min(options);
            if ~isempty(newLabel)
                L(i,j)=newLabel;
            end
        end
    end
end

time=toc

rgb = label2rgb(L,'jet',[.5 .5 .5]);
%clearvars -except L rgb

figure(2)
imshow(rgb)

hold on
plot(finalDoors(:,1),finalDoors(:,2),'d');