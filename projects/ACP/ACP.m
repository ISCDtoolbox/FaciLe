%% ACP sur les masseters
%
% Il s'agit d'identifier les formes les plus representatives des
% deplacements prix en comptes lors du morphing des masseters (a developper plus subtilement).
%

%% Importation de tous les fichiers ".depl" pour l'ACP

numFiles = 70; %Nombre de fichiers a importer (selon segmentation et morphing)
numLines = 4859;%5974; %Nombre de lignes d'un fichiers (toujours pareil puisque morphing a partir d'une ellipse)
massRot = cell(1,numFiles);
massNoRot = cell(1,numFiles);

for fileNum = 1:numFiles
    fileName = sprintf('DispRot.%d.sol',fileNum);
    massRot{fileNum} = readsol(fileName);
    massRot{fileNum} = massRot{fileNum}.';
end

for fileNum = 1:numFiles
    fileName = sprintf('DispNoRot.%d.sol',fileNum);
    massNoRot{fileNum} = readsol(fileName);
    massNoRot{fileNum} = massNoRot{fileNum}.';
end

[template,tri,edg,crn]= readmesh('templateEllipSurf.mesh');
template = template.';

%% Importation de tous les fichiers ".grad" pour le produit scalaire

gradP0Rot = cell(1,numFiles);
gradP0NoRot = cell(1,numFiles);
gradP0RotNoNorm = cell(1,numFiles);

for fileNum = 1:numFiles
    fileName = sprintf('GradP0Rot.%d.sol',fileNum); %GradP0(sur triangle)
    gradP0Rot{fileNum} = readGrad(fileName);
    fileName = sprintf('GradP0NoRot.%d.sol',fileNum); %GradP0(sur triangle)
    gradP0NoRot{fileNum} = readGrad(fileName);
    fileName = sprintf('GradP0RotNoNorm.%d.sol',fileNum); %GradP0(sur triangle)
    gradP0RotNoNorm{fileNum} = readGrad(fileName);
end

%% Calcul de la matrice des variances/covariances (partie la plus longue)
% avec "produit scalaire H1"

% Nombre de masseters a cinsiderer dans la base
numFiles = numFiles-1; % pour avoir un masseter pas dans la base

% Initialisation/Allocation de memoire (pas obligatoire mais gain de temps
% (plus rapide de remplir une matrice plutot que d'augmenter sa taille)
ARot = zeros(numFiles, numFiles);
ANoRot = zeros(numFiles, numFiles);
AGradRot = zeros(numFiles, numFiles);
AGradNoRot = zeros(numFiles, numFiles);
Fact = 0.0000; % 0 pour le produit scalaire L2 et 1 pour le produit scalaire H1

% Boucles de calcul des produits scalaire = matrice variance/covariance
for i = 1:numFiles
    for j = 1:numFiles
        aRot = sum(dot(massRot{i}, massRot{j}));
        aNoRot = sum(dot(massNoRot{i}, massNoRot{j}));
        aaRot = sum(dot(gradP0RotNoNorm{i}, gradP0RotNoNorm{j}));
        aaNoRot = sum(dot(gradP0NoRot{i}, gradP0NoRot{j}));
        ARot(i,j) = aRot + Fact*aaRot;
        ANoRot(i,j) = aNoRot + Fact*aaNoRot;
        AGradRot(i,j) = aaRot;
        AGradNoRot(i,j) = aaNoRot;
    end 
end
[VRot, DRot] = eig(ARot);
[VNoRot, DNoRot] = eig(ANoRot);
[VGradRot, DGradRot] = eig(AGradRot);
[VGradNoRot, DGradNoRot] = eig(AGradNoRot);

% A matrice des variances/covariances
% D matrice des valeurs propres rang?es par ordes croissants
% V matrice des vecteurs propres

%% Forme du masseter moyen / Composantes Principales

% Initialisation des tableaux de cellules / Allocation de memoire
bRot = cell(numFiles, numFiles);
BRot = cell(1, numFiles);
bNoRot = cell(numFiles, numFiles);
BNoRot = cell(1, numFiles);
for i = 1:numFiles
   BRot{i} =  zeros(numLines,3);
   BNoRot{i} =  zeros(numLines,3);
end

%avec seulement le gradient dans la matrice A
bGradRot = cell(numFiles, numFiles);
BGradRot = cell(1, numFiles);
bGradNoRot = cell(numFiles, numFiles);
BGradNoRot = cell(1, numFiles);
for i = 1:numFiles
   BGradRot{i} =  zeros(numLines,3);
   BGradNoRot{i} =  zeros(numLines,3);
end

% Calcul des vecteurs a utiliser dans la combinaison lineaire : les
% composantes principales j'imagine.
for i = 1:numFiles % les "70" valeur d'un vecteur propre
    for j = 1:numFiles % Numero des vecteur propres consid?r?s
        bRot{i, j} = VRot(i,j)*massRot{i};
        BRot{j} = BRot{j} + bRot{i, j};
        bNoRot{i, j} = VNoRot(i,j)*massNoRot{i};
        BNoRot{j} = BNoRot{j} + bNoRot{i, j};
        %fileName = sprintf('DispMoyRot.%d.sol',i); % Sauvegarde des etapes
        %dlmwrite(fileName, B{j} , 'delimiter', ' '); %
        bGradRot{i, j} = VGradRot(i,j)*massRot{i};
        BGradRot{j} = BGradRot{j} + bGradRot{i, j};
        bGradNoRot{i, j} = VGradNoRot(i,j)*massNoRot{i};
        BGradNoRot{j} = BGradNoRot{j} + bGradNoRot{i, j};
    end
end

% Normalisation des B{i} (composante principale)
BBRot = cell(numFiles,1);
BBNoRot = cell(numFiles,1);
for i = 1:numFiles
    BBRot{i} = BRot{i}/sqrt(sum(dot(BRot{i}, BRot{i})));% + sum(dot(gradient{i}, gradient{i}))); % normalisation en fonction du produit scalaire utilise
    BBNoRot{i} = BNoRot{i}/sqrt(sum(dot(BNoRot{i}, BNoRot{i})));
end

BBGradRot = cell(numFiles,1);
BBGradNoRot = cell(numFiles,1);
for i = 1:numFiles
    BBGradRot{i} = BGradRot{i}/sqrt(sum(dot(BGradRot{i}, BGradRot{i})));% + sum(dot(gradient{i}, gradient{i}))); % normalisation en fonction du produit scalaire utilise
    BBGradNoRot{i} = BGradNoRot{i}/sqrt(sum(dot(BGradNoRot{i}, BGradNoRot{i})));
end

% for fileNum = numFiles-10 : numFiles
%     fileName = sprintf('PrinCompo.%d.sol',fileNum);  
%     ok = writesol(fileName, BB{fileNum}.');
% end

% b seulement une etape intermediaire de calcul
% B composentes principales de notre analyse
% BB ce sont les B mais normalises


%% Composition lineaire d'un nouveau masseter a partir des composantes principales

% Coefficient pour les combinaisons lineaires
alphaRot = zeros(numFiles, numFiles+1);
alphaNoRot = zeros(numFiles, numFiles+1);
for i = 1:numFiles %numero des vecteurs propres consideres
    for j = 1:numFiles+1 %pour tous les masseters de la base + le 52 qui n'y est pas
        alphaRot(i,j) = sum(dot(massRot{j},BBRot{i}))/sum(dot(BBRot{i},BBRot{i})); %trace(transpose(A)*B) sous forme matriciel...
        alphaNoRot(i,j) = sum(dot(massNoRot{j},BBNoRot{i}))/sum(dot(BBNoRot{i},BBNoRot{i}));
    end
end
%idem Gradient
alphaGradRot = zeros(numFiles, numFiles+1);
alphaGradNoRot = zeros(numFiles, numFiles+1);
for i = 1:numFiles %numero des vecteurs propres consideres
    for j = 1:numFiles+1 %pour tous les masseters de la base + le 52 qui n'y est pas
        alphaGradRot(i,j) = sum(dot(massRot{j},BBGradRot{i}))/sum(dot(BBGradRot{i},BBGradRot{i})); %trace(transpose(A)*B) sous forme matriciel...
        alphaGradNoRot(i,j) = sum(dot(massNoRot{j},BBGradNoRot{i}))/sum(dot(BBGradNoRot{i},BBGradNoRot{i}));
    end
end


% Initialisation des tableaux de cellules / Allocation de memoire
LLRot = cell(numFiles,numFiles+1);
LRot = cell(1, numFiles+1);
LLNoRot = cell(numFiles,numFiles+1);
LNoRot = cell(1, numFiles+1);
meshRot = cell(1, numFiles+1);
meshNoRot = cell(1, numFiles+1);
solRot = cell(1, numFiles+1);
solNoRot = cell(1, numFiles+1);
for i = 1:numFiles+1
   LRot{i} =  zeros(numLines,3);
   LNoRot{i} =  zeros(numLines,3);
   meshRot{i} =  zeros(numLines,3);
   meshNoRot{i} =  zeros(numLines,3);
   solRot{i} =  zeros(numLines,3);
   solNoRot{i} =  zeros(numLines,3);
end
%idem gradient
LLGradRot = cell(numFiles,numFiles+1);
LGradRot = cell(1, numFiles+1);
LLGradNoRot = cell(numFiles,numFiles+1);
LGradNoRot = cell(1, numFiles+1);
meshGradRot = cell(1, numFiles+1);
meshGradNoRot = cell(1, numFiles+1);
solGradRot = cell(1, numFiles+1);
solGradNoRot = cell(1, numFiles+1);
for i = 1:numFiles+1
   LGradRot{i} =  zeros(numLines,3);
   LGradNoRot{i} =  zeros(numLines,3);
   meshGradRot{i} =  zeros(numLines,3);
   meshGradNoRot{i} =  zeros(numLines,3);
   solGradRot{i} =  zeros(numLines,3);
   solGradNoRot{i} =  zeros(numLines,3);
end



% Combinaison lineaire a proprement parlee pour les 70 masseters
for i = 1:numFiles %Numero des vecteur pris en compte
    for j = 1:numFiles+1 % Pour tous les masseters + le 70
        LLRot{i,j} = alphaRot(i,j)*BBRot{i};
        LRot{j} = LRot{j} + LLRot{i,j};
        LLNoRot{i,j} = alphaNoRot(i,j)*BBNoRot{i};
        LNoRot{j} = LNoRot{j} + LLNoRot{i,j};
        meshRot{j} = template + LRot{j};
        meshNoRot{j} = template + LNoRot{j};
        solRot{j} = massRot{j} - LRot{j};
        solNoRot{j} = massNoRot{j} - LNoRot{j};
    end
end
%idem Gradient
for i = 1:numFiles %Numero des vecteur pris en compte
    for j = 1:numFiles+1 % Pour tous les masseters + le 70
        LLGradRot{i,j} = alphaGradRot(i,j)*BBGradRot{i};
        LGradRot{j} = LGradRot{j} + LLGradRot{i,j};
        LLGradNoRot{i,j} = alphaGradNoRot(i,j)*BBGradNoRot{i};
        LGradNoRot{j} = LGradNoRot{j} + LLGradNoRot{i,j};
        meshGradRot{j} = template + LGradRot{j};
        meshGradNoRot{j} = template + LGradNoRot{j};
        solGradRot{j} = massRot{j} - LGradRot{j};
        solGradNoRot{j} = massNoRot{j} - LGradNoRot{j};
    end
end

% Boucle pour les sauvegardes des reconstructions des masseters a
% partir des combinaison lineaire
% for fileNum = 1:numFiles+1 
%     fileName = sprintf('ACPTotRot.%d.mesh',fileNum); 
%     okmeshRot = writemesh(fileName, meshRot{fileNum}.', tri,edg,crn);
%     fileName = sprintf('ACPTotNoRot.%d.mesh',fileNum); 
%     okmeshNoRot = writemesh(fileName, meshNoRot{fileNum}.', tri,edg,crn);
% 	fileName = sprintf('ACPTotRot.%d.sol',fileNum); 
%     oksolRot = writesol(fileName, solRot{fileNum}.');
%     fileName = sprintf('ACPTotNoRot.%d.sol',fileNum); 
%     oksolNoRot = writesol(fileName, solNoRot{fileNum}.');
% end
% idem gradient
% for fileNum = 1:numFiles+1 
%     fileName = sprintf('ACPGradTotRot.%d.mesh',fileNum); 
%     GokmeshRot = writemesh(fileName, meshGradRot{fileNum}.', tri,edg,crn);
%     fileName = sprintf('ACPGradTotNoRot.%d.mesh',fileNum); 
%     GokmeshNoRot = writemesh(fileName, meshGradNoRot{fileNum}.', tri,edg,crn);
% 	fileName = sprintf('ACPGradTotRot.%d.sol',fileNum); 
%     GoksolRot = writesol(fileName, solGradRot{fileNum}.');
%     fileName = sprintf('ACPGradTotNoRot.%d.sol',fileNum); 
%     GoksolNoRot = writesol(fileName, solGradNoRot{fileNum}.');
% end

% alpha(69, 1) coefficient de la premiere composante pour representer le
% masseter numero un.
% L{i} representation du masseter i dans la base des 69 masseters

%% Idem mais avec seulement les n = 5 premieres composantes principales

% Initialisation des tableaux de cellules / Allocation de memoire
nb = 5;%[3, 5, 10, 15, 30];% Nombre de composante principale considerees (10 on commence a pouvoir les differencier)
for n = nb
    NNRot = cell(numFiles,numFiles+1);
    NRot = cell(1, numFiles+1);
    NNNoRot = cell(numFiles,numFiles+1);
    NNoRot = cell(1, numFiles+1);
    nmeshRot = cell(1, numFiles+1);
    nmeshNoRot = cell(1, numFiles+1);
    nsolRot = cell(1, numFiles+1);
    nsolNoRot = cell(1, numFiles+1);
    for i = 1:numFiles+1
        NRot{i} =  zeros(numLines,3);
        NNoRot{i} =  zeros(numLines,3);
        nmeshRot{i} = zeros(numLines,3);
        nmeshNoRot{i} = zeros(numLines,3);
        nsolRot{i} = zeros(numLines,3);
        nsolNoRot{i} = zeros(numLines,3);
    end
    %idem gradient
%     NNGradRot = cell(numFiles,numFiles+1);
%     NGradRot = cell(1, numFiles+1);
%     NNGradNoRot = cell(numFiles,numFiles+1);
%     NGradNoRot = cell(1, numFiles+1);
%     nmeshGradRot = cell(1, numFiles+1);
%     nmeshGradNoRot = cell(1, numFiles+1);
%     nsolGradRot = cell(1, numFiles+1);
%     nsolGradNoRot = cell(1, numFiles+1);
%     for i = 1:numFiles+1
%         NGradRot{i} =  zeros(numLines,3);
%         NGradNoRot{i} =  zeros(numLines,3);
%         nmeshGradRot{i} = zeros(numLines,3);
%         nmeshGradNoRot{i} = zeros(numLines,3);
%         nsolGradRot{i} = zeros(numLines,3);
%         nsolGradNoRot{i} = zeros(numLines,3);
%     end

    % Combinaison lineaire a proprement parlee pour les 52 masseters
    for i = numFiles-n+1:numFiles %Numero des vecteur pris en compte
        for j = 1:numFiles+1 % Pour tous les masseters + le 52
            NNRot{i,j} = alphaRot(i,j)*BBRot{i};
            NRot{j} = NRot{j} + NNRot{i,j};
            NNNoRot{i,j} = alphaNoRot(i,j)*BBNoRot{i};
            NNoRot{j} = NNoRot{j} + NNNoRot{i,j};
            nmeshRot{j} = template + NRot{j};
            nmeshNoRot{j} = template + NNoRot{j};
            nsolRot{j} = LRot{j} - NRot{j};
            nsolNoRot{j} = LNoRot{j} - NNoRot{j};
            %idem gradient
%             NNGradRot{i,j} = alphaGradRot(i,j)*BBGradRot{i};
%             NGradRot{j} = NGradRot{j} + NNGradRot{i,j};
%             NNGradNoRot{i,j} = alphaGradNoRot(i,j)*BBGradNoRot{i};
%             NGradNoRot{j} = NGradNoRot{j} + NNGradNoRot{i,j};
%             nmeshGradRot{j} = template + NGradRot{j};
%             nmeshGradNoRot{j} = template + NGradNoRot{j};
%             nsolGradRot{j} = LGradRot{j} - NGradRot{j};
%             nsolGradNoRot{j} = LGradNoRot{j} - NGradNoRot{j};
        end
    end
    
    % Boucle pour les sauvegardes des reconstructions des masseters a
    % partir des combinaison lineaire des n premiere composantes principales
    for fileNum = 1:numFiles+1
%         fileName = sprintf('ACPRedBas%dRot.%d.mesh',n,fileNum);
%         okkmeshRot = writemesh(fileName, nmeshRot{fileNum}.', tri,edg,crn);
        fileName = sprintf('ACPRedBas%dNoRot.%d.mesh',n,fileNum);
        okkmeshNoRot = writemesh(fileName, nmeshNoRot{fileNum}.', tri,edg,crn);
%         fileName = sprintf('ACPRedBas%dRot.%d.sol',n,fileNum);
%         okksolRot = writesol(fileName, nsolRot{fileNum}.');
        fileName = sprintf('ACPRedBas%dNoRot.%d.sol',n,fileNum);
        okksolNoRot = writesol(fileName, nsolNoRot{fileNum}.');
    end
      % idem gradient
%     for fileNum = 1:numFiles+1
%         fileName = sprintf('ACPGradRedBas%dRot.%d.mesh',n,fileNum);
%         GokkmeshRot = writemesh(fileName, nmeshGradRot{fileNum}.', tri,edg,crn);
%         fileName = sprintf('ACPGradRedBas%dNoRot.%d.mesh',n,fileNum);
%         GokkmeshNoRot = writemesh(fileName, nmeshGradNoRot{fileNum}.', tri,edg,crn);
%         fileName = sprintf('ACPGradRedBas%dRot.%d.sol',n,fileNum);
%         GokksolRot = writesol(fileName, nsolGradRot{fileNum}.');
%         fileName = sprintf('ACPGradRedBas%dNoRot.%d.sol',n,fileNum);
%         GokksolNoRot = writesol(fileName, nsolGradNoRot{fileNum}.');
%     end
end
%% Comparaison avec et sans recalage
%
% décroissance des valeurs propres
% eigGradRecal=diag(DGradRot);
% eigRecal=diag(DRot);
% x=70:-1:1;
% hold on
% plot(x,log10(eigRecal))
% plot(x,log10(eigGradRecal))
% hold off
%
