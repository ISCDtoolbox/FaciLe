function ok = writemesh(name,coor,tri,edg,crn)
% WRITEMESH : write the mesh file name
%  example: ok = writemesh('toto.mesh',coor,tri,edg,crn)
%           coor vertex coordinates 
%           tri  triangles array 
%           edg edges array [is1,is2,ref]  
%           crn corners array 
%        

  k = strfind(name,'.mesh'); 
  if ( length(k) ~= 0 )
    in  = [ name(1:(k(1)-1)) '.mesh'];
  else
    in  = [ name '.mesh'];
  end  
  
  dim = 3; % Mettre 2 pour la 2D (j'imagine puiqu'en mettant 3 за fonctionne pour la 3D ...)
  ok = 0;

  if ( size(tri,1) ~= 3 )
    error('Wrong shaped triangles tab must be [1:3,1:NbrTri]');
  end
  
  if ( size(coor,1) ~= 3 ) %idem => selon la dimention
    error('Wrong shaped vertex tab must be [1:3,1:NbrVer]');
  end

  fid = fopen(in,'w+');
  if ( fid == -1 ) 
    error(['Cannot create mesh file ' in ]);
  else
    disp([ '% '  in ' CREATED ']);
  end


fprintf(fid,'MeshVersionFormatted 1\n');
fprintf(fid,'Dimension\n3\n');
fprintf(fid,'Vertices\n%d\n',size(coor,2));
fprintf(fid,' %f %f %f %d \n',[coor ; zeros(1,size(coor,2))]); %3D
fprintf(fid,'Triangles\n%d\n',size(tri,2));
fprintf(fid,' %d %d %d  %d\n',[tri ; zeros(1,size(tri,2))]);
fprintf(fid,'Edges\n%d\n',size(edg,2));
fprintf(fid,' %d %d %d \n',edg);
fprintf(fid,'Ridges\n%d\n',size(edg,2));
fprintf(fid,' %d \n', (1:1:size(edg,2)));
fprintf(fid,'Corners\n%d\n',length(crn));
fprintf(fid,' %d \n', crn);

fprintf(fid,'\nEnd\n');  
fclose(fid);