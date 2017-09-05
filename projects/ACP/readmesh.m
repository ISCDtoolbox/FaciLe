function [coor,tri,edg,crn] = readmesh(name)
% READMESH : open the mesh file name
%            outputs the vertices and 
%            triangles arrays
%  example: [coor,tri,edg,crn] = readmesh('toto.mesh')
%   


coor = [];
tri  = [];
crn  = [];
edg  = [];


fid = fopen(name,'r');
if ( fid == -1 ) 
 error(['Mesh file ' name ' does not exist']);
else
 disp([ '% ' name ' OPENED ']);
end


while ( ~feof(fid) )

  str = fgetl( fid );
  str = str(find(str~=' '));
  switch ( lower(str) )
  
  case 'dimension'
    dim = fscanf( fid, '%d', 1 );
    if ( dim ~= 3 ) % Mettre 2 pour la 2D (j'imagine puiqu'en mettant 3 за fonctionne pour la 3D ...)
      error(' Invalid inout mesh ');
    end  

  case 'vertices'
    NbrVer = fscanf( fid, '%d', 1 );
%     if ( dim == 2 )
%       coor = fscanf( fid, '%f %f %d', [3, NbrVer] ); % Fonctionne pas pour une obscure raison ...
%       coor = coor(1:2,:);
%     else
      coor = fscanf( fid, '%f %f %f %d', [4, NbrVer] );
      coor = coor(1:3,:);
%    end  
   
  case 'triangles'
    NbrTri = fscanf( fid, '%d', 1 );
    tri = fscanf( fid, '%d %d %d %d', [4, NbrTri] );
    tri = tri(1:3,:);
  
  case 'edges'
    NbrEdg = fscanf( fid, '%d', 1 );
    edg = fscanf( fid, '%d %d %d', [3, NbrEdg] );
  
  case 'corners'
    NbrCor = fscanf( fid, '%d', 1 );
    crn = fscanf( fid, '%d ', [1, NbrCor] );

  

  end
end

fclose(fid);
end

