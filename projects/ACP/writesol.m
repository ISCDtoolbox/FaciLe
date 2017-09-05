function ok = writesol(name,sol)
% WRITESOL : create the solution file name
%  example: ok = writesol('toto.sol',sol)
%   


dim = 3;

ok = 0;

k = strfind(name,'.sol'); 
if ( length(k) ~= 0 )
  in  = [ name(1:(k(1)-1)) '.sol'];
else
  in  = [ name '.sol'];
end  

fid = fopen(in,'w+');
if ( fid == -1 ) 
 error(['Cannot create solution file ' in ]);
else
 disp([ '% '  in ' CREATED ']);
end



fprintf(fid,'MeshVersionFormatted 1\n\n');

if ( dim == 3 )

  fprintf(fid,'Dimension\n3\n\n');
  if ( size(sol,1) ~= 4   )
    %fprintf(fid,'SolAtVertices\n%d\n1 %d\n',[size(sol,2) size(sol,1)]);
    fprintf(fid,'SolAtVertices\n%d\n1 2\n\n',size(sol,2));
  else
    fprintf(fid,'SolAtVertices\n%d\n3 1 2 1\n',size(sol,2));
  end  
    
  if ( size(sol,1) == 1 )
    fprintf(fid, ' %f \n', sol );
  elseif ( size(sol,1) == 2 )
    fprintf(fid, ' %f %f\n', sol );
  elseif ( size(sol,1) == 3 )
    fprintf(fid, ' %f %f %f\n', sol );
  elseif ( size(sol,1) == 4 )
    fprintf(fid, ' %f %f %f %f \n', sol );
  else
    error(' Wrong solution type ');  
  end

else
  error(' Invalid solution ');
  
  
end      
  
fprintf(fid,'\nEnd\n');
fclose(fid);

end

