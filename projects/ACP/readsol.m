function sol = readsol(name)
% READSOL : open the solution file name
%           outputs the solution array
%  example: sol = readsol('toto.sol')
% 

sol = [];
dim = 3;
  
fid = fopen(name,'r');
if ( fid == -1 ) 
 error(['Solution file ' name ' does not exist']);
else
 disp([ '% '  name ' OPENED ']); 
end

SizFld2 = [1,2,3];
SizFld3 = [2,3,6];


while ( ~feof(fid) )

  str = fgetl( fid );
  str = str(find(str~=' '));
  switch ( lower(str) )

  case 'solatvertices'
    NbrLin = fscanf( fid, '%d', 1 );
    NbrSol = fscanf( fid, '%d', 1 );
    SolTyp = fscanf( fid, ' %d ', NbrSol );
    if ( dim == 2 )
      SizFld = sum(SizFld2(SolTyp));
    else
      SizFld = sum(SizFld3(SolTyp));
    end
    sol = zeros(SizFld,NbrLin);
    sol = fscanf( fid, ' %f ', [SizFld, NbrLin] );
   
  end
end


fclose(fid);
end

