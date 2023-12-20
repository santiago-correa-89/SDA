function [ aoa, cl, cd, cm ] = readairfoildata_v3( polars , path )

format       = table2array(readtable( fullfile( path, polars{3} ) )) ;
[ row, col ] = size(format)          ;

for i = 1:length(polars)
    if i == 1 || i == 2   % First two columns are cylinder shapes
        aoa(:,i)   =    format(:,1) ;  
        for j = 1:row
            strucdata    = table2array( readtable( fullfile( path, polars{i} ) ) ) ;
            cl(j,i)      = strucdata(1,2);
            cd(j,i)      = strucdata(1,3);
            cm(j,i)      = strucdata(1,4);
            %fstat(:,i)      = strucdata(1,4);
            %clinv(:,i)      = strucdata(1,4);
            %clfullysep(:,i) = strucdata(1,4);
        end
    else
        strucdata       = table2array( readtable( fullfile( path, polars{i} )) );
        if length(strucdata(:,1)) ~= row
            [uniqueValues, ~, index] = unique(strucdata(:, 1));
            rowsToKeep               = accumarray(index, 1) == 1;
            filteredData             = strucdata(rowsToKeep, :) ;
            aoa(:,i)                 = format(:, 1) ;
            for j = 1:row
                cl(j,i)    =  interp1( filteredData(:,1), filteredData(:,2),  format(j, 1), 'linear', 'extrap' );
                cd(j,i)    =  interp1( filteredData(:,1), filteredData(:,3),  format(j, 1), 'linear', 'extrap' );
                cm(j,i)    =  interp1( filteredData(:,1), filteredData(:,4),  format(j, 1), 'linear', 'extrap' );
            end
        else
            aoa(:,i)        = strucdata(:,1) ;
            cl(:,i)         = strucdata(:,2) ;
            cd(:,i)         = strucdata(:,3) ;
            cm(:,i)         = strucdata(:,4) ;
            %fstat(:,i)      = strucdata(:,4);
            %clinv(:,i)      = strucdata(:,4);
            %clfullysep(:,i) = strucdata(:,4);
        end
    end
end
end