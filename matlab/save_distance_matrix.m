%  Copyright 2014 IST Austria
%
%  Contributed by: Jan Reininghaus
%
%  This file is part of DIPHA.
%
%  DIPHA is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  DIPHA is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with DIPHA.  If not, see <http://www.gnu.org/licenses/>.

function save_distance_matrix( distance_matrix, filename )
    %% open file for writing
    fid = fopen( filename, 'w' );

    %% DIPHA magic number
    fwrite( fid, 8067171840, 'int64' );

    %% file type identifier
    fwrite( fid, 7, 'int64' );

    %% total number of input points n 
    fwrite( fid, size( distance_matrix, 2 ), 'int64' );

    %% floating point values for the coordinates of the points
    fwrite( fid, distance_matrix(:), 'double' );

    %% close file
    fclose(fid);
end
