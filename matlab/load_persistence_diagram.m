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

function [dims, birth_values, death_values] = load_persistence_diagram( filename )
    %% open file for reading
    fid = fopen( filename, 'r' );

    %% make sure it is a DIPHA file
    dipha_identifier = fread( fid, 1, 'int64' );
    assert( dipha_identifier == 8067171840, 'input is not a DIPHA file' );

    %% make sure it is a persistence_diagram file
    diagram_identifier = fread( fid, 1, 'int64' );
    assert( diagram_identifier == 2, 'input is not a persistence_diagram file' );

    %% read actual data from file
    num_pairs = fread( fid, 1, 'int64' );
    dims = fread( fid, 3 * num_pairs, 'int64', 16 );
    fseek( fid, 32, 'bof' );
    birth_values = fread( fid, 3 * num_pairs, 'double', 16 );
    fseek( fid, 40, 'bof' );
    death_values = fread( fid, 3 * num_pairs, 'double', 16 );

    %% close file
    fclose( fid );    
end