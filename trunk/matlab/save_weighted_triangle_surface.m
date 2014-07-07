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

% triangles: num_triangles x 3 array of triangles with vertex indices, vertex_values: num_vertices x 1 array of vertex values
function save_weighted_triangle_surface( triangles, vertex_values, filename )
    %% convert into MATLAB format for sparse(...) etc.
    triangles = double(triangles);
    vertex_values = double(vertex_values);      

    %% connectivity
    e4n = get_e4n( triangles );
    ed4n = get_ed4n( e4n );
    ed4e = get_ed4e( triangles, ed4n );
    n4ed = get_n4ed( ed4n );
    
    %% number of simplices
    num_triangles = size(triangles, 1);
    num_edges = size(n4ed, 2);
    num_vertices = size(vertex_values,1);
    num_simplices = num_vertices + num_edges + num_triangles;

    %% dimensions
    dim = 2;
    dims = [ zeros(num_vertices, 1); ones(num_edges, 1); 2 * ones(num_triangles, 1) ];
    
    %% values
    edge_values = max( vertex_values(n4ed) )';
    triangle_values = max( edge_values(ed4e) )';
    values = [ vertex_values; edge_values; triangle_values ];
    
    %% offsets
    edge_offsets = (2 * (1:num_edges) - 2)';
    triangle_offsets = (edge_offsets(end) + 2 + 3 * (1:num_triangles) - 3)';
    offsets = [ zeros(num_vertices, 1); edge_offsets; triangle_offsets ];
    
    %% entries
    entries = [ n4ed(:); (ed4e(:) + num_vertices) ] - 1;
    num_entries = length(entries);
    
    %% write data to file
    fid = fopen( filename, 'w' );
    fwrite( fid, 8067171840, 'int64' );
    fwrite( fid, 0, 'int64' );
    fwrite( fid, 0, 'int64' );
    fwrite( fid, num_simplices, 'int64' );
    fwrite( fid, dim, 'int64' );
    fwrite( fid, dims, 'int64' );
    fwrite( fid, values, 'double' );
    fwrite( fid, offsets, 'int64' );
    fwrite( fid, num_entries, 'int64' );
    fwrite( fid, entries, 'int64' );
    fclose(fid);
end

function n4ed = get_n4ed( ed4n )
    [I, J, S] = find( triu(ed4n) );
    n4ed(:, S) = [I, J]';
end

function ed4e = get_ed4e(n4e, ed4n)  
    ed4e(1,:) = full( ed4n( sub2ind( size(ed4n), n4e(:,2), n4e(:,1) ) ) );
    ed4e(2,:) = full( ed4n( sub2ind( size(ed4n), n4e(:,3), n4e(:,2) ) ) );
    ed4e(3,:) = full( ed4n( sub2ind( size(ed4n), n4e(:,1), n4e(:,3) ) ) );
end

function ed4n = get_ed4n( e4n )
    temp = triu(e4n+e4n');
    [I,J] = find(temp);
    ed4n = sparse(I,J,1:length(I),size(e4n,1),size(e4n,1));
    ed4n = ed4n+ed4n';
end

function e4n = get_e4n( n4e )
    num_triangles = size(n4e,1);
    I = [n4e(:,1); n4e(:,2); n4e(:,3)];
    J = [n4e(:,2); n4e(:,3); n4e(:,1)];
    S = [1:num_triangles, 1:num_triangles, 1:num_triangles];
    e4n = sparse(I,J,S);
end