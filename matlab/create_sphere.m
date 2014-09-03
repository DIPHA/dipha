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

function filename = create_sphere( num_points, dimension )
    %% create filename based on parameters
    filename = ['sphere_' num2str( dimension ) '_' num2str( num_points ) '.complex'];

    %% set seed for random number generator
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',pi));
    %RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
    
    %% create actual data
    points = zeros( dimension, num_points );
    cur_num_points = 0;
    while cur_num_points < num_points
        random_point = 2 * (rand( dimension, 1 ) - 0.5 );
        if norm(random_point) < 1
            cur_num_points = cur_num_points + 1;
            points(:, cur_num_points) = random_point ./ norm( random_point );
        end
    end 
    
    %% compute distance matrix
    distance_matrix = zeros( num_points );
    for i = 1:num_points
        for j = 1:num_points
            distance_matrix(i,j) = norm( points(:, i) - points(:, j));
        end
    end
    
    %% save to disk in DIPHA format
    save_distance_matrix( distance_matrix, filename );
end