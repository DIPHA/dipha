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

function create_sphere_extrinsic_sparse_rips_data( num_points, dimension, threshold )
    %% create filename based on parameters
    filename = ['sphere_' num2str( dimension ) '_' num2str( num_points ) '_' num2str( threshold ) '.complex'];

    %% create actual data
    random_points = rand( dimension, num_points ) - 0.5;
    normalization_factor = 1 ./ sqrt( sum( random_points.^2 ) );
    points = random_points .* repmat( normalization_factor, dimension, 1 );

    %% save to disk in DIPHA format
    save_extrinsic_sparse_rips_complex( points, size( points, 1 ), threshold, filename );
end