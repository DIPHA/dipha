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

function filename = create_smooth_image_data( resolution )
    %% create filename based on parameters
    filename = ['smooth_' num2str( resolution ) '.complex'];
    
    spacing = linspace(-2,2,resolution);
    [X,Y,Z] = meshgrid(spacing,spacing,spacing);
    
    data =  1 .* sin(1 .* X) .* sin( 2 .* Y ) .* sin( 3 .* Z ) + ...
            2 .* sin(2 .* X) .* sin( 1 .* Y ) .* sin( 3 .* Z ) + ...
            3 .* sin(3 .* X) .* sin( 2 .* Y ) .* sin( 1 .* Z ) + ...
            4 .* sin(1 .* X) .* sin( 3 .* Y ) .* sin( 2 .* Z ) + ...
            5 .* sin(2 .* X) .* sin( 3 .* Y ) .* sin( 1 .* Z ) + ...
            6 .* sin(3 .* X) .* sin( 1 .* Y ) .* sin( 2 .* Z ) + ...
            1 .* cos(3 .* X) .* cos( 1 .* Y ) .* cos( 2 .* Z ) + ...
            2 .* cos(2 .* X) .* cos( 1 .* Y ) .* cos( 3 .* Z ) + ...
            3 .* cos(1 .* X) .* cos( 2 .* Y ) .* cos( 3 .* Z ) + ...
            4 .* cos(3 .* X) .* cos( 2 .* Y ) .* cos( 1 .* Z ) + ...
            5 .* cos(2 .* X) .* cos( 3 .* Y ) .* cos( 1 .* Z ) + ...
            6 .* cos(1 .* X) .* cos( 3 .* Y ) .* cos( 2 .* Z );
    
    %% save to disk in DIPHA format
    save_image_data( data, filename );
end
