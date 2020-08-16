# This code is ported to Julia from ASTRA-toolbox
# /*
# -----------------------------------------------------------------------
# Copyright: 2010-2018, imec Vision Lab, University of Antwerp
#            2014-2018, CWI, Amsterdam
# Contact: astra@astra-toolbox.com
# Website: http://www.astra-toolbox.com/
# This file is part of the ASTRA Toolbox.
# The ASTRA Toolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# The ASTRA Toolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------
# */

function geom_2vec_parallel2d(DetectorWidth::T, ProjectionAngles) where {T<:AbstractFloat}
#     ProjectionAngles = ProjectionAngles
    
    vectors = zeros(T, length(ProjectionAngles), 6);
    
    for (i, θ) in enumerate(ProjectionAngles)
        vectors[i,1] = sin(θ);
        vectors[i,2] = -cos(θ);

        vectors[i,3] = 0;
        vectors[i,4] = 0;

        vectors[i,5] = cos(θ) * DetectorWidth;
        vectors[i,6] = sin(θ) * DetectorWidth;

    end
    return vectors
end

function geom_2vec_parallel3d(DetectorSpacingX::T, DetectorSpacingY::T, ProjectionAngles) where {T<:AbstractFloat}
    vectors = zeros(Float32, length(ProjectionAngles), 12);
    
    for (i, θ) in enumerate(ProjectionAngles)
        # ray direction
        vectors[i, 1] =  sin(θ)
        vectors[i, 2] = -cos(θ)
        vectors[i, 3] = 0

        # center of detector
        vectors[i, 4] = 0
        vectors[i, 5] = 0
        vectors[i, 6] = 0

        # vector from detector pixel (0,0) to (0,1)
        vectors[i, 7] = cos(θ) * DetectorSpacingX
        vectors[i, 8] = sin(θ) * DetectorSpacingX
        vectors[i, 9] = 0

        # vector from detector pixel (0,0) to (1,0)
        vectors[i, 10] = 0
        vectors[i, 11] = 0
        vectors[i, 12] = DetectorSpacingY

    end
    return vectors
end

function geom_2vec_fan(DetectorWidth, ProjectionAngles, DistanceOriginSource, DistanceOriginDetector)
    vectors = zeros(length(ProjectionAngles), 6)
    
    for (i, θ) in enumerate(ProjectionAngles)
        vectors[i,1] =  sin(θ) * DistanceOriginSource
        vectors[i,2] = -cos(θ) * DistanceOriginSource

        vectors[i,3] = -sin(θ) * DistanceOriginDetector
        vectors[i,4] =  cos(θ) * DistanceOriginDetector

        vectors[i,5] =  cos(θ) * DetectorWidth;
        vectors[i,6] =  sin(θ) * DetectorWidth;
    end
    return vectors
end

function geom_2vec_cone(DetectorSpacingX, DetectorSpacingY, ProjectionAngles, DistanceOriginSource, DistanceOriginDetector)
    vectors = zeros(Float32, length(ProjectionAngles), 12);
    
    for (i, θ) in enumerate(ProjectionAngles)
        # ray direction
        vectors[i, 1] =  sin(θ) * DistanceOriginSource
        vectors[i, 2] = -cos(θ) * DistanceOriginSource
        vectors[i, 3] = 0

        # center of detector
        vectors[i, 4] = -sin(θ) * DistanceOriginDetector
        vectors[i, 5] =  cos(θ) * DistanceOriginDetector
        vectors[i, 6] = 0

        # vector from detector pixel (0,0) to (0,1)
        vectors[i, 7] = cos(θ) * DetectorSpacingX
        vectors[i, 8] = sin(θ) * DetectorSpacingX
        vectors[i, 9] = 0

        # vector from detector pixel (0,0) to (1,0)
        vectors[i, 10] = 0
        vectors[i, 11] = 0
        vectors[i, 12] = DetectorSpacingY

    end
    return vectors
end
