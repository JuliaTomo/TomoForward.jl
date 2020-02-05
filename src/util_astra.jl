function geom_2vec_parallel2d(DetectorWidth, ProjectionAngles)
    """ Ported from ASTRA toolbox
    https://github.com/astra-toolbox/astra-toolbox/blob/10d87f45bc9311c0408e4cacdec587eff8bc37f8/matlab/tools/astra_geom_2vec.m
    """
    
#     ProjectionAngles = ProjectionAngles
    
    vectors = zeros(length(ProjectionAngles), 6);
    
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

function geom_2vec_parallel3d(DetectorSpacingX, DetectorSpacingY, ProjectionAngles)
    """Ported form ASTRA toolbox
    https://github.com/astra-toolbox/astra-toolbox/blob/10d87f45bc9311c0408e4cacdec587eff8bc37f8/matlab/tools/astra_geom_2vec.m
    """    
    vectors = zeros(length(ProjectionAngles), 12);
    
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

function geom_2vec_cone(DetectorSpacingX, DetectorSpacingY, ProjectionAngles, DistanceOriginSource, DistanceOriginDetector)
    """Ported form ASTRA toolbox
    https://github.com/astra-toolbox/astra-toolbox/blob/10d87f45bc9311c0408e4cacdec587eff8bc37f8/matlab/tools/astra_geom_2vec.m
    """
    vectors = zeros(length(ProjectionAngles), 12);
    
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