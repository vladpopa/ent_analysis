import numpy

nt, nz, ny, nx = 4320, 121, 96, 96

dt, dz, dy, dx = 300., 40., 100., 100.

ug = 0. 
vg = 0.

model_dir = '/tera/vpopa/data_analysis' 
analysis_dir = '/tera/vpopa/bomex/analysis' 
data_dir = '/tera/vpopa/bomex/data/' 

def index_to_zyx(index):
    z = index / (ny*nx)
    index = index % (ny*nx)
    y = index / nx
    x = index % nx
    return numpy.array((z, y, x))
                
def zyx_to_index(z, y, x):
    return ny*nx*z + nx*y + x

def index_to_array_3d(index):
    temp_array = numpy.zeros((nz*ny*nx,), numpy.bool) 
    temp_array[index] = 1
    return temp_array.reshape((nz, ny, nx))

def calc_com(mask):
    pts = index_to_zyx( mask )

    z = pts[0,:].astype(float).mean()
    # Correct Center of Mass for reentrant domain
    y1 = pts[1,:].astype(float)
    x1 = pts[2,:].astype(float)
    y2 = (y1 < ny/2.)*y1 + (y1>= ny/2.)*(y1 - ny)
    x2 = (x1 < nx/2.)*x1 + (x1>= nx/2.)*(x1 - nx)
    y1m = y1.mean()
    y2m = y2.mean()
    x1m = x1.mean()
    x2m = x2.mean()
    
    if numpy.var(y2 - y2m) > numpy.var(y1 - y1m):
        y = y1m
    else:
        y = (y2m + .5)%ny - .5
        
    if numpy.var(x2 - x2m) > numpy.var(x1 - x1m):
        x = x1m
    else:
        x = (x2m + .5)%nx - .5
        
    return numpy.array((z, y, x))

def calc_distance(point1, point2):
    # Calculate distances corrected for reentrant domain
    delta_x = numpy.abs(point2[2] - point1[2])
    if delta_x >= (nx/2): delta_x = nx - delta_x
    delta_y = numpy.abs(point2[1] - point1[1])
    if delta_y >= (ny/2): delta_y = ny-delta_y
    delta_z = point2[0] - point1[0]
    return numpy.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

def expand_indexes(indexes):
    # Expand a given set of indexes to include the nearest
    # neighbour points in all directions.
    # indexes is an array of grid indexes
                    
    K_J_I = index_to_zyx( indexes )

    stack_list = [K_J_I, ]
    for item in ((-1, 0, 0), (1, 0, 0),
                 (0, -1, 0), (0, 1, 0), 
                 (0, 0, -1), (0, 0, 1)):
        stack_list.append( K_J_I + numpy.array(item)[:, numpy.newaxis] )
    
    expanded_index = numpy.hstack(stack_list)

    # re-entrant domain
    expanded_index[0, expanded_index[0, :] == nz] = nz-1
    expanded_index[0, expanded_index[0, :] < 0] = 0
    expanded_index[1, :] = expanded_index[1, :]%ny
    expanded_index[2, :] = expanded_index[2, :]%nx

    # convert back to indexes
    expanded_index = zyx_to_index(expanded_index[0, :],
                                  expanded_index[1, :],
                                  expanded_index[2, :])
                                  
    expanded_index = numpy.unique(expanded_index)
    
    return expanded_index


