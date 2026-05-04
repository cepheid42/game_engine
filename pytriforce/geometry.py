
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from ..utility import get_rotation_matrix

# set serif font
plt.rc('font', family='serif',serif = "cmr10", size=18)
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['text.usetex']=True
plt.rcParams['axes.formatter.use_mathtext']=True

# basic approach to 3D geometry:
# create 2D shape, then extrude

class Polygon:
    EPSILON = 1e-10
    def __init__(self,num_vertices=4, radius=0.1, order=1, regular=True):

        self.order = order
        # make a regular polygon with vertices at the specified radius, starting point on x-axis, in xy-plane
        angles = np.linspace(0, 2*np.pi, num_vertices+1)
        points = [np.cos(angles)*radius*(1-self.__class__.EPSILON),np.sin(angles)*radius*(1-self.__class__.EPSILON),np.zeros_like(angles)]
        self.spline = None
        self.u = None

        self.reticulate_splines(points)

    @property
    def center(self):
        return np.mean(self.points[:,:-1], axis=1)

    @property
    def max_edge_to_center_distance(self):
        return np.sqrt(np.max(np.sum((self.points - np.atleast_2d(self.center).T)**2,axis=0)))

    @property
    def perimeter(self):
        return self.side_lengths.sum()

    @property
    def points(self):
        return self.spline(self.u)

    @property
    def side_lengths(self):
        # todo: make this work for k > 1 splines
        points = self.spline(self.u)
        diffs = np.diff(points, axis=1)
        return np.sqrt(np.sum(diffs**2,axis=0))

    def reticulate_splines(self, points):
        self.spline, self.u = interpolate.make_splprep(points,k = self.order)

    def rotate(self, angle, axis = (0, 0, 1)):
        """
        Performs a rotation of the polygon about the specified axis; the polygon is re-parameterized afterward as well

        Inputs:
         - angle (double): the angle for the rotation
         - axis_vector (ndarray): a vector indicating the axis about which the rotation is to take place; normalization not necessary

        Outputs:
         - (none)
        """

        self.reticulate_splines(get_rotation_matrix(angle, axis).dot(self.points))

    def translate(self,displacement):
        self.reticulate_splines(self.points + (np.atleast_2d(displacement)*(1-self.__class__.EPSILON)).T)

class ExtrudedShape:
    def __init__(self, shapes, n_splines = 0, order=1):
        self.shapes = shapes

        self.distance_centers = np.sqrt(((np.diff(self.centers, axis=0))**2).sum(axis=1))

        self.order = order
        self.n_splines = n_splines
        self.splines = None
        self.us = None
        self.half_perimeter_locations = None
        self.half_perimeter_distances = None
        self.half_perimeter = None
        self.half_perimeter_param = None

        self.boolean_map = None
        self.max_sides = np.max([len(shape.u) for shape in self.shapes])
        self.max_perimeter = np.max([shape.perimeter for shape in self.shapes])

        if self.n_splines == 0:
            self.n_splines = 3 * (self.max_sides - 1) + 1

        self.d_top = self.shapes[-1].max_edge_to_center_distance
        self.d_bottom = self.shapes[0].max_edge_to_center_distance

        self.reticulate_splines()

    def recalculate_perimeter_values(self):
        self.half_perimeter_locations = np.array([0, self.d_bottom, *(self.d_bottom + self.distance_centers.cumsum()), self.d_bottom + self.distance_centers.cumsum()[-1] + self.d_top])
        self.half_perimeter_distances = np.diff(self.half_perimeter_locations)
        self.half_perimeter = self.half_perimeter_distances.sum()
        self.half_perimeter_param = self.half_perimeter_locations / self.half_perimeter

    @property
    def centers(self):
        return np.array([shape.center for shape in self.shapes])

    def rotate(self, angle, axis):
        for shape in self.shapes:
            shape.rotate(angle, axis)

        self.reticulate_splines()

    def translate(self, displacement):
        for shape in self.shapes:
            shape.translate(displacement)

        self.reticulate_splines()

    def reticulate_splines(self):

        self.recalculate_perimeter_values()

        self.splines = []
        self.us = []

        s = np.linspace(0, 1, self.n_splines)
        points_edge = np.array([shape.spline(s) for shape in self.shapes])

        for i in range(self.n_splines):
            points = np.concatenate([np.atleast_2d(self.centers[0,:]), points_edge[:,:,i], np.atleast_2d(self.centers[-1,:])]).T

            spline, u = interpolate.make_splprep(points, u = self.half_perimeter_param, k = self.order)
            self.splines.append(spline)
            self.us.append(u)


    def compute_resolution(self, res = 0.01):
        # side_lengths = self.shape_bottom.side_lengths
        # side_u = self.shape_bottom.u
        # s_array = []
        # for i in range(len(side_u)-1):
        #     s_array.append(np.linspace(side_u[i],side_u[i+1],int(2 * side_lengths[i] / res), endpoint=False))

        s_array = np.linspace(0, 1, int(2 * self.max_perimeter / res))

        # s_array.append([1.0])
        #
        # s_array=np.concatenate(s_array)

        t_array = [np.linspace(self.half_perimeter_param[i],self.half_perimeter_param[i+1],int(2 * self.half_perimeter_distances[i] / res), endpoint=False) for i in range(len(self.half_perimeter_param)-1)]

        t_array.append(np.array([1.0]))

        t_array=np.concatenate(t_array)

        return s_array, t_array

    def voxelize(self, grid, res = None):

        x, y, z = grid["x"], grid["y"], grid["z"]

        if res is None:
            res = grid["deltas"][0]

        s_array, t_array = self.compute_resolution(res)

        # regenerate the connecting splines
        self.n_splines = len(s_array)
        self.reticulate_splines()

        ts = t_array

        points = [spline(ts) for spline in self.splines]

        points = np.hstack(points)

        self.boolean_map=np.zeros(x.shape,dtype=bool)

        xid=np.searchsorted(x[:,0,0],points[0,:]+1e-16,'left')
        yid=np.searchsorted(y[0,:,0],points[1,:]+1e-16,'left')
        zid=np.searchsorted(z[0,0,:],points[2,:]+1e-16,'left')

        xid[xid==x.shape[0]]-=1
        yid[yid==y.shape[1]]-=1
        zid[zid==z.shape[2]]-=1

        xid-=1
        yid-=1
        zid-=1


        # fill in the edges
        self.boolean_map[xid,yid,zid]=True

        # fill in the center
        # convex shape is assumed so this is straightforward
        for k in range(len(z[0,0,:])):
            for j in range(len(y[0,:,0])):
                nonzeros=np.nonzero(self.boolean_map[:,j,k])[0]
                if len(nonzeros)>0:
                    start=nonzeros[0]
                    end=nonzeros[-1]
                    self.boolean_map[start:end, j, k]=True

class Composite3D:
    def __init__(self,shapes=(),operations=None):
        self.shapes=shapes
        self.operations=operations
        self.boolean_map=None

        if self.operations is None:
            self.operations=["or"]*(len(self.shapes)-1)

        if isinstance(self.operations, str):
            self.operations = [self.operations]

    def rotate(self,angle,axis=(0,0,1)):
        for shape in self.shapes:
            shape.rotate(angle,axis)

    def translate(self,offset):
        for shape in self.shapes:
            shape.translate(offset)

    def voxelize(self,grid):
        for shape in self.shapes:
            shape.voxelize(grid)

        self.boolean_map=self.shapes[0].boolean_map
        for k,operation in enumerate(self.operations):
            if operation=="or" or operation=="add":
                self.boolean_map = self.boolean_map | self.shapes[k+1].boolean_map
            elif operation=="and":
                self.boolean_map = self.boolean_map & self.shapes[k+1].boolean_map
            elif operation=="xor":
                self.boolean_map = self.boolean_map ^ self.shapes[k+1].boolean_map
            elif operation=="sub":
                self.boolean_map = self.boolean_map & ~self.shapes[k+1].boolean_map

class RectangularPrism(ExtrudedShape):
    def __init__(self, dimensions = (1.0,1.0,1.0), center = (0.0,0.0,0.0)):
        p1 = Polygon(num_vertices=4, radius = dimensions[0]*np.sqrt(2)/2)
        p1.rotate(np.pi/4)
        p1.translate((center[0], center[1], center[2] - dimensions[2]/2))

        p2 = Polygon(num_vertices=4, radius = dimensions[0]*np.sqrt(2)/2)
        p2.rotate(np.pi/4)
        p2.translate((center[0], center[1], center[2] + dimensions[2]/2))
        super(RectangularPrism, self).__init__([p1,p2])

class Cylinder(ExtrudedShape):
    def __init__(self, radius = 1.0, length = 1.0, center = (0.0,0.0,0.0), num_verticies = 60, centered = True):
        p1 = Polygon(num_vertices = num_verticies, radius = radius)
        p1.translate((center[0], center[1], center[2] - length/2 * centered))

        p2 = Polygon(num_vertices = num_verticies, radius = radius)
        p2.translate((center[0], center[1], center[2] + length/2 + length/2 * (not centered)))
        super(Cylinder, self).__init__([p1,p2])

class Ring(Composite3D):
    def __init__(self,outer_radius=1.0,inner_radius=0.5,length=1.0,centered=True):
        super(Ring, self).__init__(shapes=[Cylinder(radius=outer_radius,length=length,centered=centered), Cylinder(radius=inner_radius,length=length*1.5,centered=centered)], operations=["sub"])
