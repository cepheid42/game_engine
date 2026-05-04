from simsopt.field import CircularCoil
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from utilities import get_rotation_matrix


class ConstantField:
    def __init__(self,values = (0.0,0.0,0.0), field_type = "B"):
        self.values = values
        self.field_type = field_type
        self.points = None

    def set_points(self,array):
        self.points = array.copy()

    def B(self):
        return np.full(self.points.shape, self.values)

    def A(self):
        return np.full(self.points.shape, self.values)


class CoilSet:
    """
    Base class for magnetic coil sets using CircularCoil

    Contains rotation, translation, and assembly functions
    """

    def __init__(self):
        self.coils = []
        self.spline = None
        self.coils_s = None

    def compute_lineout_spline(self):
        """
        Computes a parametrized piecewise cubic spline through the centers of the coils

        """
        if len(self.coils) == 0:
            return

        points=np.zeros((3,len(self.coils)))

        for i, coil in enumerate(self.coils):
            points[:,i]=coil.center

        self.spline, self.coils_s = interpolate.make_splprep(points, s=0)

    def rotate(self, angle, axis = (0.0, 1.0, 0.0)):
        """
        Rotates all coils in the coil set around the specified axis

        Note: recreates each coil on the fly, as the 'normal' parameter doesn't seem to be editable after CircularCoil creation

        :param angle: The rotation angle in radians
        :param axis: The axis around which the rotation is done; (0.0, 1.0, 0.0), y-axis by default
        :return: None
        """
        rotation_matrix = get_rotation_matrix(angle, axis)
        for i, coil in enumerate(self.coils):
            self.coils[i] = CircularCoil(I=coil.I, r0=coil.r0, center=rotation_matrix @ coil.center, normal=coil.normal + np.array([0, angle]))

        self.compute_lineout_spline()

    def translate(self, displacement):
        """
        Translates all coils in the coilset the specified displacement

        :param displacement: The distance to shift in meters
        :return: None
        """
        displacement = np.array(displacement)
        for coil in self.coils:
            coil.center += displacement

        self.compute_lineout_spline()

    def assemble(self):
        """
        Adds up all the coils in the coilset

        :return: MagneticFieldSum object
        """
        final_sum = self.coils[0]
        for i, coil in enumerate(self.coils):
            if i == 0: continue
            final_sum += coil
        return final_sum


class GDMTCoils(CoilSet):
    """
    Gas Dynamic Mulit-Mirror Trap coil set class

    :param currents_main: nominal currents in main mirror segment coils
    :param currents_nozzle: nominal currents in nozzle coils
    :param currents_end_cell: nominal currents in end cell segment coils (last coil is local nozzle coil)
    :param r_main: nominal coil radii in main mirror segment coils
    :param r_nozzle: nominal coil radii in nozzle coils
    :param r_end_cell: nominal coil radii in end cell segment coils (last coil is local nozzle coil)
    :param l_main: nozzle-to-nozzle length of main mirror segment
    :param l_end_cell: nozzle-to-nozzle length of each end cell segment
    :param n_end_cells: number of end cells on each end of the main mirror
    """

    def __init__(self, currents_main = tuple([1000.0] * 4), currents_nozzle = (4000.0, 4000.0), currents_end_cell = (250.0, 3500.0), r_main = tuple([0.12] * 4), r_nozzle = (0.04, 0.04), r_end_cell = (0.12, 0.05), l_main=0.3, l_end_cell = 0.075, n_end_cells=2):
        super().__init__()
        self.currents_main = currents_main
        self.currents_nozzle = currents_nozzle
        self.currents_end_cell = currents_end_cell
        self.r_main = r_main
        self.r_nozzle = r_nozzle
        self.r_end_cell = r_end_cell
        self.n_end_cells = n_end_cells
        self.l_main = l_main
        self.l_end_cell = l_end_cell

        # a little sanity checking
        if len(self.currents_main) != len(self.r_main):
            raise ValueError(f"main mirror current and radius arrays should be the same length: # of currents is {len(self.currents_main)} but # of radii is {len(self.r_main)}!")
        if len(self.currents_nozzle) != len(self.r_nozzle):
            raise ValueError(f"nozzle current and radius arrays should be the same length: # of currents is {len(self.currents_nozzle)} but # of radii is {len(self.r_nozzle)}!")
        if len(self.currents_end_cell) != len(self.r_end_cell):
            raise ValueError(f"end cell current and radius arrays should be the same length: # of currents is {len(self.currents_end_cell)} but # of radii is {len(self.r_end_cell)}!")

        # main mirror, oriented along z, centered around z=0
        # add nozzles first
        for i, (current, radius) in enumerate(zip(self.currents_nozzle, self.r_nozzle)):
            self.coils.append(CircularCoil(I= current, r0=radius, center=[0, 0, -self.l_main / 2 + i * self.l_main]))

        # add main coils
        for i, (current, radius) in enumerate(zip(self.currents_main, self.r_main)):
            self.coils.append(CircularCoil(I= current, r0=radius, center=[0, 0, -self.l_main / 2 + (i + 1) * self.l_main / (len(self.currents_main) + 1)]))

        # add compensator coils
        # self.coils.append(CircularCoil(I= 500, r0=0.13, center=[-0.0, 0, -self.l_main / 2 + (len(self.currents_main) +0.5) * self.l_main / (len(self.currents_main) + 1)], normal=[0,-0/180 * np.pi]))
        # self.coils.append(CircularCoil(I= 2000, r0=0.13, center=[0.02, 0, -self.l_main / 2 + 0.5 * self.l_main / (len(self.currents_main) + 1)], normal=[0,20/180 * np.pi]))

        # tandem end cells to make multi-mirror
        # each end cell consists of a main coil and a nozzle coil
        for i in range(n_end_cells):
            for j, (current, radius) in enumerate(zip(self.currents_end_cell, self.r_end_cell)):
                self.coils.append(CircularCoil(I= current, r0=radius, center=[0, 0, -self.l_main / 2 - i * self.l_end_cell - (j + 1) * self.l_end_cell / (len(self.currents_end_cell))]))
                self.coils.append(CircularCoil(I= current, r0=radius, center=[0, 0, self.l_main / 2 + i * self.l_end_cell + (j + 1) * self.l_end_cell / (len(self.currents_end_cell))]))

        # sort the coils along z
        new_coils = []
        zs = np.zeros(len(self.coils))
        for i, coil in enumerate(self.coils):
            zs[i] = coil.center[2]

        for i, new_idx in enumerate(np.argsort(zs)):
            new_coils.append(self.coils[new_idx])

        self.coils = new_coils

        # compute the main mirror ratio
        main_b = self.assemble().set_points(np.array([[0,0,0],[0,0,self.l_main/2]])).B()
        self.actual_mirror_ratio = main_b[1,2]/main_b[0,2]

        self.actual_end_cell_mirror_ratio = 0.0
        if n_end_cells > 0:
            end_b = self.assemble().set_points(np.array([[0,0,self.l_main / 2 + (len(self.currents_end_cell) - 1) * self.l_end_cell / len(self.currents_end_cell)],[0,0,self.l_main / 2 + self.l_end_cell]])).B()
            self.actual_end_cell_mirror_ratio = end_b[1,2]/end_b[0,2]

        self.compute_lineout_spline()

        print(self.actual_mirror_ratio, self.actual_end_cell_mirror_ratio)


    @property
    def length(self):
        """
        :return: The overall length of the GDMT
        """
        return self.l_main + self.n_end_cells * self.l_end_cell * 2


class ELMOCoils(CoilSet):
    """
    ELMO bumpy torus segment coil set class

    :param currents: nominal coil currents
    :param r: nominal coil radii
    :param circle_radius: radius of centerline of end coils (for matching to other elements)
    :param arc: angle in radians that the segment subtends
    :param spline_points: control points for a cubic spline for the centers of the coils (at least 4 points needed)
    :param ends: boolean indicating whether end coils should be included
    """

    def __init__(self, currents = tuple([2000.0] * 9), r=tuple([0.08] * 9), circle_radius = 0.3, arc=np.pi / 2, spline_points = None, ends=False):
        super().__init__()
        self.currents = currents
        self.r = r
        self.circle_radius = circle_radius
        self.arc = arc
        self.spline_points = spline_points

        # a little sanity checking
        if len(self.currents) != len(self.r):
            raise ValueError(f"current and radius arrays should be the same length: # of currents is {len(self.currents)} but # of radii is {len(self.r)}!")

        if self.spline_points is None:
            # construct a cubic spline through control points on a circle
            self.spline_points = np.array([
                [0.0, 0.0, self.circle_radius],
                [self.circle_radius * np.sin(self.arc/3), 0, self.circle_radius * np.cos(self.arc/3)],
                [self.circle_radius * np.sin(self.arc*2/3), 0, self.circle_radius * np.cos(self.arc*2/3)],
                [self.circle_radius * np.sin(self.arc), 0, self.circle_radius * np.cos(self.arc)]]).T

        spline, u = interpolate.make_splprep(self.spline_points, k=3)

        # interface coil test
        s = 0.5 / (len(self.currents) + 1)
        interface_coil = CircularCoil(I = 500, r0 = 0.16, center = spline(s)+ np.array([.019,0.0,0.015]), normal=[0.0, -4.3/180 * np.pi + np.pi/2])
        # self.coils.append(interface_coil)

        for i, (current, radius) in enumerate(zip(self.currents, self.r)):
            s = (i + 1) / (len(self.currents) + 1)
            dot = np.dot(spline.derivative()(s), (1.0, 0.0, 0.0)) / np.linalg.norm(spline.derivative()(s))
            # detect when we need to flip the sign of the angle with the cross product
            cross = np.cross(spline.derivative()(s), (1.0, 0.0, 0.0)) / np.linalg.norm(spline.derivative()(s))
            # todo: generalize this for things not in the zx plane
            angle = - np.arccos(dot) * np.sign(cross[1])
            self.coils.append(CircularCoil(I=current, r0=radius, center=spline(s), normal=[0.0, angle + np.pi/2]))

        self.compute_lineout_spline()


class PCLMSection:
    """
    Polygonal Closed-Loop Mirror coil set class

    :param currents_gdmt_main: nominal currents in main mirror segment coils
    :param currents_gdmt_nozzle: nominal currents in nozzle coils
    :param currents_gdmt_end_cell: nominal currents in end cell segment coils (last coil is local nozzle coil)
    :param currents_elmo_main: nominal currents in main elmo segment coils
    :param r_gdmt_main: nominal coil radii in main mirror segment coils
    :param r_gdmt_nozzle: nominal coil radii in nozzle coils
    :param r_gdmt_end_cell: nominal coil radii in end cell segment coils (last coil is local nozzle coil)
    :param r_gdmt_main: nominal coil radii in main elmo segment coils
    :param l_gdmt_main: nozzle-to-nozzle length of each main mirror segment
    :param l_gdmt_end_cell: nozzle-to-nozzle length of each end cell segment
    :param n_gdmt_end_cells: number of end cells on each end of each main mirror
    :param angle: angle of ELMO segment between two GDMT segments in degrees
    :param radius: radius of centerline of ELMO segment in meters
    :param n_gdmt_end_cells: number of end cells on either side of main mirror in GDMT segments
    :param angle: angle in radians that the elmo segment subtends
    :param radius: radius of centerline of end coils (for matching to other elements)
    """

    def __init__(self,
                 currents_gdmt_main = tuple([1000.0] * 4),
                 currents_gdmt_nozzle = (4000.0, 4000.0),
                 currents_gdmt_end_cell = (250.0, 3500.0),
                 currents_elmo_main = tuple([2000.0] * 9),
                 r_gdmt_main = tuple([0.12] * 4),
                 r_gdmt_nozzle = (0.04, 0.04),
                 r_gdmt_end_cell = (0.12, 0.05),
                 r_elmo_main = tuple([0.08] * 9),
                 l_gdmt_main=0.3, l_gdmt_end_cell = 0.075, n_gdmt_end_cells=2,
                 angle = 90, radius=0.5, spline_points = None, testing = False):

        self.angle = angle * np.pi/180
        self.radius = radius
        self.mirror_length = l_gdmt_main

        # form the GDMT and ELMO segments
        self.gdmt1 = GDMTCoils(currents_gdmt_main, currents_gdmt_nozzle, currents_gdmt_end_cell, r_gdmt_main, r_gdmt_nozzle, r_gdmt_end_cell, l_gdmt_main, l_gdmt_end_cell, n_gdmt_end_cells)
        self.gdmt2 = GDMTCoils(currents_gdmt_main, currents_gdmt_nozzle, currents_gdmt_end_cell, r_gdmt_main, r_gdmt_nozzle, r_gdmt_end_cell, l_gdmt_main, l_gdmt_end_cell, n_gdmt_end_cells)

        if angle == 0:
            self.elmo = None
            self.radius = 0
        else:
            self.elmo = ELMOCoils(currents_elmo_main, r_elmo_main, circle_radius = radius, arc = self.angle, spline_points = spline_points)
            self.elmo.rotate(-np.pi/2)

        # rotate and translate the various segments into place
        self.gdmt1_translate = np.array([-self.radius, 0, -self.gdmt1.length/2 * (testing == False)])
        self.gdmt1.translate(self.gdmt1_translate)

        self.gdmt2.rotate(self.angle)
        self.gdmt2_translate = np.array([-self.radius * np.cos(self.angle) + self.gdmt2.length/2 * np.sin(self.angle), 0, self.gdmt2.length/2 * np.cos(self.angle) + self.radius * np.sin(self.angle) + 0.01 * (angle == 0)])
        self.gdmt2.translate(self.gdmt2_translate)

        # create axial spline
        if angle == 0:
            self.all_coils = self.gdmt1.coils + self.gdmt2.coils
        else:
            self.all_coils = self.gdmt1.coils + self.elmo.coils + self.gdmt2.coils

        self.coil_centers=np.zeros((3,len(self.all_coils)))
        self.coil_normals=np.zeros((2,len(self.all_coils)))

        for i, coil in enumerate(self.all_coils):
            self.coil_centers[:,i]=coil.center
            self.coil_normals[:,i]=coil.normal

        self.spline, self.coils_s = interpolate.make_splprep(self.coil_centers, s=0)

        if angle == 0:
            if testing:
                self.field = self.gdmt1.assemble()
            else:
                self.field = self.gdmt1.assemble() + self.gdmt2.assemble()
        else:
            self.field = self.gdmt1.assemble() + self.elmo.assemble() + self.gdmt2.assemble()

    def get_lineout_center(self,n_points = 500):
        """
        Evaluates the center spline at the requested number of points
        :param n_points: number of points to interpolate the spline at (default 500)
        :return: interpolated spline, transposed to have shape (n_points, 3)
        """
        return self.spline(np.linspace(0,1,n_points)).T

    def get_spline_outer(self,n_points = 500, psi=0.0, r_frac = 0.5, radius = 0.0):
        """
        Evaluates an offset spline at the requested number of points
        :param radius: constant radius to use if nonzero
        :param r_frac: fraction of each coil radius to use if radius is zero
        :param psi: angle from x-axis to use
        :param n_points: number of points to interpolate the spline at (default 500)
        :return: interpolated spline, transposed to have shape (n_points, 3)
        """

        # todo: make this work in full 3D (rather than just planar zx)

        coil_outers = np.copy(self.coil_centers)
        for i, coil in enumerate(self.all_coils):
            r = radius
            if r == 0:
                r = r_frac * coil.r0
            # print(np.cos(self.coil_normals[1,i]),np.sin(self.coil_normals[1,i]))
            coil_outers[0,i] += np.sin(psi) * np.sin(self.coil_normals[1,i] + np.pi/2) * r
            coil_outers[1,i] += np.cos(psi) * r
            coil_outers[2,i] += np.sin(psi) * np.cos(self.coil_normals[1,i] + np.pi/2) * r

        outer_spline, coils_s_outers = interpolate.make_splprep(coil_outers, s=0, k=1)

        return outer_spline, coils_s_outers

        # result = outer_spline(np.linspace(0,1,n_points)).T
        #
        # return result

    def dump_coil_data(self, filename="pclm_coil_data.txt"):
        """
        Write coil data (centers, rotation, current, radius) to a text file
        Coils appear in the same order as the center interpolating spline
        units are SI
        :return:
        """
        coil_data = np.zeros((len(self.all_coils),7))

        for i,coil in enumerate(self.all_coils):
            coil_data[i,:] = [*coil.center, *coil.normal, coil.I, coil.r0]

        np.savetxt(filename, coil_data, fmt="% 10.7e", delimiter=", ", header="center_x       center_y       center_z         phi             theta           current         radius")

    def get_lineout(self, plot = True, show = False):
        points = self.get_lineout_center()

        distance = np.sqrt((np.diff(points,axis=0)**2).sum(axis=1)).cumsum()

        field = self.field
        lineout_field=field.set_points(points).B()

        if plot:
            plt.figure(figsize=(10,6))
            plt.plot(distance,np.sqrt((lineout_field**2).sum(axis=1))[:-1])
            plt.xlabel("distance along centerline (m)")
            plt.ylabel("B (T)")
            plt.title("Centerline lineout of magnetic field strength")
            plt.xlim([distance[0],distance[-1]])
            plt.grid(ls="--")
            plt.savefig("elmo-lineout.png", bbox_inches="tight", pad_inches=0.05)
            if show:
                plt.show()

        return distance, lineout_field

    def get_magnetic_field_deviation_along_axis(self, radius = 0.02, num_lineouts = 4, plot = True, show = False):
        s = np.linspace(0,1,500)
        points = self.spline(s).T

        tangents = self.spline.derivative()(s).T

        offaxis_points = np.zeros((num_lineouts, len(points), 3))
        r = radius

        offaxis_lineouts_b=[]

        # todo: vectorize this
        for i, psi in enumerate(np.linspace(0, 2*np.pi, num_lineouts, endpoint = False)):
            for j, tangent in enumerate(tangents):
                dot = np.dot(tangent, np.array([0.0,0.0,1.0])) / np.linalg.norm(tangent)
                # detect when we need to flip the sign of the angle with the cross product
                cross = np.cross(tangent, np.array([0.0,0.0,1.0])) / np.linalg.norm(tangent)
                # todo: generalize this for things not in the zx plane
                angle = - np.arccos(dot) * np.sign(cross[1])
                # print(angle+np.pi/2)
                # angle += np.pi/2
                # angle = -np.pi/2 + np.pi/2
                offaxis_points[i, j, :] = points[j,:]
                offaxis_points[i, j, 0] += r * np.sin(psi) * np.cos(angle)
                offaxis_points[i, j, 1] += r * np.cos(psi) #* r
                offaxis_points[i, j, 2] += r * np.sin(psi) * (-np.sin(angle))

            print(offaxis_points[i,:,2])

            self.field.set_points(offaxis_points[i,:,:])
            offaxis_lineouts_b.append(np.sqrt(((self.field.B())**2).sum(axis=1)))

        print(offaxis_lineouts_b)

        #
        #
        # offaxis_splines=[]
        # offaxis_splines_u=[]
        # offaxis_lineouts_b=[]
        #
        # # print(num_lineouts)
        # # r = radius
        # #
        # # for i, psi in enumerate(np.linspace(0, 2*np.pi, num_lineouts, endpoint = False)):
        # #     print(psi)
        # #     for j, point in enumerate(points):
        # #         offaxis_points[i, j, :] = point
        # #         offaxis_points[i, j, 0] += np.sin(psi) * np.sin(self.coil_normals[1,j] + np.pi/2) * r
        # #         offaxis_points[i, j, 1] += np.cos(psi) * r
        # #         offaxis_points[i, j, 2] += np.sin(psi) * np.cos(self.coil_normals[1,j] + np.pi/2) * r
        # #         print(self.coil_normals[1,j],np.sin(psi) * np.cos(self.coil_normals[1,j] + np.pi/2) * r)
        # #
        # # print(offaxis_points)
        # # exit()
        #
        # for i in range(num_lineouts):
        #     spline, u = self.get_spline_outer(psi=2*np.pi/num_lineouts*i, radius=radius)
        #     offaxis_splines.append(spline)
        #     offaxis_splines_u.append(u)
        #     self.field.set_points(offaxis_splines[i](offaxis_splines_u[i]).T)
        #     offaxis_lineouts_b.append(np.sqrt(((self.field.B())**2).sum(axis=1)))

        offaxis_lineouts_b=np.array(offaxis_lineouts_b)

        # print(offaxis_lineouts[0])
        # print(offaxis_lineouts_b[0])
        # print(offaxis_lineouts[1])
        # print(offaxis_lineouts_b[1])
        # offaxis_lineouts=np.array(offaxis_lineouts)

        average = offaxis_lineouts_b.mean(axis=0).squeeze()

        deviation = ((offaxis_lineouts_b - average[np.newaxis,:])/average[np.newaxis,:])

        print(deviation)

        if plot:
            plt.figure(figsize=(10,6))
            for i in range(deviation.shape[0]):
                plt.plot(deviation[i,:-1], label=rf"${360/num_lineouts*i:0.1f}^{{\circ}}$")
                # break
            plt.xlabel("distance along centerline (m)")
            plt.ylabel("fraction deviation in B")
            plt.title(f"Deviation from average in B (r = {radius} m")
            # plt.xlim([distance[0],distance[-1]])
            # plt.xlim([2.0,3.0])
            # plt.ylim([-.0001,.0001])
            # plt.xlim([0,22])
            plt.grid(ls="--")
            plt.legend(title="Angle")
            plt.savefig("elmo-lineout-b-deviation.png", bbox_inches="tight", pad_inches=0.05)
            if show:
                plt.show()
            plt.close()

        # return points, offaxis_lineouts, offaxis_lineouts_b
