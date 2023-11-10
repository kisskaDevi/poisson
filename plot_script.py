from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
import argparse


class point:
    def __init__(self, x: float, y: float) -> None:
        self.x = x
        self.y = y

    def __init__(self, in_pair: str) -> None:
        pair = in_pair.replace('(','').replace(')','').split(',')
        self.x = float(pair[0])
        self.y = float(pair[1])


class grid(list[list[point]]): 
    def __init__(self, filename: str) -> None:
        super().__init__([point(pair) for pair in line.split('\t') if pair != '\n'] for line in open(filename, "r"))
        self.extent = (min(min([pair.x for pair in line]) for line in self), 
                       max(max([pair.x for pair in line]) for line in self),
                       min(min([pair.y for pair in line]) for line in self), 
                       max(max([pair.y for pair in line]) for line in self))

    def get_np_x(self) -> np.array:
        return np.array([[pair.x for pair in line] for line in self])
    
    def get_np_y(self) -> np.array:
        return np.array([[pair.y for pair in line] for line in self])


class function:
    def __init__(self, u: float, gx: float, gy: float) -> None:
        self.u = u
        self.gx = gx
        self.gy = gy
        self.gabs = sqrt(self.gx * self.gx + self.gy * self.gy)

    def __init__(self, in_function: str) -> None:
        func = in_function.replace('(','').replace(')','').split(',')
        self.u = float(func[0])
        self.gx = float(func[1])
        self.gy = float(func[2])
        self.gabs = sqrt(self.gx * self.gx + self.gy * self.gy)


class field(list[list[function]]):
    def __init__(self, filename: str) -> None:
        super().__init__([function(func) for func in line.split('\t') if func != '\n'] for line in open(filename, "r"))
    
    def get_np_u(self) -> np.array:
        return np.array([[f.u for f in line] for line in self])

    def get_np_gx(self) -> np.array:
        return np.array([[f.gx for f in line] for line in self])
    
    def get_np_gy(self) -> np.array:
        return np.array([[f.gy for f in line] for line in self])
    
    def get_np_gabs(self) -> np.array:
        return np.array([[f.gabs for f in line] for line in self])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--path', type=str, help='path to data', default="./build/Release/res")
    parser.add_argument('--type', type=str, help='plot type')
    args = parser.parse_args()
    
    xy_grid = grid(args.path + "/xy_grid.txt")
    u_field = field(args.path + "/u_field.txt")

    # exact = np.array([[float(val) for val in line.split('\t') if val != '\n'] for line in open(args.path + "/ex_func.txt", "r")])
    value = u_field.get_np_u()

    fig = plt.figure(figsize = (12,10))
    match args.type:
        case "2d":
            plt.xlabel('x')
            plt.ylabel('y')
            plt.imshow(np.transpose(value), interpolation='bilinear', cmap=plt.cm.jet, origin='lower', extent=xy_grid.extent, vmin=np.min(value), vmax=np.max(value)) 
            plt.colorbar()
            plt.quiver(xy_grid.get_np_x(), xy_grid.get_np_y(), - u_field.get_np_gx(), - u_field.get_np_gy(), units='width')
        case "3d":
            ax = plt.axes(projection='3d')
            surf = ax.plot_surface(xy_grid.get_np_x(), xy_grid.get_np_y(), value, cmap = plt.cm.cividis)
    plt.show()