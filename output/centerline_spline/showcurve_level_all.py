import vtk
import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
import time
from scipy.spatial.transform import Rotation as R


def create_custom_colormap_azure():
    """创建自定义颜色映射：黑色 -> 蓝色 -> 天蓝 -> 白色"""
    cdict = {
        'red': [
            [0.0, 0 / 255, 0 / 255],  # 黑色
            [0.25, 0 / 255, 0 / 255],  # 蓝色
            [0.5, 135 / 255, 135 / 255],  # 天蓝
            [1.0, 255 / 255, 255 / 255]  # 白色
        ],
        'green': [
            [0.0, 0 / 255, 0 / 255],  # 黑色
            [0.25, 0 / 255, 0 / 255],  # 蓝色
            [0.5, 206 / 255, 206 / 255],  # 天蓝
            [1.0, 255 / 255, 255 / 255]  # 白色
        ],
        'blue': [
            [0.0, 0 / 255, 0 / 255],  # 黑色
            [0.25, 255 / 255, 255 / 255],  # 蓝色
            [0.5, 250 / 255, 250 / 255],  # 天蓝
            [1.0, 255 / 255, 255 / 255]  # 白色
        ]
    }
    return LinearSegmentedColormap('CustomCmap', cdict)


def create_custom_colormap_pure():
    """创建自定义颜色映射：统一颜色"""
    RR = 0 / 255
    GG = 191 / 255
    BB = 255 / 255
    cdict = {
        'red': [
            [0.0, RR, RR],
            [1.0, RR, RR]
        ],
        'green': [
            [0.0, GG, GG],
            [1.0, GG, GG]
        ],
        'blue': [
            [0.0, BB, BB],
            [1.0, BB, BB]
        ]
    }
    return LinearSegmentedColormap('CustomCmap', cdict)


def create_oriented_disc(center, normal, radius, resolution=30):
    """创建具有特定方向的圆盘，使用Plane并旋转以兼容旧版PyVista"""
    plane = pv.Plane(center=(0, 0, 0), direction=(0, 0, 1), i_size=2 * radius, j_size=2 * radius,
                     i_resolution=resolution, j_resolution=resolution)

    default_normal = np.array([0, 0, 1])
    target_normal = normal / np.linalg.norm(normal)

    cross = np.cross(default_normal, target_normal)
    if np.linalg.norm(cross) > 1e-6:
        axis = cross / np.linalg.norm(cross)
        angle = np.arccos(np.dot(default_normal, target_normal)) * 180 / np.pi
        rot = R.from_rotvec(axis * angle * np.pi / 180)
        plane.points = rot.apply(plane.points)

    plane.points += center
    distances = np.linalg.norm(plane.points - center, axis=1)
    plane = plane.extract_points(distances <= radius)
    return plane


def add_yellow_discs(plotter, points, bounds=None, radius=0.08, resolution=360, thickness=0.01):
    """在每个点上添加带厚度的黄色圆盘，法向量为指向下一个点的方向"""
    if bounds is None or len(points) < 2:
        return

    for i in range(len(points) - 1):
        current_point = points[i]
        next_point = points[i + 1]

        if not is_inside_bounds(current_point, bounds):
            continue

        direction = next_point - current_point
        if np.linalg.norm(direction) < 1e-6:
            continue

        normal = direction / np.linalg.norm(direction)
        cylinder = pv.Cylinder(center=(0, 0, 0), direction=(0, 0, 1), radius=radius, height=thickness,
                               resolution=resolution)
        z_axis = np.array([0, 0, 1])
        if not np.allclose(normal, z_axis) and not np.allclose(normal, -z_axis):
            rot_axis = np.cross(z_axis, normal)
            rot_axis = rot_axis / np.linalg.norm(rot_axis)
            rot_angle = np.degrees(np.arccos(np.dot(z_axis, normal)))
            cylinder.rotate_vector(rot_axis, angle=rot_angle, inplace=True)
        cylinder.translate(current_point, inplace=True)
        plotter.add_mesh(cylinder, color='yellow', opacity=0.9)

def add_yellow_discs_T(plotter, points, vectors, bounds=None, radius=0.08, resolution=360, thickness=0.01):
    """在每个点上添加带厚度的黄色圆盘，法向量为指向下一个点的方向"""
    if bounds is None or len(points) < 2:
        return

    for i in range(len(points) - 1):
        current_point = points[i]
        next_point = points[i + 1]

        if not is_inside_bounds(current_point, bounds):
            continue

        direction = vectors[i]
        if np.linalg.norm(direction) < 1e-6:
            continue

        normal = direction / np.linalg.norm(direction)
        cylinder = pv.Cylinder(center=(0, 0, 0), direction=(0, 0, 1), radius=radius, height=thickness,
                               resolution=resolution)
        z_axis = np.array([0, 0, 1])
        if not np.allclose(normal, z_axis) and not np.allclose(normal, -z_axis):
            rot_axis = np.cross(z_axis, normal)
            rot_axis = rot_axis / np.linalg.norm(rot_axis)
            rot_angle = np.degrees(np.arccos(np.dot(z_axis, normal)))
            cylinder.rotate_vector(rot_axis, angle=rot_angle, inplace=True)
        cylinder.translate(current_point, inplace=True)
        plotter.add_mesh(cylinder, color='yellow', opacity=0.9)

def add_cylinder_around_curve(plotter, points, radius=0.1, num_points=20, color='cyan', opacity=0.3, bounds=None,
                              extend_length=0.2):
    """
    对每一对相邻点组成的小段曲线，单独生成带前后拓展的圆柱
    :param plotter: pyvista plotter对象
    :param points: 曲线点集 (N, 3)
    :param radius: 圆柱半径
    :param num_points: 每个小段的插值点数
    :param color: 圆柱颜色
    :param opacity: 透明度
    :param bounds: 边界范围
    :param extend_length: 每个小段首尾拓展的长度（相对小段自身长度的比例）
    """
    if len(points) < 2 or bounds is None:
        return

    # 遍历每一对相邻点，作为独立小段
    for i in range(len(points) - 1):
        p1 = points[i]
        p2 = points[i + 1]
        segment_vec = p2 - p1
        segment_len = np.linalg.norm(segment_vec)
        if segment_len < 1e-6:
            continue

        # 小段的单位方向向量
        dir_vec = segment_vec / segment_len

        # 计算小段的拓展点：起点向前拓展，终点向后拓展
        extend_dist = extend_length * segment_len
        p1_extend = p1 - dir_vec * extend_dist  # 小段起点向前拓展
        p2_extend = p2 + dir_vec * extend_dist  # 小段终点向后拓展

        # 构建拓展后的小段点集
        extended_segment = np.vstack([p1_extend, p1, p2, p2_extend])

        # 对拓展后的小段进行插值，生成平滑曲线
        t = np.linspace(0, 1, len(extended_segment))
        t_fine = np.linspace(0, 1, num_points)
        x = np.interp(t_fine, t, extended_segment[:, 0])
        y = np.interp(t_fine, t, extended_segment[:, 1])
        z = np.interp(t_fine, t, extended_segment[:, 2])
        fine_points = np.column_stack((x, y, z))

        # 创建小段的样条曲线
        spline = pv.PolyData(fine_points)
        lines = np.hstack([np.ones((num_points - 1, 1), dtype=int) * 2,
                           np.arange(num_points - 1)[:, np.newaxis],
                           np.arange(1, num_points)[:, np.newaxis]])
        spline.lines = lines.ravel()

        # 生成小段的圆柱
        cylinder = spline.tube(radius=radius, capping=False)

        # 裁剪边界外的部分
        cyl_points = cylinder.points
        inside = np.array([is_inside_bounds(point, bounds) for point in cyl_points])
        if np.any(inside):
            cylinder = cylinder.extract_points(inside)
            plotter.add_mesh(cylinder, color=color, opacity=opacity)


def read_curve_data(filename):
    """读取曲线数据文件"""
    with open(filename, 'r') as file:
        Nv = int(file.readline().strip())
        Nps = [int(file.readline().strip()) for _ in range(Nv)]

    Bread = np.loadtxt(filename, skiprows=1 + Nv)
    return Nv, Nps, Bread


def calculate_curvature(BB_p):
    """计算曲率（向量化优化）"""
    Np = BB_p.shape[1]
    curvature = np.zeros(Np)

    r1 = BB_p[:, :-2]
    r2 = BB_p[:, 1:-1]
    r3 = BB_p[:, 2:]

    dr1 = r2 - r1
    dr2 = r3 - r2
    T1 = dr1 / np.linalg.norm(dr1, axis=0)
    T2 = dr2 / np.linalg.norm(dr2, axis=0)
    dT = T2 - T1
    ds = (np.linalg.norm(dr1, axis=0) + np.linalg.norm(dr2, axis=0)) / 2
    curvature[1:-1] = np.linalg.norm(dT, axis=0) / ds

    curvature[0] = curvature[1]
    curvature[-1] = curvature[-2]
    return curvature


def is_inside_bounds(point, bounds):
    """检查点是否严格在边界内（不包括边界）"""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    return (x_min < point[0] < x_max and
            y_min < point[1] < y_max and
            z_min < point[2] < z_max)


def create_pyvista_lines(BB_p, log_curvature, cmap, bounds):
    """创建 pyvista 线段，仅保留完全在边界内的线段"""
    points = []
    lines = []
    colors = []

    for j in range(BB_p.shape[1] - 1):
        point1 = BB_p[:, j]
        point2 = BB_p[:, j + 1]

        if is_inside_bounds(point1, bounds) and is_inside_bounds(point2, bounds):
            dist = np.linalg.norm(point1 - point2)
            if dist > np.pi:
                continue

            avg_log_curvature = (log_curvature[j] + log_curvature[j + 1]) / 2
            color = cmap((avg_log_curvature - np.min(log_curvature)) / (np.max(log_curvature) - np.min(log_curvature)))
            points.append(point1)
            points.append(point2)
            lines.append([2, len(points) - 2, len(points) - 1])
            colors.append(color[:3])

    poly_data = pv.PolyData()
    if points:
        poly_data.points = np.array(points)
        poly_data.lines = np.array(lines)
        poly_data["colors"] = np.array(colors)
    return poly_data


def add_spheres(color, plotter, points, radius=0.05, bounds=None):
    """在每个严格位于边界内的点位置添加球体"""
    if bounds is None:
        return

    unique_points = np.unique(points, axis=0)
    for point in unique_points:
        if is_inside_bounds(point, bounds):
            sphere = pv.Sphere(radius=radius, center=point)
            plotter.add_mesh(sphere, color=color, opacity=1.0)


def add_cubes_between_points(plotter, points, side_length=0.16, opacity=0.5, bounds=None, color='yellow'):
    """在每两个相邻且完全位于边界内的点之间添加长方体"""
    if bounds is None:
        return

    expand = side_length / 2
    for i in range(len(points) - 1):
        current_point = points[i]
        next_point = points[i + 1]

        if not (is_inside_bounds(current_point, bounds) and is_inside_bounds(next_point, bounds)):
            continue

        x_min, x_max = min(current_point[0], next_point[0]), max(current_point[0], next_point[0])
        y_min, y_max = min(current_point[1], next_point[1]), max(current_point[1], next_point[1])
        z_min, z_max = min(current_point[2], next_point[2]), max(current_point[2], next_point[2])

        x_min -= expand
        x_max += expand
        y_min -= expand
        y_max += expand
        z_min -= expand
        z_max += expand

        center = [(x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2]
        x_length = x_max - x_min
        y_length = y_max - y_min
        z_length = z_max - z_min

        cube = pv.Cube(center=center, x_length=x_length, y_length=y_length, z_length=z_length)
        plotter.add_mesh(cube, color=color, opacity=opacity)
        plotter.add_mesh(
            cube.extract_all_edges(),
            color='black',
            opacity=1.0,
            line_width=6,
            render_lines_as_tubes=True
        )


def add_bounding_box(plotter, bounds, color='black', line_width=4):
    """添加黑色边界框"""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    edges = [
        np.array([[x_min, y_min, z_min], [x_max, y_min, z_min]]),
        np.array([[x_min, y_max, z_min], [x_max, y_max, z_min]]),
        np.array([[x_min, y_min, z_min], [x_min, y_max, z_min]]),
        np.array([[x_max, y_min, z_min], [x_max, y_max, z_min]]),
        np.array([[x_min, y_min, z_max], [x_max, y_min, z_max]]),
        np.array([[x_min, y_max, z_max], [x_max, y_max, z_max]]),
        np.array([[x_min, y_min, z_max], [x_min, y_max, z_max]]),
        np.array([[x_max, y_min, z_max], [x_max, y_max, z_max]]),
        np.array([[x_min, y_min, z_min], [x_min, y_min, z_max]]),
        np.array([[x_max, y_min, z_min], [x_max, y_min, z_max]]),
        np.array([[x_min, y_max, z_min], [x_min, y_max, z_max]]),
        np.array([[x_max, y_max, z_min], [x_max, y_max, z_max]])
    ]

    for edge in edges:
        plotter.add_mesh(pv.Line(edge[0], edge[1]), color=color, line_width=line_width)


def add_bounding_box_p(plotter, bounds, color='black', line_width=4):
    """添加带偏移的黑色边界框"""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    shift = 0.2
    x_min -= shift
    y_min -= shift
    z_min -= shift
    x_max += shift
    y_max += shift
    z_max += shift
    edges = [
        np.array([[x_min, y_min, z_min], [x_max, y_min, z_min]]),
        np.array([[x_min, y_max, z_min], [x_max, y_max, z_min]]),
        np.array([[x_min, y_min, z_min], [x_min, y_max, z_min]]),
        np.array([[x_max, y_min, z_min], [x_max, y_max, z_min]]),
        np.array([[x_min, y_min, z_max], [x_max, y_min, z_max]]),
        np.array([[x_min, y_max, z_max], [x_max, y_max, z_max]]),
        np.array([[x_min, y_min, z_max], [x_min, y_max, z_max]]),
        np.array([[x_max, y_min, z_max], [x_max, y_max, z_max]]),
        np.array([[x_min, y_min, z_min], [x_min, y_min, z_max]]),
        np.array([[x_max, y_min, z_min], [x_max, y_min, z_max]]),
        np.array([[x_min, y_max, z_min], [x_min, y_max, z_max]]),
        np.array([[x_max, y_max, z_min], [x_max, y_max, z_max]])
    ]

    for edge in edges:
        plotter.add_mesh(pv.Line(edge[0], edge[1]), color=color, line_width=line_width)


def add_black_background(plotter, bounds):
    """为后面的三个面添加背景"""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    shift = 0.2
    x_min -= shift
    y_min -= shift
    z_min -= shift
    x_max += shift
    y_max += shift
    z_max += shift

    plane1 = pv.Plane(center=[x_min, (y_min + y_max) / 2, (z_min + z_max) / 2],
                      direction=[1, 0, 0],
                      i_size=(z_max - z_min),
                      j_size=(y_max - y_min))
    plane2 = pv.Plane(center=[(x_min + x_max) / 2, y_max, (z_min + z_max) / 2],
                      direction=[0, 1, 0],
                      i_size=(z_max - z_min),
                      j_size=(x_max - x_min))
    plane3 = pv.Plane(center=[(x_min + x_max) / 2, (y_min + y_max) / 2, z_min],
                      direction=[0, 0, 1],
                      i_size=(x_max - x_min),
                      j_size=(y_max - y_min))

    color = '#f0f0f0'
    opa = 1.0
    plotter.add_mesh(plane1, color=color, opacity=opa)
    plotter.add_mesh(plane2, color=color, opacity=opa)
    plotter.add_mesh(plane3, color=color, opacity=opa)

def add_arrow_for_N_vector(plotter, points, N_vectors, bounds=None, scale=0.1, color='green', opacity=1.0):
    """在每个点上添加表示N向量的箭头"""
    if bounds is None or len(points) != len(N_vectors):
        return

    for i in range(len(points)):
        current_point = points[i]
        N_vector = N_vectors[i]

        if not is_inside_bounds(current_point, bounds):
            continue

        # 创建箭头
        arrow = pv.Arrow(start=current_point, direction=N_vector, scale=scale, tip_length=0.1, tip_radius=0.05, shaft_radius=0.02)
        plotter.add_mesh(arrow, color=color, opacity=opacity)

def main():
    custom_cmap = create_custom_colormap_pure()
    file_prefix = 'centerline_spline'
    file_start = 1
    file_end = 1

    size = 1024
    plotter = pv.Plotter(window_size=[size, size], notebook=False, off_screen=False)
    all_points = []
    all_curvatures = []

    x_s = -np.pi
    x_e = np.pi
    length = x_e - x_s
    bounds = [x_s, x_e, x_s, x_e, x_s, x_e]

    for i_f in range(file_start, file_end + 1):
        filename = 'centerline_spline_00001.dat'
        print(f"Processing file: {filename}")
        Nv, Nps, Bread = read_curve_data(filename)
        filename = 'centerline_spline_T_00001.dat'
        Nv, Nps, T_list = read_curve_data(filename)
        # 读取N向量数据
        filename = 'centerline_spline_N_00001.dat'
        print(f"Processing file: {filename}")
        Nv, Nps, N_list = read_curve_data(filename)
        all_points.append(Bread[:, :3])
        all_points.append(T_list[:, :3])
        all_points.append(N_list[:, :3])

        for i in range(Nv):
            Np = Nps[i]
            BB_p = np.zeros((3, Np+1))
            TT_p = np.zeros((3, Np+1))
            NN_p = np.zeros((3, Np+1))
            head = sum(Nps[:i]) if i > 0 else 0

            for j in range(Np+1):
                BB_p[:, j] = Bread[head + j%Np, :3]
                TT_p[:, j] = T_list[head + j%Np, :3]
                NN_p[:, j] = N_list[head + j%Np, :3]  # 读取N向量数据

            current_points = BB_p.T
            T_points = TT_p.T
            N_points = NN_p.T  # N向量数据
            add_spheres('#b2182b', plotter, current_points, radius=0.04, bounds=bounds)
            #add_yellow_discs(plotter, current_points, bounds=bounds, radius=0.5)
            #add_yellow_discs_T(plotter, current_points, T_points, bounds=bounds, radius=0.36)
            # add_cubes_between_points(plotter, current_points, side_length=0.48, opacity=0.1, bounds=bounds,
            #                          color='#de77ae')

            # 调用修改后的函数，设置每个小段的拓展长度比例
            # add_cylinder_around_curve(plotter, current_points, radius=0.24, num_points=20,
            #                           color='cyan', opacity=0.2, bounds=bounds, extend_length=0.24)

            # 添加N向量箭头
            add_arrow_for_N_vector(plotter, current_points, N_points, bounds=bounds, scale=0.3, color='#fdae6b')

            curvature = calculate_curvature(BB_p)
            log_curvature = np.log(curvature + np.finfo(float).eps)
            all_curvatures.append(log_curvature)
            poly_data = create_pyvista_lines(BB_p, log_curvature, custom_cmap, bounds)
            plotter.add_mesh(poly_data, scalars="colors", rgb=True, line_width=40)

        filename = '../../centerline_input.dat'
        print(f"Processing file: {filename}")
        Nv, Nps, Bread = read_curve_data(filename)
        all_points.append(Bread[:, :3])

        for i in range(Nv):
            Np = Nps[i]
            BB_p = np.zeros((3, Np))
            head = sum(Nps[:i]) if i > 0 else 0

            for j in range(Np):
                BB_p[:, j] = Bread[head + j, :3]

            current_points = BB_p.T
            add_spheres('#4155FF', plotter, current_points, radius=0.08, bounds=bounds)

    if all_curvatures:
        all_curvatures = np.concatenate(all_curvatures)
        min_curv = np.min(all_curvatures)
        max_curv = np.max(all_curvatures)
        dummy = pv.PolyData()
        dummy.points = np.array([[0, 0, 0]])
        dummy["scalars"] = [0]

    plotter.enable_anti_aliasing()
    light = pv.Light()
    light.intensity = 0.5
    light.position = (5, 0, 0)
    plotter.add_light(light)

    plotter.view_isometric()
    plotter.camera.azimuth = -45
    plotter.camera.elevation = -35
    #plotter.camera.view_angle = 30
    plotter.camera.parallel_scale = 2.0
    ratio = length / (2 * np.pi)
    shift = [0.76 * ratio, 0, -0.0 * ratio]
    new_position = [plotter.camera.position[i] + shift[i] for i in range(3)]
    new_focal_point = [plotter.camera.focal_point[i] + shift[i] for i in range(3)]
    plotter.camera.position = new_position
    plotter.camera.focal_point = new_focal_point
    plotter.camera.SetParallelProjection(True)

    output_image_path = 'show_level_all.jpeg'
    # plotter.screenshot(output_image_path)
    plotter.show()


if __name__ == "__main__":
    main()