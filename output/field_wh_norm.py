import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm
import time

def create_custom_colormap():
    cdict = {
        'red':   [[0.0,  16/255, 16/255],
                  [0.25, 109/255, 109/255],
                  [0.5,  1.0, 1.0],
                  [0.75, 220/255, 220/255],
                  [1.0,  109/255, 109/255]],
        'green': [[0.0,  70/255, 70/255],
                  [0.25, 173/255, 173/255],
                  [0.5,  1.0, 1.0],
                  [0.75, 109/255, 109/255],
                  [1.0,  1/255, 1/255]],
        'blue':  [[0.0,  128/255, 128/255],
                  [0.25, 209/255, 209/255],
                  [0.5,  1.0, 1.0],
                  [0.75, 87/255, 87/255],
                  [1.0,  31/255, 31/255]]
    }
    return LinearSegmentedColormap('CustomCmap', cdict)

def read_binary_data(file_path, header_size, nvar):
    with open(file_path, 'rb') as f:
        f.seek(header_size)
        data = np.fromfile(f, dtype=np.float32, count=-1)
        N = int(round((len(data)/nvar)**(1/3)))
        if len(data) != N**3 * nvar:
            raise ValueError(f"Expected {N**3 * nvar} elements but read {len(data)} elements.")
        reshaped_data = data.reshape((N, N, N, nvar), order='F')
    return reshaped_data, N

def create_pyvista_grid(x_coords, y_coords, z_coords, vor2, helicity):
    X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
    grid = pv.StructuredGrid(X, Y, Z)
    grid.point_data['vor'] = np.sqrt(vor2).flatten(order='F')
    grid.point_data['helicity'] = helicity.flatten(order='F')
    return grid

def smooth_contours(contours, total_iterations=30, n_steps=3):
    smoothed_contours = contours.copy()
    step_iterations = total_iterations // n_steps
    for _ in tqdm(range(n_steps), desc="Smoothing mesh"):
        smoothed_contours = smoothed_contours.smooth(n_iter=step_iterations, relaxation_factor=0.05)
    return smoothed_contours

def add_background_plane(plotter, position, normal, size, color='black'):
    plane = pv.Plane(center=position, direction=normal, i_size=size[0], j_size=size[1])
    plotter.add_mesh(plane, color=color, opacity=1)

def main():
    custom_cmap = create_custom_colormap()
    #file_path = 'F:/projects/DNS/oldcode/Steady/DNS_E1k2[p0=-0.666]_Re360_Steady_HITnew_N256/output/w_h006250.dat'
    folder_name = "./"
    file_path = folder_name + 'field_wh.dat'
    header_size = 0 #48*4
    nvar = 2

    reshaped_data, N = read_binary_data(file_path, header_size, nvar)
    vor = reshaped_data[:, :, :, 0]
    vor2 = vor**2
    helicity = reshaped_data[:, :, :, 1]

    #
    rms_vor = np.sqrt(np.mean(vor2))
    rms_helicity = np.sqrt(np.mean(np.square(helicity)))
    iso_w = rms_vor*3/np.sqrt(2)
    helicity_max = 3*rms_helicity


    x_coords = np.linspace(0, 2*np.pi, N)
    y_coords = np.linspace(0, 2*np.pi, N)
    z_coords = np.linspace(0, 2*np.pi, N)
    grid = create_pyvista_grid(x_coords, y_coords, z_coords, vor2, helicity)

    plotter = pv.Plotter(window_size=[2048, 2048], notebook=False, off_screen=True)
    contours = grid.contour(isosurfaces=[iso_w], scalars='vor')

    start_time = time.time()
    contours = smooth_contours(contours)
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.2f} seconds")

    bounds = grid.bounds
    plotter.add_mesh(contours, scalars='helicity', cmap=custom_cmap, clim=[-helicity_max, helicity_max], opacity=1.0)
    plotter.add_mesh(pv.Box(bounds=bounds), color='black', style='wireframe', line_width=4)
    plotter.remove_scalar_bar()

    background_planes = [
        {'position': (0, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2), 'normal': (1, 0, 0), 'size': ((bounds[1] - bounds[0]), (bounds[3] - bounds[2]), (bounds[5] - bounds[4]))},
        {'position': ((bounds[0] + bounds[1]) / 2, 0, (bounds[4] + bounds[5]) / 2), 'normal': (0, -1, 0), 'size': ((bounds[1] - bounds[0]), (bounds[3] - bounds[2]), (bounds[5] - bounds[4]))},
        {'position': ((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, 0), 'normal': (0, 0, -1), 'size': ((bounds[1] - bounds[0]), (bounds[3] - bounds[2]), 0)}
    ]
    for plane in background_planes:
        add_background_plane(plotter, **plane)

    plotter.enable_anti_aliasing()
    light = pv.Light()
    light.intensity = 0.5
    plotter.add_light(light)

    plotter.view_isometric()
    plotter.camera.azimuth = -21
    plotter.camera.elevation = -15
    plotter.camera.view_angle = 23.1
    shift = [0.76, 0, -0.35]
    new_position = [plotter.camera.position[i] + shift[i] for i in range(3)]
    new_focal_point = [plotter.camera.focal_point[i] + shift[i] for i in range(3)]
    plotter.camera.position = new_position
    plotter.camera.focal_point = new_focal_point

    # figure_name = f"w{round(iso_w)}_h{round(helicity_max)}.png"
    figure_name = "field_wh_norm.png"
    output_image_path = folder_name + figure_name
    plotter.screenshot(output_image_path)
    plotter.show()

if __name__ == "__main__":
    main()