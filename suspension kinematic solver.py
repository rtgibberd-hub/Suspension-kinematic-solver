import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


# --- Geometry dictionary (placeholder values) ---
geometry = {
    "lower_wishbone": {
        "pivot_1": [500.928, -166.143, -132.921],
        "pivot_2": [750.830, -164.172, -132.921],
        "outboard": [645.103, -555.018, -125.19]
    },
    "upper_wishbone": {
        "pivot_1": [510.756, -261, 25.28],
        "pivot_2": [755.668, -261,25.28],
        "outboard": [645.103, -520.624, 53.465]
    },
    "damper": {
        "chassis_mount": [645, -265, 260],
        "wishbone_mount": [644.081, -511.323, -91.76]
    },
    "tie_rod": {
        "rack_mount": [579.443, -273.296, -66.932],
        "knuckle_link": [574.947, -544.567, -67.189]
    },
    "knuckle": {
        "lower_wishbone_connection": [644.481, -555.018, -125.19],
        "upper_wishbone_connection": [645.103, -520.624, 53.465],
        "tie_rod_connection":  [574.947, -544.567, -67.189]
    },
    "wheel_axle": {
        "axle_knuckle": [645.103, -511.323, -29.721],
        "axle_wheel": [645.103, -667, -29.721]
    },
    "wheel": {
        "center": [645.103, -649.736, -29.721], 
        "top": [645.103, -649.736, 173.479], 
        "bottom": [645.103, -649.736, -232.921]
        
    }, 
    "contact_patch":{
        "front_inboard": [0, 0, 0], 
        "front_outbaord": [], 
        "back_inbaord":[],
        "back_outbaord": [] 
       
    }, 
    "rocker": {
        "damper_mount":[644.531, -359.541, 123.250], 
        "chasis_mount": [645.103, -265, 69], 
        "ARB_mount": [644.531, -312.945, 79.888], 
        "pushrod_mount":[644.531, -360.438, 91.263]
    }
}

# --- Constants ---
boundary = 5
initial_guess = 0
max_iterations = 1000

# --- Extract geometry vectors ---
axle_knuckle = geometry["wheel_axle"]["axle_knuckle"]
axle_wheel = geometry["wheel_axle"]["axle_wheel"]
wheel_top  = geometry ["wheel"]["top"]
wheel_bottom = geometry ["wheel"]["bottom"]
wheel_center = geometry["wheel"]["center"]

damp_lwb_mt = geometry["damper"]["wishbone_mount"]
damp_chassis_mt = geometry["damper"]["chassis_mount"]
damper_mount = geometry["rocker"]["pushrod_mount"]
P1 = geometry["lower_wishbone"]["pivot_1"]
P2 = geometry["lower_wishbone"]["pivot_2"]
P3 = geometry["upper_wishbone"]["pivot_1"]
P4 = geometry["upper_wishbone"]["pivot_2"]

lwb_outboard_coord = geometry["lower_wishbone"]["outboard"]
uwb_outboard_coord = geometry["upper_wishbone"]["outboard"]

lwb_knuckle = geometry["knuckle"]["lower_wishbone_connection"]
uwb_knuckle = geometry["knuckle"]["upper_wishbone_connection"]
knuckle_v = np.array(uwb_knuckle) - np.array(lwb_knuckle)

lw_wishbone_v = np.array(lwb_outboard_coord) - np.array(P1)
uw_wishbone_v = np.array(uwb_outboard_coord) - np.array(P3)

uwb_rot_axis = np.array(P3) - np.array(P4)
uwb_pivot = P3

lwb_rot_axis = np.array(P2) - np.array(P1)
lwb_pivot = P1

outboard_damp_mt = np.array(damp_lwb_mt) - np.array(P1)
axle_v = np.array(axle_wheel) - np.array(axle_knuckle)
len_axle_uwb = np.array(uwb_outboard_coord) - np.array(axle_knuckle)
wheel_v = np.array(wheel_top) - np.array(wheel_bottom)



# --- Main simulation ---
def main():

    # --- Lengths ---
    global damper_length
    damper_length = length_calc(damp_chassis_mt, lwb_outboard_coord)
    global knuckle_length
    knuckle_length = length_calc(lwb_knuckle, uwb_knuckle)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    print (f"original damper length is: {damper_length}")

    changed_damper_length = input("What is the length of the damper you want to simulate?\n")
    clean_input_damp_len = float(changed_damper_length.strip())

    new_lwb_outboard = move_damp_mt(clean_input_damp_len, lwb_outboard_coord, P1, P2)

    new_uwb_outboard = move_damp_mt(clean_input_damp_len, uwb_outboard_coord, P3, P4)


    new_damper_mount = move_damp_mt(clean_input_damp_len, damp_lwb_mt, P1, P2)

    angle_change_damper = angle_finder(np.array(damp_lwb_mt) - np.array(P1), np.array(new_damper_mount) - np.array(P1))

    upd_knuckle_v = np.array(new_uwb_outboard) - np.array(new_lwb_outboard)
    delta_kingpin_angle = angle_finder(upd_knuckle_v, knuckle_v)
    print(f"Change in camber is: {delta_kingpin_angle * (180/np.pi)} degrees")

    knuckle_bottom_to_axle = length_calc(lwb_knuckle, axle_knuckle)
    dir_upd_knuckle_v = upd_knuckle_v / np.linalg.norm(upd_knuckle_v)
    new_knuckle_axle = new_lwb_outboard + dir_upd_knuckle_v * knuckle_bottom_to_axle


    axle_change_v = new_knuckle_axle - np.array(axle_knuckle)
    knuckle_dir_old = np.array(uwb_knuckle) - np.array(lwb_knuckle)
    knuckle_dir_new = np.array(new_uwb_outboard) - np.array(new_lwb_outboard)
    rotation_axis = np.cross(knuckle_dir_old, knuckle_dir_new)    
    rotation_axis /= np.linalg.norm(rotation_axis)

    theta = angle_finder(knuckle_dir_old, knuckle_dir_new)

    axle_wheel_local = np.array(axle_wheel) - np.array(axle_knuckle)
    rotated_axle_wheel_local = rodriguez_rotation(axle_wheel_local, rotation_axis, theta)
    upd_axle_wheel = new_knuckle_axle + rotated_axle_wheel_local
    upd_axle_v = upd_axle_wheel - new_knuckle_axle
    upd_unit_wheel_v_perp = upd_axle_v / np.linalg.norm(upd_axle_v)
    wheel_width =  length_calc(wheel_center, axle_wheel)
    upd_wheel_centre = upd_axle_wheel + upd_unit_wheel_v_perp * wheel_width

    wheel_top_local = np.array(wheel_top) - np.array(axle_wheel)
    wheel_bottom_local = np.array(wheel_bottom) - np.array(axle_wheel)

    rotated_top_local = rodriguez_rotation(wheel_top_local, rotation_axis, theta)
    rotated_bottom_local = rodriguez_rotation(wheel_bottom_local, rotation_axis, theta)

    upd_wheel_top = upd_axle_wheel + rotated_top_local
    upd_wheel_bottom = upd_axle_wheel + rotated_bottom_local
    plot_points_and_dampers(new_lwb_outboard, ax)
    plot_initial_vectors(ax)
    plot_final_vectors(ax, new_lwb_outboard, new_uwb_outboard, new_damper_mount, new_knuckle_axle, upd_axle_wheel, upd_wheel_centre)
    plt.show()






    return new_lwb_outboard
   
def simulate_geometry(changed_damper_length):
    new_lwb_outboard = move_damp_mt(changed_damper_length, lwb_outboard_coord, P1, P2)
    new_uwb_outboard = move_damp_mt(changed_damper_length, uwb_outboard_coord, P3, P4)
    new_damper_mount = move_damp_mt(changed_damper_length, damp_lwb_mt, P1, P2)

    upd_knuckle_v = np.array(new_uwb_outboard) - np.array(new_lwb_outboard)
    knuckle_bottom_to_axle = length_calc(lwb_knuckle, axle_knuckle)
    dir_upd_knuckle_v = upd_knuckle_v / np.linalg.norm(upd_knuckle_v)
    new_knuckle_axle = new_lwb_outboard + dir_upd_knuckle_v * knuckle_bottom_to_axle

    knuckle_dir_old = np.array(uwb_knuckle) - np.array(lwb_knuckle)
    knuckle_dir_new = np.array(new_uwb_outboard) - np.array(new_lwb_outboard)
    rotation_axis = np.cross(knuckle_dir_old, knuckle_dir_new)
    rotation_axis /= np.linalg.norm(rotation_axis)
    theta = angle_finder(knuckle_dir_old, knuckle_dir_new)

    axle_wheel_local = np.array(axle_wheel) - np.array(axle_knuckle)
    rotated_axle_wheel_local = rodriguez_rotation(axle_wheel_local, rotation_axis, theta)
    upd_axle_wheel = new_knuckle_axle + rotated_axle_wheel_local
    upd_axle_v = upd_axle_wheel - new_knuckle_axle
    upd_unit_wheel_v_perp = upd_axle_v / np.linalg.norm(upd_axle_v)
    wheel_width = length_calc(wheel_center, axle_wheel)
    upd_wheel_centre = upd_axle_wheel + upd_unit_wheel_v_perp * wheel_width

    delta_kingpin_angle = angle_finder(np.array(uwb_knuckle) - np.array(lwb_knuckle),
                                   np.array(new_uwb_outboard) - np.array(new_lwb_outboard))
    camber_deg = np.degrees(delta_kingpin_angle)
    return new_lwb_outboard, new_uwb_outboard, new_damper_mount, new_knuckle_axle, upd_axle_wheel, upd_wheel_centre, camber_deg


# --- Rotation functions ---

def vector_rotator(vector1, new_vector1, pivot_pt, end_pt):
    a = vector1 / np.linalg.norm(vector1)
    b = new_vector1 / np.linalg.norm(new_vector1)
    cross = np.cross(a, b)
    norm_cross = np.linalg.norm(cross)

    if norm_cross < 1e-8:
        if np.dot(a, b) > 0:
            R = np.eye(3)
        else:
            perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
            axis = np.cross(a, perp)
            axis /= np.linalg.norm(axis)
            R = rodriguez_rotation_matrix(axis, np.pi)
    else:
        axis = cross / norm_cross
        theta = np.arctan2(norm_cross, np.dot(a, b))
        kx = np.array([[0, -axis[2], axis[1]],
                       [axis[2], 0, -axis[0]],
                       [-axis[1], axis[0], 0]])
        R = np.eye(3) + np.sin(theta) * kx + (1 - np.cos(theta)) * (kx @ kx)

    return pivot_pt + R @ (end_pt - pivot_pt)


def rodriguez_rotation(vector, rot_axis, theta_rad):
    k = rot_axis / np.linalg.norm(rot_axis)
    v = vector
    return v * np.cos(theta_rad) + np.cross(k, v) * np.sin(theta_rad) + k * np.dot(k, v) * (1 - np.cos(theta_rad))


def rodriguez_rotation_matrix(axis, theta):
    kx = np.array([[0, -axis[2], axis[1]],
                   [axis[2], 0, -axis[0]],
                   [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(theta) * kx + (1 - np.cos(theta)) * (kx @ kx)

# --- Utility functions ---
def length_calc(a, b):
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    while len(a) < 3:
        a = np.append(a, 0.0)
    while len(b) < 3:
        b = np.append(b, 0.0)
    return np.linalg.norm(b - a)

def angle_finder(v1, v2):
    return np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))



def plot_cylinder_along_vector(ax, start_point, end_point, radius, resolution=50, color='black', alpha=1):
    # Vector from start to end
    v = np.array(end_point) - np.array(start_point)
    height = np.linalg.norm(v)
    v_unit = v / height

    # Create orthonormal basis
    not_v = np.array([1, 0, 0]) if abs(v_unit[0]) < 0.9 else np.array([0, 1, 0])
    n1 = np.cross(v_unit, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v_unit, n1)

    # Create circle in XY plane
    theta = np.linspace(0, 2 * np.pi, resolution)
    circle = np.array([radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)])

    # Sweep circle along vector
    z = np.linspace(0, height, resolution)
    cylinder = np.zeros((3, resolution, resolution))
    for i in range(resolution):
        point = start_point + v_unit * z[i]
        cylinder[:, i, :] = point[:, None] + n1[:, None] * circle[0] + n2[:, None] * circle[1]

    # Plot surface
    ax.plot_surface(cylinder[0], cylinder[1], cylinder[2], color=color, alpha=alpha, edgecolor='red')

def shortest_length_calc(point, vec_pt_1, vec_pt_2):
    Point = np.array(point)
    A = np.array(vec_pt_1)
    B = np.array(vec_pt_2)
    pivot_vector = B - A
    pivot_unit = pivot_vector / np.linalg.norm(pivot_vector)
    AP = Point - A
    proj_length = np.dot(AP, pivot_unit)
    proj_point = A + proj_length * pivot_unit
    radius_vector = proj_point - Point
    radius_length = np.linalg.norm(radius_vector)
    return radius_vector, radius_length


def move_damp_mt(changed_damper_length, damp_lwb_mount, pivot_pt1, pivot_pt2):
    # Define pivot axis
    A = np.array(pivot_pt1)
    B = np.array(pivot_pt2)
    axis = B - A
    axis_unit = axis / np.linalg.norm(axis)

    # Vector from axis to outboard point
    P = np.array(damp_lwb_mount)
    AP = P - A
    proj_length = np.dot(AP, axis_unit)
    center = A + proj_length * axis_unit  # closest point on axis
    radius_vector = P - center

    radius_length = np.linalg.norm(radius_vector)
    arc_length = changed_damper_length - damper_length
    theta = arc_length / radius_length

    rotated_vector = rodriguez_rotation(radius_vector, axis, theta)
    new_coord = center + rotated_vector

    return new_coord





def plot_points_and_dampers(new_lwb_outboard, ax):
    # Plot pivots
    for pt, label in zip([P1, P2, P3, P4], ['P1', 'P2', 'P3', 'P4']):
        ax.scatter(*pt, color='blue')
        ax.text(*pt, label, fontsize=10)

    # Original damper
    ax.plot([damp_chassis_mt[0], lwb_outboard_coord[0]],
            [damp_chassis_mt[1], lwb_outboard_coord[1]],
            [damp_chassis_mt[2], lwb_outboard_coord[2]],
            color='green', label='Original Damper')

    # Changed damper
    ax.plot([damp_chassis_mt[0], new_lwb_outboard[0]],
            [damp_chassis_mt[1], new_lwb_outboard[1]],
            [damp_chassis_mt[2], new_lwb_outboard[2]],
            color='red', linestyle='--', label='Changed Damper')


def plot_initial_vectors(ax):
    ax.set_title("Initial Geometry Vectors")

    def draw_vector(start, end, label, color='blue'):
        start = np.array(start)
        end = np.array(end)
        vec = end - start
        ax.quiver(*start, *vec, color=color, arrow_length_ratio=0.1)
        ax.text(*end, label, fontsize=9, color=color)

    # Lower wishbone
    draw_vector(P1, lwb_outboard_coord, "Lower WB", 'red')
    draw_vector(P2, lwb_outboard_coord, "Lower WB", 'red')

    # Upper wishbone
    draw_vector(P3, uwb_outboard_coord, "Upper WB", 'green')
    draw_vector(P4, uwb_outboard_coord, "Upper WB", 'green')

    # Damper
    draw_vector(damp_chassis_mt, damp_lwb_mt, "Damper", 'purple')

    # Knuckle
    draw_vector(lwb_knuckle, uwb_knuckle, "Knuckle", 'orange')

    # Axle
    draw_vector(axle_knuckle, axle_wheel, "Axle", 'cyan')

    wheel_radius = length_calc(wheel_top, wheel_bottom) / 2
    plot_cylinder_along_vector(ax, axle_wheel, wheel_center, wheel_radius)

def plot_final_vectors(ax, new_lwb_outboard, new_uwb_outboard, new_damper_mount, new_knuckle_axle, upd_axle_wheel, upd_wheel_centre):    
    ax.set_title("Final Geometry Vectors")

    def draw_vector(start, end, label, color='blue'):
        start = np.array(start)
        end = np.array(end)
        vec = end - start
        ax.plot([start[0], end[0]],
                [start[1], end[1]],
                [start[2], end[2]],
                linestyle='--', color=color)
        ax.text(*end, label, fontsize=9, color=color)

    # Lower wishbone
    draw_vector(P1, new_lwb_outboard, "Lower WB", 'red')
    draw_vector(P2, new_lwb_outboard, "Lower WB", 'red')

    # Upper wishbone
    draw_vector(P3, new_uwb_outboard, "Upper WB", 'green')
    draw_vector(P4, new_uwb_outboard, "Upper WB", 'green')

    # Damper
    draw_vector(damp_chassis_mt, new_damper_mount, "Damper", 'purple')

    # Knuckle
    draw_vector(new_lwb_outboard, new_uwb_outboard, "Knuckle", 'orange')

    # Axle
    draw_vector(new_knuckle_axle, upd_axle_wheel, "Axle", 'cyan')

    wheel_radius = length_calc(wheel_top, wheel_bottom) / 2
    plot_cylinder_along_vector(ax, upd_axle_wheel, upd_wheel_centre, wheel_radius)
    
def animate_suspension():
    fig = plt.figure(figsize=(16, 8))
    fig.set_constrained_layout(True)
    ax3d = fig.add_subplot(121, projection='3d')
    ax2d = fig.add_subplot(122)
    fig.suptitle("Suspension Animation and Camber Angle")

    wheel_radius = length_calc(wheel_top, wheel_bottom) / 2

    # Full sweep: compression and rebound
    single_sweep = np.concatenate([
        np.linspace(damper_length - 140, damper_length + 10, 30),
        np.linspace(damper_length + 10, damper_length - 140, 30)
    ])

    # Repeat the sweep 3 times
    damper_range = np.tile(single_sweep, 3)


    # Preallocate camber history for consistent indexing
    camber_history = np.zeros(len(damper_range))

    # Set fixed axis limits once
    ax2d.set_xlim(np.min(damper_range), np.max(damper_range))
    ax2d.set_ylim(-5, 5)
    ax2d.set_xlabel("Damper Length (mm)")
    ax2d.set_ylabel("Camber Angle (°)")
    ax2d.set_title("Camber vs Damper Length")

    def update(frame):
        ax3d.cla()
        ax3d.view_init(elev=5, azim=25)
        ax3d.set_xlim(400, 900)
        ax3d.set_ylim(-800, -100)
        ax3d.set_zlim(-400, 400)
        ax3d.set_title(f"Damper Length: {damper_range[frame]:.1f} mm")

        # Ground plane
        ground_z = -232.921
        x = np.linspace(400, 900, 2)
        y = np.linspace(-800, -100, 2)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, ground_z)
        ax3d.plot_surface(X, Y, Z, color='gray', alpha=0.5)

        # Geometry update
        new_lwb_outboard, new_uwb_outboard, new_damper_mount, new_knuckle_axle, upd_axle_wheel, upd_wheel_centre, camber_deg = simulate_geometry(damper_range[frame])
        plot_final_vectors(ax3d, new_lwb_outboard, new_uwb_outboard, new_damper_mount, new_knuckle_axle, upd_axle_wheel, upd_wheel_centre)

        # Update camber history and plot
        # 2D camber plot
        camber_history[frame] = camber_deg
        ax2d.cla()
        ax2d.plot(damper_range[:frame+1], camber_history[:frame+1], color='black')
        # Add moving dot at current frame
        dot = ax2d.plot(damper_range[frame], camber_history[frame], 'ro', markersize=8)

        # Fixed axis limits — adjust these as needed
        ax2d.set_xlim(350, 600)   # Example range for damper length
        ax2d.set_ylim(-2, 8)      # Example range for camber angle

    ani = animation.FuncAnimation(fig, update, frames=len(damper_range), interval=100)
    plt.show()
    plt.close(fig)

fig_static = plt.figure(figsize=(10, 8))
ax_static = fig_static.add_subplot(111, projection='3d')
plot_points_and_dampers(lwb_outboard_coord, ax_static)
plot_initial_vectors(ax_static)
ax_static.set_title("Original Suspension Geometry")
plt.show()



main()
animate_suspension()


 
