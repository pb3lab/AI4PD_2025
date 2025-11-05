import os
import argparse
import numpy as np
from Bio import PDB

def sample_random_points_within_sphere(radius, sampling_size):
    """
    Generate random points uniformly distributed within a sphere's volume.
    """
    points = []
    for _ in range(sampling_size):
        # Random radius, ensuring uniform sampling within the sphere's volume
        r = radius * np.cbrt(np.random.uniform(0, 1))  # Correct volume scaling with cube root

        # Random spherical coordinates
        theta = np.arccos(2 * np.random.uniform(0, 1) - 1)  # Correctly scaled polar angle
        phi = 2 * np.pi * np.random.uniform(0, 1)           # Azimuthal angle

        # Convert spherical coordinates to Cartesian
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        points.append(np.array([x, y, z]))
    return points

import numpy as np

def sample_homogeneous_points(radius, sampling_size):
    """
    Generate points uniformly distributed within a sphere by dynamically allocating layers
    and populating each sub-sphere surface using the Fibonacci algorithm, with layer rotation.

    Parameters:
        radius (float): The radius of the sphere.
        sampling_size (int): The target number of points (excluding the center).

    Returns:
        List[np.array]: A list of Cartesian coordinates for the sampled points.
    """
    points = []

    # Initialize debug information storage
    rotation_angles = []  # Store rotation angles for each layer
    layer_avg_distances = []  # Store average intra-layer distances
    interlayer_avg_distances = []  # Store average inter-layer distances
    layer_points_coords_all = []  # Store points for each layer


    # Add the center point (radius = 0)
    points.append(np.array([0.0, 0.0, 0.0]))

    # Compute the allowable range for sampling size
    min_sampling_size = int(sampling_size * 0.6)
    max_sampling_size = int(sampling_size * 1.4)

    # Compute volume of the entire sphere
    total_volume = (4 / 3) * np.pi * radius**3

    # Determine optimal number of layers
    num_layers = int(np.cbrt(sampling_size))  # Approximate layer count
    if num_layers < 2:
        num_layers = 2  # Ensure at least two layers
    layer_volumes = []  # Store the volumes of each shell

    # Calculate layer volumes and radii
    layer_radii = np.linspace(0, radius, num_layers + 1)[1:]  # Exclude the center
    for i, r in enumerate(layer_radii, start=1):
        inner_volume = (4 / 3) * np.pi * ((r - radius / num_layers)**3) if i > 1 else 0
        shell_volume = (4 / 3) * np.pi * r**3 - inner_volume
        layer_volumes.append(shell_volume)

    # Allocate points per layer proportionally to shell volume
    total_allocated_points = 0
    layer_points = []
    for volume in layer_volumes:
        points_in_layer = max(1, int((volume / total_volume) * max_sampling_size))
        layer_points.append(points_in_layer)
        total_allocated_points += points_in_layer

    # Ensure at least 4 points in the first layer
    if layer_points[0] < 4:
        layer_points[0] = 4
        total_allocated_points += (4 - layer_points[0])

    # Truncate or expand sampling size if necessary
    if total_allocated_points > max_sampling_size:
        excess = total_allocated_points - max_sampling_size
        for i in range(len(layer_points)):
            if excess <= 0:
                break
            if layer_points[i] > 1:
                layer_points[i] -= 1
                excess -= 1
    elif total_allocated_points < min_sampling_size:
        deficit = min_sampling_size - total_allocated_points
        for i in range(len(layer_points)):
            if deficit <= 0:
                break
            layer_points[i] += 1
            deficit -= 1

    # Suggest ideal sampling size if adjustments were significant
    ideal_sampling_size = sum(layer_points)
    if ideal_sampling_size != sampling_size:
        print(f"### Recommended Sampling Size ###")
        print(f"For optimal coverage, use a sampling size of approximately {ideal_sampling_size}.")

    # Print debug information
    print(f"### Sampling Information ###")
    print(f"Input Radius: {radius}")
    print(f"Requested Sampling Size: {sampling_size}")
    print(f"Adjusted Sampling Size (including center point): {ideal_sampling_size + 1}")
    print(f"Number of Layers: {num_layers}")
    print(f"Layer Radii: {layer_radii}")
    print(f"Points Per Layer: {layer_points}")
    print(f"###########################")

    # Populate each layer using the Fibonacci sphere algorithm
    for layer_idx, (layer_radius, num_points) in enumerate(zip(layer_radii, layer_points), start=1):
        phi = (1 + np.sqrt(5)) / 2  # Golden ratio
        rotation_angle = (layer_idx * 137.5) % 360  # Rotate each layer slightly for better coverage
        rotation_matrix = np.array([
            [np.cos(np.radians(rotation_angle)), -np.sin(np.radians(rotation_angle)), 0],
            [np.sin(np.radians(rotation_angle)), np.cos(np.radians(rotation_angle)), 0],
            [0, 0, 1]
        ])
        for i in range(num_points):
            z = 1 - (2 * i + 1) / num_points  # Evenly spaced along the z-axis
            theta = 2 * np.pi * (i / phi % 1)  # Spread azimuthally
            r_xy = np.sqrt(1 - z**2)
            x = layer_radius * r_xy * np.cos(theta)
            y = layer_radius * r_xy * np.sin(theta)
            z = layer_radius * z
            point = np.dot(rotation_matrix, np.array([x, y, z]))
            points.append(point)

        # Store rotation angle for this layer
        rotation_angles.append(rotation_angle)

        # Calculate intra-layer average distance for current layer
        layer_points_coords = []
        for i in range(num_points):
            z = 1 - (2 * i + 1) / num_points  # Evenly spaced along the z-axis
            theta = 2 * np.pi * (i / phi % 1)  # Spread azimuthally
            r_xy = np.sqrt(1 - z**2)
            x = layer_radius * r_xy * np.cos(theta)
            y = layer_radius * r_xy * np.sin(theta)
            z = layer_radius * z
            layer_points_coords.append(np.dot(rotation_matrix, np.array([x, y, z])))

        # Compute distances within the current layer
        if len(layer_points_coords) > 1:
            distances = [
                np.linalg.norm(np.array(layer_points_coords[i]) - np.array(layer_points_coords[j]))
                for i in range(len(layer_points_coords)) for j in range(i + 1, len(layer_points_coords))
            ]
            avg_intra_layer_distance = np.mean(distances)
        else:
            avg_intra_layer_distance = 0  # Single point layer

        layer_avg_distances.append(avg_intra_layer_distance)

        # Compute inter-layer distances for this layer with all adjacent layers
        interlayer_distances_to_adjacent = []
        for adj_layer_idx in range(layer_idx - 1, layer_idx + 2):
            if adj_layer_idx > 0 and adj_layer_idx <= num_layers and adj_layer_idx != layer_idx:
                prev_or_next_layer_points = points[-len(layer_points_coords):]
                interlayer_distances = [
                    np.linalg.norm(np.array(pt) - np.array(prev_or_next_layer_points[j]))
                    for pt in layer_points_coords
                    for j in range(len(prev_or_next_layer_points))
                ]
                avg_interlayer_distance = np.mean(interlayer_distances)
                interlayer_distances_to_adjacent.append((adj_layer_idx, avg_interlayer_distance))

        interlayer_avg_distances.append(interlayer_distances_to_adjacent)

    # Print rotation and distance information
    print(f"### Rotation and Distance Information ###")
    for idx, (angle, intra_dist, inter_distances) in enumerate(
        zip(rotation_angles, layer_avg_distances, interlayer_avg_distances), start=1
    ):
        print(f"Layer {idx}:")
        print(f"  Rotation Angle: {angle:.2f} degrees")
        print(f"  Avg Intra-Layer Distance: {intra_dist:.4f} units")
        for adj_layer, inter_dist in inter_distances:
            print(f"  Avg Inter-Layer Distance to Layer {adj_layer}: {inter_dist:.4f} units")
    print(f"###########################################")

    return points

### INSERT LATTICE FLAG ###
def sample_homogeneous_points_lattice(radius, sampling_size):
    """
    Generate points uniformly distributed within a sphere using an optimized lattice structure.

    Parameters:
        radius (float): The radius of the sphere.
        sampling_size (int): The target number of points.

    Returns:
        List[np.array]: A list of Cartesian coordinates for the sampled points.
    """
    # Compute lattice spacing for desired point density
    lattice_spacing = (4/3 * np.pi * radius**3 / sampling_size)**(1/3)

    # Generate points in a cube larger than the sphere
    grid_range = np.arange(-radius, radius + lattice_spacing, lattice_spacing)
    points = []

    for x in grid_range:
        for y in grid_range:
            for z in grid_range:
                # Check if the point lies within the sphere
                if x**2 + y**2 + z**2 <= radius**2:
                    points.append(np.array([x, y, z]))

    # If we have more points than needed, reduce to exact sampling_size
    if len(points) > sampling_size:
        # Use spatially aware sampling to trim excess points
        indices = np.random.choice(len(points), sampling_size, replace=False)
        points = [points[i] for i in indices]

    # If we have fewer points, fill the remaining by random sampling within the sphere
    while len(points) < sampling_size:
        r = radius * np.cbrt(np.random.uniform(0, 1))
        theta = np.arccos(2 * np.random.uniform(0, 1) - 1)
        phi = 2 * np.pi * np.random.uniform(0, 1)

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        points.append(np.array([x, y, z]))

    return points

def sample_points_on_sphere_surface(radius, sampling_size):
    """
    Generate homogeneous points on the sphere surface using the Fibonacci sphere method.
    """
    points = []
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    for i in range(sampling_size):
        z = 1 - (2 * i + 1) / sampling_size  # Evenly spaced along z-axis
        theta = 2 * np.pi * (i / phi % 1)    # Spread azimuthally using golden ratio
        r = np.sqrt(1 - z**2)
        points.append(radius * np.array([r * np.cos(theta), r * np.sin(theta), z]))
    return points


def process_pdb(input_pdb, output_pdb):
    """
    Process the output PDB file to keep only the ORI line, append it to the input PDB, and add a TER line.
    """
    with open(output_pdb, 'r') as file:
        lines = file.readlines()

    # Find the ORI line
    ori_line = None
    for line in lines:
        if "ORI ORI X" in line:
            ori_line = line.strip()
            break

    if not ori_line:
        print(f"Error: ORI line not found in {output_pdb}.")
        return

    # Read the input PDB content
    with open(input_pdb, 'r') as file:
        input_lines = file.readlines()

    # Write the processed PDB content back to the output file
    with open(output_pdb, 'w') as file:
        file.writelines(input_lines)  # Write input PDB content
        file.write(ori_line + "\n")  # Add the ORI line at the bottom
        file.write("TER\n")  # Add the TER line


def random_ORI_move(input_pdb, output_dir, radius=5, sampling_size=100, initial_ORI_coord=None, random=True, homogenous_surface_only=False):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)

    # Add ORI atom if initial coordinates are provided
    if initial_ORI_coord:
        # Process the input string into coordinates
        coords = initial_ORI_coord.split()
        if len(coords) != 3:
            print("Error: --initial_ORI_coord must contain exactly three values (x y z).")
            import sys
            sys.exit(1)
        x, y, z = map(float, coords)

        last_residue = list(structure[0].get_residues())[-1]
        ori_res = PDB.Residue.Residue(('H', last_residue.id[1] + 1, ' '), "ORI", ' ')
        ori_atom = PDB.Atom.Atom("ORI", (x, y, z), 1.0, 1.0, ' ', 'ORI', ori_res, element='ORI')
        ori_res.add(ori_atom)

        # Create chain X if it doesn't exist and add the ORI residue
        if "X" not in structure[0]:
            structure[0].add(PDB.Chain.Chain("X"))
        structure[0]["X"].add(ori_res)
    else:
        ori_atom = None
        for atom in structure.get_atoms():
            if atom.name == 'ORI' and atom.parent.parent.id == 'X':
                ori_atom = atom
                break
        if not ori_atom:
            print("No ORI atom found on chain X in the PDB file.")
            print("Use --initial_ORI_coord for the PDB file without ORI atom.")
            return

    ori_position = ori_atom.coord

    # Generate sampling points
    if homogenous_surface_only:
        points = sample_points_on_sphere_surface(radius, sampling_size)
    elif random:
        points = sample_random_points_within_sphere(radius, sampling_size)
    else:
        points = sample_homogeneous_points(radius, sampling_size)

    # Generate output PDB files starting from 01
    for i, displacement in enumerate(points, start=1):
        new_position = ori_position + displacement
        ori_atom.coord = new_position

        output_pdb = os.path.join(output_dir, f"{os.path.basename(input_pdb)[:-4]}_ORI_{i:02d}.pdb")
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)

        # Process the PDB file to include only ORI line, input content, and TER
        process_pdb(input_pdb, output_pdb)

        print(f"Processed ORI atom positioned to {output_pdb} "
              f"({'randomly' if random else 'homogeneously (volume + surface)'}"
              f"{' on surface only' if homogenous_surface_only else ''})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomize or homogeneously sample the ORI token coordination.")
    parser.add_argument('--input_pdb', type=str, required=True, help='Input PDB file including an initial ORI token.')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for PDB files.')
    parser.add_argument('--radius', type=float, default=5, help='Radius to randomize ORI token.')
    parser.add_argument('--sampling_size', type=int, default=100, help='Number of ORI tokens to generate.')
    parser.add_argument('--initial_ORI_coord', type=str, default=None, help='Coordination for an initial ORI token. Example: "0 0 0"')
    parser.add_argument('--random', action='store_true', help='If specified, use random sampling. Default is homogeneous.')
    parser.add_argument('--homogenous_distribution_on_sphere_surface_only', action='store_true',
                        help='If specified, generate homogeneous points on the sphere surface only.')

    args = parser.parse_args()

    # Check for mutually exclusive flags
    if args.random and args.homogenous_distribution_on_sphere_surface_only:
        print("\n### WARNING: TERMINATING SCRIPT ###")
        print("\nCannot specify both --random and --homogenous_distribution_on_sphere_surface_only.")
        import sys
        sys.exit(1)

    # Just pass the arguments to the random_ORI_move function.
    # Do not attempt to modify 'structure' or create the ORI atom here.
    random_ORI_move(
        args.input_pdb,
        args.output_dir,
        radius=args.radius,
        sampling_size=args.sampling_size,
        initial_ORI_coord=args.initial_ORI_coord,
        random=args.random,
        homogenous_surface_only=args.homogenous_distribution_on_sphere_surface_only
    )


