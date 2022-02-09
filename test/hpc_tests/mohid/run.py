from PYOPATRA import *
import os
import h5py
import numpy as np

datadir = os.path.expandvars("$WORK/oil-spill-data")

if __name__ == "__main__":

	# read in lat/long for MOHID grid	
	if rank == 0:
		with h5py.File(datadir+"/Lagrangian_test_file.hdf5", "r") as l:
			latitude = np.array(l["Grid"]["Latitude"][:])
			longitude = np.array(l["Grid"]["Longitude"][:])

	else:
		latitude = longitude = None

	# read in velocities
	fp = MOHIDStyleFileParser()
	fp.read(datadir+"/currents_small.hdf5", latitude, longitude)

	# set up mesh
	if rank == 0: print("Setting up mesh")
	tm2d = TriangularMesh2D()
	tm2d.setup_mesh(fp, 2)

	if rank == 0:
		print("Mesh setup complete. Validating mesh . . .")
	
	locations = tm2d.get_vertex_locations()
	velocities = tm2d.get_velocities(1)
	print(rank, locations.shape)
	print(rank, velocities.shape) 
	sample_inds = np.array([1250938, 1250935, 1249537, 1250936, 1250937, 1249532, 1249536,
       1249535, 1249533, 1249534])
	sample_locations = locations[sample_inds]
	expected_locations = np.array([[18.9565745, 39.402099 ],
		[18.9565745, 39.372099 ],
		[18.9465695, 39.392099 ],
		[18.9565745, 39.382099 ],
		[18.9565745, 39.392099 ],
		[18.9465695, 39.342099 ],
		[18.9465695, 39.382099 ],
		[18.9465695, 39.372099 ],
		[18.9465695, 39.352099 ],
		[18.9465695, 39.362099 ]])
	assert np.allclose(sample_locations, expected_locations), f"Locations {sample_locations} on rank {rank} don't match {expected_locations}!"

	sample_velocities = velocities[sample_inds]
	expected_velocities = np.array([[0.03049565, 0.95504117],
		[0.04714984, 0.95514542],
		[0.03991551, 0.95566958],
		[0.04154022, 0.95638305],
		[0.03591438, 0.95645416],
		[0.06250262, 0.95646954],
		[0.04337472, 0.95688128],
		[0.04733628, 0.9574005 ],
		[0.0567146 , 0.95741165],
		[0.05173129, 0.95758718]])
	
	assert np.allclose(sample_velocities, expected_velocities), f"Sample velocities {sample_velocities} don't match {expected_velocities} for rank {rank}!"

	if rank == 0:
		print("Mesh validation complete.")
