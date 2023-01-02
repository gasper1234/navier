import numpy as np


# diffus rho
def diffus(grid, n, k):
	N_y = len(grid)
	N_x = len(grid[0])
	#določimo začetno vrednost, ne vem še katera je najboljša
	old_vals = np.copy(grid)
	# initial value, will be updated
	new_vals = np.ones_like(grid)*np.average(grid)
	for reps in range(n):
		for y in range(N_y):
			for x in range(N_x):
				# index selection
				ind_ym, ind_yp = y-1, (y+1) % N_y
				ind_xm, ind_xp = x-1, (x+1) % N_x
				ave_val = new_vals[ind_yp, x] + new_vals[ind_ym, x] + new_vals[y, ind_xp] + new_vals[y, ind_xm]
				grid[y, x] = (old_vals[y, x] + k * ave_val / 4) / (1+k)
		new_vals = np.copy(grid)
	return grid

def diffus_speed(vx_grid, vy_grid, n, k):
	N_y = len(vx_grid)
	N_x = len(vx_grid[0])
	#določimo začetno vrednost, ne vem še katera je najboljša
	old_vals_x = np.copy(vx_grid)
	old_vals_y = np.copy(vy_grid)
	# initial value, will be updated
	new_vals_x = np.ones_like(vx_grid)*np.average(vx_grid)
	new_vals_y = np.ones_like(vy_grid)*np.average(vy_grid)
	for reps in range(n):
		for y in range(N_y):
			for x in range(N_x):
				# index selection
				ind_ym, ind_yp = y-1, (y+1) % N_y
				ind_xm, ind_xp = x-1, (x+1) % N_x
				ave_val_x = new_vals_x[ind_yp, x] + new_vals_x[ind_ym, x] + new_vals_x[y, ind_xp] + new_vals_x[y, ind_xm]
				ave_val_y = new_vals_y[ind_yp, x] + new_vals_y[ind_ym, x] + new_vals_y[y, ind_xp] + new_vals_y[y, ind_xm]
				vx_grid[y, x] = (old_vals_x[y, x] + k * ave_val_x / 4) / (1+k)
				vy_grid[y, x] = (old_vals_y[y, x] + k * ave_val_y / 4) / (1+k)
		new_vals_x = np.copy(vx_grid)
		new_vals_y = np.copy(vy_grid)
	return vx_grid, vy_grid

# advect rho
def advect(rho_grid, vx_grid, vy_grid, dt):
	N_y = len(vx_grid)
	N_x = len(vx_grid[0])
	grid = np.zeros_like(rho_grid)
	for y in range(N_y):
		for x in range(N_x):
			x_val = x-vx_grid[y, x]*dt
			x_ind = int(np.floor(x_val)) % N_x
			x_above = (x_val) % N_x - np.floor(x_val) % N_x
			y_val = y-vy_grid[y, x]*dt
			y_ind = int(np.floor(y_val)) % N_y
			y_above = (y_val) % N_y - np.floor(y_val) % N_y
			aa = rho_grid[y_ind, x_ind] * (1-x_above) * (1-y_above)
			y_new_ind = (y_ind+1) % N_y
			ba = rho_grid[y_new_ind, x_ind] * (1-x_above) * y_above
			x_new_ind = (x_ind+1) % N_x
			ab = rho_grid[y_ind, x_new_ind] * x_above * (1-y_above)
			bb = rho_grid[y_new_ind, x_new_ind] * x_above * y_above
			grid[y, x] = aa + ba + ab + bb
	return grid

def advect_speed(vx_grid, vy_grid, dt):
	N_y = len(vx_grid)
	N_x = len(vx_grid[0])
	new_vx_grid = np.zeros_like(vx_grid)
	new_vy_grid = np.zeros_like(vy_grid)
	for y in range(N_y):
		for x in range(N_x):
			x_val = x+vx_grid[y, x]*dt
			x_ind = int(np.floor(x_val)) % N_x
			x_above = (x_val) % N_x - np.floor(x_val) % N_x
			y_val = y+vy_grid[y, x]*dt
			y_ind = int(np.floor(y_val)) % N_y
			y_above = (y_val) % N_y - np.floor(y_val) % N_y
			new_vx_grid[y_ind, x_ind] += vx_grid[y, x] * (1-x_above) * (1-y_above)
			new_vy_grid[y_ind, x_ind] += vy_grid[y, x] * (1-x_above) * (1-y_above)
			y_new_ind = (y_ind+1) % N_y
			new_vx_grid[y_new_ind, x_ind] += vx_grid[y, x] * (1-x_above) * y_above
			new_vy_grid[y_new_ind, x_ind] += vy_grid[y, x] * (1-x_above) * y_above
			x_new_ind = (x_ind+1) % N_x
			new_vx_grid[y_ind, x_new_ind] += vx_grid[y, x] * x_above * (1-y_above)
			new_vy_grid[y_ind, x_new_ind] += vy_grid[y, x] * x_above * (1-y_above)
			new_vx_grid[y_new_ind, x_new_ind] += vx_grid[y, x] * x_above * y_above
			new_vy_grid[y_new_ind, x_new_ind] += vy_grid[y, x] * x_above * y_above
	# hitrost se širi naprej gostota pa "nazaj"
	return new_vx_grid, new_vy_grid

def correct_diverg(vx_grid, vy_grid, n):
	N_y = len(vx_grid)
	N_x = len(vx_grid[0])
	vx_div = np.zeros_like(vx_grid)
	vy_div = np.zeros_like(vy_grid)
	p_grid =  np.zeros_like(vx_grid)
	p_grid_new = np.copy(p_grid)
	grad = np.zeros_like(vx_grid)
	for reps in range(n):
		for y in range(N_y):
			for x in range(N_x):
				ind_ym, ind_yp = y-1, (y+1) % N_y
				ind_xm, ind_xp = x-1, (x+1) % N_x
				grad[y, x] = vy_grid[ind_yp, x] - vy_grid[ind_ym, x] + vx_grid[y, ind_xp] - vx_grid[y, ind_xm]
				p_grid_new[y, x] = (p_grid[ind_yp, x] + p_grid[ind_ym, x] + p_grid[y, ind_xp] + p_grid[y, ind_xm] - grad[y, x]) / 4
		p_grid = p_grid_new
	for y in range(N_y):
		for x in range(N_x):
			ind_ym, ind_yp = y-1, (y+1) % N_y
			ind_xm, ind_xp = x-1, (x+1) % N_x
			vx_div[y, x] = (p_grid[y, ind_xp] - p_grid[y, ind_xm]) / 2
			vy_div[y, x] = (p_grid[ind_yp, x] - p_grid[ind_ym, x]) / 2
	return vx_grid - vx_div, vy_grid - vy_div




