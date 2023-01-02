from RT_fluid_fun_bound import *

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

# vse implementirano za PRP

# definiramo gosto, hitrostno polje (x, y), vsi 2D arrayi
N = 40
rho_grid = np.zeros((N, N))
vx_grid = np.zeros((N, N))
vy_grid = np.zeros((N, N))

#meje
bound = np.full((N, N), False)
for j in range(N):
	bound[j, 0] = True
	bound[j, -1] = True

for j in range(5, 10):
	for i in range(18, 23):
		bound[j, i] = True
print(bound)	


# začetno stanje
for j in range(0, 1):
	for i in range(N//5, 4*N//5):
		if i % 4 == 0:
			rho_grid[j, i] = 1
for j in range(N):
	for i in range(1, N-1):
		if not bound[j, i]:
			vy_grid[j, i] = 1
print(vy_grid)






# podatki
dt = 0.25
N_sim = 300
# saverji
grid_s = np.zeros((N_sim, len(rho_grid), len(rho_grid[0])))
vx_grid_s = np.zeros((N_sim, len(rho_grid), len(rho_grid[0])))
vy_grid_s = np.zeros((N_sim, len(rho_grid), len(rho_grid[0])))

print(vy_grid[:, 5])

for i in range(N_sim):
	# trenutno brez obeh difuzij
	#rho_grid = diffus(rho_grid, 10, dt/4)
	rho_grid = advect(rho_grid, vx_grid, vy_grid, dt)
	# spušča črnilo v vodo
	for k in range(0, 1):
		for l in range(N//5, 4*N//5):
			if l % 4 == 0:
				rho_grid[k, l] = 1
	#vx_grid, vy_grid = diffus_speed(vx_grid, vy_grid, 10, dt*2)
	vx_grid, vy_grid = advect_speed(vx_grid, vy_grid, dt)
	vx_grid, vy_grid = correct_diverg(vx_grid, vy_grid, bound, 10)
	vx_grid, vy_grid = clear_speed(vx_grid, vy_grid, bound)	
	print(np.sum(vy_grid))
	grid_s[i] = rho_grid
	vx_grid_s[i] = vx_grid
	vy_grid_s[i] = vy_grid
	# pri PRP se lahko druge robne pogoje aplicira vmes (morda)
	for j in range(2):
		for i in range(2*N//5, 3*N//5):
			vy_grid[j, i] = 1

#animacija

fig, ax = plt.subplots()

# definiramo za plotanje
x, y = np.meshgrid(np.linspace(0, len(vx_grid)-1, len(vx_grid)), np.linspace(0, len(vy_grid[0])-1, len(vy_grid[0])))
speed = np.sqrt(vx_grid**2 + vy_grid**2)
lw = speed/speed.max()*1.5

rect = patches.Rectangle((18, 5), 4, 5, linewidth=2, edgecolor='k', facecolor='none', zorder=10)
def animate(n):
	ax.clear()
	speed = np.sqrt(vx_grid_s[n]**2 + vy_grid_s[n]**2)
	lw = speed/speed.max()*1.5
	ax.streamplot(x, y, vx_grid_s[n], vy_grid_s[n], linewidth=lw)
	ax.imshow(grid_s[n], cmap='Greys', vmin=0, vmax=1)
	ax.add_patch(rect)
	plt.show()

animate(3)

ani = animation.FuncAnimation(fig, animate, frames = N_sim-1, repeat = True)
ani.save('im_vecna_napredno.gif')
plt.show()
