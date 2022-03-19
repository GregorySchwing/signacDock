# init.py
import signac

project = signac.init_project('ensemble-docking-vina')

for p in [5.63, 7.1]:
    for N in [1, 10, 100]:
        for n in range(1, N+1, 1):
            for d in range(0, 960, 1):
                sp = {'pH': p, 'N': N, 'n': n, "d" : d}
                job = project.open_job(sp)
                job.init()    
