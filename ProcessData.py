import numpy as np
import csv
import sys
import os
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator


class ProblemData:

    def __init__(self):
        self.colors = [
                        (0.5, 0.0, 0.0),  # 'red',
                        (1.0, 0.6, 0.0),  # 'orange',
                        (0.0, 0.5, 1.0),  # 'cyan',
                        (0.0, 0.0, 0.5),  # 'blue',
                        (0.0, 0.0, 0.0)   # 'black'
                      ]   
    
    def fill(self, nn, nel, ned, length, depth, I, freq, dt, nsteps, nrec, nshots, nlvls, save_solution, sample_l, sample_d):
        #General config attributes
        self.nn = int(nn)
        self.nel = int(nel)
        self.ned = int(ned)
        self.ne = self.nel*self.ned
        self.length = float(length)
        self.depth = float(depth)
        self.I = float(I)
        self.freq = float(freq)
        self.dt = float(dt)
        self.nsteps = int(nsteps)
        self.nrec = int(nrec)
        self.nshots = int(nshots)
        self.nlvls = int(nlvls)
        self.save_solution = bool(int(save_solution))
        self.sample_l = int(sample_l)
        self.sample_d = int(sample_d)
        
        #Model attributes
        self.pulse = np.zeros((self.nsteps))
        self.control = np.zeros((self.nn))
        self.levels = np.zeros((self.nlvls))
        self.raw_solution = np.zeros((self.sample_l*self.sample_d*self.nsteps))
        self.raw_connectivity = np.zeros((4*self.ne))
        self.receivers_x = np.zeros((self.nrec))
        self.receivers_y = np.zeros((self.nrec))
        self.sources_x = np.zeros((self.nshots))
        self.sources_y = np.zeros((self.nshots))
        
    def update(self):
        self.solution = np.zeros((self.nn,self.nsteps))
        index = 0
        for t in range(self.nsteps):
            for n in range(self.sample_l*self.sample_d):
                self.solution[n,t] = self.raw_solution[index]
                index += 1
                    
        self.connectivity = np.zeros((self.ne,4),dtype=int) #NOT USING IT YET
        index = 0
        for e in range(self.ne):
            for n in range(4):
                self.connectivity[e,n] = int(self.raw_connectivity[index])
                index += 1
        
        #Receivers and sources in a list of coordinates [(xi,yi)]
        self.receivers = [(self.receivers_x[i],self.receivers_y[i]) for i in range(self.nrec)]
        self.sources = [(self.sources_x[i],self.sources_y[i]) for i in range(self.nshots)]
        
    def log(self):
        print(f"Number of nodes: {self.nn}")
        print(f"Number of elements in the length direction: {self.nel}")
        print(f"Number of elements in the depth direction: {self.ned}")
        print(f"Length of the domain: {self.length}")
        print(f"Depth of the domain: {self.depth}")
        print(f"Ricker pulse intensity: {self.I}")
        print(f"Ricker pulse frequency: {self.freq}")
        print(f"Delta time of newmark: {self.dt}")
        print(f"Number of time steps for the dynamic solution: {self.nsteps}")
        print(f"Number of receivers: {self.nrec}")
        print(f"Number of shots: {self.nshots}")
        print(f"Number of control function levels: {self.nlvls}")
        print(f"Levels: {[self.levels[i] for i in range(self.nlvls)]}")
        print(f"Saving solution to render propagation: {self.save_solution}")
        
    def save(self, folder_name):
        if os.path.exists(f"./{folder_name}") is False:
            os.mkdir(folder_name)
            
        self.plt_pulse(folder_name)
        self.plt_contour(folder_name)
        
        if self.save_solution:
            self.render_propagation(folder_name)
        
    def render_propagation(self, folder_name, shot=0, sample_size=10, save=False, show=True):
        if save is False:
            show = True

        ims = []
        steps = self.nsteps//sample_size

        aux = np.zeros((self.sample_l*self.sample_d))
        axField = np.zeros((self.sample_d,self.sample_l))

        fig, ax = plt.subplots(dpi=300)

        for t in range(steps):
            time = sample_size*t
            if time < self.nsteps:
                for j in range(self.sample_d):
                    for i in range(self.sample_l):
                        axField[(self.sample_d - 1) - j, i] = self.solution[i + j * self.sample_l, time]

                im = ax.imshow(axField, cmap='binary', animated=True)
                ims.append([im])

        ani = animation.ArtistAnimation(fig,ims,interval=5,blit=True,repeat_delay=800)
        
        if show:
            plt.show()
        
        
    def plt_contour(self, folder_name, fig_number=0, name="contour", fill=True, pltSR=True, show=False):
        vp = np.zeros((self.nel+1, self.ned+1))
        for j in range(self.ned+1):
            for i in range(self.nel+1):
                vp[i,j] = self.control[i+j*(self.nel+1)]
                
        plt.figure(fig_number, dpi=300)
        
        if fill:
            func = plt.contourf
        else:
            func = plt.contour

        plot = func(vp,
                    cmap=None,
                    extent=(0.0, self.length, 0.0, self.depth),
                    vmin=self.levels[0], vmax=self.levels[-1],
                    colors=self.colors,
                    levels=self.levels
                   )

        plt.xlim(0.0, self.length)
        plt.ylim(0.0, self.depth)
        plt.xlabel('X (km)')
        plt.ylabel('Z (km)')
        
        if pltSR:
            s = np.zeros((self.nshots, 2))
            r = np.zeros((self.nrec, 2))
            for source in range (self.nshots):
                for i in range(2):
                    s[source,i] = self.sources[source][i]
            for receiver in range (self.nrec):
                for i in range(2):
                    r[receiver,i] = self.receivers[receiver][i]
            plt.scatter(s[:,0], s[:,1], color=(0.0, 0.8, 1.0), marker='o', s=20)
            plt.scatter(r[:,0], r[:,1], color='cyan', marker='x', s=14)

        plt.gca().set_aspect('equal')
        plt.savefig(f"./images/{name}_{fig_number}.png")
        if show:
            plt.show()
        plt.close()
        
    def plt_pulse(self, folder_name, show=False):
        fig, ax = plt.subplots(dpi=300)
        ax.plot([i for i in range(self.nsteps)], self.pulse, 'b')
        plt.xlim(0.0, self.nsteps)
        plt.xlabel('Steps in time')
        plt.ylabel('Pulse amplitude')
        plt.savefig("./images/pulse.png")
        if show:
            plt.show()
        plt.close()
        
    @staticmethod
    def parse_csv_file(file_path, model):
        '''
        Index:
        
        Row 0 - Problem configuration
        Row 1 - Ricker pulse
        Row 2 - Control function
        Row 3 - Control levels
        Row 4 - Receiver xcoord
        Row 5 - Receiver ycoord
        Row 6 - Sources xcoord
        Row 7 - Sources ycoord
        Row 8 - Mesh connectivity
        Row 9 - Solution field (NOT FULL)
        
        '''
        row_counter = 0
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if row_counter == 0:
                    model.fill(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14])
                if row_counter == 1:
                    for i in range(len(row)):
                        model.pulse[i] = float(row[i])
                if row_counter == 2:
                    for i in range(len(row)):
                        model.control[i] = float(row[i])
                if row_counter == 3:
                    for i in range(len(row)):
                        model.levels[i] = float(row[i])
                if row_counter == 4:
                    for i in range(len(row)):
                        model.receivers_x[i] = float(row[i])
                if row_counter == 5:
                    for i in range(len(row)):
                        model.receivers_y[i] = float(row[i])
                if row_counter == 6:
                    for i in range(len(row)):
                        model.sources_x[i] = float(row[i])
                if row_counter == 7:
                    for i in range(len(row)):
                        model.sources_y[i] = float(row[i])
                if row_counter == 8:
                    for i in range(4*model.ne):
                        model.raw_connectivity[i] = float(row[i])
                if row_counter == 9:
                    for i in range(len(row)-1):
                        model.raw_solution[i] = float(row[i])
                        
                else:
                    pass
                row_counter += 1
                
        model.update()


def main():
    print("\nPOST-PROCESSING DATA\n    ...   \n")
    problem = ProblemData()
    csv.field_size_limit(sys.maxsize)
    ProblemData.parse_csv_file(file_path="./Output/data.csv", model=problem)
    problem.log()
    problem.save(folder_name="images")
    
    
    
    print("\n    ...   \n\nDone. ")

if __name__ == "__main__":
    main()
    
