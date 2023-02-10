import numpy as np
import csv
import sys

'''
Validation considering following problem parameters:

%%%% KVALID , MVALID, SVALID %%%%
Domain = 2x2 km
Mesh = 50x50 elements
T = 2.6s
dt = 0.002s
I = 1.0
freq = 2Hz
vel = 1.0km/s (only background)
Single shot = (x=1.5, y=1.5)
No absorbing boundaries

%%%% KVALID2 , MVALID2, SVALID2 %%%% TODO
Domain = 2x2 km
Mesh = 80x80 elements
T = 2.6s
dt = 0.002s
I = 1.0
freq = 2Hz
vel1 = 1.0km/s 
vel2 = 3.0km/s (single centered square 0.3x0.3 km)
Four shots = (x=0.0, y=0.0), (x=2.0, y=0.0), (x=0.0, y=2.0), (x=2.0, y=2.0)
No absorbing boundaries

'''

class Valid:

    def __init__(self):
        self.tol = 10**(-7)

    def fill(self, nn, nsteps, nshots):
        self.nn = int(nn)
        self.nsteps = int(nsteps)
        self.nshots = int(nshots)
        self.raw_stiffness = np.zeros((self.nn*self.nn))
        self.raw_mass = np.zeros((self.nn*self.nn))
        self.raw_solution = np.zeros((self.nn*self.nsteps*self.nshots))
        
    def update(self):
        self.solution = np.zeros((self.nn,self.nsteps,self.nshots))
        index = 0
        for s in range(self.nshots):
            for t in range(self.nsteps):
                for n in range(self.nn):
                    self.solution[n,t,s] = self.raw_solution[index]
                    index += 1
                    
        self.stiffness = np.zeros((self.nn,self.nn))
        index = 0
        for j in range(self.nn):
            for i in range(self.nn):
                self.stiffness[i,j] = self.raw_stiffness[index]
                index += 1
                
        self.mass = np.zeros((self.nn,self.nn))
        index = 0
        for j in range(self.nn):
            for i in range(self.nn):
                self.mass[i,j] = self.raw_mass[index]
                index += 1
                
    def run(self, kvalid: str, mvalid: str, svalid: str):
        self.val_stiffness(example=kvalid)
        self.val_mass(example=mvalid)
        self.val_solution(example=svalid)
    
    def val_stiffness(self, example: str):
        with open(f'./validationData/{example}.npy', 'rb') as f:
            kvalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(kvalid, self.stiffness)
            if result < self.tol:
                print("The stiffness matrix was validated with success.\n")
            else:
                print(f"The stiffness matrix is not correct. Error value: {result}\n")
        return 
        
    def val_mass(self, example: str):
        with open(f'./validationData/{example}.npy', 'rb') as f:
            mvalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(mvalid, self.mass)
            if result < self.tol:
                print("The mass matrix was validated with success.\n")
            else:
                print(f"The mass matrix is not correct. Error value: {result}\n")
        return 
    
    def val_solution(self, example: str):
        with open(f'./validationData/{example}.npy', 'rb') as f:
            svalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(svalid, self.solution, time=True)
            if result < self.tol*100*self.nn*self.nshots*self.nsteps:
                print("The solution was validated with success.")
            else:
                print(f"The solution is not correct. Error value: {result}")
        return 
            
    @staticmethod
    def least_squares(matA, matB, time=False):
        value = 0
        if time is False:
            for j in range(matA.shape[1]):
                for i in range(matA.shape[0]):
                    value += (matA[i,j]-matB[i,j])**2
            return value
        if time is True:
            for k in range(matA.shape[2]):
                for j in range(matA.shape[1]):
                    for i in range(matA.shape[0]):
                        value += (matA[i,j,k]-matB[i,j,k])**2
            return value
        
    @staticmethod
    def parse_csv_file(file_path, model):
        '''
        Index:
        Row 0 - Problem configuration
        Row 1 - Solution field (FULL)
        Row 2 - Stiffness matrix
        Row 3 - Mass matrix
        
        '''
        row_counter = 0
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if row_counter == 0:
                    model.fill(row[0],row[8],row[10])
                if row_counter == 1:
                    for i in range(model.nshots*model.nsteps*model.nn):
                        model.raw_solution[i] = float(row[i])
                if row_counter == 2:
                    for i in range(model.nn*model.nn):
                        model.raw_stiffness[i] = float(row[i])
                if row_counter == 3:
                    for i in range(model.nn*model.nn):
                        model.raw_mass[i] = float(row[i])
                else:
                    pass
                row_counter += 1
                
        model.update()



def main():
    print("\nVALIDATING DATA\n    ...   \n")
    valid = Valid()
    csv.field_size_limit(sys.maxsize)
    Valid.parse_csv_file(file_path="./Output/valid.csv", model=valid)
    valid.run("kvalid","mvalid","svalid")
    

    print("\n    ...   \n\nDone. ")
    
    
if __name__ == "__main__":
    main()