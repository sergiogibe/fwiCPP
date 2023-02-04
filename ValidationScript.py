import numpy as np
import csv
import sys

'''
Validation considering bellow problem parameters:

Domain = 2x2 km
Mesh = 50x50 elements
T = 2.6s
dt = 0.002s
I = 1.0
freq = 2Hz
vel = 1.0km/s

'''

class Valid:

    def __init__(self):
        self.tol = 10**(-10)

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
                
    def run(self):
        self.val_stiffness()
        self.val_mass()
        self.val_solution()
    
    def val_stiffness(self):
        with open('./validationData/kvalid.npy', 'rb') as f:
            kvalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(kvalid, self.stiffness)
            if result < self.tol:
                print("The stiffness matrix was validated with success.\n")
            else:
                print(f"The stiffness matrix is not correct. Error: {result}\n")
        return 
        
    def val_mass(self):
        with open('./validationData/mvalid.npy', 'rb') as f:
            mvalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(mvalid, self.mass)
            if result < self.tol:
                print("The mass matrix was validated with success.\n")
            else:
                print(f"The mass matrix is not correct. Error: {result}\n")
        return 
    
    def val_solution(self):
        with open('./validationData/svalid.npy', 'rb') as f:
            svalid = np.load(f,allow_pickle=True)
            result = Valid.least_squares(svalid, self.solution, time=True)
            if result < self.tol:
                print("The solution was validated with success.")
            else:
                print(f"The solution is not correct. Error: {result}")
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
        Row 9 - Solution field
        *Row 10 - Validation script only (Stiffness)
        *Row 11 - Validation script only (Mass)
        
        '''
        row_counter = 0
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if row_counter == 0:
                    model.fill(row[0],row[8],row[10])
                if row_counter == 9:
                    for i in range(model.nshots*model.nsteps*model.nn):
                        model.raw_solution[i] = float(row[i])
                if row_counter == 10:
                    for i in range(model.nn*model.nn):
                        model.raw_stiffness[i] = float(row[i])
                if row_counter == 11:
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
    Valid.parse_csv_file(file_path="./dataraw/data.csv", model=valid)
    valid.run()
    

    print("\n    ...   \n\nDone. ")
    
    
if __name__ == "__main__":
    main()