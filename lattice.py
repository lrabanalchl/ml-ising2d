class Ising2D:
    """A versatile class for instantianting either the Ising Model or a Lattice Gauge Theory
    in two-dimensional space.
    
    Args:
        lattice_size (int): Linear size of the 2D lattice.
        coupling     (int, optional): Strength of interaction.
        model        (str, optional): Model we want to instantiate.
    
    Attributes:
        *args:            Arguments passed when an instance's class is created.
        spins:            Spin configuration.
        energy:           Total energy of the current configuration.
        magnetization:    Magnetization of the current configuration.
        crit_temperature: Critical temperature.
        neighbors:        Indexes of nearest neighbors.
    
    Methods:
        high_T_spins    --->  Generates a random configuration of spins. Simulates high-T phase.
        get_energy      --->  Computes the total energy of the current configuration.
        get_magnetization ->  Computes the magnetization of the current configuration.
        sweep           --->  Perform one Monte Carlo step. Flip spins.
    """
    def __init__(self, lattice_size, coupling=1, model='ising'):
        self.lattice_size = lattice_size
        self.n_sites = lattice_size ** 2
        
        self.model = model
        if self.model == 'ising':
            # One spin on each site
            self.n_spins = lattice_size ** 2
        elif self.model == 'gauge':
            # One spin on each link
            self.n_spins = 2 * (lattice_size ** 2)
            
        self.coupling = coupling
        self.crit_temperature = 2.0 / (np.log(1.0 + np.sqrt(2)) * self.coupling)
        
        """Identify the indexes of the four nearest neighbors of a spin 
        in a 2D lattice with periodic boundary conditions.
        
        Since we have `n_spins` spins in the lattice, we will need a
        (n_spins, 4) array to store the neighbors.
        
        An attribute `self.neighbours` is created.
        """
        self.neighbours = np.zeros((self.n_sites, 4), dtype=np.int)
        
        for index in range(self.n_sites):
            # Neighbour to the right:
            self.neighbours[index, 0] = index + 1
            if index % self.lattice_size == (self.lattice_size - 1):
                self.neighbours[index, 0] = index + 1 - self.lattice_size
      
            # Downwards neighbour:
            self.neighbours[index, 1] = index + self.lattice_size
            if index >= (self.n_sites - self.lattice_size):
                self.neighbours[index, 1] = index + self.lattice_size - self.n_sites
          
            # Neighbour to the left:
            self.neighbours[index, 2] = index - 1
            if index % self.lattice_size == 0:
                self.neighbours[index, 2] = index - 1 + self.lattice_size
  
            # Upwards neighbour:
            self.neighbours[index, 3] = index - self.lattice_size
            if index <= (self.lattice_size - 1):
                self.neighbours[index, 3] = index - self.lattice_size + self.n_sites
    
    def high_T_spins(self, random_state=0):
        """This function generates a random configuration of spins.
        We simulate a high-T phase by creating a random 
        sample of either +1 or -1 spins.
        
        An attribute `self.spins` is created.
        """
        random.seed(random_state)
        
        self.spins = np.zeros(self.n_spins, dtype=np.int)
        for index in range(self.n_spins):
            self.spins[index] = 2*random.randint(0, 1) - 1
        
    def get_energy(self):
        """Computes the total energy of the current spin configuration.
        This computation avoids double counting.
        
        An attribute `self.energy` is created.
        """
        self.energy = 0
        
        if self.model == 'ising':
            for index in range(self.n_sites):
                self.energy += -self.coupling * (self.spins[index] * self.spins[self.neighbours[index, 0]] + 
                                                 self.spins[index] * self.spins[self.neighbours[index, 1]] )
        elif self.model == 'gauge':
            for index in range(self.n_sites):
                self.energy += -self.coupling * (self.spins[2 * index] * self.spins[(2 * index) + 1] * 
                                                 self.spins[2 * self.neighbours[index, 1]] * 
                                                 self.spins[(2 * self.neighbours[index, 0]) + 1])
    
    def get_magnetization(self):
        """Computes the magnetization of the current spin configuration.
        
        An attribute `self.magnetization` is created.
        """
        self.magnetization = np.sum(self.spins)

    def sweep(self, temperature):
        """One sweep through the lattice is to attempt to flip all spins.
        
        This is one Monte Carlo step. After one Monte Carlo step we have
        obtained a new configuration from the old.

        We accept a spin flip with a probability equal to the Metropolis
        function. The probability of accepting the move depends on the
        energy difference.
        
        Args:
            temperature (float): Value of temperature.
        """
        for _ in range(self.n_spins):
            # Randomly choose which spin to consider flipping
            spin_location = random.randint(0, self.n_spins - 1)
            
            # Change in energy of the proposed move by considering:
            if self.model == 'ising': # 1. only the nearest neighbours
                deltaE = 0
                for neigh in range(4):
                    deltaE += 2*self.coupling*self.spins[site]*self.spins[self.neighbours[site, neigh]]
            elif self.model == 'gauge': # 2. the two plaquettes it will affect
                first_plaquette = spin_location // 2
                
                # Get second_plaquette based on whether the spin is on a horizontal or vertical link
                if (spin_location % 2) == 0:
                    second_plaquette = self.neighbours[first_plaquette, 3]
                else:
                    second_plaquette = self.neighbours[first_plaquette, 2]
    
                deltaE = 2 * self.coupling * ((self.spins[2 * first_plaquette] * self.spins[(2 * first_plaquette) + 1] * 
                                               self.spins[2 * self.neighbours[first_plaquette, 1]] * 
                                               self.spins[(2 * self.neighbours[first_plaquette, 0]) + 1]) +
                                              (self.spins[2 * second_plaquette] * self.spins[(2 * second_plaquette) + 1] * 
                                               self.spins[2 * self.neighbours[second_plaquette, 1]] * 
                                               self.spins[(2 * self.neighbours[second_plaquette, 0]) + 1]))        
            # Flipping the spin
            if (deltaE <= 0) or (random.random() < np.exp(-deltaE / temperature)):
                self.spins[spin_location] = -self.spins[spin_location]
