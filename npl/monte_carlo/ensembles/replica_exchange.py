# monte_carlo/Ensembles/replica_exchange.py

import numpy as np
from ..mc_moves import InsertionMove, DeletionMove, DisplacementMove  # Import specific move classes

class Replica:
    def __init__(self, atoms, mu, T, U, N):
        """
        Initialize a single replica with its properties.
        Args:
            atoms (Atoms): The atom configuration for the replica.
            mu (float): The chemical potential.
            T (float): The temperature.
            U (float): The potential energy of the configuration.
            N (int): The number of particles in the replica.
        """
        self.atoms = atoms
        self.mu = mu
        self.T = T
        self.U = U
        self.N = N


class ReplicaExchangeGCMC:
    def __init__(self,
                 replicas,
                 exchange_interval,
                 beta,
                 rng_acceptance,
                 volume,
                 masses,
                 min_distance,
                 max_distance):
        """
        Initialize the ReplicaExchangeGCMC object.

        Args:
            replicas (list of Replica): The list of replicas.
            exchange_interval (int): The interval between replica exchanges.
            beta (float): The inverse temperature (1/kT).
            rng_acceptance (RandomNumberGenerator): A random number generator for move acceptance.
            volume (float): The system volume.
            masses (dict): A dictionary of particle masses.
            min_distance (float): Minimum particle distance for insertion moves.
            max_distance (float): Maximum particle distance for insertion moves.
        """
        self.replicas = replicas
        self.exchange_interval = exchange_interval
        self._beta = beta
        self.rng_acceptance = rng_acceptance
        self.volume = volume
        self.masses = masses
        self.min_distance = min_distance
        self.max_distance = max_distance

    def perform_exchange(self, replica_i, replica_j):
        """
        Perform the exchange of configurations between two replicas.

        Args:
            replica_i (Replica): The first replica.
            replica_j (Replica): The second replica.
        """
        # Calculate the energy difference and decide whether to exchange replicas
        energy_diff = replica_i.U - replica_j.U
        exchange_prob = np.exp(self._beta * energy_diff)

        if exchange_prob > self.rng_acceptance.get_uniform():
            # Perform the exchange
            replica_i.atoms, replica_j.atoms = replica_j.atoms, replica_i.atoms
            replica_i.U, replica_j.U = replica_j.U, replica_i.U
            replica_i.N, replica_j.N = replica_j.N, replica_i.N
            # You can update other properties like mu or T if needed

    def simulate_replica(self, replica, steps):
        """
        Simulate one replica for a given number of steps.
        This will be executed in parallel for each replica.

        Args:
            replica (Replica): The replica to simulate.
            steps (int): The number of Monte Carlo steps to simulate.

        Returns:
            Replica: The final state of the replica.
        """
        for step in range(steps):
            move_type = self.rng_acceptance.get_uniform()

            if move_type < 0.33:
                move = InsertionMove(replica)  # Create an insertion move object
                move_result = move.perform_move()
                if self._acceptance_condition(*move_result):
                    move.apply_move()
            elif move_type < 0.66:
                move = DeletionMove(replica)  # Create a deletion move object
                move_result = move.perform_move()
                if self._acceptance_condition(*move_result):
                    move.apply_move()
            else:
                move = DisplacementMove(replica)  # Create a translation move object
                move_result = move.perform_move()
                if self._acceptance_condition(*move_result):
                    move.apply_move()

        return replica  # Return the final state of the replica after the simulation

    def run_simulation(self, total_MC_steps):
        """
        Run the full simulation, simulating replicas in parallel and performing replica exchanges.
        
        Args:
            total_MC_steps (int): The total number of Monte Carlo steps for the simulation.
        """
        steps_per_replica = total_MC_steps // len(self.replicas)
        
        # Parallelize the simulation of replicas
        with Pool(processes=len(self.replicas)) as pool:
            results = pool.starmap(self.simulate_replica, [(replica, steps_per_replica) for replica in self.replicas])

        # Perform replica exchanges after each batch of steps
        self.exchange_replicas(results)

    def exchange_replicas(self, replicas):
        """
        Perform replica exchange between replicas at regular intervals.
        
        Args:
            replicas (list of Replica): The list of replicas to perform the exchange.
        """
        for i in range(0, len(replicas) - 1, 2):
            replica_i = replicas[i]
            replica_j = replicas[i + 1]
            self.perform_exchange(replica_i, replica_j)
