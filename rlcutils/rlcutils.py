import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.optimize import minimize

class Monomer:
    """
    Fits the RLC normal vectors to monomers.
    """
    def __init__(self, a1: list, a2: list, a3: list, monomer_width: float = 0.403,
    ) -> None:
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self._length = monomer_width

    def __repr__(self):
        return "Monomer at {}".format(self.center)

    @staticmethod
    def _proj(v1, v2):
        """
        Function that calculates the projection of v2 onto v1. Used for Gram-Schmidt orthogonalization.
        :param v1: The first vector.
        :param v2: The second vector.
        :return:
        """
        return np.dot(v1, v2) / np.dot(v1, v1) * v1

    @property
    def center(self):
        """
        Function to calculate the center of the monomer.
        """
        return (self._a1.position + self._a2.position) / 2

    @property
    def t1(self):
        """
        Function to calculate the t1 vector.
        """
        # cross product guarantees orthogonality to t2 and t3
        t1 = np.cross(self.t2, self.t3)
        return t1 / np.linalg.norm(t1)

    @property
    def t2(self):
        """
        Function to calculate the t2 vector.
        """
        t2 = self._a3.position - self.center
        t2 -= self._proj(self.t3, t2)  #  orthogonalize to t3
        return t2 / np.linalg.norm(t2)

    @property
    def t3(self):
        """
        Function to calculate the t3 vector.
        """
        t3 = self._a1.position - self._a2.position
        return t3 / np.linalg.norm(t3)

class Polymer:
    """
    Constructs a polymer from a collection of monomers.
    """
    def __init__(self, atoms, monomers, universe):
        self.atoms = atoms
        self.monomers = monomers
        self.universe = universe
        self.A1 = None
        self.A2 = None
        self.A3 = None

    @property
    def t1(self):
        return [i.t1 for i in self.monomers]

    @property
    def t2(self):
        return [i.t2 for i in self.monomers]

    @property
    def t3(self):
        return [i.t3 for i in self.monomers]

    @property
    def dt1ds(self):
        t = np.array(self.t1)
        return (t[1:] - t[:-1]) / self.length

    @property
    def dt2ds(self):
        t = np.array(self.t2)
        return (t[1:] - t[:-1]) / self.length

    @property
    def dt3ds(self):
        t = np.array(self.t3)
        return (t[1:] - t[:-1]) / self.length
    
    def update(self, result):
        # Updates polymer properties from fitting results. for now, just A values
        self.A1 = result["A"][0]
        self.A2 = result["A"][1]
        self.A3 = result["A"][2]

class Universe:
    """
    The universe class. Everything starts here.
    """
    def __init__(self, tpr_file: str, xtc_file: str
    ) -> None:
        """
        Initializes the universe class. This universe class will be built on top of the MDAnalysis universe class.
        :param tpr_file: The path to the topology file (.tpr) of the GROMACS simulation.
        :param xtc_file: The path to the trajectory file (.xtc) of the GROMACS simulation.
        :return:
        """
        self._mdu = mda.Universe(tpr_file, xtc_file)

class RibbonLikeChain(Universe):
    """
    The ribbon-like chain class. Inherits from the Universe class.
    """
    def __init__(self, tpr_file: str, xtc_file: str,
    ) -> None:
        """
        Initializes the ribbon-like chain class.
        :param tpr_file: The path to the topology file (.tpr) of the GROMACS simulation.
        :param xtc_file: The path to the trajectory file (.xtc) of the GROMACS simulation.
        :param cutoff: The cutoff length in nanometers. This is used to calculate the persistence length.
        """
        super().__init__(tpr_file, xtc_file) # Initialize the Universe class.

        self._polymers = []
        self._monomer_width = None
        self._t1t1 = None
        self._t2t2 = None
        self._t3t3 = None
        self.result = {}
    
    @staticmethod
    def _ribbon_autocorrelation(series
    ) -> list:
        """
        Function to calculate the autocorrelation function for the ribbon-like chain. Keeps the tangent vectors pointing in the same general direction,
        which hides cis-trans rotations in favor of preserving the ribbon shape.
        :param series: The series of vectors to calculate the autocorrelation function for.
        :return:
        """
        # if dot product with previous element (s_prev) is negative, flip the vector and keep going.
        result = []
        s_prev = series[0]
        for i, s in enumerate(series):
            neighbor_dot = np.dot(s, s_prev)
            base_dot = np.dot(s, series[0])
            if neighbor_dot > 0:
                result.append(base_dot)
                s_prev = s
            else:
                result.append(-base_dot)
                s_prev = -s
        return result
    
    @staticmethod
    def _tau_correlation(series
    ) -> list:
        """
        Function to calculate a smoothed autocorrelation function.
        :param series: The series of vectors to calculate the autocorrelation function for.
        :return:
        """
        result = []
        for tau in range(len(series)):
            tau_result = []
            for i in range(len(series) - tau):
                base_dot = np.dot(series[i], series[i + tau])
                tau_result.append(base_dot)
            result.append(np.mean(tau_result))
        return result        

    def gen_polymers(self, residues: int | list[int], atom_selection: tuple, monomer_width: float = 0.403
    ) -> None:
        """
        Function to generate a list of polymers in the simulations.
        :param residues: If int, the number of polymer residues in the simulation. If list, should be a list of residue ids.
        :param atom_selection: Tuple of MDAnalysis atom selection language. There should only be three entries.
                               These three atoms are used to define t1, t2, and t3. t3 will be pointing from atom 1 to atom 2.
                               t2 will be orthogonal to t3 and pointing from atom 3 to the center of atoms 1 and 2.
                               t1 will be the cross product of t2 and t3. This definition was used in the simulation RLC
                               paper for polythiophene-based polymers.
        :param monomer_width: The width of the monomer in nanometers.
        :return:
        """
        self._monomer_width = monomer_width

        if isinstance(residues, int):
            residues = list(range(1,residues+1))
        elif isinstance(residues, list):
            pass
        else:
            raise TypeError("residues must be either an int or a list of ints.")

        for n in residues:
            atoms = self._mdu.select_atoms("resid {}".format(n), updating=True)
            a1list = self._mdu.select_atoms("resid {} and {}".format(n, atom_selection[0]), updating=True,
            )
            a2list = self._mdu.select_atoms("resid {} and {}".format(n, atom_selection[1]), updating=True,
            )
            a3list = self._mdu.select_atoms("resid {} and {}".format(n, atom_selection[2]), updating=True,
            )    

            monomers = [
                Monomer(a1, a2, a3, monomer_width) for a1, a2, a3 in zip(a1list, a2list, a3list)
            ]
            self._polymers.append(Polymer(atoms, monomers, self._mdu))

        return self
    
    def calc_ocf(self, n_timesteps: int,
    ) -> None:
        """
        Function to calculate the orientation correlation function for the triple of RLC normal vectors.
        :param n_timesteps: The number of timesteps to use when calculating the autocorrelation.
                            Calculates the OCF using the last n_timesteps of the simulation trajectory.
        :return:
        """
        assert len(self._polymers) > 0, "No polymers have been generated. Run gen_polymers() first."

        n_monomers = len(self._polymers[0].monomers)
        t1t1 = np.zeros((n_timesteps * len(self._polymers), n_monomers))
        t2t2 = np.zeros((n_timesteps * len(self._polymers), n_monomers))
        t3t3 = np.zeros((n_timesteps * len(self._polymers), n_monomers))
        x = 0

        last_frame = self._mdu.trajectory.n_frames
        for frame in self._mdu.trajectory[int(last_frame - n_timesteps):]:
            for p in self._polymers:
                t1t1[x, :] = self._ribbon_autocorrelation(p.t1)
                t2t2[x, :] = self._ribbon_autocorrelation(p.t2)
                t3t3[x, :] = self._tau_correlation(p.t3)
                x += 1
        
        self._t1t1 = np.mean(t1t1[:], axis=0)
        self._t2t2 = np.mean(t2t2[:], axis=0)
        self._t3t3 = np.mean(t3t3[:], axis=0)

        return self
    
    def fit(self, cutoff: int | float = None,
    ) -> None:
        """
        Function to fit the orientation correlation functions to the RLC model.
        :param: cutoff: The cutoff length in nanometers.
        :return:
        """
        assert len(self._t1t1) > 0, "Orientatin correlation functions have not been generated. Run calc_ocf() first."

        if cutoff:
            assert isinstance(cutoff, (int, float))
        
        indices = np.array(range(len(self._t1t1)))
        L_full = np.array(indices * self._monomer_width)

        # account for cutoff if desired
        idxs = list(range(len(L_full))) if cutoff is None else np.where(L_full < cutoff)
        L = L_full[idxs]
        t1_vec = self._t1t1[idxs]
        t2_vec = self._t2t2[idxs]
        t3_vec = self._t3t3[idxs]

        def objective(A):
            A1, A2, A3 = A[0], A[1], A[2]
            c1 = 1 / (2 * A2) + 1 / (2 * A3)
            c2 = 1 / (2 * A1) + 1 / (2 * A3)
            c3 = 1 / (2 * A1) + 1 / (2 * A2)

            t1_err = np.linalg.norm(np.exp(-c1 * L) - t1_vec)
            t2_err = np.linalg.norm(np.exp(-c2 * L) - t2_vec)
            t3_err = np.linalg.norm(np.exp(-c3 * L) - t3_vec)
            return t1_err + t2_err + t3_err
        
        A0 = np.array([1, 1, 1])
        opt_results = minimize(
            objective, A0, method="nelder-mead", options={"xatol": 1e-8, "disp": False}
        )

        lp = np.zeros(3)
        lp[0] = (0.5 * ((1 / opt_results.x[1]) + (1 / opt_results.x[2]))) ** (-1)
        lp[1] = (0.5 * ((1 / opt_results.x[0]) + (1 / opt_results.x[2]))) ** (-1)
        lp[2] = (0.5 * ((1 / opt_results.x[0]) + (1 / opt_results.x[1]))) ** (-1)

        self.result["opt"] = opt_results
        self.result["A"] = opt_results.x

        self.result["lp"] = lp

        self.result["L_full"] = L_full
        self.result["t1t1_full"] = self._t1t1
        self.result["t2t2_full"] = self._t2t2
        self.result["t3t3_full"] = self._t3t3

        self.result["t1t1_data"] = t1_vec
        self.result["t2t2_data"] = t2_vec
        self.result["t3t3_data"] = t3_vec

        A1, A2, A3 = opt_results.x
        c1 = 1 / (2 * A2) + 1 / (2 * A3)
        c2 = 1 / (2 * A1) + 1 / (2 * A3)
        c3 = 1 / (2 * A1) + 1 / (2 * A2)

        self.result["t1t1_fit"] = np.exp(-c1 * L)
        self.result["t2t2_fit"] = np.exp(-c2 * L)
        self.result["t3t3_fit"] = np.exp(-c3 * L)

        self.result["L"] = L

        return self

    def get_sim_data(self
    ) -> pd.DataFrame:
        """
        Function to return the orientation correlation functions calculated from the simulation data.
        :return:
        """
        return pd.DataFrame({"L_full": self.result["L_full"], "t1t1_full": self.result["t1t1_full"], "t2t2_full": self.result["t2t2_full"], "t3t3_full": self.result["t3t3_full"]})

    def get_fit_data(self
    ) -> pd.DataFrame:
        """
        Function to return the orientation correlation functions fitted to the RLC model.
        :return:
        """
        return pd.DataFrame({"L_fit": self.result["L"], "t1t1_fit": self.result["t1t1_fit"], "t2t2_fit": self.result["t2t2_fit"], "t3t3_fit": self.result["t3t3_fit"]})
