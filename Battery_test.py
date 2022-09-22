import pybamm
import pybamm as pb

class batterytest:
    @staticmethod
    def test():

        model = pb.lithium_ion.DFN()
        sim = pb.Simulation(model)
        sim.solve([0, 3600])
        sim.plot(linestyles=["-."])

        return 0


