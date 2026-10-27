"""Independent matrix construction and Fourier checks of paper conventions."""
import unittest
import numpy as np
from qwalk.model import Walk, rotation, gain, bulk_discriminant, bulk_regime


class Operators(unittest.TestCase):
    def test_product_and_departure_phase(self):
        for model,gamma in [("U1",.2),("U2",.7)]:
            w=Walk(model,gamma,6,phase=.37*np.pi,defect=2)
            length=len(w.sites)
            shift=np.zeros((2*length,2*length),dtype=complex)
            direction=1 if model=="U1" else -1
            for j in range(length):
                for spin,d in [(0,direction),(1,-direction)]:
                    shift[2*((j+d)%length)+spin,2*j+spin]=np.exp(1j*w.phase*(w.sites[j]==w.defect))
            coin=lambda m:np.kron(np.eye(length),m)
            g=coin(gain(model,gamma))
            if model=="U1":
                u=shift@np.linalg.inv(g)@coin(rotation(model,w.theta2))@shift@g@coin(rotation(model,w.theta1))
            else:
                r1,r2=coin(rotation(model,w.theta1/2)),coin(rotation(model,w.theta2/2))
                u=r1@shift@r2@g@r2@shift@r1
            np.testing.assert_allclose(w.matrix(),u,atol=1e-14)

    def test_unitary_limits(self):
        for model,gamma in [("U1",0),("U2",1)]:
            u=Walk(model,gamma,10,phase=.4*np.pi).matrix()
            np.testing.assert_allclose(u.conj().T@u,np.eye(len(u)),atol=2e-14)

    def test_bulk_trace_matches_real_space_plane_wave(self):
        for model,gamma in [("U1",.2),("U2",.7)]:
            w=Walk(model,gamma,10,phase=0)
            for m in [0,2,5,9]:
                k=2*np.pi*m/len(w.sites)
                plane=np.exp(1j*k*w.sites)
                u=np.column_stack([w.step(plane[:,None]*np.eye(2)[j])[w.half_size] for j in range(2)])
                np.testing.assert_allclose(np.trace(u)/2,bulk_discriminant(model,gamma,w.theta1,w.theta2,k),atol=1e-14)
                np.testing.assert_allclose(u[0,0]*u[1,1]-u[0,1]*u[1,0],1,atol=1e-13)

    def test_phase_regimes_and_invalid_inputs(self):
        for model,gammas in [("U1",[.01,.2,.3]),("U2",[.95,.7,.1])]:
            for gamma,expected in zip(gammas,["unbroken","partially broken","completely broken"]):
                w=Walk(model,gamma)
                self.assertEqual(bulk_regime(model,gamma,w.theta1,w.theta2)["regime"],expected)
        with self.assertRaises(ValueError): Walk("U2",0)
        with self.assertRaises(ValueError): Walk("other",.1)


if __name__=="__main__": unittest.main()
