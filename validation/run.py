"""Run as python -m validation.run; writes a measured, portable JSON report."""
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linear_sum_assignment, brentq
from qwalk import Walk, bulk_regime, bulk_discriminant
from .paper_anchors import TABLES


def even_modes(w):
    values, profiles, residuals = w.localized_modes(radius=20, minimum_mass=.95)
    # A full step preserves parity on the infinite line. The paper's tables
    # concern the sector reached from n=0; odd-ring boundary mixing must be tiny.
    even = profiles[w.sites % 2 == 0].sum(axis=0) > .9999
    return values[even], profiles[:, even], residuals[even]


def match(w, reference):
    values, profiles, residuals = even_modes(w)
    if len(values) != 4:
        raise AssertionError(f"expected four even localized modes, found {len(values)}")
    sampled = profiles[np.isin(w.sites, [-4,-2,0,2,4])].T
    rows, cols = linear_sum_assignment(np.max(abs(reference[:,None,:]-sampled[None,:,:]),axis=2))
    return values[cols], profiles[:,cols], residuals[cols]


def main():
    report = {"implementation": "2026 paper-based reimplementation; numerical checks only",
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "table_absolute_tolerance": 1.1e-6,
              "eigen_half_sizes": [80,120], "evolution_half_sizes": [650,800], "tables": {}}
    for label, (model, gamma, anchors) in TABLES.items():
        reference = np.asarray(anchors)
        small, large = Walk(model,gamma,80), Walk(model,gamma,120)
        values, profiles, residuals = match(small,reference)
        v2,p2,r2 = match(large,reference)
        selected = np.isin(small.sites,[-4,-2,0,2,4])
        error = float(np.max(abs(profiles[selected].T-reference)))
        profile_difference = float(np.max(abs(profiles-p2[40:-40])))
        odd_mass = float(profiles[small.sites%2!=0].sum(axis=0).max())
        edge_mass = float(profiles[abs(small.sites)>70].sum(axis=0).max())
        a = Walk(model,gamma,650).evolve(300)
        b = Walk(model,gamma,800).evolve(300)[:,150:-150]
        final_peak = int(np.argmax(a[-1])-650)
        expected_peak = 2 if model == "U1" else -2
        regime = bulk_regime(model,gamma,small.theta1,small.theta2)
        assert error < 1.1e-6, (label,error)
        assert max(residuals) < 1e-10
        assert odd_mass < 1e-10
        assert profile_difference < 1e-8
        assert np.max(abs(a-b)) < 1e-12
        assert np.max(abs(a.sum(axis=1)-1)) < 1e-12
        assert final_peak == expected_peak
        report["tables"][label] = {"model":model,"gamma":gamma,"bulk":regime,
            "max_abs_printed_probability_error":error,"max_eigen_residual":float(max(residuals)),
            "max_odd_sector_mass":odd_mass,"eigen_profile_size_difference":profile_difference,
            "max_eigenvalue_modulus_size_difference":float(max(abs(abs(values)-abs(v2)))),
            "max_eigen_boundary_mass_abs_n_gt_70":edge_mass,
            "moduli_in_table_order":abs(values).tolist(),"computed_even_probabilities":profiles[selected].T.tolist(),
            "peak_sites_T0_T30_T100_T300":(np.argmax(a[[0,30,100,300]],axis=1)-650).tolist(),
            "probability_at_final_peak":float(a[-1,650+final_peak]),
            "evolution_size_difference":float(np.max(abs(a-b))),
            "max_probability_sum_error":float(np.max(abs(a.sum(axis=1)-1))),
            "final_mass_on_minus2_0_plus2":float(a[-1,[648,650,652]].sum())}
    theta1,theta2=small.theta1,small.theta2
    thresholds={}
    for model,brackets in [("U1",[(0,.1),(.1,.4)]),("U2",[(.8,1),(.4,.8)])]:
        roots=[]
        for index,bracket in enumerate(brackets):
            field="max_abs_d0" if index==0 else "min_abs_d0"
            roots.append(brentq(lambda g:bulk_regime(model,g,theta1,theta2)[field]-1,*bracket))
        thresholds[model]=roots
    report["bulk_thresholds_unbroken_partial_then_partial_complete"]=thresholds
    # Appendix C changes theta2, with no use of Appendix D anchors.
    report["appendix_C_qualitative"] = []
    for model, gammas in [("U1",[.1,.3,.5]),("U2",[.8,.5,.2])]:
        for gamma in gammas:
            w=Walk(model,gamma,650,theta2=2*np.pi/3)
            p=w.evolve()
            peak=int(w.sites[np.argmax(p[-1])])
            assert peak == (2 if model=="U1" else -2)
            report["appendix_C_qualitative"].append({"model":model,"gamma":gamma,"theta2":float(w.theta2),
                "bulk":bulk_regime(model,gamma,w.theta1,w.theta2),"T300_peak":peak})
    # Four Fig. 5 time-evolution cases (caption uses U1 gamma=.5 for panel c).
    report["figure5_qualitative"] = []
    for model,gamma,phase in [("U1",.2,.3*np.pi),("U1",.5,.5*np.pi),("U2",.7,.3*np.pi),("U2",.1,.5*np.pi)]:
        w=Walk(model,gamma,650,phase=phase)
        p=w.evolve()
        central=p[:,[648,650,652]].sum(axis=1)
        # Qualitative criterion: more than 80% lies away from these three
        # sites at T=300. No digitized figure or precise curve is claimed.
        assert central[-1] < .2
        report["figure5_qualitative"].append({"model":model,"gamma":gamma,"phase_over_pi":phase/np.pi,
            "mass_minus2_0_plus2_T0_T30_T100_T300":central[[0,30,100,300]].tolist(),
            "T300_peak":int(w.sites[np.argmax(p[-1])])})
    report["boundary_note"] = "For T<=300, support travels at most 600 sites. N=650 and N=800 cannot wrap from initial n=0; odd ring does not change the occupied parity in these evolutions. Eigensystems use smaller rings and explicit size/edge diagnostics."
    report["all_checks_passed"]=True
    out=Path(__file__).with_name("report.json")
    out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"all_checks_passed":True,"report":str(out),"table_max_errors":{k:v["max_abs_printed_probability_error"] for k,v in report["tables"].items()},"bulk_thresholds":thresholds},indent=2))


if __name__ == "__main__":
    main()
