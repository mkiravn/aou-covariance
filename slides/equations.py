"""Estimator equations as PNGs for the deck (matplotlib mathtext)."""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm"})
EQ = {
    "grm": r"$A_{ij}=\frac{1}{M}\sum_{m}\frac{(x_{im}-2p_m)(x_{jm}-2p_m)}{2p_m(1-p_m)}$",
    "pcrelate": r"$\hat\phi_{ij}=\frac{\sum_m (x_{im}-2\mu_{im})(x_{jm}-2\mu_{jm})}{4\sum_m\sqrt{\mu_{im}(1-\mu_{im})\,\mu_{jm}(1-\mu_{jm})}}$",
    "pcrelate_mu": r"$2\mu_{im}=\mathrm{E}[x_{im}\mid \mathrm{PC}_{1..4}(i)]$",
    "wg": r"$\hat\beta_{ij}=\frac{M_{ij}-\bar M}{1-\bar M}$",
    "wg_m": r"$M_{ij}=\frac{1}{M}\sum_m \frac{1+(x_{im}-1)(x_{jm}-1)}{2},\quad \bar M=\mathrm{mean}_{i\neq j}\,M_{ij}$",
    "rel": r"$r_{ij}=A_{ij}=2\hat\phi_{ij}=2\hat\beta_{ij}$",
    "cov": r"$y_i\,y_j=a+b\,r_{ij}+\varepsilon_{ij}$",
    "jk": r"$\mathrm{SE}=\sqrt{\frac{B-1}{B}\sum_{b}(\hat\theta_{(b)}-\bar\theta)^2},\quad B=50$",
    "F": r"$F=1-\frac{H_{\mathrm{obs}}}{2p(1-p)}=\frac{\mathrm{Var}(x)}{2p(1-p)}-1$",
}
os.makedirs(os.path.join(os.path.dirname(__file__), "eq"), exist_ok=True)
for k, e in EQ.items():
    fig = plt.figure(figsize=(0.1, 0.1))
    fig.text(0, 0, e, fontsize=30, color="#1A1A1A")
    fig.savefig(os.path.join(os.path.dirname(__file__), "eq", f"{k}.png"), dpi=300, bbox_inches="tight", pad_inches=0.05, transparent=True)
    plt.close(fig)
