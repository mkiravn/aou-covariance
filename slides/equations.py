"""Estimator equations as PNGs for the deck (matplotlib mathtext)."""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm"})
EQ = {
    "grm": r"$A_{ij}=\frac{1}{M}\sum_{m}\frac{(x_{im}-2p_m)(x_{jm}-2p_m)}{2p_m(1-p_m)}$",
    "pcrelate": r"$\hat\phi_{ij}=\frac{\sum_m (x_{im}-2\mu_{im})(x_{jm}-2\mu_{jm})}{4\sum_m\sqrt{\mu_{im}(1-\mu_{im})\,\mu_{jm}(1-\mu_{jm})}}$",
    "pcrelate_mu": r"$2\mu_{im}=\mathrm{E}[x_{im}\mid \mathrm{PC}_{1..4}(i)]$",
    "wg": r"$\hat k_{jj'S^*}=\frac{\tilde A_{jj'}-\tilde A_{S^*}}{1-\tilde A_{S^*}}$",
    "wg_m": r"$\tilde A_{jj'}=\frac{1}{L}\sum_l \frac{1+(X_{jl}-1)(X_{j'l}-1)}{2},\quad \tilde A_{S^*}=\mathrm{mean}_{j\neq j'}\,\tilde A_{jj'}$",
    "rel": r"$r_{jj'}=A_{jj'}=2\hat\phi_{jj'}=2\hat k_{jj'S^*}$",
    "cov": r"$y_i\,y_j=a+b\,x_{ij}+\varepsilon_{ij}$",
    "h2": r"${\hat{h}}^{2}=\hat{b}$",
    "jk": r"$\mathrm{SE}=\sqrt{\frac{B-1}{B}\sum_{b}(\hat\theta_{(b)}-\bar\theta)^2},\quad B=50$",
    "F": r"$F=1-\frac{H_{\mathrm{obs}}}{2p(1-p)}=\frac{\mathrm{Var}(x)}{2p(1-p)}-1$",
}
os.makedirs(os.path.join(os.path.dirname(__file__), "eq"), exist_ok=True)
for k, e in EQ.items():
    fig = plt.figure(figsize=(0.1, 0.1))
    fig.text(0, 0, e, fontsize=30, color="#1A1A1A")
    fig.savefig(os.path.join(os.path.dirname(__file__), "eq", f"{k}.png"), dpi=300, bbox_inches="tight", pad_inches=0.05, transparent=True)
    plt.close(fig)
