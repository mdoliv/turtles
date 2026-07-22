import Models_2D
import dadi
from glob import glob
import csv

# Model functions
func_nomig = Models_2D.no_mig
func_ex_nomig = dadi.Numerics.make_extrap_log_func(func_nomig)

func_symmig = Models_2D.sym_mig
func_ex_symmig = dadi.Numerics.make_extrap_log_func(func_symmig)

func_scsymmig = Models_2D.sec_contact_sym_mig
func_ex_scsymmig = dadi.Numerics.make_extrap_log_func(func_scsymmig)

func_asymmig = Models_2D.asym_mig
func_ex_asymmig = dadi.Numerics.make_extrap_log_func(func_asymmig)

func_scasymmig = Models_2D.sec_contact_asym_mig
func_ex_scasymmig = dadi.Numerics.make_extrap_log_func(func_scasymmig)

# Species comparison
spp = "Ei x Lo"
pair = "ei_lo"

# FS file
fs = dadi.Spectrum.from_file(f"../dadi/All.{pair}.folded.downproj.fs")
ns = fs.sample_sizes

# Grid sizes
pts_init = max(fs.sample_sizes)
pts = [int(pts_init * 1.1) + 2, int(pts_init * 1.2) + 4, int(pts_init * 1.3) + 6]

# Bootstrap runs
all_boot = [
    dadi.Spectrum.from_file(boot) for boot in glob(f"../dadi/{pair}_bootstrap/*.fs")
]

print(f"Read {len(all_boot)} bootstrap replicates")

# Optimised parameters
no_mig_params = []
sym_mig_params = []
sec_contact_params = []
asym_mig_params = []
sec_contact_asym_mig_params = []
results_file = "../dadi/best_model_comparison.csv"

with open(results_file, "r") as results_fopen:
    reader = csv.reader(results_fopen)
    for row in reader:
        if row[0] != spp:
            continue
        if row[1] == "no_mig":
            no_mig_params = [float(row[3]), float(row[4]), float(row[5])]
        elif row[1] == "split_mig":
            sym_mig_params = [
                float(row[3]),
                float(row[4]),
                float(row[7]),
                float(row[5]),
            ]
        elif row[1] == "sec_contact_sym_mig":
            sec_contact_params = [
                float(row[3]),
                float(row[4]),
                float(row[7]),
                float(row[5]),
                float(row[6]),
            ]
        elif row[1] == "asym_mig":
            asym_mig_params = [
                float(row[3]),
                float(row[4]),
                float(row[8]),
                float(row[9]),
                float(row[5]),
            ]
        elif row[1] == "sec_contact_asym_mig":
            sec_contact_asym_mig_params = [
                float(row[3]),
                float(row[4]),
                float(row[8]),
                float(row[9]),
                float(row[5]),
                float(row[6]),
            ]

# Evaluate models with optimised params
ll_nomig = dadi.Inference.ll_multinom(func_ex_nomig(no_mig_params, ns, pts), fs)
ll_symmig = dadi.Inference.ll_multinom(func_ex_symmig(sym_mig_params, ns, pts), fs)
ll_scsymmig = dadi.Inference.ll_multinom(
    func_ex_scsymmig(sec_contact_params, ns, pts), fs
)
ll_asymmig = dadi.Inference.ll_multinom(func_ex_asymmig(asym_mig_params, ns, pts), fs)
ll_scasymmig = dadi.Inference.ll_multinom(
    func_ex_scasymmig(sec_contact_asym_mig_params, ns, pts), fs
)

with open(f"{pair}_sym_mig.tsv", "w+") as sym_mig_file:
    print("no mig vs. sym mig")
    print(f"likelihood of no mig: {ll_nomig}")
    print(f"likelihood of sym mig: {ll_symmig}")

    p_lrt_sym = list(no_mig_params)
    p_lrt_sym.insert(2, 0)

    adj_sym = dadi.Godambe.LRT_adjust(
        func_ex_symmig, pts, all_boot, p_lrt_sym, fs, nested_indices=[2], multinom=True
    )
    D_sym = 2 * (ll_symmig - ll_nomig)
    D_adj_sym = D_sym * adj_sym
    pval_adj_sym = dadi.Godambe.sum_chi2_ppf(D_adj_sym, weights=(0.5, 0.5))

    print(
        f"p-value for rejecting no_mig with sym_mig (GIM adjusted) {pval_adj_sym:.4f}"
    )
    sym_mig_file.write("pair\tll_nomig\tll_symmig\tpval_GIM\n")
    sym_mig_file.write(f"{pair}\t{ll_nomig}\t{ll_symmig}\t{pval_adj_sym}")

with open(f"{pair}_sc_sym_mig.tsv", "w+") as sc_sym_mig_file:
    print("sym mig vs. sec contact")
    print(f"likelihood of sym mig: {ll_symmig}")
    print(f"likelihood of sec contact: {ll_scsymmig}")

    p_lrt_sc = list(sym_mig_params)
    p_lrt_sc.insert(3, 0)

    adj_sc = dadi.Godambe.LRT_adjust(
        func_ex_scsymmig, pts, all_boot, p_lrt_sc, fs, nested_indices=[3], multinom=True
    )
    D_sc = 2 * (ll_scsymmig - ll_symmig)
    D_adj_sc = D_sc * adj_sc
    pval_adj_sc = dadi.Godambe.sum_chi2_ppf(D_adj_sc, weights=(0.5, 0.5))

    print(
        f"p-value for rejecting sym_mig with sec_contact_mig (GIM adjusted) {pval_adj_sc:.4f}"
    )
    sc_sym_mig_file.write("pair\tll_symmig\tll_scsymmig\tpval_GIM\n")
    sc_sym_mig_file.write(f"{pair}\t{ll_symmig}\t{ll_scsymmig}\t{pval_adj_sc}")

with open(f"{pair}_asym_mig.tsv", "w+") as asym_mig_file:
    print("asym mig vs. sym mig")
    print(f"likelihood of sym mig: {ll_symmig}")
    print(f"likelihood of asym mig: {ll_asymmig}")

    p_lrt_asym = list(sym_mig_params)
    p_lrt_asym.insert(3, sym_mig_params[2])

    adj_asym = dadi.Godambe.LRT_adjust(
        func_ex_asymmig,
        pts,
        all_boot,
        p_lrt_asym,
        fs,
        nested_indices=[3],
        multinom=True,
    )
    D_asym = 2 * (ll_asymmig - ll_symmig)
    D_adj_asym = D_asym * adj_asym
    pval_adj_asym = dadi.Godambe.sum_chi2_ppf(D_adj_asym, weights=(0, 1))

    print(
        f"p-value for rejecting sym_mig with asym_mig (GIM adjusted) {pval_adj_asym:.4f}"
    )
    asym_mig_file.write("pair\tll_symmig\tll_asymmig\tpval_GIM\n")
    asym_mig_file.write(f"{pair}\t{ll_symmig}\t{ll_asymmig}\t{pval_adj_asym}")

with open(f"{pair}_sc_asym_mig.tsv", "w+") as sc_asym_mig_file:
    print("asym mig vs. sec contact asym mig")
    print(f"likelihood of asym mig: {ll_asymmig}")
    print(f"likelihood of sec contact asym mig: {ll_scasymmig}")

    p_lrt_sca = list(asym_mig_params)
    p_lrt_sca.insert(5, 0)

    adj_sca = dadi.Godambe.LRT_adjust(
        func_ex_scasymmig,
        pts,
        all_boot,
        p_lrt_sca,
        fs,
        nested_indices=[5],
        multinom=True,
    )
    D_sca = 2 * (ll_scasymmig - ll_asymmig)
    D_adj_sca = D_sca * adj_sca
    pval_adj_sca = dadi.Godambe.sum_chi2_ppf(D_adj_sca, weights=(0.5, 0.5))

    print(
        f"p-value for rejecting asym_mig with sec_contact_asym_mig (GIM adjusted) {pval_adj_sca:.4f}"
    )
    sc_asym_mig_file.write("pair\tll_asymmig\tll_scasymmig\tpval_GIM\n")
    sc_asym_mig_file.write(f"{pair}\t{ll_asymmig}\t{ll_scasymmig}\t{pval_adj_sca}")
