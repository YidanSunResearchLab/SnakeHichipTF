import os
import numpy as np
import cooler
import argparse
import csv

def compute_quantiles_sparse(clr, chrs_to_use, quantiles):
    all_nonzero = []

    total_entries = 0
    for chrom in chrs_to_use:
        if chrom not in clr.chromnames:
            continue
        print(f"Fetching {chrom} as sparse")
        mat = clr.matrix(balance=False, sparse=True).fetch(chrom)
        all_nonzero.append(mat.data)
        shape = mat.shape
        total_entries += shape[0] * shape[1]

    nonzero_values = np.concatenate(all_nonzero)
    nonzero_q = np.quantile(nonzero_values, quantiles)
    all_q = np.quantile(nonzero_values[nonzero_values>1], quantiles)
    
    # Compute quantiles for all values, accounting for zeros
    #num_zeros = total_entries - nonzero_values.size
    #all_q = np.zeros(len(quantiles))
    #frac_zeros = num_zeros / total_entries if total_entries > 0 else 0
    #for i, q in enumerate(quantiles):
    #    if q <= frac_zeros:
    #        all_q[i] = 0
    #    else:
    #        q_adjusted = (q - frac_zeros) / (1 - frac_zeros)
    #        all_q[i] = np.quantile(nonzero_values, q_adjusted)    

    return all_q, nonzero_q

def write_csv(filename, quantiles, values):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Quantile", "Value", "Log1p_Value"])
        for q, v in zip(quantiles, values):
            writer.writerow([f"{int(q*100)}th", v, np.log1p(v)])

def main():
    parser = argparse.ArgumentParser(description="Fast quantile computation from Hi-C data using sparse matrices")
    parser.add_argument('--input', required=True)
    parser.add_argument('--resolution', type=int, default=5000)
    parser.add_argument('--output_dir', default='quantile_output')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.input.endswith('.mcool'):
        clr = cooler.Cooler(f"{args.input}::/resolutions/{args.resolution}")
    else:
        clr = cooler.Cooler(args.input)

    use_chr_prefix = any(c.startswith('chr') for c in clr.chromnames)
    chrs_to_use = [f"chr{i}" if use_chr_prefix else str(i) for i in range(1, 23)] + [f"chrX" if use_chr_prefix else "X"]

    quantiles = [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]
    all_q, nonzero_q = compute_quantiles_sparse(clr, chrs_to_use, quantiles)

    write_csv(os.path.join(args.output_dir, "quantiles_all.csv"), quantiles, all_q)
    write_csv(os.path.join(args.output_dir, "quantiles_nonzero.csv"), quantiles, nonzero_q)

    print("✅ Quantiles saved.")

if __name__ == "__main__":
    main()
