import matplotlib.pyplot as plt
import numpy as np
import math as m

def extract_fit_results(filename, observable_names):
    """
    Extract fit results, errors, and Chi²/ndf values from a text file.
    """
    results = []
    errors = []
    errors_stat = []
    errors_templ = []
    chi2_ndf_values = []
    slope = []
    counter = 0
    try:
        with open(filename, 'r') as file:
            for line in file:
                if not line.strip():
                    continue
                
                # Extract fit results (using split)
                if "DoLiTeFit" in line and "Result:" in line:
                    print(line)
                    parts = line.split()
                    print(parts[12])
                    try:
                        result_idx = parts.index("Result:")
                        result_val = float(parts[result_idx + 1])
                        error_val = float(parts[6])
                        error_templ_val = float(parts[12])
                        # Find error after "+-"
                        #for i, part in enumerate(parts):
                        #    if part == "+/-":
                        #        error_val = float(parts[i + 1])
                        #        break
                        #else:
                        #    continue  # Skip if no "+-"
                        
                        results.append(result_val)
                        errors.append(error_val)
                        errors_templ.append(error_templ_val)
                    except (ValueError, IndexError):
                        continue
                
                # Extract Chi²/ndf
                elif "Chi^2/ndf" in line:
                    print(line)
                    parts = line.split()
                    try:
                        chi2_idx = parts.index("Chi^2/ndf")
                        chi2_val = float(parts[chi2_idx + 1])
                        chi2_ndf_values.append(chi2_val)
                    except (ValueError, IndexError):
                        continue

                elif "STAT_DATA" in line:
                    print(line)
                    parts = line.split()
                    try:
                        error_stat_val = float(parts[1])
                        errors_stat.append(error_stat_val)
                    except (ValueError, IndexError):
                        continue
                elif "Tmpl_169.00" in line:
                    data_line = next(file).strip().split()
                    slope.append(abs(float(data_line[10]) - float(data_line[4]))/float(data_line[0]))
                    print("Counter ",counter)
                    counter += 1
                    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return [], [], [], []
    except Exception as e:
        print(f"Error reading file: {e}")
        return [], [], [], []
    
    # Trim all lists to the same length
    min_len = min(len(results), len(errors), len(chi2_ndf_values), len(observable_names))
    results = results[:min_len]
    errors = errors[:min_len]
    errors_templ = errors_templ[:min_len]
    errors_stat = errors_stat[:min_len]
    chi2_ndf_values = chi2_ndf_values[:min_len]
    observable_names = observable_names[:min_len]
    slope = slope[:min_len]
    
    return results, errors, errors_templ, errors_stat, chi2_ndf_values, observable_names, slope

def create_two_panel_plot(observable_names, results, errors, errors_templ, errors_stat, chi2_ndf_values, sensitivity, title="Fit Results and Chi²/ndf"):
    """
    Create two-panel plot sharing only the x-axis (not y-axis).
    - Upper: Fit results with error bars
    - Lower: Chi²/ndf values
    """
    #fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [2, 1, 1, 1]})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    #fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10), sharex=True)

    # Upper panel: Fit results
    total_error = np.sqrt(np.array(errors_templ)**2+np.array(errors)**2)
    weights = 1/total_error**2
    mean = np.sum(np.array(results)*weights) / np.sum(weights)
    print("Mean ",mean)
    ax1.fill_between(observable_names, np.array(results)-total_error, np.array(results)+total_error, step='mid',  color='#8C000F', label="total error")
    ax1.fill_between(observable_names, np.array(results)-np.array(errors), np.array(results)+np.array(errors), step='mid', color='#7BC8F6', label="exp. error")
    ax1.errorbar(observable_names, results, yerr=errors_stat, fmt='o', capsize=5, color='black',
                 ecolor='black', elinewidth=2, capthick=2, markersize=5, label="stat. error")
    ax1.set_ylabel("Fit Result", fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_ylim(169.1, 175.9)
    ax1.axhline(y=172.5, color='red', linestyle='--', linewidth=2, alpha=0.7, label=r'$m_{t}^{true}=172.5$')

    ax1.text(-0.5, 170.2, f"$m_{{t}}^{{fit,avg}} = {mean:.2f}$",
         fontsize=12, bbox=dict(boxstyle="round,pad=0.3", alpha=0.7, facecolor="white"))
    #ax1.text(-0.5, 171.35, r"9 $m_{bl}$ and 4 $m_{bW}$ bins",
    #     fontsize=12, bbox=dict(boxstyle="round,pad=0.3", alpha=0.7, facecolor="white"))
    ax1.legend(loc='upper left')
    
    ax2.clear()
    for i, val in enumerate(chi2_ndf_values):
        marker = '^' if val > 4 else 'o'
        if val > 4:
            val = 3.95
        #color = 'red' if val > 4 else 'green'
        #size = 10 if val > 4 else 80
        ax2.plot(observable_names[i], val, marker=marker, color='blue', markersize=8)
    
    #ax2.plot(observable_names, chi2_ndf_values, 'o', color='blue', markersize=8, linewidth=2)
    ax2.set_ylabel("Chi²/ndf", fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0,4)
    ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Target: 172.5')
    #
    #ax3.clear()
    #for i, val in enumerate(sensitivity):
    #    x = val/(errors_templ[i]/results[i])
    #    marker = '^' if x > 4 else 'o'
    #    if x > 4:
    #        x = 3.95
    #    ax3.plot(observable_names[i], x, marker=marker, color='blue', markersize=8)
    #    #print(f"Slope {val} and rel error {errors_templ[i]/results[i]} for bin {i}")
    #ax3.set_ylabel(r"$\frac{|d\sigma_{174}-d\sigma_{171}|}{d\sigma_{172.5}}/\frac{\sigma_{templ}}{m_{t,fit}}$", fontsize=12)
    #ax3.grid(True, alpha=0.3)
    #ax3.set_ylim(0,4)
    #ax3.axhline(y=1.0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    #
    #ax4.clear()
    #for i, val in enumerate(sensitivity):
    #    x = val/(total_error[i]/results[i])
    #    marker = '^' if x > 4 else 'o'
    #    if x > 4:
    #        x = 3.95
    #    ax4.plot(observable_names[i], x, marker=marker, color='blue', markersize=8)
    #        
    #ax4.set_ylabel(r"$\frac{|d\sigma_{174}-d\sigma_{171}|}{d\sigma_{172.5}}/\frac{\sigma_{tot.}}{m_{t,fit}}$", fontsize=12)
    #ax4.grid(True, alpha=0.3)
    #ax4.set_ylim(0,4)
    #ax4.axhline(y=1.0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    plt.savefig("fit_results_with_chi2.pdf", format="pdf", dpi=300, bbox_inches='tight')
    

def create_latex_table(observable_names, results, errors, errors_templ, errors_stat,
                       output_file="table.tex"):
    """
    Creates a LaTeX table from the input lists and writes it to a text file.

    Args:
        observable_names: list of observable name strings
        results:          list of result values
        errors:           list of total error values
        errors_templ:     list of template error values
        errors_stat:      list of statistical error values
        output_file:      path to the output .tex file
    """
    if not (len(observable_names) == len(results) == len(errors)
            == len(errors_templ) == len(errors_stat)):
        raise ValueError("All input lists must have the same length.")

    header = (
        r"\begin{table}[h]" + "\n"
        r"  \centering" + "\n"
        r"  \begin{tabular}{lcccc}" + "\n"
        r"    \hline\hline" + "\n"
        r"    Observable & Result & Exp. Error (incl stat) & Templ. Error & Stat. Error \\" + "\n"
        r"    \hline" + "\n"
    )

    rows = []
    counter = 0
    for name, res, err, err_t, err_s in zip(
            observable_names, results, errors, errors_templ, errors_stat):
        row = f"    {name} & {res} & {err} & {err_t} & {err_s} \\\\"
        rows.append(row)
        counter += 1
        if ( counter % 14 == 0 ):
            rows.append(r"    \hline" + "\n")
    footer = (
        r"    \hline\hline" + "\n"
        r"  \end{tabular}" + "\n"
        r"  \caption{Results}" + "\n"
        r"  \label{tab:results}" + "\n"
        r"\end{table}"
    )

    latex = header + "\n".join(rows) + "\n" + footer

    with open(output_file, "w") as f:
        f.write(latex)

    print(f"LaTeX table written to '{output_file}'")


    
def main():
    # Configuration
    input_filename = "log.fit.all.txt"
    #input_filename = "log.fit.mbl9.mbw4.drbl19.X.data.txt"

    #observable_names = []
    #for i in range(8):
    #    observable_names.append(rf"$m(bl) (0-{i})$")
    #    #for i in range(12, 7, -1):       # 12 down to 8
    #    for j in range(9, 4, -1):   # 8 down to 3
    #        #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j})$")
    #        for k in range(24, 14, -1): # 10 down to 3
                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}). ptw (0-{k})$")
                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}), mwbbl (0-{k})$")
                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}), eta(l) (0-{k})$")

                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}). rbq (0-{k})$")
                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}), rbj (0-{k})$")
                #observable_names.append(rf"$m(bl) (0-{i}), m(bW) (0-{j}), y(Whad) ({24-k}-{k})$")



    observable_names = [
        r"$\Delta R(bl)$",
        r"$\Delta R(bW)$",
        r"$\eta(l)$",
        r"$m_{T}(l)$",
        r"$m(bl)$",
        r"$m(bW)$",
        r"$m^{minimax}(Wbbl)$",
        r"$m(Wbbl)$",
        r"$m^{min}(bl)$",
        r"$m_{T}^{min}(bl)$",
        r"$p_{T}(Wbbl)$",
        r"$p_{T}(Wbbl\nu)$",
        r"$p_{T}(b_1)$",
        r"$p_{T}(b_2)$",
        r"$p_{T}(bl)$",
        r"$p_{T}(bl\nu)$",
        r"$p_{T}(bW_{had})$",
        r"$p_{T}(l)$",
        r"$E_{T}^{miss}$",
        r"$p_{T}(owj_1)$",
        r"$p_{T}(owj_2)$",
        r"$p_{T}(W_{had})$",
        r"$R_{bj}$",
        r"$R_{bq}$",
        r"$y(W_{had})$",
        r"$m(bl) (1-11)$",
        r"$m(bW) (1-8)$",
        r"$m(Wbbl) (1-8)$",
	r"$m(bl) + m(bW)$",
	r"$m(bl) + m(bW) + y(W)$",
        r"$m(bl) + m(bW) + y(W) + m(Wbbl)$",
    ]
        
    # Extract data
    results, errors, errors_templ, errors_stat, chi2_ndf_values, names, slopes = extract_fit_results(input_filename, observable_names)
    
    # Check if we have data
    if not results:
        print("No fit results or Chi²/ndf values found.")
        return
    
    # Create the two-panel plot
    create_two_panel_plot(
        observable_names=names,
        results=results,
        errors=errors,
        errors_templ=errors_templ,
        errors_stat=errors_stat,
        chi2_ndf_values=chi2_ndf_values,
        sensitivity=slopes,
        title="Fit Results and Chi²/ndf"
    )

    create_latex_table(
        observable_names=names,
        results=results,
        errors=errors,
        errors_templ=errors_templ,
        errors_stat=errors_stat
        )
    
if __name__ == "__main__":
    main()

