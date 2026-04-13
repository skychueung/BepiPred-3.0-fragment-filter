### IMPORTS AND STATIC PATHS ###
from bp3 import bepipred3
from bp3.fragment_filter import run_fragment_filter
from bp3.physicochemical_filter import run_physicochemical_filter, run_physicochemical_filter_layered
from pathlib import Path
import argparse
import zipfile

WORK_DIR = Path(Path(__file__).parent.resolve())

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Make B-cell epitope predictions from fasta file.")
parser.add_argument("-i", required=True, action="store", dest="fasta_file", type=Path, help="Fasta file contianing antigens")
parser.add_argument("-o", required=True, action="store", dest="out_dir", type=Path, help="Output directory to store B-cell epitope predictions.")
parser.add_argument("-pred", action="store", choices=["mjv_pred", "vt_pred"], required=True, dest="pred", help="Majorty vote ensemble prediction or variable threshold predicition on average ensemble posistive probabilities. ")
parser.add_argument("-add_seq_len", action="store_true", dest="add_seq_len", help="Add sequence lengths to esm-encodings. Default is false. On the web server this option is set to true.")
parser.add_argument("-esm_dir", action="store", default=WORK_DIR / "esm_encodings", dest="esm_dir", type=Path, help="Directory to save esm encodings to. Default is current working directory.")
parser.add_argument("-t", action="store", default=0.1512, type=float, dest="var_threshold", help="Threshold to use, when making predictions on average ensemble positive probability outputs. Default is 0.1512")
parser.add_argument("-top", action="store", default=0.2, type=float, dest="top_cands", help="Top percentage of epitope residues. Default is 0.2")
parser.add_argument("-rolling_window_size", default=9, type=int, dest="rolling_window_size", help="Window size to use for rolling average on B-cell epitope probability scores. Default is 9")
parser.add_argument("-plot_linear_epitope_scores", action="store_true", dest="plot_linear_epitope_scores", help="Use linear B-cell epitope probability scores for plot. Default is false.")
parser.add_argument("-z", action="store_true", dest="zip_results", help="Specify option to create zip the bepipred-3.0 results (except the interactive .html figure). Default is false.")

# fragment filter patch
parser.add_argument("--fragment_filter", action="store_true", dest="fragment_filter",
                    help="Run optional post-processing fragment filter on predicted peptide table.")
parser.add_argument("--frag_input_table", action="store", dest="frag_input_table", type=Path, default=None,
                    help="Optional explicit input peptide table (.tsv/.csv). If omitted, the script will try to auto-detect a table in out_dir.")
parser.add_argument("--frag_min_len", action="store", dest="frag_min_len", type=int, default=5,
                    help="Minimum fragment length. Default is 5.")
parser.add_argument("--frag_max_len", action="store", dest="frag_max_len", type=int, default=15,
                    help="Maximum fragment length. Default is 15.")
parser.add_argument("--frag_mode", action="store", dest="frag_mode", choices=["keep", "split"], default="keep",
                    help="Fragment mode. 'keep' keeps only parent peptides within length range. 'split' generates sliding-window fragments.")
parser.add_argument("--frag_step", action="store", dest="frag_step", type=int, default=1,
                    help="Sliding window step size for split mode. Default is 1.")
parser.add_argument("--frag_best_per_parent", action="store_true", dest="frag_best_per_parent",
                    help="Keep only the best-ranked fragment per parent peptide.")
parser.add_argument("--frag_top_n", action="store", dest="frag_top_n", type=int, default=None,
                    help="Keep top N ranked fragments per parent peptide in split mode.")
parser.add_argument("--frag_output_name", action="store", dest="frag_output_name", default=None,
                    help="Optional output file name for fragment filter result, e.g. fragment_candidates.tsv")

# physicochemical filter patch
parser.add_argument("--physicochemical_filter", action="store_true", dest="physicochemical_filter",
                    help="Run optional physicochemical filtering on fragment-filter output.")
parser.add_argument("--phys_min_charge", action="store", dest="phys_min_charge", type=float, default=0.0,
                    help="Minimum net charge threshold at pH 7.4. Default is > 0.")
parser.add_argument("--phys_max_gravy", action="store", dest="phys_max_gravy", type=float, default=0.0,
                    help="Maximum GRAVY threshold. Default is < 0.")
parser.add_argument("--phys_min_pi", action="store", dest="phys_min_pi", type=float, default=None,
                    help="Minimum pI threshold.")
parser.add_argument("--phys_max_pi", action="store", dest="phys_max_pi", type=float, default=None,
                    help="Maximum pI threshold.")
parser.add_argument("--phys_max_cys", action="store", dest="phys_max_cys", type=int, default=None,
                    help="Maximum allowed cysteine count.")
parser.add_argument("--phys_exclude_high_disulfide_risk", action="store_true", dest="phys_exclude_high_disulfide_risk",
                    help="Exclude peptides with higher/high disulfide risk.")
parser.add_argument("--phys_output_mode", action="store", dest="phys_output_mode",
                    choices=["single", "layered"], default="single",
                    help="Output mode: 'single' (one combined file, default) or 'layered' (separate files per filter step).")

args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
var_threshold = args.var_threshold
pred = args.pred
add_seq_len = args.add_seq_len
esm_dir = args.esm_dir
top_cands = args.top_cands
rolling_window_size = args.rolling_window_size
plot_linear_epitope_scores = args.plot_linear_epitope_scores
zip_results = args.zip_results

fragment_filter = args.fragment_filter
frag_input_table = args.frag_input_table
frag_min_len = args.frag_min_len
frag_max_len = args.frag_max_len
frag_mode = args.frag_mode
frag_step = args.frag_step
frag_best_per_parent = args.frag_best_per_parent
frag_top_n = args.frag_top_n
frag_output_name = args.frag_output_name

physicochemical_filter = args.physicochemical_filter
phys_min_charge = args.phys_min_charge
phys_max_gravy = args.phys_max_gravy
phys_min_pi = args.phys_min_pi
phys_max_pi = args.phys_max_pi
phys_max_cys = args.phys_max_cys
phys_exclude_high_disulfide_risk = args.phys_exclude_high_disulfide_risk
phys_output_mode = args.phys_output_mode

### FUNCTIONS ###
def zip_function(result_files, outfile):
    zipf = zipfile.ZipFile(outfile, 'w')
    for result_file in result_files:
        zipf.write(result_file, arcname=result_file.name)
    zipf.close()

### MAIN ###

## Load antigen input and create ESM-2 encodings ##

# if you have the esm2 model stored locally, you can use this command.
# To work you need both esm2_t33_650M_UR50D.pt and the esm2_t33_650M_UR50D-contact-regression.pt
# stored in same directory.
# MyAntigens = bepipred3.Antigens(
#     fasta_file,
#     esm_dir,
#     add_seq_len=add_seq_len,
#     run_esm_model_local=str(WORK_DIR / "models" / "esm2_t33_650M_UR50D.pt")
# )

MyAntigens = bepipred3.Antigens(fasta_file, esm_dir, add_seq_len=add_seq_len)
MyBP3EnsemblePredict = bepipred3.BP3EnsemblePredict(
    MyAntigens,
    rolling_window_size=rolling_window_size,
    top_pred_pct=top_cands
)
MyBP3EnsemblePredict.run_bp3_ensemble()

out_dir.mkdir(parents=True, exist_ok=True)

MyBP3EnsemblePredict.create_toppct_files(out_dir)
MyBP3EnsemblePredict.create_csvfile(out_dir)

fragment_output_path = None

# fragment filter patch
if fragment_filter:
    try:
        fragment_output_path = run_fragment_filter(
            out_dir=out_dir,
            min_len=frag_min_len,
            max_len=frag_max_len,
            mode=frag_mode,
            step=frag_step,
            best_only=frag_best_per_parent,
            top_n=frag_top_n,
            input_table=frag_input_table,
            output_name=frag_output_name,
        )
    except Exception as e:
        print(f"[fragment_filter] WARNING: fragment filter failed: {e}")

# physicochemical filter patch
if physicochemical_filter:
    try:
        if fragment_output_path is None:
            raise ValueError("physicochemical_filter requires fragment_filter output. Please enable --fragment_filter first.")

        if phys_output_mode == "layered":
            stem = fragment_output_path.stem
            run_physicochemical_filter_layered(
                input_table=fragment_output_path,
                output_dir=out_dir,
                stem=stem,
                sequence_col="fragment",
                min_charge=phys_min_charge,
                max_gravy=phys_max_gravy,
                min_pi=phys_min_pi,
                max_pi=phys_max_pi,
                max_cys=phys_max_cys,
                exclude_high_disulfide_risk=phys_exclude_high_disulfide_risk,
            )
        else:
            phys_output_name = fragment_output_path.stem + "_physicochemical.tsv"
            phys_output_path = out_dir / phys_output_name

            run_physicochemical_filter(
                input_table=fragment_output_path,
                output_table=phys_output_path,
                sequence_col="fragment",
                length_col="fragment_length",
                start_col="fragment_start",
                end_col="fragment_end",
                min_charge=phys_min_charge,
                max_gravy=phys_max_gravy,
                min_pi=phys_min_pi,
                max_pi=phys_max_pi,
                max_cys=phys_max_cys,
                exclude_high_disulfide_risk=phys_exclude_high_disulfide_risk,
            )
    except Exception as e:
        print(f"[physicochemical_filter] WARNING: physicochemical filter failed: {e}")

## B-cell epitope predictions ##
if pred == "mjv_pred":
    MyBP3EnsemblePredict.bp3_pred_majority_vote(out_dir)
elif pred == "vt_pred":
    MyBP3EnsemblePredict.bp3_pred_variable_threshold(out_dir, var_threshold=var_threshold)

# generate plots (generating graphs for a maximum of 40 proteins)
MyBP3EnsemblePredict.bp3_generate_plots(
    out_dir,
    num_interactive_figs=50,
    use_rolling_mean=plot_linear_epitope_scores
)

# zip results
if zip_results:
    print("Zipping results")
    result_files = [f for f in out_dir.glob("*") if f.suffix != ".html"]
    zip_function(result_files, out_dir / "bepipred3_results.zip")
    for result_file in result_files:
        result_file.unlink()