# Counts of Template-Primer Combinations

This script uses data from the study *"Quantitating primer-template interactions using deconstructed PCR"* to generate a counts file of template-primer combinations for statistical analysis.

---

### Requirements
Download the following files from the repository and upload them to your home directory:
- `Counts_Script.py`
- `metadata.csv`
- `maps.csv`
- `run.csv`

---

### Usage

**Step 1:** Create a text file listing your desired templates (space-delimited) from *Table 1: Experimental pools of primers and templates* in the study. An example file `templates_exampleB27.txt` is provided in the repository.

**Step 2:** Repeat the process for primers.

**Step 3:** Run the script, replacing `--exp`, `--t`, and `--p` with your experiment name and file names:
```bash
python3 Counts_Script.py --run run.csv --meta metadata.csv --exp B27 --map map.csv --t templates.txt --p primers.txt
```

---

### Output
A file named `final.csv` will be saved to your current working directory containing the template-primer counts.
