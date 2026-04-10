# Counts of Template-Primer Combinations
  # CURRENTLY NOT WORKING PROPERLY 

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
| Rpl | Run | Primer_Template_Combo | 5' | Middle | 3' | Annealing Temp | PCR Method | Count |
|-----|-----|-----------------------|----|--------|----|----------------|------------|-------|
| 1 | SRR8399877 | ST09V34 | 0 | 1 | 0 | 55C | TAS | 50 |
| 1 | SRR8399795 | ST09V30 | 1 | 0 | 1 | 45C | DePCR | 100 |
| . | . | . | . | . | . | . | . | . |
| . | . | . | . | . | . | . | . | . |
| 8 | SRR8399793 | ST09V33 | 0 | 0 | 0 | 55C | DePCR | 0 |
