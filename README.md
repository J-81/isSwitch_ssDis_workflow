## Snakemake based workflow for isSwitch Assignment

---

## isSwitch Description
- isSwitch describes the experimentally known variance in secondary structure across highly sequence conserved proteins.
- These variations may play a role in describing allosteric mechanisms and intrinsically disordered regions.

---

## Environment setup

- A user should have Conda installed
- Tested with Conda version 4.8.3
- Tested on Ubuntu 18.04

> conda env create -f envs/main.yml
> conda activate snakemake_workflows

---
## Usage (local machine only)

> snakemake --use-conda
