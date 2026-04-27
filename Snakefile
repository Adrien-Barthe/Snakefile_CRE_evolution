import os
import glob
import re
import itertools

# Rules executing on the head node
localrules: all, create_models, merge_mutations_counts, merge_m3_ave_time_fix, create_R_scripts

# ==============================================================================
# CURRENT SIMULATION PARAMETERS
# Defines the exact models that SLiM needs to run during this execution.
# ==============================================================================
DOMINANCES = ["0.1", "0.7"] # Dominance coeff
SELECTIONS = ["-0.001", "-0.01", "-0.0425", "-0.1", "-0.3"] # Selection coeff for coding mutations (sM2)
SELECTIONS_M3 = ["-0.0001"] # Selection coeff for regulatory mutations (sM3)
PD_L = ["0.5"] # Physical distance (kb)
CRE_L = ["1"] # CRE size (kb)

# Fixed parameters for all models
# You can modify the number of replicates for a test, but delete manually the models data afterward ( in SLIM/ and 05-results/)
REPLICATS = range(1, 1001) # Generates 1000 simulations per model
GENERATION = 10000 # Number of simulated generations before data collection

# ==============================================================================
# PLOT LOCKING SYSTEM (TARGET PLOT PARAMETERS)
# Defines the combinations data that needs to be present to run the PDFs rules.
# ==============================================================================
TARGET_H = ["0.1", "0.25", "0.5", "0.7"]  #  All possible values for h
TARGET_REC = ["0.5", "1", "5", "10", "500"] #  All possible values for pd
TARGET_SM3 = ["-0", "-0.0001", "-0.001"] #  All possible values for sM3
TARGET_S = ["-0.001", "-0.01", "-0.0425", "-0.1", "-0.3"] #  All possible values for sM2

# String consisting concatenated pd values ( is used in fixedMutAcrossRECOMB script )
PD_STR_COMMA = ",".join(TARGET_REC)
PD_STR_USCORE = "_".join(TARGET_REC)

# ==============================================================================
# DETECTION OF EXISTING MODELS
# Search the results folder to see what models data has already been computed.
# ==============================================================================
# Generates the combinations that will be computed now
CURRENT_COMBOS = set(itertools.product(DOMINANCES, PD_L, CRE_L, SELECTIONS_M3, SELECTIONS))

# Looks for combinations already completed on the disk
EXISTING_COMBOS = set()
if os.path.exists("05-results"):
    for folder in glob.glob("05-results/fyonetal-h*_rec*"):
        # Value extraction via regular expressions
        m = re.search(r'fyonetal-h([0-9.]+)_rec([0-9.]+)kb_cre([0-9.]+)kb_sM3([0-9.-]+)_s([0-9.-]+)', folder)
        if m:
            EXISTING_COMBOS.add(m.groups()) # Adds a tuple (h, rec, cre, sm3, s)

# Merges what exists with what will be done to unlock the plotting rules
ALL_AVAILABLE = CURRENT_COMBOS.union(EXISTING_COMBOS)

# Priority order: always compute epistasis before the null model
ruleorder: run_slim_epist_all > run_slim_null_all

# ==============================================================================
# Determines exactly which files Snakemake must produce at the end.
# ==============================================================================
# A. Base targets: Simulation files and raw statistics
ALL_TARGETS = [
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist_g{g}_{mut}_{type}Mut.out", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, g=GENERATION, cl=CRE_L, pl=PD_L, mut=["m1", "m2", "m3"], type=["s", "f"]),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_g{g}_{mut}_{type}Mut.out", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, g=GENERATION, cl=CRE_L, pl=PD_L, mut=["m1", "m2", "m3"], type=["s", "f"]),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist_g{g}_m3_aveTimeFix.out", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, g=GENERATION, cl=CRE_L, pl=PD_L),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_g{g}_m3_aveTimeFix.out", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, g=GENERATION, cl=CRE_L, pl=PD_L),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist/pi_fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist.txt", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, cl=CRE_L, pl=PD_L),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist/sumStats_fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist_haps.txt", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, cl=CRE_L, pl=PD_L),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}/pi_fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}.txt", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, cl=CRE_L, pl=PD_L),
    expand("05-results/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}/sumStats_fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_haps.txt", h=DOMINANCES, s=SELECTIONS, sm3=SELECTIONS_M3, cl=CRE_L, pl=PD_L)
]

# B. Dynamic targets: Adding PDFs only if the required data is available
# Lock for H plot
for pl in PD_L:
    for cl in CRE_L:
        for sm3 in SELECTIONS_M3:
            for s in SELECTIONS:
                if all((h_val, pl, cl, sm3, s) in ALL_AVAILABLE for h_val in TARGET_H):
                    ALL_TARGETS.append(f"05-results/pdfs_h_analysis/SinglePlot_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}.pdf")

# Lock for Recombination plot
for h in DOMINANCES:
    for cl in CRE_L:
        for sm3 in SELECTIONS_M3:
            if all((h, pl_val, cl, sm3, s_val) in ALL_AVAILABLE for pl_val in TARGET_REC for s_val in TARGET_S):
                ALL_TARGETS.append(f"05-results/pdfs_rec_analysis/fyonetal-h{h}_cre{cl}kb_sM3{sm3}_PD_{PD_STR_USCORE}_acrossRecomb.done")

# Lock for SM3 plot
for h in DOMINANCES:
    for pl in PD_L:
        for cl in CRE_L:
            for s in SELECTIONS:
                if all((h, pl, cl, sm3_val, s) in ALL_AVAILABLE for sm3_val in TARGET_SM3):
                    ALL_TARGETS.append(f"05-results/pdfs_sm3_analysis/SinglePlot_h{h}_rec{pl}kb_cre{cl}kb_s{s}.pdf")

# Lock for SM2 global plots
for h in DOMINANCES:
    for pl in PD_L:
        for cl in CRE_L:
            for sm3 in SELECTIONS_M3:
                if all((h, pl, cl, sm3, s_val) in ALL_AVAILABLE for s_val in TARGET_S):
                    ALL_TARGETS.append(f"05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_SFSstats_violins.pdf")
                    ALL_TARGETS.append(f"05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_HapsPlots_violins.pdf")
                    ALL_TARGETS.append(f"05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_FixedStats_violins.pdf")

rule all:
    input:
        ALL_TARGETS


# --- Generate SLiM models (Eidos scripts) ---
# It creates two files per parameter combination: one with epistasis (ASE) and one without (null).
rule create_models:
    output:
        epist = "Models/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}_epist.txt",
        no_epist = "Models/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}.txt",
    params:
        fin_sim = int(GENERATION) + 1 # Calculates the exact generation to end the simulation
    shell:
        """
        mkdir -p Models

        # model epist
        cat << 'EOF' > {output.epist}
// set up a simple simulation with neutral and slightly deleterious variants
initialize()
{{
    // define global constants 
    defineConstant("GL", 50000);
    defineConstant("PD", asInteger({wildcards.pl})*1000);
    defineConstant("CL", asInteger({wildcards.cl})*1000);
    defineConstant("LT", GL + PD + CL);
    defineConstant("dASE", 0.7);
    defineConstant("domM2", {wildcards.h});

    initializeMutationRate(1.25e-7);
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", domM2, "f", {wildcards.s});
    
    initializeMutationType("m3", 0.5, "f", {wildcards.sm3});

    initializeGenomicElementType("g1", c(m1,m2), c(1,2));
    initializeGenomicElementType("g2", m3, 1);
    initializeGenomicElement(g1, 0, GL - 1);
    initializeGenomicElement(g2, GL + PD, LT - 1);
    initializeRecombinationRate(1.25e-7);
}}

1 early() {{
    defineConstant("simID", getSeed());
    sim.addSubpop("p1", 1000);
    hOne = (1 - dASE) ^ -(log(domM2)/log(2));
    hTwo = dASE ^ -(log(domM2)/log(2));
    sim.setValue("h1", hOne);
    sim.setValue("h2", hTwo);
}}

{GENERATION} late() {{ sim.outputFull(); }}
{GENERATION} late() {{ sim.outputFixedMutations(); }}
{params.fin_sim} early() {{ sim.simulationFinished(); }}

mutationEffect(m2) {{
    mutID = mut.id;
    mut2G1 = individual.genome1.mutationsOfType(m2);
    mut2G2 = individual.genome2.mutationsOfType(m2);
    countSm3G1 = individual.genome1.countOfMutationsOfType(m3);
    countSm3G2 = individual.genome2.countOfMutationsOfType(m3);
    if (homozygous)
        return 1.0 + mut.selectionCoeff;
    else
        if ((any(mut2G1.id==mutID)) & (countSm3G1 > countSm3G2)) 
            return 1.0 + sim.getValue("h2") * mut.selectionCoeff;
        else if (any(mut2G1.id==mutID) & (countSm3G1 < countSm3G2)) 
            return 1.0 + sim.getValue("h1") * mut.selectionCoeff;
        else if (any(mut2G2.id==mutID) & (countSm3G1 < countSm3G2)) 
            return 1.0 + sim.getValue("h2") * mut.selectionCoeff;
        else if (any(mut2G2.id==mutID) & (countSm3G1 > countSm3G2)) 
            return 1.0 + sim.getValue("h1") * mut.selectionCoeff;
        else
            return 1.0 + (mut.mutationType.dominanceCoeff * mut.selectionCoeff);
}}
EOF

        # null model
        cat << 'EOF' > {output.no_epist}
// set up a simple simulation with neutral and slightly deleterious variants
initialize()
{{
    defineConstant("GL", 50000 );
    defineConstant("PD", asInteger({wildcards.pl})*1000);
    defineConstant("CL", asInteger({wildcards.cl})*1000);
    defineConstant("LT", GL + PD + CL);
    defineConstant("dASE", 0.7);
    defineConstant("domM2", {wildcards.h});

    initializeMutationRate(1.25e-7);
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", domM2, "f", {wildcards.s});
    initializeMutationType("m3", 0.5, "f", {wildcards.sm3});

    initializeGenomicElementType("g1", c(m1,m2), c(1,2));
    initializeGenomicElementType("g2", m3, 1);
    initializeGenomicElement(g1, 0, GL - 1);
    initializeGenomicElement(g2, GL + PD, LT - 1);
    initializeRecombinationRate(1.25e-7);
}}

1 early() {{
    defineConstant("simID", getSeed());
    sim.addSubpop("p1", 1000);
}}

{GENERATION} late() {{ sim.outputFull(); }}
{GENERATION} late() {{ sim.outputFixedMutations(); }}
{params.fin_sim} early() {{ sim.simulationFinished(); }}
EOF
        """



# --- Run Epistatic SLiM Simulations ---
# Executes the SLiM simulations for the epistatic models in parallel using xargs
rule run_slim_epist_all:
    input:
        model = "Models/{model_name}.txt"
    output:
        # checking a single file as output is faster than 1000, which make the DAG faster
        done_flag = "SLIM/{model_name}/.slim_simulations.done"
    threads: 20
    resources:
        slurm_partition = "fast",
        mem_mb = 10000,
    wildcard_constraints:
        # Ensures this rule only applies to models containing "epist" in their name
        model_name = ".*epist.*" 
    priority: 50
    params:
        # Converts the range of replicates into a space-separated string for Bash
        reps = " ".join(map(str, REPLICATS)),
        gen = GENERATION
    shell:
        """
        # Create directories, in theory snakemake does it by itself
        mkdir -p SLIM/{wildcards.model_name}
        mkdir -p SLIM/{wildcards.model_name}/counts
        mkdir -p "05-results/{wildcards.model_name}"

        # Define the path for the temporary execution script
        SCRIPT="SLIM/{wildcards.model_name}/run_model_rep.sh"

        cat << 'EOF' > $SCRIPT
#!/bin/bash
REP=$1 # Replicate ID passed by xargs

# Define output filenames dynamically based on the replicate ID and model name
OUT_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}.out"
SEG_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.segMutations"
FIX_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.fixedMutations"
GEN_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.genomes"
PREFIX1="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_step1_"
PREFIX2="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_step2_"

# Run the SLiM simulation
slim {input.model} > $OUT_FILE

# Parse the output using csplit to extract Segregating and Fixed mutations
csplit -z --prefix $PREFIX1 $OUT_FILE /'#OUT: {params.gen} {params.gen}'/ '{{*}}'
mv ${{PREFIX1}}01 ${{SEG_FILE}}_temp
mv ${{PREFIX1}}02 $FIX_FILE
rm -f ${{PREFIX1}}*
        
# Further split the temporary file to isolate the Genomes block
csplit -z --prefix $PREFIX2 ${{SEG_FILE}}_temp /Genomes/ '{{*}}'
mv ${{PREFIX2}}00 $SEG_FILE
mv ${{PREFIX2}}01 $GEN_FILE
rm -f ${{PREFIX2}}* ${{SEG_FILE}}_temp

# Count the number of mutations (m1, m2, m3) and save to .count files
# (The '|| true' prevents grep from returning an error code if 0 mutations are found)
(grep "m1" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m1_s.count
(grep "m1" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m1_f.count
(grep "m2" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m2_s.count
(grep "m2" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m2_f.count
(grep "m3" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_s.count
(grep "m3" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_f.count

# Calculate the average time to fixation specifically for regulatory mutations (m3)
m3fixed=$( (grep "m3" $FIX_FILE || true) | wc -l )
if [ "$m3fixed" -gt 0 ]; then
    # Extracts columns 8 and 9 (generations), subtracts them, and computes the average via awk
    (grep "m3" $FIX_FILE  || true) | awk '{{print $9-$8}}' | awk '{{ sum += $1; n++ }} END {{ if (n > 0) print sum / n; }}' > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_repAveTimeFix.out
else
    # If no m3 fixed, output 0 
    echo "0" > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_repAveTimeFix.out
fi
EOF
        # ==========================================
        
        # Make the script executable
        chmod +x $SCRIPT
        
        # Execute the script in parallel using xargs
        # -n 1: Passes one replicate ID per script execution
        # -P {threads}: Runs up to the defined number of threads simultaneously
        echo "{params.reps}" | xargs -n 1 -P {threads} $SCRIPT
        
        # Clean up the temporary script and mark the rule as successfully completed
        rm -f $SCRIPT
        touch {output.done_flag}
        """





# --- Run Null SLiM Simulations ---
# Unlike the epistatic rule, this rule processes replicates sequentially 
# using a standard Bash 'for' loop and only requests 1 thread.
rule run_slim_null_all:
    input:
        model = "Models/{model_name}.txt"
    output:
        # same than in epistasis model
        done_flag = "SLIM/{model_name}/.slim_simulations.done"
    priority: 50
    threads: 1
    resources:
        slurm_partition = "fast",
        mem_mb = 1000,
    params:
        # Converts the range of replicates into a space-separated string for Bash
        reps = " ".join(map(str, REPLICATS)),
        gen = GENERATION
    shell:
        """
        # Create directories, in theory snakemake does it by itself
        mkdir -p SLIM/{wildcards.model_name}
        mkdir -p SLIM/{wildcards.model_name}/counts
        mkdir -p "05-results/{wildcards.model_name}"

        # Iterate sequentially over all replicates
        for REP in {params.reps}; do
          
            # Define output filenames dynamically based on the replicate ID
            OUT_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}.out"
            SEG_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.segMutations"
            FIX_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.fixedMutations"
            GEN_FILE="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_g{params.gen}.genomes"
            PREFIX1="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_step1_"
            PREFIX2="SLIM/{wildcards.model_name}/{wildcards.model_name}_${{REP}}_step2_"

            # Run the SLiM simulation
            slim {input.model} > $OUT_FILE

            # Parse the output using csplit to extract Segregating and Fixed mutations
            csplit -z --prefix $PREFIX1 $OUT_FILE /'#OUT: {params.gen} {params.gen}'/ '{{*}}'
            mv ${{PREFIX1}}01 ${{SEG_FILE}}_temp
            mv ${{PREFIX1}}02 $FIX_FILE
            rm -f ${{PREFIX1}}*
            
            # Further split the temporary file to isolate the Genomes block
            csplit -z --prefix $PREFIX2 ${{SEG_FILE}}_temp /Genomes/ '{{*}}'
            mv ${{PREFIX2}}00 $SEG_FILE
            mv ${{PREFIX2}}01 $GEN_FILE
            rm -f ${{PREFIX2}}* ${{SEG_FILE}}_temp

            # Count the number of mutations (m1, m2, m3) and save to .count files
            # (The '|| true' prevents grep from crashing the script if 0 mutations are found)
            (grep "m1" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m1_s.count
            (grep "m1" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m1_f.count
            (grep "m2" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m2_s.count
            (grep "m2" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m2_f.count
            (grep "m3" $SEG_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_s.count
            (grep "m3" $FIX_FILE || true) | wc -l > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_f.count

            # Calculate the average time to fixation specifically for regulatory mutations (m3)
            m3fixed=$( (grep "m3" $FIX_FILE || true) | wc -l )
            if [ "$m3fixed" -gt 0 ]; then
                # Extracts columns 8 and 9 (generations), subtracts them, and computes the average via awk
                (grep "m3" $FIX_FILE  || true) | awk '{{print $9-$8}}' | awk '{{ sum += $1; n++ }} END {{ if (n > 0) print sum / n; }}' > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_repAveTimeFix.out
            else
                # If no m3 fixed, output 0
                echo "0" > SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_${{REP}}_m3_repAveTimeFix.out
            fi

        done
        
        # Mark the rule as successfully completed
        touch {output.done_flag}
        """



# --- Merge Mutation Counts ---
# Aggregates the individual '.count' files (generated per replicate) into a 
# single summary '.out' file for a specific mutation type (m1, m2, or m3) 
# and state (segregating 's' or fixed 'f').
rule merge_mutations_counts:
    input:
        # Waits for the flag ensuring all SLiM replicates for this model are complete
        slim_done = "SLIM/{model_name}/.slim_simulations.done",
    output:
        # The final concatenated summary file placed in the 05-results/ folder
        count_summary = "05-results/{model_name}/{model_name}_g"+ str(GENERATION) +"_{mut}_{type}Mut.out",
    shell:
        """
        # Concatenates all individual replicate counts for the specified mutation type and state
        cat SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_*_{wildcards.mut}_{wildcards.type}.count > {output.count_summary}
        """

# --- Merge M3 Average Time to Fixation ---
# Aggregates the individual files containing the average time to fixation 
# for regulatory mutations (m3) across all replicates into a single summary file.
rule merge_m3_ave_time_fix:
    input:
        # Waits for the flag ensuring all SLiM replicates for this model are complete
        slim_done = "SLIM/{model_name}/.slim_simulations.done",
    output:
        # The final concatenated summary file placed in the 05-results/ folder
        m3_fix_time_summary = "05-results/{model_name}/{model_name}_g"+ str(GENERATION) +"_m3_aveTimeFix.out"
    shell:
        """
        # Concatenates all individual replicate fixation times for m3
        cat SLIM/{wildcards.model_name}/counts/{wildcards.model_name}_*_m3_repAveTimeFix.out > {output.m3_fix_time_summary}
        """



# --- Compute Summary Statistics ---
# Uses an R script to calculate genetic statistics (pi, per site stats, and haplotypes) 
# from the raw SLiM outputs. 
# Implements a dynamic batching system with 'xargs' to process replicates in parallel 
rule compute_sum_stats:
    input:
        # Ensures all SLiM simulations are fully completed before starting the R analysis
        slim_done = "SLIM/{model_name}/.slim_simulations.done",
        script = "02-scripts/computingSumStats_SLIM.R",
        funcs = "02-scripts/sumStats_functions.R"
    output:
        # Final aggregated statistics files
        pi = "05-results/{model_name}/pi_{model_name}.txt",
        per_site = "05-results/{model_name}/sumStats_{model_name}_perSite.txt",
        haps = "05-results/{model_name}/sumStats_{model_name}_haps.txt"
    params:
        max_rep = max(REPLICATS),
        gen = GENERATION
    threads: 10
    resources:
        mem_mb = 5000,   
        slurm_partition = "fast"
    shell:
        """
        # Clean up any leftover temporary batch files from previous interrupted runs
        rm -f SLIM/{wildcards.model_name}/pi_{wildcards.model_name}_*.txt
        rm -f SLIM/{wildcards.model_name}/sumStats_{wildcards.model_name}_*.txt
        
        MAX_REP={params.max_rep}
        THREADS={threads}
        
        # Dynamic Batch Size Calculation
        # Divides the total number of replicates by the number of allocated threads.
        # (The '+ THREADS - 1' ensures the division rounds up to cover all replicates).
        BATCH_SIZE=$(( (MAX_REP + THREADS - 1) / THREADS ))
        
        # Generate Start-End Pairs for each batch
        PAIRS=""
        for (( start=1; start<=MAX_REP; start+=BATCH_SIZE )); do
            end=$(( start + BATCH_SIZE - 1 ))
            if [ $end -gt $MAX_REP ]; then end=$MAX_REP; fi
            PAIRS="$PAIRS $start $end"
        done
        
        # Create the execution wrapper for R
        SCRIPT="SLIM/{wildcards.model_name}/run_R_batch.sh"
        cat << 'EOF' > $SCRIPT
#!/bin/bash
START=$1
END=$2
# Calls the R script with the specific start and end replicate indices
Rscript --vanilla {input.script} SLIM {wildcards.model_name} $START $END {params.gen}
EOF
        chmod +x $SCRIPT
   
        # Parallel Execution via xargs
        # -n 2: Takes exactly two arguments from PAIRS (start and end) per job
        # -P {threads}: Runs up to the maximum allocated threads simultaneously
        echo $PAIRS | xargs -n 2 -P {threads} $SCRIPT
        rm -f $SCRIPT

        # ==========================================
        # 6. Stitching the Batch Outputs Together
        # Iterates through the generated chunks, keeps the header from the first file,
        # and appends only the data rows from subsequent files.
        # ==========================================
        
        # Merge Pi statistics
        first=true
        for f in SLIM/{wildcards.model_name}/pi_{wildcards.model_name}_*.txt; do
            if [ "$first" = true ]; then
                head -n 1 "$f" > {output.pi} # Keep header
                first=false
            fi
            tail -n +2 "$f" >> {output.pi}   # Append data without header
        done

        # Merge Per-Site statistics
        first=true
        for f in SLIM/{wildcards.model_name}/sumStats_{wildcards.model_name}_perSite_*.txt; do
            if [ "$first" = true ]; then
                head -n 1 "$f" > {output.per_site}
                first=false
            fi
            tail -n +2 "$f" >> {output.per_site}
        done

        # Merge Haplotype statistics
        first=true
        for f in SLIM/{wildcards.model_name}/sumStats_{wildcards.model_name}_haps_*.txt; do
            if [ "$first" = true ]; then
                head -n 1 "$f" > {output.haps}
                first=false
            fi
            tail -n +2 "$f" >> {output.haps}
        done

        # 7. Final Cleanup: remove the temporary chunk files
        rm -f SLIM/{wildcards.model_name}/*_100.txt
        """

# --- Generate SFS, Haplotypes, and Fixed Mutations Plots ---
# Creates PDF box plots comparing different selection coefficients (sM2) for different statistics.
rule plot_sum_stats:
    input:
        script_sfs = "02-scripts/plotting_SFS_Stats.R",
        script_hap = "02-scripts/plottingHapStats.R",
        script_fixed = "02-scripts/plotting_Fix_Stats.R",
        haps_epist = expand("05-results/fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_epist/sumStats_fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_epist_haps.txt", s=TARGET_S),
        haps_null = expand("05-results/fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}/sumStats_fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_haps.txt", s=TARGET_S)
    output:
        pdf_sfs = "05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_SFSstats_violins.pdf",
        pdf_hap = "05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_HapsPlots_violins.pdf",
        pdf_fix = "05-results/pdfs_sm2_analysis/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_FixedStats_violins.pdf",
        mean_hap = "05-results/sumStats/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_mean_hapStats.txt",
        median_hap = "05-results/sumStats/fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}_median_hapStats.txt"
    params:
        targetModel = "fyonetal-h{h}_rec{pl}kb_cre{cl}kb_sM3{sm3}",
        reps = max(REPLICATS)
    threads: 1
    resources:
        mem_mb = 1000,
        partition = "fast"
    shell:
        """

        Rscript --vanilla {input.script_sfs} $PWD/05-results {params.targetModel} {params.reps}
        Rscript --vanilla {input.script_hap} $PWD/05-results {params.targetModel} {params.reps}
        Rscript --vanilla {input.script_fixed} $PWD/05-results {params.targetModel} {params.reps}
        """

# --- Recombination Analysis Plots ---
# Analyzes and plots fixed mutations as a function of the physical distance (recombination rate).
# Uses a dummy '.done' file because the R script generates multiple PDFs in the output directory.
rule plot_fixed_mut_across_recomb:
    input:
        haps_epist = expand("05-results/fyonetal-h{{h}}_rec{pl}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_epist/sumStats_fyonetal-h{{h}}_rec{pl}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_epist_haps.txt", s=TARGET_S, pl=TARGET_REC),
        haps_null = expand("05-results/fyonetal-h{{h}}_rec{pl}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}/sumStats_fyonetal-h{{h}}_rec{pl}kb_cre{{cl}}kb_sM3{{sm3}}_s{s}_haps.txt", s=TARGET_S, pl=TARGET_REC),        
        script = "02-scripts/fixedMutAcrossRECOMB.R"
    output:
        dummy = "05-results/pdfs_rec_analysis/fyonetal-h{h}_cre{cl}kb_sM3{sm3}_PD_" + PD_STR_USCORE + "_acrossRecomb.done"
    params:
        suffix = "fyonetal-h{h}_rec",
        cre_size = "{cl}",
        pd_comma = PD_STR_COMMA,
        pd_uscore = PD_STR_USCORE,
        sm3 = "{sm3}"
    threads: 1
    resources:
        mem_mb = 1000,
        slurm_partition = "fast"
    shell:
        """
        Rscript --vanilla {input.script} $PWD/05-results {params.suffix} {params.cre_size} {params.pd_comma} {params.pd_uscore} {params.sm3}
        touch {output.dummy}
        """


# --- Dominance Coefficient (h) Analysis Plots ---
# Generates comparison plots for a specific model across different dominance coefficients (h values).
rule plot_h_stats:
    input:
        script = "02-scripts/plotting_h_Stats.R",
        haps_epist = expand("05-results/fyonetal-h{h}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{{s}}_epist/sumStats_fyonetal-h{h}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{{s}}_epist_haps.txt", h=TARGET_H),
        haps_null = expand("05-results/fyonetal-h{h}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{{s}}/sumStats_fyonetal-h{h}_rec{{pl}}kb_cre{{cl}}kb_sM3{{sm3}}_s{{s}}_haps.txt", h=TARGET_H)    
    output:
        pdf = "05-results/pdfs_h_analysis/SinglePlot_rec{pl}kb_cre{cl}kb_sM3{sm3}_s{s}.pdf"
    params:
        h_comma = ",".join(TARGET_H)
    threads: 1
    resources:
        mem_mb = 1000,
        slurm_partition = "fast"
    shell:
        """
        Rscript --vanilla {input.script} $PWD/05-results {wildcards.pl} {wildcards.cl} {wildcards.sm3} {wildcards.s} {params.h_comma}
        """

# --- Regulatory Selection (sM3) Analysis Plots ---
# Generates comparison plots for a specific model across different regulatory selection coefficients (sM3 values).
rule plot_sm3_stats:
    input:
        script = "02-scripts/plotting_sm3_Stats.R",
        haps_epist = expand("05-results/fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{sm3}_s{{s}}_epist/sumStats_fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{sm3}_s{{s}}_epist_haps.txt", sm3=TARGET_SM3),
        haps_null = expand("05-results/fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{sm3}_s{{s}}/sumStats_fyonetal-h{{h}}_rec{{pl}}kb_cre{{cl}}kb_sM3{sm3}_s{{s}}_haps.txt", sm3=TARGET_SM3)
    output:
        pdf = "05-results/pdfs_sm3_analysis/SinglePlot_h{h}_rec{pl}kb_cre{cl}kb_s{s}.pdf"
    params:
        sm3_comma = ",".join(TARGET_SM3)
    threads: 1
    resources:
        mem_mb = 1000,
        slurm_partition = "fast"
    shell:
        """
        Rscript --vanilla {input.script} $PWD/05-results {wildcards.pl} {wildcards.cl} {wildcards.h} {wildcards.s} {params.sm3_comma}
        """

# Create all the R scripts needed for prior rules
rule create_R_scripts:
    output:
        sum_stat = "02-scripts/computingSumStats_SLIM.R",
        sum_stat_f = "02-scripts/sumStats_functions.R",
        sfs = "02-scripts/plotting_SFS_Stats.R",
        hap = "02-scripts/plottingHapStats.R",
        fix_reco = "02-scripts/fixedMutAcrossRECOMB.R",
        per_h_stats = "02-scripts/plotting_h_Stats.R",
        fixed_mut = "02-scripts/plotting_Fix_Stats.R",
        per_sm3_stats = "02-scripts/plotting_sm3_Stats.R" 
    priority: 51
    shell:
        """

        cat << 'EOF' > {output.sum_stat}
#!/usr/bin/env Rscript
#############################
##
## This script takes the output of SLiM and builds a 
## genome image of the mutations: varPos on the columns, genomes on the rows
## it takes as input *.segMutations and *.genomes 
##
##
## by Isabel Alves - Fev 2023
##
#############################


#Version récente
args = commandArgs(trailingOnly=TRUE)
library(ggplot2) #R 4.3
library(tidyverse) 
source("02-scripts/sumStats_functions.R")

wrkDir <- args[1]
sceName <- args[2]
simStart <- as.numeric(args[3])
simEnd <- as.numeric(args[4]) 
setwd(paste0(wrkDir, "/", sceName))

# examples : c("twoLoci_rec_cre100kb_s0.001","twoLoci_rec_cre100kb_s0.001_epist", "twoLoci_rec_cre100kb_s0.0425", "twoLoci_rec_cre100kb_s0.0425_epist","twoLoci_rec_cre100kb_s0.1", "twoLoci_rec_cre100kb_s0.1_epist", "twoLoci_rec_cre100kb_s0.3", "twoLoci_rec_cre100kb_s0.3_epist")
#scenario name

genNb <-  as.numeric(args[5])# from 0 4999
popSize <- 1000 

# arguments to pass - TEST 
# sceName <- "twoLoci_rec10kb_cre5kb_s-0.0425"
# genNb <- 10000 # from 0 4999
# nbOfSims <- 100
# wrkDir <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/CRE_evolution/05-Results/sims_antoine_all/experience_haplotypes"
# popSize <- 1000
# simID <- 100
# setwd(wrkDir)
#test pi ans S computation
# set.seed(189989)
# seqTest <- sample(c("A","T", "G", "C"), 50, replace = T)
# df_seq <- data.frame(rbind(seqTest,seqTest,seqTest,seqTest,seqTest,seqTest,seqTest,
#                            seqTest,seqTest,seqTest))
# df_seq[c(2,5,6,9),4] <- "G"
# df_seq[c(7),17] <- "T"
# df_seq[c(2,4,5,6,7,8,9,10),44] <- "A"
# df_seq[c(4,8),50] <- "C"
# 
# test_uniqHaps <- df_seq %>% group_by_all %>% count
# test_HapFreq <- test_uniqHaps$n/sum(test_uniqHaps$n)
# pitest <- computing_Pi(test_uniqHaps[,-ncol(test_uniqHaps)], test_HapFreq, 10)
# segTest <- 4/sum(1/rep(1:(10-1)))


#creating output files with headers
sumStats_perSite_fName <- paste0("sumStats_", sceName, "_perSite_", simStart, "_", simEnd, ".txt")
sumStats_pi_fName <- paste0("pi_", sceName, "_", simStart, "_", simEnd, ".txt")
sumStats_haps_fName <- paste0("sumStats_", sceName, "_haps_", simStart, "_", simEnd, ".txt")

cat(paste("Nb_m1", "Nb_m2", "Nb_m3", "AveHet", "AveHet_m1", "AveHet_m2", "AveHet_m3", sep = "\\t"), sep = "\\n", file = sumStats_perSite_fName)
cat(paste("Overall_Nucleotide_diversity", "CS_Nucleotide_diversity", "CRE_Nucleotide_diversity", 
          "Overall_SegSites", "Overall_ThetaS", "CS_SegSites", "CS_ThetaS", 
          "CRE_SegSites", "CRE_ThetaS", "OverallTajimasD", "TajimasDCS", "TajimasDCRE", sep = "\\t"), sep = "\\n", file = sumStats_pi_fName)
cat(paste("avePerHap_nbMut2_ontoM3haps", "avePerHap_nbMut2_ontoNOM3haps", 
          "NbHaps_w_m3_w_m2", "NbHaps_w_m3_wo_m2", "NbHaps_wo_m3_w_m2", "NbHaps_wo_m3_wo_m2", "Odds_m3_m2","p-value_m3_m2",
          "avePerHap_nbMut1_ontoM3haps", "avePerHap_nbMut1_ontoNOM3haps",
          "NbHaps_w_m3_w_m1", "NbHaps_w_m3_wo_m1", "NbHaps_wo_m3_w_m1", "NbHaps_wo_m3_wo_m1", "Odds_m3_m1","p-value_m3_m1", "nbUniqHap", 
          "CorrM2M3acrossHaps", "pValCorrM2M3acrossHaps", 
          "meanNbm2ontoM3plus","meanNbm2ontoM3minus","medianNbm2ontoM3plus","medianNbm2ontoM3minus",
          sep = "\\t"), 
    sep = "\\n", file = sumStats_haps_fName)


for(simID in simStart:simEnd) {{

  #open the mutation section of SLIM's output
  #twoLoci_rec5kb_cre1kb_s-0.001_epist_54_g10000.segMutations 
  #setwd(paste0(wrkDir, "/sims/", sceName))
  #setwd("/Users/ialves/Dropbox/SLIM/sims/twoLoci_rec_cre1kb_s0.3/")
  #simID <- 98
  segSitesTable <- openSegSites(sceName, simID, genNb)
  #dim(segSitesTable)
  #checking whether there're variants present more than once in the mut table
  # These sites are excluded from all the following analysis 
  cleanSegSitesTable <- detectMultiallelicSites(segSitesTable, TRUE)
  #dim(cleanSegSitesTable)
  nbOfMultiallelicSites <- nrow(segSitesTable)-nrow(cleanSegSitesTable)
  rm(segSitesTable)
  
  #heterozygosity 
  hetTable <- computingHeterozygosity(cleanSegSitesTable, popSize)
  mut_tblName <- paste0(sceName, "_mutProperties_", simID, ".txt")
  
  #cat(paste("Nb_m1", "Nb_m2", "Nb_m3", "AveHet", "AveHet_m1", "AveHet_m2", "AveHet_m3", sep = "\\t"), sep = "\\n", file = sumStats_perSite_fName)
  cat(paste(sum(hetTable$muType == "m1"), sum(hetTable$muType == "m2"), sum(hetTable$muType == "m3"), mean(hetTable$Heterozyg),
            mean(hetTable$Heterozyg[hetTable$muType == "m1"]),
            mean(hetTable$Heterozyg[hetTable$muType == "m2"]),
            mean(hetTable$Heterozyg[hetTable$muType == "m3"]), sep = "\\t"), sep = "\\n", file = sumStats_perSite_fName, append = T)
  # reading the genomes 
  genMatrix <- getGenomesFromSLIMoutput(sceName, simID, genNb, popSize)
  #dim(genMatrix)
  df_genomes <- data.frame(genMatrix)
  #retireving unique haplotypes 
  uniq_haps <- df_genomes %>% group_by_all %>% count
  #computing haplotype frequencies
  haps_freq <- uniq_haps$n/sum(uniq_haps$n)
  #computing nb of m2 and m3 onto the same haplotypes
  # the idea is to check whether there si a negative correlation between the
  # nb of m2 and m3. 
  m2m3AcrossHapsDF <- data.frame(m2=apply(uniq_haps, 1, function(x) {{sum(x=="m2")}}), m3=apply(uniq_haps, 1, function(x) {{sum(x=="m3")}}))
  # m2m3AcrossHapsDF %>% ggplot(aes(x=m2, y=m3)) + geom_point() + theme_classic() + geom_smooth(method='lm', formula= y~x)
  cor_m2_m3_haps <- cor.test(m2m3AcrossHapsDF$m2, m2m3AcrossHapsDF$m3, method="spearman", continuity = T)$estimate
  corPval_m2_m3_haps <- cor.test(m2m3AcrossHapsDF$m2, m2m3AcrossHapsDF$m3, method="spearman", continuity = T)$p.value
  
  # getting the number of m2 onto the ind genomes with lowest (underexpressed) vs highest nb of m3 (overexpressed)
  kkk <- computing_Nbm2_Nbm3_acrossGenomes(df_genomes, popSize)
  # kkk %>% pivot_longer(col=1:2, names_to = "backgroundM3", values_to="m2") %>% ggplot(aes(x=backgroundM3, y=m2)) + geom_boxplot() + stat_compare_means(method="wilcox", paired=T) + 
  # xlab("Enhancer background") + ylab("Number of deleterious mutations")
  # mean and median below are across individuals - to print haps file
  meanNbm2 <- apply(kkk, 2, mean,na.rm=T) 
  medianNbm2 <- apply(kkk, 2, median,na.rm=T)
  
  #####--------------------
  #### Compute m3 background per site per individual
  #####--------------------
  # mut_tb will gain two new columns with the info
  # on how many ind are het at m2/m1 and the m2 or m1 allele is onto
  # the m3-- background
  posM2 <- sort(unlist(hetTable %>% filter(muType == "m2") %>% select(Pos)))
  posM1 <- sort(unlist(hetTable %>% filter(muType == "m1") %>% select(Pos)))

  #apply(gen_m[gen_m$IndID %in% indHet,], 1, function(x) {{ sum(x == "m3")}})
  hetTable <- hetTable %>% mutate(`propIndM2_m3-`= rep(NA, nrow(hetTable)))
  hetTable <- hetTable %>% mutate(`propIndM1_m3-`= rep(NA, nrow(hetTable)), ObsHet=rep(NA, nrow(hetTable)))
  
  
  tmpM2 <- countingMut1or2_m3minus_combinations_perSite("m2", popSize)
  if(length(tmpM2$df) == length(posM2)) {{
    names(tmpM2$df) <- posM2
  }}
  tmpM1 <- countingMut1or2_m3minus_combinations_perSite("m1", popSize)
  if(length(tmpM1$df) == length(posM1)) {{
    names(tmpM1$df) <- posM1
  }}
  hetTable$`propIndM2_m3-`[match(posM2, hetTable$Pos)] <- tmpM2$df
  hetTable$`propIndM1_m3-`[match(posM1, hetTable$Pos)] <- tmpM1$df
  hetTable$ObsHet[match(posM2, hetTable$Pos)] <- tmpM2$nbhetInd
  hetTable$ObsHet[match(posM1, hetTable$Pos)] <- tmpM1$nbhetInd
  write.table(hetTable, file=mut_tblName, quote = F, row.names = F, col.names = T, sep = "\\t")
  #####--------------------
  #### end
  #####--------------------
  
  
  #computing pi
  nucleotideDiv <- computing_Pi(uniq_haps[,-ncol(uniq_haps)], haps_freq, popSize)
  # computing Theta S
  thetaS <- nrow(cleanSegSitesTable)/sum(1/rep(1:((popSize*2)-1)))
  
  ### Within CRE and CS 
  #computing pi 
  nbofM3Mut <- nrow(cleanSegSitesTable %>% filter(V3 == "m3"))
  if(nbofM3Mut > 0 ) {{
    
    CREsites <- MutPositions("m3")
    CREgenoMatrix <- df_genomes %>% select(colnames(df_genomes)[colnames(df_genomes) %in% paste0("X",CREsites)])
    CSgenoMatrix <- df_genomes %>% select(colnames(df_genomes)[!colnames(df_genomes) %in% paste0("X",CREsites)])
    uniqHapsCRE <- CREgenoMatrix %>% group_by_all %>% count
    uniqHapsCS <- CSgenoMatrix %>% group_by_all %>% count
    hapsFreqCRE <- uniqHapsCRE$n/sum(uniqHapsCRE$n)
    hapsFreqCS <- uniqHapsCS$n/sum(uniqHapsCS$n)
    
    piCRE <- computing_Pi(uniqHapsCRE[,-ncol(uniqHapsCRE)], hapsFreqCRE, popSize)
    piCS <- computing_Pi(uniqHapsCS[,-ncol(uniqHapsCS)], hapsFreqCS, popSize)
    
  }} else {{
    CSgenoMatrix <- df_genomes 
    uniqHapsCS <- CSgenoMatrix %>% group_by_all %>% count
    hapsFreqCS <- uniqHapsCS$n/sum(uniqHapsCS$n)
    # popSize here is in nb of diploid individuals
    piCS <- computing_Pi(uniqHapsCS[,-ncol(uniqHapsCS)], hapsFreqCS, popSize)
    piCRE <- 0
    
  }}

  #computing S 
  CRESegSites <- nrow(cleanSegSitesTable %>% filter(V3 == "m3"))
  CREThetaSegSites <- CRESegSites/sum(1/rep(1:((popSize*2)-1)))
  CSSegSites <- nrow(cleanSegSitesTable %>% filter(V3 == "m1" | V3 == "m2"))
  CSThetaSegSites <- CSSegSites/sum(1/rep(1:((popSize*2)-1)))
  
  #computing Tajimas'D 
  tajDCRE <- piCRE - CREThetaSegSites
  tajDCS <- piCS - CSThetaSegSites
  ### ------------
  overTajimasD <- nucleotideDiv - thetaS
  cat(paste(nucleotideDiv, piCS, piCRE, nrow(cleanSegSitesTable), thetaS, CSSegSites, CSThetaSegSites,
            CRESegSites, CREThetaSegSites, overTajimasD, tajDCS, tajDCRE, sep = "\\t"), sep = "\\n", file = sumStats_pi_fName, append = T)
  
  #computing the number of mut of type X per haplotype
  nb_m3_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m3")
  nb_m2_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m2")
  nb_m1_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m1")
  
  #the vector below contains the ave nb of m2 (or m1) mut (>0) within haps containing or not m3 mut.
  m2_onm3 <- nb_mX_onto_m3_background(df_genomes, nb_m2_per_hap, nb_m3_per_hap, "m2")
  names(m2_onm3) <- c("mean_nb_mut2_ontoM3haps", "mean_nb_mut2_ontoNoM3haps", "nb_hap_w_m3_w_m2",
                      "nb_hap_w_m3_wo_m2","nb_hap_wo_m3_w_m2","nb_hap_wo_m3_wo_m2")
  m1_onm3 <- nb_mX_onto_m3_background(df_genomes, nb_m1_per_hap, nb_m3_per_hap, "m1")
  names(m1_onm3) <- c("mean_nb_mut1_ontoM3haps", "mean_nb_mut1_ontoNoM3haps", "nb_hap_w_m3_w_m1",
                      "nb_hap_w_m3_wo_m1","nb_hap_wo_m3_w_m1","nb_hap_wo_m3_wo_m1")
  
  #are there more haplotypes with m3 and m2 than m3 and no m2 ? 
  contingencyTable_m2_m3 <- matrix(m2_onm3[-c(1:2)], ncol = 2)
  # fisher.test(contingencyTable_m2_m3)$p.value
  # fisher.test(contingencyTable_m2_m3)$odd
  contingencyTable_m1_m3 <- matrix(m1_onm3[-c(1:2)], ncol = 2)
  cat(paste(paste(m2_onm3, collapse = "\\t"), fisher.test(contingencyTable_m2_m3)$estimate,fisher.test(contingencyTable_m2_m3)$p.value,
            paste(m1_onm3, collapse = "\\t"), fisher.test(contingencyTable_m1_m3)$estimate,fisher.test(contingencyTable_m1_m3)$p.value,
            dim(uniq_haps)[1], unlist(cor_m2_m3_haps), unlist(corPval_m2_m3_haps), paste(meanNbm2, collapse = "\\t"), paste(medianNbm2, collapse = "\\t"),
            sep = "\\t"), sep = "\\n", file = sumStats_haps_fName, append = T)
  # LAST PART: print nb uniq haps, corr nb m2 and nb m3, respective pvalue, mean nb m2 onto m3++, mean nb m2 onto m3--
  # median nb m2 onto m3++, median nb m2 onto m3--

}} 
  

  
  #####
  ##### to be finished - fev 7, 2023
  #####
  #m3 SFS entries - do a table of the variables below
  # hap_counts_per_m2 <- apply(df_genomes, 2, function(x) {{sum(x=="m2")}})
  # hap_counts_per_m2_condm3 <- apply(df_genomes[which(nb_haps_w_m3>0),], 2, function(x) {{sum(x=="m2")}})
  # plot(table(hap_counts_per_m2)/sum(table(hap_counts_per_m2)), type = "b")
  # lines(table(hap_counts_per_m2_condm3)/sum(table(hap_counts_per_m2_condm3)), type = "b")
  # 
  # hap_counts_m3 <- apply(df_genomes, 2, function(x) {{sum(x=="m3")}})
  # hap_counts_m1 <- apply(df_genomes, 2, function(x) {{sum(x=="m1")}})
    
    
EOF

        cat << 'EOF' > {output.sum_stat_f}

#############################
##
## Functions
##
#############################

#open the mutation table 
openSegSites <- function(k, y, z) {{
  # k = sceName
  # y = simID (replicate)
  # z = genNb
  x <- file(paste0(k, "_", y, "_g", z, ".segMutations"), "r")
  line <- readLines(x, 1)
  count <- 1
  while(!grepl("^Individuals:+", line)) {{
    line <- readLines(x, 1)
    count <- count + 1
  }}
  #! note that skip = 5 requires knowing that there are originally 5 lines of header in this file
  mut_tb <- read.table(paste0(k, "_", y, "_g", z, ".segMutations"),
                       skip = 5, nrows = count-6 , sep = " ")
  return(mut_tb)
}}
# handle multiallelic sites
detectMultiallelicSites <- function(segM, clean) {{
  # segM = segSitesTable
  #clean = T/F; if T then multiallelic sites will be removed
  #otherwise renamed by giving the physical position followed by X.1, X.2 etc
  reNameMut <- table(segM$V4)[table(segM$V4) > 1]
  if(length(reNameMut) > 0) {{
    if(clean) {{
      index_to_replace <- which(segM$V4 %in% as.numeric(names(reNameMut)))
      segM <- segM[-index_to_replace,]
    }} else {{
      for(i in 1:length(reNameMut)) {{
        index_to_replace <- which(segM$V4 == as.numeric(names(reNameMut[i])))
        segM$V4[index_to_replace] <- paste0(names(reNameMut[i]), ".", 1:reNameMut[i])
      }} #end of for along the multiallelic sites
    }} #end of the if clean T/F
  }} #end of the if there are multiallelic sites
  return(segM)
}} 

# computing heterozygosity and adding it to hte clean table
#computing heterozygosity  
computingHeterozygosity <- function(segTbl, ps) {{
  # segTbl = cleanSegSitesTable
  # ps = popSize
  
  p <- segTbl$V9/(ps*2) 
  q <- 1-p
  het <- 2*p*q
  newTbl <- segTbl %>% mutate(Heteroz = het)
  names(newTbl) <- c("SiteIDI", "SiteIDII", "muType", "Pos", "selCoeff", "domCoeff", "subPop",
                     "genOrigin", "absFreq", "Heterozyg")
  return(newTbl)
}}
###------------

# read the genome outputed by slim
getGenomesFromSLIMoutput <- function(k, y, z, ps) {{
  # k = sceName
  # y = simID (replicate)
  # z = genNb
  # ps = popSize
  
  gen_l <- list(ps*2)
  for(genome in 0:(ps*2-1)) {{
    
    #genome <- 0
    tmp_v <- as.numeric(scan(paste0(k, "_", y, "_g", z, ".genomes"), skip = 1+genome, nlines = 1, sep = " ",                                      what = character())[-c(1,2)])
    cleanTmp_v <- tmp_v[tmp_v %in% hetTable$SiteIDI]
    
    genomeAnc <- rep("A", nrow(hetTable))
    genomeAnc[hetTable$SiteIDI %in% cleanTmp_v] <- hetTable$muType[hetTable$SiteIDI %in% cleanTmp_v]
    gen_l[[genome+1]] <- genomeAnc
    rm(tmp_v, cleanTmp_v)
  }}
  
  # tranform list to matrix
  gen_m <- do.call(rbind,gen_l)
  gen_m <- gen_m[,order(hetTable$Pos, decreasing = F)]
  genNames <- paste0("G", 1:(ps*2))
  colnames(gen_m) <- as.character(hetTable$Pos[order(hetTable$Pos, decreasing = F)])
  rownames(gen_m) <- genNames
  return(gen_m)
}}  
##--------------------

## Computing nucleotide diversity 
####
##
computing_Pi <- function(haps_m, haps_frq_v, ps) {{
  
  # haps_m <- uniq_haps[,-ncol(uniq_haps)]
  # haps_frq_v <- haps_freq
  # ps <- popSize
  
  dissimilarity.matrix <- apply(t(haps_m),2,function(x)colSums(x!=t(haps_m)))
  diag(dissimilarity.matrix)<-0
  
  pi_tmp <- matrix(, ncol = ncol(dissimilarity.matrix), nrow = nrow(dissimilarity.matrix))
  for(h1 in 1:(nrow(haps_m)-1)) {{
    for(h2 in (h1+1):nrow(haps_m)) {{
      pi_tmp[h1, h2] <- haps_frq_v[h1]*haps_frq_v[h2]*dissimilarity.matrix[h1,h2]
    }}
  }}
  pi <- ((ps*2)/((ps*2)-1))*2*sum(apply(pi_tmp, 1, sum, na.rm=T))
  return(pi)
}}  
##--------------------
#computing the number of mut of type X per haplotype 
computing_nb_mutX_per_haplotype <- function(genomes,muType) {{
  
  nb_haps_w_mutX <- apply(genomes, 1, function(x){{ sum(x == muType)}})
  return(nb_haps_w_mutX)
}}

#test whether the average nb of m2 mut differs in haplotypes w or wo m3 mut
nb_mX_onto_m3_background <- function(genomes, haps_w_mX, haps_w_m3, muType) {{
  
  index_m3 <- which(haps_w_m3 > 0)
  index_no_m3 <- which(haps_w_m3 == 0)
  
  index_m2 <- which(haps_w_mX > 0)
  index_no_m2 <- which(haps_w_mX == 0)
  
  #intersect haps with m3 and m2
  nb_haps_w_m3_w_m2 <- intersect(index_m3, index_m2)
  nb_haps_w_m3_wo_m2 <- intersect(index_m3, index_no_m2)
  nb_haps_wo_m3_w_m2 <- intersect(index_no_m3, index_m2)
  nb_haps_wo_m3_wo_m2 <- intersect(index_no_m3, index_no_m2)
  
  mean_nb_m2_within_m3_haps <- mean(computing_nb_mutX_per_haplotype(genomes[nb_haps_w_m3_w_m2,], muType))
  mean_nb_m2_within_nom3_haps <- mean(computing_nb_mutX_per_haplotype(genomes[nb_haps_wo_m3_w_m2,], muType))
  return(c(mean_nb_m2_within_m3_haps, mean_nb_m2_within_nom3_haps, 
           length(nb_haps_w_m3_w_m2), length(nb_haps_w_m3_wo_m2), length(nb_haps_wo_m3_w_m2), 
           length(nb_haps_wo_m3_wo_m2)))
}}

# get allele counts from the mutation table 
alleleCounts <- function(mutType) {{
  
  x <- as.numeric(unlist(cleanSegSitesTable %>% filter(V3 == mutType)  %>% select(V9)))
  
  return(x)
}}
# get allele physical positions from the mutation table 
MutPositions <- function(mutType) {{
  
  x <- as.numeric(unlist(cleanSegSitesTable %>% filter(V3 == mutType)  %>% select(V4)))
  
  return(x)
}}

# Compute LD statistics for each pairwise combination of CREalleles and CS (mut2) alleles
computing_LD <- function(PhyPosCREallele) {{
  
  #PhyPosCREallele <- 100552
  for(i in CS_delSites) {{
    
    #i <- 47226
    df.tmp <- df_genomes %>% select(colnames(df_genomes)[colnames(df_genomes) %in% c(paste0("X",PhyPosCREallele), paste0("X",i))])
    countHapsLD <- df.tmp %>% group_by_all %>% count
    
    mObsFreq <- matrix(rep(0,4), ncol = 2)
    colnames(mObsFreq) <- c("A", "m3")
    rownames(mObsFreq) <- c("A", "m2")
    mObsFreq[1,1] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "A" & countHapsLD[,2] == "A"),3]/(popSize*2))
    mObsFreq[1,2] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "A" & countHapsLD[,2] == "m3"),3]/(popSize*2))
    mObsFreq[2,1] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "m2" & countHapsLD[,2] == "A"),3]/(popSize*2))
    mObsFreq[2,2] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "m2" & countHapsLD[,2] == "m3"),3]/(popSize*2))
    
    mExpFreqIndep <- matrix(rep(0,4), ncol = 2)
    colnames(mExpFreqIndep) <- c("A", "m3")
    rownames(mExpFreqIndep) <- c("A", "m2")
    freqM3 <- as.numeric(cleanSegSitesTable %>% filter(V4 == PhyPosCREallele) %>% select(V9))/(popSize*2)
    freqM2 <- as.numeric(cleanSegSitesTable %>% filter(V4 == i) %>% select(V9))/(popSize*2)
    
    mExpFreqIndep[1,1] <- (1-freqM2)*(1-freqM3)
    mExpFreqIndep[1,2] <- (1-freqM2)*(freqM3)
    mExpFreqIndep[2,1] <- (freqM2)*(1-freqM3)
    mExpFreqIndep[2,2] <- (freqM2)*(freqM3)
    
    #mExpFreqIndep
    mObsFreq[is.na(mObsFreq)] <- 0
    XiSqPvalue <- chisq.test(mObsFreq*(popSize*2), mExpFreqIndep*(popSize*2))$p.value
    D <- mObsFreq[1,1]*mObsFreq[2,2]-mObsFreq[1,2]*mObsFreq[2,1]
    rSq <- D^2/((1-freqM2)*freqM2*(1-freqM3)*freqM3)
    cat(paste(PhyPosCREallele, i,D, rSq, XiSqPvalue, sep = "\\t"), sep = "\\n", file = sumStats_LD_fName, append = T)
    
  }}
}}

# get the nb of heterozygous ind PER SITE and how many ind have a m2 or m1 onto a m3++ background
countingMut1or2_m3minus_combinations_perSite <- function(mutationType, nbInd){{
  
  # geral variables : df_genomes & hetTable 
  # debugging 
  # mutationType <- "m2" or "m1"
  # posVector <- posM2
  # nbInd <- 1000
  
  posVector <- sort(unlist(hetTable %>% filter(muType == mutationType) %>% select(Pos)))
  # posM1 <- sort(unlist(mut_tb %>% filter(muType == "m1") %>% select(Pos)))
  posM3 <- sort(unlist(hetTable %>% filter(muType == "m3") %>% select(Pos)))
  gen_m <- df_genomes %>% mutate(IndID=rep(1:nbInd, each=2))
  
  propIndforMutTypeM_m3minus <- c()
  
  for(m in 1:length(posVector)) {{
    #m <- 55
    indCarryingM2 <- gen_m$IndID[which(gen_m[,colnames(gen_m) == paste0("X",posVector[m])] == mutationType)]
    indHet <- as.numeric(names(table(indCarryingM2)[table(indCarryingM2) == 1]))
    if(length(indHet) > 0) {{ # there are het ind for the specific m2
      
      subDF <- gen_m[gen_m$IndID %in% indHet, c(which(colnames(gen_m) == paste0("X", posVector[m])),c(which(colnames(gen_m) %in% paste0("X",posM3)),which(colnames(gen_m) =="IndID")))]
      m3G1G2 <- matrix(apply(subDF, 1, function(x) {{ sum(x == "m3")}}), ncol = 2, byrow = T)
      
      genOne <- 1
      genTwo <- 2
      nbIndM2_m3Minus <- c()
      
      for(ind in 1:length(indHet)) {{
        
        #ind <- 1
        idxm2 <- which(subDF[genOne:genTwo,which(colnames(subDF) == paste0("X",posVector[m]))] == mutationType)
        nbIndM2_m3Minus[ind] <- sum(m3G1G2[ind,idxm2] < m3G1G2[ind,-idxm2])
        
        genOne <- genOne + 2 
        genTwo <- genTwo + 2 
      }}
      propIndM2_m3minus <- sum(nbIndM2_m3Minus)/length(indHet)
      names(propIndM2_m3minus) <- posVector[m]
      propIndforMutTypeM_m3minus[m] <- propIndM2_m3minus
      
    }} else {{ # no m2 het ind 
      propIndM2_m3minus <- NA
      names(propIndM2_m3minus) <- posVector[m]
      propIndforMutTypeM_m3minus[m] <- propIndM2_m3minus
      indHet <- 0
    }}

  }}
  return(list(df=propIndforMutTypeM_m3minus, nbhetInd=indHet))
}}
####-----------------------------

# get the nb of m2 mutations onto m3++ and m3-- backgrounds across individuals
computing_Nbm2_Nbm3_acrossGenomes <- function(mHapFull, sampleSize) {{
  genOne <- 1
  genTwo <- 2
  m2OverExp <- c()
  m2UnderExp <- c()
  
  for(ind in 1:sampleSize) {{
    
    #ind <- 2
    
    if(ind == 1) {{
      m2m3PerInd <- matrix(c(apply(mHapFull[genOne:genTwo,], 1, function(x) {{sum(x=="m2")}}), apply(mHapFull[genOne:genTwo,], 1, function(x) {{sum(x=="m3")}})), ncol = 2, byrow = F)
      m2m3AcrossInd <- m2m3PerInd
      
      # getting nb of m2 onto overexpressed background per ind (+m3). 
      tmpM2OverExp <- m2m3PerInd[which(m2m3PerInd[,2] == max(m2m3PerInd[,2])),1]
      # getting nb of m2 onto underexpressed background per ind (-m3).
      tmpM2UnderExp <- m2m3PerInd[which(m2m3PerInd[,2] == min(m2m3PerInd[,2])),1]
      if(length(tmpM2OverExp) == 2) {{
        
        m2OverExp <- NA
        m2UnderExp <- NA
        
      }} else {{
        
        m2OverExp <- m2m3PerInd[which(m2m3PerInd[,2] == max(m2m3PerInd[,2])),1]
        m2UnderExp <- m2m3PerInd[which(m2m3PerInd[,2] == min(m2m3PerInd[,2])),1]
        
      }}
      
    }} else {{
      m2m3PerInd <- matrix(c(apply(mHapFull[genOne:genTwo,], 1, function(x) {{sum(x=="m2")}}), apply(mHapFull[genOne:genTwo,], 1, function(x) {{sum(x=="m3")}})), ncol = 2, byrow = F)
      m2m3AcrossInd <- rbind(m2m3AcrossInd,m2m3PerInd)  
      
      # getting nb of m2 onto overexpressed background per ind (+m3). 
      tmpM2OverExp <- m2m3PerInd[which(m2m3PerInd[,2] == max(m2m3PerInd[,2])),1]
      # getting nb of m2 onto underexpressed background per ind (-m3).
      tmpM2UnderExp <- m2m3PerInd[which(m2m3PerInd[,2] == min(m2m3PerInd[,2])),1]
      
      if(length(tmpM2OverExp) == 2) {{
        m2OverExp <- c(m2OverExp, NA)
        m2UnderExp <- c(m2UnderExp, NA)
      }} else {{
        m2OverExp <- c(m2OverExp,tmpM2OverExp)
        m2UnderExp <- c(m2UnderExp,tmpM2UnderExp)
        
      }}
    }}
    
    genOne <- genOne + 2
    genTwo <- genTwo + 2
  }}
  # here it gives the number of m2 and m3 onto the same genome (out of 2000)
  m2m3AcrossIndDF <- data.frame(m2m3AcrossInd)
  colnames(m2m3AcrossIndDF) <- c("m2", "m3")
  rownames(m2m3AcrossIndDF) <- rownames(mHapFull)
  
  # here it gives the number of m2 onto the ind genome with the highest nb of m3
  # int this case only individuals carrying a dif nb of m3 are taken into account 
  # in the computation
  m2OverUnderExp <- data.frame(cbind(m2OverExp,m2UnderExp))
  colnames(m2OverUnderExp) <- c("m2ontoM3+", "m2ontoM3-")
  return(m2OverUnderExp)
}}
#############################
################-- End of functions
######
EOF
        cat << 'EOF' > {output.sfs}
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
#############################
###
###
###
### FUNCTIONS
###
###
#############################
plottingHaplotypeSumStats <- function(listINPUT, nameCol, pName, yAxisLabel) {{
  # global variables: namesStats, selCoeff, scenarios, majorScenarios, aseScenario
  # listINPUT = perSiteStats
  # nameCol = sumStats
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  selCoeffLabels <- paste0("s=", selCoeff)
  aveNbMutPerHap_m <- matrix(as.numeric(unlist(lapply(listINPUT, "[", which(namesStats == nameCol)))), 
                             ncol=length(selCoeff))
  colnames(aveNbMutPerHap_m) <- names(listINPUT)
  df_haps <- data.frame(aveNbMutPerHap_m)
  df_haps <- df_haps %>% select(names(df_haps)) %>%
    pivot_longer(., cols = names(df_haps), names_to = "Scenario", values_to = "sumStat") %>% 
    mutate(selectionCoeff=rep(selCoeff, replicates), aseModel=rep(aseScenario,replicates), model=rep(majorScenarios, replicates)) 
  df_haps$sumStat[is.na(df_haps$sumStat)] <- 0
  
  # this part performed a violin plot which is confusing to read
  # p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
  #   labs(title = plotTitle, x = "Selection model", y = yLabel) + 
  #   geom_violin(aes(color=Scenario), trim=FALSE, position = position_dodge(0.9)) + 
  #   geom_boxplot(aes(fill=Scenario), width=0.5, position = position_dodge(0.9)) + 
  #   scale_x_discrete(labels=unique(selCoeffLabels)) + 
  #   scale_fill_manual(values=colorsClasses) + scale_color_manual(values=colorsClasses) +
  #   guides(color = "none") + theme_classic()
   
  # added subtitle that says the model parameters for each graphs
  modelSubtitle <- gsub("_", " | ", gsub("fyonetal-", "", targetModel))
  
  p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
    labs(title = pName, subtitle = modelSubtitle, x = "Selection model", y = yAxisLabel) + 
    geom_boxplot(aes(color=Scenario)) + 
    scale_x_discrete(labels=unique(selCoeffLabels)) + 
    scale_color_manual(values=colorsClasses) + 
    guides(color = "none") + theme_classic(base_size = 15)


  # p-values
  pValuesList <- list()
  for(s in unique(selCoeff)) {{ 
    #print(s)
    pValuesList[[s]] <- wilcox.test(unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "noASE") %>% select("sumStat")), 
                                    unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "ASE") %>% select("sumStat")), alternative = "two")$p.value
    
  }}
  
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F, rownames=unique(selCoeff))
  fichier_pvalues <- paste0(outputPval, "/", targetModel, "_", nameCol, "_pValues.txt")
  write.table(pValues_v, file = fichier_pvalues,
              quote = F, row.names = T, col.names = F, sep = "\\t")
  print(p + stat_compare_means(aes(group = aseModel), label = "p.signif"))
  
}}
########################--------------------
#########
### End of functions


#############################
###
###
###       MAIN
###
###
#############################

args = commandArgs(trailingOnly=TRUE)
#setting working directory 
wrkDir <- args[1]
setwd(wrkDir)

outputPDF <- paste0(wrkDir, "/pdfs_sm2_analysis")
if(!dir.exists(outputPDF)) {{
  dir.create(outputPDF)
}}

outputPval <- paste0(wrkDir, "/p_values")
if(!dir.exists(outputPval)) {{
  dir.create(outputPval)
}}

#####-----------------
##
## VARIABLES
##
#####-----------------
# nb of replicates / model
replicates <- as.numeric(args[3])

# setting target model 
targetModel <- args[2]

# selection coefficients
selCoeff <- rep(c("-0.001", 
                  "-0.01", 
                  "-0.0425",
                  "-0.1", 
                  "-0.3"), each=2)
# scenarios
scenarios <- paste0(targetModel, "_s", paste0(selCoeff, c("", "_epist")))
majorScenarios <- paste0(targetModel, "_s=", selCoeff)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")
#####-----------------
#####-------- End of variables
#####----


# opening the pdf 
pdf(file=paste0(outputPDF, "/", targetModel, "_SFSstats_violins.pdf"), height = 6, width = 6)

# reading summary stats SLiM
perSiteStats <- list()
pi <- list()

for (s in scenarios) {{
  print(s);
  
  perSiteStats[[s]] <- read.table(paste0(wrkDir, "/", s, "/sumStats_", s, "_perSite.txt"), header = T, sep = "\\t", 
                                  na.strings = "NaN")
  pi[[s]] <- read.table(file=paste0(wrkDir, "/", s, "/pi_", s, ".txt"), header = T, sep = "\\t", na.strings = NA)
}}

# processing SFS stats 
namesStats <- names(perSiteStats[[1]])
namesPi <- names(pi[[1]])

#######################
##
### Plotting & computing p-Values
##
#######################
## plot Nb of deleterious mutations - coding sequence 
sumStats <- "Nb_m2"
plotTitle <- "Deleterious coding"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot nb of neutral coding
sumStats <- "Nb_m1"
plotTitle <- "Neutral coding"
yLabel <- "Number of mutations"
#===================================
plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)

## plot nb of neutral cis-regulatory
sumStats <- "Nb_m3"
plotTitle <- "regulatory"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot overall heterozygosity 
sumStats <- "AveHet"
plotTitle <- "Overall heterozygosity"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m1 heterozygosity 
sumStats <- "AveHet_m1"
plotTitle <- "Heterozygosity neutral coding"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m2 heterozygosity 
sumStats <- "AveHet_m2"
plotTitle <- "Heterozygosity deleterious coding"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m3 heterozygosity 
sumStats <- "AveHet_m3"
plotTitle <- "Heterozygosity regulatory"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
###------------------------
###           Pi 
###------------------------
# plot overall pi
namesStats <- names(pi[[1]])
 
sumStats <- "Overall_Nucleotide_diversity"
plotTitle <- expression("Overall nucleotide diversity ("*pi*")")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Coding region pi
sumStats <- "CS_Nucleotide_diversity"
plotTitle <- expression("Nucleotide diversity ("*pi*") - Coding seq.")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Regulatory region pi
sumStats <- "CRE_Nucleotide_diversity"
plotTitle <- expression("Nucleotide diversity ("*pi*") - CRE")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Overall_SegSites
sumStats <- "Overall_SegSites"
plotTitle <- "Overall segrating sites (S)"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Overall_theta S
sumStats <- "Overall_ThetaS"
plotTitle <- expression("Overall "*theta*"(S)")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CS_SegSites 
sumStats <- "CS_SegSites"
plotTitle <- "Segrating sites (S) - Coding seq."
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CS_ThetaS
sumStats <- "CS_ThetaS"
plotTitle <-  expression(theta*"(S) - Coding Seq.")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CRE_SegSites
sumStats <- "CRE_SegSites"
plotTitle <-  "Segrating sites (S) - CRE"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CRE_ThetaS
sumStats <- "CRE_ThetaS"
plotTitle <-  expression(theta*"(S) - CRE")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## OverallTajimasD
sumStats <- "OverallTajimasD"
plotTitle <-  "Overall Tajima's D"
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## TajimasD - CS
sumStats <- "TajimasDCS"
plotTitle <-  "Tajima's D - Coding seq."
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
##TajimasD - CRE
sumStats <- "TajimasDCRE"
plotTitle <-  "Tajima's D - CRE"
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
dev.off()


EOF
        cat << 'EOF' > {output.hap}
library(ggplot2)
library(tidyverse)
library(wesanderson)
library(ggpubr)# AD: ajout de cette librairie pour stat_compare_means
#############################
###
###
###
### FUNCTIONS
###
###
#############################
plottingHaplotypeSumStats <- function(listINPUT, nameCol, pName, yAxisLabel,natScale=T) {{
  # global variables: namesStats, selCoeff, scenarios, majorScenarios, aseScenario
  # listINPUT = pi
  # nameCol = sumStats
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  selCoeffLabels <- paste0("s=", selCoeff)
  aveNbMutPerHap_m <- matrix(as.numeric(unlist(lapply(listINPUT, "[", which(namesStats == nameCol)))), 
                             ncol=length(selCoeff))
  colnames(aveNbMutPerHap_m) <- names(listINPUT)
  df_haps <- data.frame(aveNbMutPerHap_m)
  df_haps <- df_haps %>% select(names(df_haps)) %>%
    pivot_longer(., cols = names(df_haps), names_to = "Scenario", values_to = "sumStat") %>% 
    mutate(selectionCoeff=rep(selCoeff, replicates), aseModel=rep(aseScenario,replicates), model=rep(majorScenarios, replicates)) 
  df_haps$sumStat[is.na(df_haps$sumStat)] <- 0

# added subtitle that says the model parameters for each graphs
  modelSubtitle <- gsub("_", " | ", gsub("fyonetal-", "", targetModel))
  
  if(natScale) {{
    p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
      # On ajoute subtitle = cleanSubtitle ici
      labs(title = pName, subtitle = modelSubtitle, x = "Selection model", y = yAxisLabel) + 
      geom_boxplot(aes(color=Scenario)) + 
      scale_x_discrete(labels=unique(selCoeffLabels)) + 
      scale_color_manual(values=colorsClasses) +
      guides(color = "none") + theme_classic()
    
  }} else {{
    p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
      # Et on l'ajoute aussi ici pour l'échelle log
      labs(title = pName, subtitle = modelSubtitle, x = "Selection model", y = yAxisLabel) + 
      coord_transform(y = "log10") +
      geom_boxplot(aes(color=Scenario)) + 
      scale_x_discrete(labels=unique(selCoeffLabels)) + 
      scale_color_manual(values=colorsClasses) +
      guides(color = "none") + theme_classic()
  }}
  
  # p-values
  pValuesList <- list()

  for(s in unique(selCoeff)) {{ 
    #print(s)
    pValuesList[[s]] <- wilcox.test(unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "noASE") %>% select("sumStat")), 
                                    unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "ASE") %>% select("sumStat")), alternative = "two")$p.value
  }}
  
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F, rownames=unique(selCoeff))
  fichier_pvalues <- paste0(outputPval, "/", targetModel, "_", nameCol, "_pValues.txt")
  write.table(pValues_v, file = fichier_pvalues,
              quote = F, row.names = T, col.names = F, sep = "\\t")
  print(p + stat_compare_means(aes(group = aseModel), label = "p.signif"))
  
}}
## ---------------------

# computing mean, median and fisher exact test on correlation m2-m3 p-values
computingSumStatsHapStats <- function(listHaps) {{
  
  #listHaps <- hapStats
  # computing mean and median
  sStatsMean <- lapply(listHaps, function(x) {{
    apply(na.omit(x), 2, mean, na.rm=T)
  }})
  
  sStatsMedian <- lapply(listHaps, function(x) {{
    apply(x, 2, median, na.rm=T)
  }})
  
  # table format 
  meanDF <- data.frame(do.call(rbind, sStatsMean))
  medianDF <- data.frame(do.call(rbind, sStatsMedian))
  write.table(meanDF, file = paste0(outputSS, "/", targetModel, "_mean_hapStats.txt"), quote = F, col.names = T, sep = "\\t")
  write.table(medianDF, file = paste0(outputSS, "/", targetModel, "_median_hapStats.txt"), quote = F, col.names = T, sep = "\\t")
  
  # computing nb of non-sign and sign correlations in noEpist and epist models
  noSignCorr <- unlist(lapply(listHaps, function(x) {{
    x %>% filter(pValCorrM2M3acrossHaps > 0.05) %>% nrow()
  }}))
  
  signCorr <-unlist(lapply(listHaps, function(x) {{
    x %>% filter(pValCorrM2M3acrossHaps <= 0.05) %>% nrow()
  }}))
  
  COUNT <- 1 
  fisherExTest <- list()
  while(COUNT < length(majorScenarios)) {{
    contgTable <- matrix(rbind(noSignCorr[COUNT:(COUNT+1)], 
                               signCorr[COUNT:(COUNT+1)]), ncol = 2)
    fisherExTest[[majorScenarios[COUNT]]] <- c(fisher.test(contgTable)$estimate, fisher.test(contgTable)$p.value,
                                               fisher.test(contgTable)$conf.int[1], fisher.test(contgTable)$conf.int[2])
    names(fisherExTest[[majorScenarios[COUNT]]]) <- c("Odds.Ratio", "p-value", "2.5%CI", "97.5%CI")
    COUNT <- COUNT + 2
  }}
  # transformign table format
  exactFTestCorr <- do.call(rbind, fisherExTest)
  # printing 
  write.table(exactFTestCorr, file = paste0(outputSS,"/", targetModel, "_FisherExTestCorr_hapStats.txt"), quote = F, col.names = T, row.names = T, sep = "\\t")
  
}}
## ---------------------
########################--------------------
#########
### End of functions

#############################
###
###
###       MAIN
###
###
#############################

#setting working directory
args = commandArgs(trailingOnly=TRUE) 
wrkDir <- args[1]
setwd(wrkDir)

# creating a pdf folder in case it doesn't exist to save pdf files with results
outputPDF <- paste0(wrkDir, "/pdfs_sm2_analysis")
if(!dir.exists(outputPDF)) {{
  dir.create(outputPDF)
}}

outputPval <- paste0(wrkDir, "/p_values")
if(!dir.exists(outputPval)) {{
  dir.create(outputPval)
}}

# creating a sumStats folder in case it doesn't exist to save pdf files with results
if(!dir.exists(paste0(wrkDir, "/sumStats"))) {{
  dir.create(paste0(wrkDir, "/sumStats"))
  outputSS <- paste0(wrkDir, "/sumStats")
}} else {{
  outputSS <- paste0(wrkDir, "/sumStats")
}}


#####-----------------
##
## VARIABLES
##
#####-----------------
# nb of replicates / model
replicates <- as.numeric(args[3])

# setting target model 
targetModel <- args[2]

# selection coefficients
selCoeff <- rep(c("-0.001", 
              "-0.01", 
              "-0.0425",
              "-0.1", 
              "-0.3"), each=2)
# scenarios
scenarios <- paste0(targetModel, "_s", paste0(selCoeff, c("", "_epist")))
majorScenarios <- paste0(targetModel, "_s=", selCoeff)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")
#####-----------------
#####-------- End of variables
#####----

# opening the pdf 
pdf(file=paste0(outputPDF, "/", targetModel, "_HapsPlots_violins.pdf"), height = 6, width = 6)

# reading summary stats SLiM
hapStats <- list()
pi <- list()

for (s in scenarios) {{
  print(s);
  
  hapStats[[s]] <- read.table(paste0(s, "/sumStats_", s, "_haps.txt"), header = T, sep = "\\t", 
                              na.strings = "NA")
  
}}

# processing hap stats 
listHaps <- list()
namesStats <- names(hapStats[[1]])
listHaps <- lapply(hapStats, function(x) {{x %>% mutate(NbofMut3Haps=apply(x[,c(3,4)], 1, sum))}})
listHaps <- lapply(listHaps, function(x) {{x %>% mutate(relNbHaps_m3_m2=NbHaps_w_m3_w_m2/NbofMut3Haps, 
                                                       relNbHaps_NOm3_m2=NbHaps_wo_m3_w_m2/(2000-NbofMut3Haps))}})
#######################
##
### Plotting & computing p-Values
##
#######################
#Deleterious mutations onto overexpressed haps
sumStats <- "avePerHap_nbMut2_ontoM3haps"
plotTitle <- "Ave Nb of deleterious mutations onto overexpressed haplotypes"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#Deleterious mutations onto underexpressed haps
sumStats <- "avePerHap_nbMut2_ontoNOM3haps"
plotTitle <- "Ave Nb of deleterious mutations onto underexpressed haplotypes"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#nb overexpressed haps with deleterious mutations
sumStats <- "NbHaps_w_m3_w_m2"
plotTitle <- "Nb of overexpressed haps with deleterious muts"
yLabel <- "Number of haplotypes"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#nb underexpressed haps with deleterious mutations
sumStats <- "NbHaps_wo_m3_w_m2"
plotTitle <- "Nb of underexpressed haps with deleterious muts"
yLabel <- "Number of haplotypes"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#nb of unique haps
sumStats <- "nbUniqHap"
plotTitle <- "Nb of unique haplotypes"
yLabel <- "Number of haplotypes"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

# correlation between m2 and m3. 
sumStats <- "CorrM2M3acrossHaps"
plotTitle <- "Correlation m2 mut - m3 mut within haps"
yLabel <- "Correlation"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

# p-values correlation between m2 and m3. 
sumStats <- "pValCorrM2M3acrossHaps"
plotTitle <- "Cor. pvalue m2 mut - m3 mut within haps"
yLabel <- "P-value"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

# ave nv of m2 onto m3+. 
sumStats <- "meanNbm2ontoM3plus"
plotTitle <- "Ave. nb. m2 onto m3+"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel, natScale = T)

# ave nv of m2 onto m3+. 
sumStats <- "meanNbm2ontoM3minus"
plotTitle <- "Ave. nb. m2 onto m3-"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel, natScale = T)
dev.off()

#######################
##
### computing mean and median
##
#######################
computingSumStatsHapStats(hapStats)


EOF

        cat << 'EOF' > {output.fix_reco}
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(wesanderson)
library(stringr)
library(formattable)

####-------------------
## functions 2
####-------------------



# facet_wrap ~selectionCoeff each plot has two lines ASE vs noASE.

plottingAveSumStatsAcrossRecCategoriesGroupSelection <- function(sumStatLabel,titlePlot,statsYlabel,hapYN,piYN,logScale){{
  x <- sumStatLabel
  list_sumStat <- lapply(scenarios, function(m) {{
    if(perHap) {{
      tmpSS <- read.table(paste0(m, "/sumStats_", m, "_haps.txt"), 
                          header = T, sep = "\\t", na.strings = NA)
    }} else {{
      if(pi){{
        tmpSS <- read.table(paste0(m, "/pi_", m, ".txt"), 
                            header = T, sep = "\\t", na.strings = NA)
      }} else {{
        tmpSS <- read.table(paste0(m, "/sumStats_", m, "_perSite.txt"), 
                            header = T, sep = "\\t", na.strings = NA)
      }}
    }}
    
    return(tmpSS %>% select({{{{x}}}}))
  }})
  names(list_sumStat) <- scenarios
  list_sumStat <- lapply(list_sumStat, function(k) {{k[which(is.na(k)),1] <- NA; return(k)}})
  df_sumStat <- do.call(cbind, list_sumStat)
  colnames(df_sumStat) <- scenarios
  mean_sumStat <- data.frame(apply(df_sumStat, 2, mean, na.rm=T))
  colnames(mean_sumStat) <- "meanSS"
  mean_sumStat <- mean_sumStat %>% mutate(FullModel=rownames(mean_sumStat), 
                                          PhysicalDist = phyDistLabel,
                                          GenetDist = rep(recRate, each=10),
                                          selectionCoeff = rep(selCoeff, 5),
                                          ASEModel = rep(aseScenario, 5),
                                          distanceLabels=scientific(GenetDist, digits = 2))
  
  mean_sumStat$distanceLabels <- factor(mean_sumStat$distanceLabels, levels=unique(mean_sumStat$distanceLabels))
  mean_sumStat$selectionCoeff <- factor(mean_sumStat$selectionCoeff,
                                        levels=unique(selCoeff))
  
  if(logScale) {{
    ggplot(mean_sumStat, aes(x=log10(GenetDist), y=log10(meanSS), group = selectionCoeff, colour = selectionCoeff)) + 
      geom_line() + 
      geom_point() +
      facet_wrap(~ASEModel) + 
      labs(title=titlePlot,
           subtitle = creTitle,
           color = "Selection coeff.") +
      ylab(statsYlabel) +
      xlab("Recombination rate - log10") +
      theme_bw() + 
      scale_color_manual(values = colorsClasses[seq(1,10,by=2)]) +
      #scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
      xlab(xLabel)
    ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_", sumStatLabel,"_acrossGenetDist", pdTag, ".pdf"), width = 20, height = 10, units = "cm")
    
  }} else {{
    ggplot(mean_sumStat, aes(x=log10(GenetDist), y=meanSS, group = selectionCoeff, colour = selectionCoeff)) + 
      geom_line() + 
      geom_point() +
      facet_wrap(~ASEModel) + 
      labs(title=titlePlot,
           subtitle = creTitle,
           color = "Selection coeff.") +
      ylab(statsYlabel) +
      xlab("Recombination rate - log10") +
      theme_bw() + 
      scale_color_manual(values = colorsClasses[seq(1,10,by=2)]) + #
      #scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
      xlab(xLabel)
    ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_", sumStatLabel,"_acrossGenetDist", pdTag, ".pdf"), width = 20, height = 10, units = "cm")  
  }}
}}

# facet_wrap per selection. ASE vs noASE within each plot 
plottingAveSumStatsAcrossRecCategoriesGroupASE <- function(sumStatLabel,titlePlot,statsYlabel,hapYN,piYN,logScale){{
  x <- sumStatLabel
  list_sumStat <- lapply(scenarios, function(m) {{
    if(perHap) {{
      tmpSS <- read.table(paste0(m, "/sumStats_", m, "_haps.txt"), 
                          header = T, sep = "\\t", na.strings = NA)
    }} else {{
      if(pi){{
        tmpSS <- read.table(paste0(m, "/pi_", m, ".txt"), 
                            header = T, sep = "\\t", na.strings = NA)
      }} else {{
        tmpSS <- read.table(paste0(m, "/sumStats_", m, "_perSite.txt"), 
                            header = T, sep = "\\t", na.strings = NA)
      }}
    }}
    
    return(tmpSS %>% select({{{{x}}}}))
  }})
  names(list_sumStat) <- scenarios
  list_sumStat <- lapply(list_sumStat, function(k) {{k[which(is.na(k)),1] <- NA; return(k)}})
  df_sumStat <- do.call(cbind, list_sumStat)
  colnames(df_sumStat) <- scenarios
  mean_sumStat <- data.frame(apply(df_sumStat, 2, mean, na.rm=T))
  colnames(mean_sumStat) <- "meanSS"
  mean_sumStat <- mean_sumStat %>% mutate(FullModel=rownames(mean_sumStat), 
                                          PhysicalDist = phyDistLabel,
                                          GenetDist = rep(recRate, each=10),
                                          selectionCoeff = rep(selCoeff, 5),
                                          ASEModel = rep(aseScenario, 5),
                                          distanceLabels=scientific(GenetDist, digits = 2))
  
  mean_sumStat$distanceLabels <- factor(mean_sumStat$distanceLabels, levels=unique(mean_sumStat$distanceLabels))
  mean_sumStat$selectionCoeff <- factor(mean_sumStat$selectionCoeff,
                                        levels=unique(selCoeff))
  
  if(logScale) {{
    ggplot(mean_sumStat, aes(x=log10(GenetDist), y=log10(meanSS), group = ASEModel, colour = ASEModel)) + 
      geom_line() + 
      geom_point() +
      facet_wrap(~selectionCoeff, scales = "free") + 
      labs(title=titlePlot,
           subtitle = creTitle,
           color = "Selection coeff.") +
      ylab(statsYlabel) +
      xlab("Recombination rate - log10") +
      theme_bw() + 
      scale_color_manual(values = colorsClasses) + #[seq(1,10,by=2)]
      #scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
      xlab(xLabel)
    ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_", sumStatLabel,"_acrossGenetDistGroupASE", pdTag, ".pdf"), width = 20, height = 10, units = "cm")

  }} else {{
    ggplot(mean_sumStat, aes(x=log10(GenetDist), y=meanSS, group = ASEModel, colour = ASEModel)) + 
      geom_line() + 
      geom_point() +
      facet_wrap(~selectionCoeff, scales = "free") + 
      labs(title=titlePlot,
           subtitle = creTitle,
           color = "Selection coeff.") +
      ylab(statsYlabel) +
      xlab("Recombination rate - log10") +
      theme_bw() + 
      scale_color_manual(values = colorsClasses) + #[seq(1,10,by=2)]
      #scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
      xlab(xLabel)
      ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_", sumStatLabel,"_acrossGenetDistGroupASE", pdTag, ".pdf"), width = 20, height = 10, units = "cm")
  }}
}}

#############################
###
###
###
### FUNCTIONS
###
###
#############################
collectingMutationsAcrossRuns <-  function(modelTags, mut_fixed) {{
  
  fullTableM2fixed <- list()
  fullTableM3fixed <- list()
  
  for(targetModel in modelTags) {{
    
    # debugging
    #targetModel <- scenarios[2]
    
    hasASE <- unlist(gregexpr(pattern ='epist', targetModel))
    whatRecomb <- unlist(strsplit(unlist(strsplit(targetModel, split = "_"))[2], split = "rec"))[2]

    ##### +1 pour les cas avec les mutations déléteres sur la cre
    globalTag <- paste0(unlist(strsplit(targetModel, split = "_"))[2:4], collapse = "_") 
    selTag <- unlist(strsplit(targetModel, split = "_"))[5]
    ######
    
    if(hasASE > 0) {{
      ASEscenario <- "ASE"
    }} else {{
      ASEscenario <- "noASE"
    }}
    
    
    #fixed mutations
    # openTableM2Fixed <- read.table(paste0(targetModel, "/", targetModel, "_g10000_m2_fMut.out"), sep = "\\t", header = F)
    # openTableM2Fixed <- openTableM2Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM2Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM2Fixed)))
    
    #fullTableM2fixed <- rbind(fullTableM2fixed, openTableM2Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM2Fixed))))
    
    if(mut_fixed) {{
      openTableM3Fixed <- read.table(paste0(targetModel, "/", targetModel, "_g10000_m3_fMut.out"), sep = "\\t", header = F)
    }} else {{
      openTableM3Fixed <- read.table(paste0(targetModel, "/", targetModel, "_g10000_m3_sMut.out"), sep = "\\t", header = F)
    }}
    
    openTableM3Fixed <- openTableM3Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM3Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM3Fixed)), 
                                                    selectionCoeff=rep(selTag, nrow(openTableM3Fixed)))
    
    openTableM3Fixed <- openTableM3Fixed %>% 
      mutate(Distance=rep(whatRecomb, nrow(openTableM3Fixed))) %>% 
      mutate(GenetDist=case_when(Distance == "0.5kb" ~ recRate[names(recRate) == "0.5kb"], Distance == "1kb" ~ recRate[names(recRate) == "1kb"], 
                                 Distance == "5kb" ~ recRate[names(recRate) == "5kb"], 
                                 Distance == "10kb" ~ recRate[names(recRate) == "10kb"], 
                                 Distance == "500kb" ~ recRate[names(recRate) == "500kb"]))
    if(hasASE > 0) {{
      fullTableM3fixed[[globalTag]][[selTag]] <- rbind(fullTableM3fixed[[globalTag]][[selTag]],openTableM3Fixed) 
      
    }} else {{
      fullTableM3fixed[[globalTag]][[selTag]] <- openTableM3Fixed
    }}
    #fullTableM3fixed[[globalTag]][[selTag]] <- print(fullTableM3fixed[[globalTag]][[selTag]], row.names=F)
    rm(openTableM3Fixed)
    #fullTableM3fixed[[globalTag]][selTag] <- rbind(fullTableM3fixed, openTableM3Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM3Fixed))))
    
  }}
  return(fullTableM3fixed)
}}
###---------------
collectingTimeToFixationAcrossRuns <-  function(modelTags, mut_fixed) {{
  
  fullTableM3fixed <- list()
  
  for(targetModel in modelTags) {{
    
    # debugging
    #targetModel <- scenarios[2]
    
    hasASE <- unlist(gregexpr(pattern ='epist', targetModel))
    whatRecomb <- unlist(strsplit(unlist(strsplit(targetModel, split = "_"))[2], split = "rec"))[2]


    ###### +1 pour les mutations déléteres sur la CRE
    globalTag <- paste0(unlist(strsplit(targetModel, split = "_"))[2:4], collapse = "_")
    selTag <- unlist(strsplit(targetModel, split = "_"))[5]
    ######
    
    if(hasASE > 0) {{
      ASEscenario <- "ASE"
    }} else {{
      ASEscenario <- "noASE"
    }}
    
    
    #fixed mutations
    # openTableM2Fixed <- read.table(paste0(targetModel, "/", targetModel, "_g10000_m2_fMut.out"), sep = "\\t", header = F)
    # openTableM2Fixed <- openTableM2Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM2Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM2Fixed)))
    
    #fullTableM2fixed <- rbind(fullTableM2fixed, openTableM2Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM2Fixed))))
    

    openTableM3Fixed <- read.table(paste0(targetModel, "/", targetModel, "_g10000_m3_aveTimeFix.out"), sep = "\\t", header = F)
    
    openTableM3Fixed <- openTableM3Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM3Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM3Fixed)), 
                                                    selectionCoeff=rep(selTag, nrow(openTableM3Fixed)))
    
    openTableM3Fixed <- openTableM3Fixed %>% 
      mutate(Distance=rep(whatRecomb, nrow(openTableM3Fixed))) %>% 
      mutate(GenetDist=case_when(Distance == "0.5kb" ~ recRate[names(recRate) == "0.5kb"], Distance == "1kb" ~ recRate[names(recRate) == "1kb"], 
                                 Distance == "5kb" ~ recRate[names(recRate) == "5kb"], 
                                 Distance == "10kb" ~ recRate[names(recRate) == "10kb"], 
                                 Distance == "500kb" ~ recRate[names(recRate) == "500kb"]))
    if(hasASE > 0) {{
      fullTableM3fixed[[globalTag]][[selTag]] <- rbind(fullTableM3fixed[[globalTag]][[selTag]],openTableM3Fixed) 
      
    }} else {{
      fullTableM3fixed[[globalTag]][[selTag]] <- openTableM3Fixed
    }}
    #fullTableM3fixed[[globalTag]][[selTag]] <- print(fullTableM3fixed[[globalTag]][[selTag]], row.names=F)
    rm(openTableM3Fixed)
    #fullTableM3fixed[[globalTag]][selTag] <- rbind(fullTableM3fixed, openTableM3Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM3Fixed))))
    
  }}
  return(fullTableM3fixed)
}}
###---------------
plottingFixedMutations <- function(nameModel, pName, yAxisLabel) {{
  
  # global variables: listFixedMut, sValues
  # nameModel = x (parsed)
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  # nameModel = "rec5kb_cre1kb"
  # pName = "Neutral regulatory: fixed mutations"
  # yAxisLabel = "Number of mutations"
  
  tmpDF <- do.call(rbind, listFixedMut[[nameModel]])
  selCoeffLabels <- paste0("s=", sValues)
  majorScenarios <- paste0(nameModel, "_",selCoeffLabels)
  colorsClasses <- wes_palette("Zissou1", length(majorScenarios)*2, type = "continuous")
  
  
  p <- ggplot(tmpDF, aes(x=selectionCoeff, y=V1)) + 
    labs(title = plotTitle, subtitle = nameModel,x = "Selection model", y = yLabel) + 
    geom_boxplot(aes(color=Model)) + 
    scale_x_discrete(labels=selCoeffLabels) + 
    scale_color_manual(values=colorsClasses) + 
    guides(color = "none") + theme_classic()
  #return(print(p+stat_compare_means(aes(group = ASEModel), label = "p.signif")))
  return(p + stat_compare_means(aes(group = ASEModel), label = "p.signif"))
}}
# -----------
plottingTimeToFix <- function(nameModel, pName, yAxisLabel) {{
  
  # global variables: listFixedMut, sValues
  # nameModel = x (parsed)
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  # nameModel = "rec5kb_cre1kb"
  # pName = "Neutral regulatory: fixed mutations"
  # yAxisLabel = "Number of mutations"
  
  tmpDF <- do.call(rbind, listTimeToFix[[nameModel]])
  tmpDF$V1[tmpDF$V1 == 0] <- NA
  selCoeffLabels <- paste0("s=", sValues)
  majorScenarios <- paste0(nameModel, "_",selCoeffLabels)
  colorsClasses <- wes_palette("Zissou1", length(majorScenarios)*2, type = "continuous")
  
  
  p <- ggplot(tmpDF, aes(x=selectionCoeff, y=V1)) + 
    labs(title = plotTitle, subtitle = nameModel,x = "Selection model", y = yLabel) + 
    geom_boxplot(aes(color=Model)) + 
    scale_x_discrete(labels=selCoeffLabels) + 
    scale_color_manual(values=colorsClasses) + 
    guides(color = "none") + theme_classic()
  #return(print(p+stat_compare_means(aes(group = ASEModel), label = "p.signif", method = "wilcox.test")))
  return(p + stat_compare_means(aes(group = ASEModel), label = "p.signif", method = "wilcox.test"))
}}

# -----------
plotNbFixedMutFunctionRecombination <- function(coefSelection,xName) {{
  
  # global: listFixedMut
  # coefSelection="-0.0425"
  #   xName="Recomb. rate (rho)"
  
  sTag <- paste0("s", selCoeffTag)
  l <- lapply(listFixedMut, function(x) {{x[[sTag]]}})
  ll <- do.call(rbind, l)
  ll$Distance <- factor(ll$Distance, levels=unique(ll$Distance))
  ll$GenetDist <- scientific(ll$GenetDist, digits = 2)
  ll$GenetDist <- factor(ll$GenetDist, levels=unique(ll$GenetDist))
  p <- ll %>% ggplot(aes(x=Distance, y=V1, group = GenetDist)) + 
    geom_boxplot() + 
    facet_wrap(~ASEModel) + 
    labs(title=paste0("Regulatory mutations: ", creTitle),
         subtitle = paste0("s=", coefSelection)) +
    ylab("Number of fixed mutations") +
    theme_bw() + 
    scale_x_discrete(labels=unique(ll$GenetDist)) + 
    xlab(xName)
  return(p + geom_smooth(method='lm', se=F,aes(group=1)))
}}
###------------------
computingMeanSDTimeToFixation <- function(nameModel, l, time=F) {{
  
  # global variables: listFixedMut, sValues
  # nameModel = x (parsed)
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  # nameModel = "rec0.5kb_cre5kb"
  # pName = "Neutral regulatory: fixed mutations"
  # yAxisLabel = "Number of mutations"
  
  tmpDF <- do.call(rbind, l[[nameModel]])
  if(time) {{
    tmpDF$V1[tmpDF$V1 == 0] <- NA # remove sims where no reg mut went to fixation   
  }}
  df_meanTime <- tmpDF %>% group_by(GenetDist, selectionCoeff, ASEModel) %>% summarise(meanTime=mean(V1, na.rm = T), sdTime=sd(V1,na.rm=T))
  
  return(df_meanTime)
}}
########################--------------------
#########
### End of functions


#############################
###
###
###       MAIN
###
###
#############################

args = commandArgs(trailingOnly=TRUE)
wrkDir <- args[1]
setwd(wrkDir)


# creating a pdf folder in case it doesn't exist to save pdf files with results
outputPDF <- paste0(wrkDir, "/pdfs_rec_analysis")
if(!dir.exists(outputPDF)) {{
  dir.create(outputPDF)
}}


#####-----------------
##
## VARIABLES
##
#####-----------------
#### added by july 2025

suffixModel <- args[2]   # ex: "fyonetal-h0.25_rec"
creSizeRaw <- args[3]    # ex: "1"
pd_string <- args[4]     # ex: "0.5,1,5"
pd_uscore <- args[5]     # ex: "0.5_1_5"
sM3_val <- args[6]       # ex: "-0.001"

# create an uniq name for each model to name the pdfs
h_val <- sub(".*-h(.*)_rec", "\\1", suffixModel)
pdf_prefix <- paste0("fyonetal-h", h_val, "_cre", creSizeRaw, "kb_sM3", sM3_val)


pdTag <- paste0("_PD_", pd_uscore)
creLabel <- paste0("kb_cre", creSizeRaw, "kb_sM3", sM3_val, "_s")
creTitle <- paste0("cre - ", creSizeRaw, "kb")
sValues <-  c(-0.001, -0.01, -0.0425, -0.1, -0.3)

# recombination rates
phyDist <- as.numeric(unlist(strsplit(pd_string, ",")))
recRate <- 1.25e-7*phyDist*1000
names(recRate) <- paste0(as.character(phyDist), "kb")
phyDistLabel <- rep(phyDist, each=2*length(sValues))

# selection coefficients
selCoeff <- rep(sValues, each=2)
aseSuffix <- rep(c("", "_epist"), length(selCoeff)/2)

# scenarios
scenarios <- paste0(suffixModel, phyDistLabel, creLabel, selCoeff, aseSuffix)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)
####

# scenarios
majorScenarios <- paste0(suffixModel, phyDistLabel[1], creLabel, "=", selCoeff)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")



#######################
##
### Plotting & computing p-Values
##
#######################
## plot Nb of fixed - neutral regulatory
fixed <- T
listFixedMut <- collectingMutationsAcrossRuns(scenarios, fixed)

xLabel <- "Recombination rate"
selCoeffTag <- "-0.3"


if(fixed) {{
  plotTitle <- "Regulatory: Fixed mutations"
  boxplotPerCRE <- paste0(outputPDF,"/", pdf_prefix, pdTag, "_fixed.pdf")
  boxplotAcrossRec <- paste0("pdfs_rec_analysis/", pdf_prefix, "_acrossGenetDist_s=", selCoeffTag, pdTag, "_fixed.pdf")
}} else {{
  plotTitle <- "Regulatory: Seg. mutations"
  boxplotPerCRE <- paste0(outputPDF,"/", pdf_prefix, pdTag, "_seg.pdf")
  boxplotAcrossRec <- paste0("pdfs_rec_analysis/", pdf_prefix, "_acrossGenetDist_s=", selCoeffTag, pdTag, "_seg.pdf")
}}

yLabel <- "Number of mutations"
fixPlots <- lapply(names(listFixedMut), function(x) {{
  scenarioRecCRE <- x; 
  plottingFixedMutations(scenarioRecCRE,plotTitle, yLabel)
  }})

ggarrange(plotlist = fixPlots, ncol = 3, nrow = 2, 
          labels = LETTERS[1:length(fixPlots)]) %>% ggexport(filename = boxplotPerCRE, height = 8, width = 12)
#===================================         
## plot nb. of fixed reg mutations as a function of recombination 
plotNbFixedMutFunctionRecombination(selCoeffTag, xLabel)
ggsave(file=boxplotAcrossRec, width = 20, height = 10, units = "cm")
#===================================
#===================================         
## plot average time to fixation of reg mutations as a function of recombination 
listTimeToFix <- collectingTimeToFixationAcrossRuns(scenarios, T)
yLabel <- "Time to fixation"
fixPlots <- lapply(names(listTimeToFix), function(x) {{
  scenarioRecCRE <- x;

  plottingTimeToFix(scenarioRecCRE,"Regulatory: fixed mutations", yLabel)
}})
# computing mean time to fix
meanTtoFixation <- lapply(names(listTimeToFix), function(x) {{
  scenarioRecCRE <- x; 
  computingMeanSDTimeToFixation(scenarioRecCRE, listTimeToFix)
}})
names(meanTtoFixation) <- names(listTimeToFix)

meanTtoFixation
m <- do.call(rbind, meanTtoFixation) 

m <- m %>% mutate(distanceLabels=scientific(GenetDist, digits = 2))
m$distanceLabels <- factor(m$distanceLabels, levels=unique(m$distanceLabels))

#colorsClasses <- wes_palette("Zissou1", length(majorScenarios)*2, type = "continuous")
ggplot(m, aes(x=log10(GenetDist), y=meanTime, group = selectionCoeff, colour = selectionCoeff)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~ASEModel) + 
  labs(title="Regulatory mutations: time to fixation",
       subtitle = creTitle,
       color = "Selection coeff.") +
  ylab("Average time to fixation (g)") +
  xlab("Recombination rate - log10") +
  theme_bw() + 
  scale_color_manual(values = colorsClasses[seq(1,10,by=2)])
  #scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
  xlab(xLabel)
ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_timeToFixation_acrossGenetDist", pdTag, ".pdf"), width = 20, height = 10, units = "cm")
#===================================
#===================================         
## plot average number of mutations that went to fixation

# computing mean time to fix
meanM3Number <- lapply(names(listFixedMut), function(x) {{
  scenarioRecCRE <- x; 
  computingMeanSDTimeToFixation(scenarioRecCRE,listFixedMut,time=F)
}})
names(meanM3Number) <- names(listFixedMut)

m <- do.call(rbind, meanM3Number) 

m <- m %>% mutate(distanceLabels=scientific(GenetDist, digits = 2))
m$distanceLabels <- factor(m$distanceLabels, levels=unique(m$distanceLabels))

#colorsClasses <- wes_palette("Zissou1", length(majorScenarios)*2, type = "continuous")
ggplot(m, aes(x=log10(GenetDist), y=meanTime, group = selectionCoeff, colour = selectionCoeff)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~ASEModel) + 
  labs(title="Regulatory mutations: nb of fixed",
       subtitle = creTitle,
       color = "Selection coeff.") +
  ylab("Average number of fixed mutations") +
  xlab("Recombination rate - log10") +
  theme_bw() + 
  scale_color_manual(values = colorsClasses[seq(1,10,by=2)]) +
#scale_x_discrete(labels=log10(unique(m$distanceLabels))) + 
xlab(xLabel)
ggsave(file=paste0("pdfs_rec_analysis/", pdf_prefix, "_number_acrossGenetDist", pdTag, ".pdf"), width = 20, height = 10, units = "cm")

##-----------------
##
##
##-----------------
sumStatTag <- "CRE_Nucleotide_diversity"
titlePlotString <- expression("Nucleotide diversity ("*pi*")")
statsYString <- expression("Nucleotide diversity ("*pi*")")
perHap <- F
pi <- T
isYLog <- F

plottingAveSumStatsAcrossRecCategoriesGroupSelection(sumStatTag, titlePlotString,statsYString,perHap, pi, isYLog)

# plotting heterozygosity 
sumStatTag <- "AveHet_m3"
titlePlotString <- "Regulatory heterozygosity"
statsYString <- "Ave. heterozygosity"
perHap <- F
pi <- F
isYLog <- F

plottingAveSumStatsAcrossRecCategoriesGroupSelection(sumStatTag, titlePlotString,statsYString,perHap, pi, isYLog)



# av m2 onto --m3
sumStatTag <- "meanNbm2ontoM3minus"
titlePlotString <- "Mean nb deleterious mutations onto underexpressed background"
statsYString <- "Mean nb deleterious mutations"
perHap <- T
pi <- T
isYLog <- F

plottingAveSumStatsAcrossRecCategoriesGroupASE(sumStatTag, titlePlotString,statsYString,perHap, pi, isYLog)

# av m2 onto --m3
sumStatTag <- "CorrM2M3acrossHaps"
titlePlotString <- "Cor. deleterious mutations vs. reg mutations"
statsYString <- "Correlation"
perHap <- T
pi <- T
isYLog <- F

plottingAveSumStatsAcrossRecCategoriesGroupSelection(sumStatTag, titlePlotString,statsYString,perHap, pi, isYLog)



EOF


        cat << 'EOF' > {output.per_h_stats}

# ==============================================================================
# plotting different statistics, for a specific model, per h values
# ==============================================================================


library(ggplot2)
library(tidyverse)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
wrkDir <- args[1]
setwd(wrkDir)

# Create output directories if they do not exist
if(!dir.exists("pdfs_h_analysis")) {{
  dir.create("pdfs_h_analysis")
}}

if(!dir.exists("p_values")) {{
  dir.create("p_values")
}}

#models parameters
rec_val  <- args[2] # Physical distance in kb
cre_val  <- args[3] # CRE length in kb
sM3_val  <- args[4] # Selection coeff of Regulatory mutations
sm2_val  <- args[5] # selection Coeff on deleterious coding mutations
h_string <- args[6] # string of dominance coefficients separated by comma

# Unsplit the string
h_values   <- unlist(strsplit(h_string, ","))


ase_models <- c("noASE", "ASE") 

#generate all the models names possible with parameters
scenarios_df <- expand.grid(
  h = h_values,
  ase = ase_models,
  stringsAsFactors = FALSE
)

# Create model name for each combinations
scenarios_df <- scenarios_df %>%
  mutate(
    epist_suffix = ifelse(ase == "ASE", "_epist", ""),
    model_name = paste0("fyonetal-h", h, "_rec", rec_val, "kb_cre", cre_val, "kb_sM3", sM3_val, "_s", sm2_val, epist_suffix)
  )


#gathering data from different files
df_timeFix_list <- list()
df_fixedM3_list <- list()
df_segM3_list   <- list()
df_piCRE_list   <- list()

for (i in 1:nrow(scenarios_df)) {{
  row <- scenarios_df[i, ]
  m <- row$model_name
  
  #time to fixation
  f_ttf <- paste0(m, "/", m, "_g10000_m3_aveTimeFix.out")
  if (file.exists(f_ttf)) {{
    tmp <- read.table(f_ttf, header=FALSE)
    tmp$V1[tmp$V1 == 0] <- NA
    df_timeFix_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #number of fixed mutations
  f_fix <- paste0(m, "/", m, "_g10000_m3_fMut.out")
  if (file.exists(f_fix)) {{
    tmp <- read.table(f_fix, header=FALSE)
    df_fixedM3_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #number of segregating mutations
  f_seg <- paste0(m, "/", m, "_g10000_m3_sMut.out")
  if (file.exists(f_seg)) {{
    tmp <- read.table(f_seg, header=FALSE)
    df_segM3_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #nucleotide diversity
  f_pi <- paste0(m, "/pi_", m, ".txt")
  if (file.exists(f_pi)) {{
    tmp <- read.table(f_pi, header=TRUE, sep="\\t")
    df_piCRE_list[[i]] <- data.frame(row, value = tmp$CRE_Nucleotide_diversity)
  }}
}}

# Combine the list into a single dataframe
df_timeFix <- bind_rows(df_timeFix_list)
df_fixedM3 <- bind_rows(df_fixedM3_list)
df_segM3   <- bind_rows(df_segM3_list)
df_piCRE   <- bind_rows(df_piCRE_list)

#plotting function
plot_stat_vs_h <- function(data, title, y_label, subtitle) {{
  
  #if no data returns an empty pdf
  if(nrow(data) == 0) return(ggplot() + ggtitle("No data"))
  
  # Ensure 'h' values are treated as ordered factors
  data$h <- factor(data$h, levels = sort(as.numeric(unique(data$h))))
  

  # Create the base boxplot comparing ASE and noASE models
  p <- ggplot(data, aes(x = h, y = value, fill = ase)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), alpha=0.8) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Dominance coefficient (h)",
      y = y_label,
      fill = "Model"
    ) +
    scale_fill_manual(values = c("noASE" = "#E1AF00", "ASE" = "#3B9AB2")) + 
    theme(legend.position = "bottom")
  
  # Initialize a list to store p-values for each 'h' level
  pValuesList <- list()
  h_levels <- levels(data$h)
  
  # Perform Wilcoxon tests between noASE and ASE for each 'h' value
  for(h_val_current in h_levels) {{
    val_noASE <- data %>% filter(h == h_val_current & ase == "noASE") %>% pull(value)
    val_ASE   <- data %>% filter(h == h_val_current & ase == "ASE") %>% pull(value)
    
    if(sum(!is.na(val_noASE)) > 0 && sum(!is.na(val_ASE)) > 0) {{
      pValuesList[[as.character(h_val_current)]] <- wilcox.test(val_noASE, val_ASE, alternative = "two")$p.value
    }} else {{
      pValuesList[[as.character(h_val_current)]] <- NA
    }}
  }}
  
  # Format the p-values matrix for exporting
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F)
  rownames(pValues_v) <- h_levels
  
  # Format the title for a safe filename and save the computed p-values to a text file
  safe_title <- gsub("[^A-Za-z0-9]", "_", title)
  base_model_name <- paste0("rec", rec_val, "kb_cre", cre_val, "kb_sM3", sM3_val, "_s", sm2_val)
  fichier_pvalues <- paste0("p_values/Hstats_", base_model_name, "_", safe_title, "_pValues.txt")
  
  write.table(pValues_v, file = fichier_pvalues, quote = F, row.names = T, col.names = F, sep = "\\t")
  print(paste("p values saved in :", fichier_pvalues))
  
  # Add p-value significance stars to the plot using ggpubr
  p_final <- p + stat_compare_means(aes(group = ase), label = "p.signif")
  
  return(p_final)
}}


# Define a subtitle with current model parameters
modelSubtitle <- paste0("CRE: ", cre_val, "kb | REC: ", rec_val, "kb | sM3: ", sM3_val, " | sM2: ", sm2_val)
pdf_name <- paste0("pdfs_h_analysis/SinglePlot_rec", rec_val, "kb_cre", cre_val, "kb_sM3", sM3_val, "_s", sm2_val, ".pdf")

pdf(file = pdf_name, width = 8, height = 6)

# Generate and append each plot to the PDF
print(plot_stat_vs_h(df_timeFix, "Average Time to Fixation in CRE", "Generations", modelSubtitle))
print(plot_stat_vs_h(df_fixedM3, "Fixed Mutations in CRE", "Number of mutations", modelSubtitle))
print(plot_stat_vs_h(df_segM3, "Segregating Mutations in CRE", "Number of mutations", modelSubtitle))
print(plot_stat_vs_h(df_piCRE, expression("Nucleotide Diversity ("*pi*") in CRE"), expression(pi), modelSubtitle))

dev.off()
print(paste("PDF generated :", pdf_name))

EOF

        cat << 'EOF' > {output.fixed_mut}
library(ggplot2)
library(tidyverse)
library(wesanderson)
library(ggpubr)

#############################
###
###       MAIN
###
#############################


args = commandArgs(trailingOnly=TRUE)
wrkDir <- args[1]
setwd(wrkDir)

# folder for PDFs
outputPDF <- paste0(wrkDir, "/pdfs_sm2_analysis")
if(!dir.exists(outputPDF)) {{
  dir.create(outputPDF)
}}

# folder for p values
outputPval <- paste0(wrkDir, "/p_values")
if(!dir.exists(outputPval)) {{
  dir.create(outputPval)
}}

# nb of replicates / model
replicates <- as.numeric(args[3])
targetModel <- args[2]

# Define base selection coefficients for coding mutations (sM2)
selCoeff_base <- c("-0.001", "-0.01", "-0.0425", "-0.1", "-0.3")
selCoeff <- rep(selCoeff_base, each=2)
# Generate combinations of scenarios (with and without epistasis/ASE)
scenarios <- paste0(targetModel, "_s", paste0(selCoeff, c("", "_epist")))
majorScenarios <- paste0(targetModel, "_s=", selCoeff)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)


# Configure plot aesthetics (colors and labels)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")
selCoeffLabels <- paste0("s=", selCoeff_base)

# Create a clean subtitle dynamically derived from the model name
modelSubtitle <- gsub("_", " | ", gsub("fyonetal-", "", targetModel))


#############################
###
###       FONCTION
###
#############################


# Function to read data, plot boxplots, and compute Wilcoxon tests
plot_fixed_stats <- function(file_suffix, pName, yAxisLabel, remove_zeros = FALSE) {{
  
  list_df <- list()

  # Loop through all scenarios to read the corresponding data files
  for (i in 1:length(scenarios)) {{
    s <- scenarios[i]
    file_path <- paste0(wrkDir, "/", s, "/", s, file_suffix)
    
    # Only process the file if it exists
    if(file.exists(file_path)) {{
      df <- read.table(file_path, sep = "\\t", header = F)
      # Ensure the file is not empty
      if(nrow(df) > 0) {{
        
        # Optionally remove zeroes (useful for Time to Fixation where 0 means no fixation occurred)
        if(remove_zeros) {{
          df$V1[df$V1 == 0] <- NA 
        }}
        
        # Annotate the dataframe for ggplot
        df$Scenario <- s
        df$selectionCoeff <- selCoeff[i]
        df$aseModel <- aseScenario[i]
        df$model <- majorScenarios[i]
        
        list_df[[s]] <- df
      }}
    }}
  }}
  # Check if data was successfully gathered
  # Exit if missing (e.g., missing m1 data in older h=0.25 or h=0.5 runs)
  print(wrkDir)
  if(length(list_df) == 0) {{
    print(paste("No data for :", file_suffix))
    return(NULL)
  }}
  
  # Combine all lists into a single dataframe
  df_haps <- do.call(rbind, list_df)
  df_haps <- df_haps %>% rename(sumStat = V1)
  
  # Set factor levels to ensure correct ordering on the x-axis
  df_haps$model <- factor(df_haps$model, levels = unique(majorScenarios))
  df_haps$Scenario <- factor(df_haps$Scenario, levels = scenarios)
  
  # Build the base ggplot
  p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
    labs(title = pName, subtitle = modelSubtitle, x = "Selection model", y = yAxisLabel) + 
    geom_boxplot(aes(color=Scenario)) + 
    scale_x_discrete(labels=unique(selCoeffLabels)) + 
    scale_color_manual(values=colorsClasses) + 
    guides(color = "none") + theme_classic(base_size = 15)
  
  # Initialize list for p-values
  pValuesList <- list()
  
  # Perform Wilcoxon rank-sum tests for each selection coefficient
  for(s in unique(selCoeff_base)) {{ 
    val_noASE <- df_haps %>% filter(selectionCoeff == s & aseModel == "noASE") %>% pull(sumStat)
    val_ASE   <- df_haps %>% filter(selectionCoeff == s & aseModel == "ASE") %>% pull(sumStat)
    

    if(sum(!is.na(val_noASE)) > 0 && sum(!is.na(val_ASE)) > 0) {{
      pValuesList[[as.character(s)]] <- wilcox.test(val_noASE, val_ASE, alternative = "two")$p.value
    }} else {{
      pValuesList[[as.character(s)]] <- NA
    }}
  }}
  # Format the p-value results and export to a text file
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F)
  rownames(pValues_v) <- unique(selCoeff_base)
  

  safe_name <- gsub("[^A-Za-z0-9]", "_", pName) 
  fichier_pvalues <- paste0(outputPval, "/", targetModel, "_", safe_name, "_pValues.txt")
  
  write.table(pValues_v, file = fichier_pvalues, quote = F, row.names = T, col.names = F, sep = "\\t")
  print(paste("P-values saved in :", fichier_pvalues))

  # Add significance stars to the plot
  p_final <- p + stat_compare_means(aes(group = aseModel), label = "p.signif")
  
  
  print(p_final)
}}

#############################
###
###       PDF GENERATIOn
###
#############################

pdf_name <- paste0(outputPDF, "/", targetModel, "_FixedStats_violins.pdf")
pdf(file = pdf_name, width = 8, height = 6)

#time of fixation of m3
plot_fixed_stats(
  file_suffix = "_g10000_m3_aveTimeFix.out", 
  pName = "Regulatory mutations: time to fixation", 
  yAxisLabel = "Time to fixation (generations)",
  remove_zeros = TRUE 
)

#number of fixed m3
plot_fixed_stats(
  file_suffix = "_g10000_m3_fMut.out", 
  pName = "Regulatory mutations: Fixed mutations (m3)", 
  yAxisLabel = "Number of mutations",
  remove_zeros = FALSE
)

#number of fixed m2
plot_fixed_stats(
  file_suffix = "_g10000_m2_fMut.out", 
  pName = "Deleterious coding: Fixed mutations (m2)", 
  yAxisLabel = "Number of mutations",
  remove_zeros = FALSE
)

#for h= 0.25 or 0.5, which were made by older scripts than the snakefile, we don't have m1 data, so there won't be any m1 graph, unless you rerun the simulations for 0.25 and 0.5
#number of fixed m1
plot_fixed_stats(
  file_suffix = "_g10000_m1_fMut.out", 
  pName = "Neutral coding: Fixed mutations (m1)", 
  yAxisLabel = "Number of mutations",
  remove_zeros = FALSE
)

dev.off()

print(paste("PDF generated :", pdf_name))

EOF
        cat << 'EOF' > {output.per_sm3_stats}

library(ggplot2)
library(tidyverse)
library(ggpubr)


# Parse command-line arguments provided by Snakemake
args = commandArgs(trailingOnly=TRUE)
wrkDir <- args[1]
setwd(wrkDir)

# Create output directory if they don't exist
outputPDF <- paste0(wrkDir, "/pdfs_sm3_analysis")
if(!dir.exists(outputPDF)) {{
  dir.create(outputPDF)
}}
outputPval <- paste0(wrkDir, "/p_values")
if(!dir.exists(outputPval)) {{
  dir.create(outputPval)
}}

#Models parameters
rec_val <- args[2]       # Physical distance (kb)
cre_val <- args[3]       # CRE size (kb)
h_val <- args[4]         # Dominance coefficient
sm2_val <- args[5]       # Selection on coding mutations
sm3_string <- args[6]    # Comma-separated selection coefficients for regulatory mutations (m3)

# Define the base target model string
targetModel <- paste0("fyonetal-h", h_val, "_rec", rec_val, "kb_cre", cre_val, "kb")

# Split the string of 'sm3' values
sm3_values   <- unlist(strsplit(sm3_string, ","))
ase_models <- c("noASE", "ASE") 

# Generate all the models names possible with parameters
scenarios_df <- expand.grid(
  sm3 = sm3_values,
  ase = ase_models,
  stringsAsFactors = FALSE
)
 # Construct the exact model names corresponding to each combination
scenarios_df <- scenarios_df %>%
  mutate(
    epist_suffix = ifelse(ase == "ASE", "_epist", ""),
    model_name = paste0("fyonetal-h", h_val, "_rec", rec_val, "kb_cre", cre_val, "kb_sM3", sm3, "_s", sm2_val, epist_suffix)
  )

# Initialize empty lists to store data extracted from the simulation files
df_timeFix_list <- list()
df_fixedM3_list <- list()
df_segM3_list   <- list()
df_piCRE_list   <- list()

# Loop through each generated model name to extract relevant statistics
for (i in 1:nrow(scenarios_df)) {{
  row <- scenarios_df[i, ]
  m <- row$model_name
  
  #time to fixation
  f_ttf <- paste0(m, "/", m, "_g10000_m3_aveTimeFix.out")
  if (file.exists(f_ttf)) {{
    tmp <- read.table(f_ttf, header=FALSE)
    tmp$V1[tmp$V1 == 0] <- NA
    df_timeFix_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #number of fixed mutations
  f_fix <- paste0(m, "/", m, "_g10000_m3_fMut.out")
  if (file.exists(f_fix)) {{
    tmp <- read.table(f_fix, header=FALSE)
    df_fixedM3_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #number of segregating mutations
  f_seg <- paste0(m, "/", m, "_g10000_m3_sMut.out")
  if (file.exists(f_seg)) {{
    tmp <- read.table(f_seg, header=FALSE)
    df_segM3_list[[i]] <- data.frame(row, value = tmp$V1)
  }}
  
  #nucleotide diversity
  f_pi <- paste0(m, "/pi_", m, ".txt")
  if (file.exists(f_pi)) {{
    tmp <- read.table(f_pi, header=TRUE, sep="\\t")
    df_piCRE_list[[i]] <- data.frame(row, value = tmp$CRE_Nucleotide_diversity)
  }}
}}

# Combine the lists into single dataframes for ggplot
df_timeFix <- bind_rows(df_timeFix_list)
df_fixedM3 <- bind_rows(df_fixedM3_list)
df_segM3   <- bind_rows(df_segM3_list)
df_piCRE   <- bind_rows(df_piCRE_list)

# Function to plot a given statistic against 'sm3' values and compute Wilcoxon p-values
plot_stat_vs_sm3 <- function(data, title, y_label, subtitle) {{
  
  #if no data returns a pdf page with no plots
  if(nrow(data) == 0) return(ggplot() + ggtitle("No data"))
  
  # Ensure 'sm3' values are treated as ordered factors (decreasing order) for correct x-axis plotting
  valeurs_uniques <- unique(data$sm3)
  niveaux_tries <- valeurs_uniques[order(as.numeric(valeurs_uniques), decreasing = TRUE)]
  data$sm3 <- factor(data$sm3, levels = niveaux_tries)

  # Create the base boxplot comparing ASE and noASE models
  p <- ggplot(data, aes(x = sm3, y = value, fill = ase)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), alpha=0.8) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Selection coefficient for regulatory mutations (sm3)",
      y = y_label,
      fill = "Model"
    ) +
    scale_fill_manual(values = c("noASE" = "#E1AC01", "ASE" = "#3B9AB2")) + 
    theme(legend.position = "bottom")
  
  # Initialize a list to store p-values for each 'sm3' level
  pValuesList <- list()
  
  # Perform Wilcoxon rank-sum tests between noASE and ASE for each 'sm3' value
  for(val in niveaux_tries) {{ 
    
    val_noASE <- data %>% filter(sm3 == val & ase == "noASE") %>% pull(value)
    val_ASE   <- data %>% filter(sm3 == val & ase == "ASE") %>% pull(value)
    
    if(sum(!is.na(val_noASE)) > 0 && sum(!is.na(val_ASE)) > 0) {{
      pValuesList[[as.character(val)]] <- wilcox.test(val_noASE, val_ASE, alternative = "two")$p.value
    }} else {{
      pValuesList[[as.character(val)]] <- NA
    }}
  }}
  
  # Format the p-values matrix for exporting
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F)
  rownames(pValues_v) <- niveaux_tries
  
  # Format the title for a safe filename and save the computed p-values to a text file
  safe_title <- gsub("[^A-Za-z0-9]", "_", title)
  fichier_pvalues <- paste0(outputPval, "/", targetModel, "_s", sm2_val, "_", safe_title, "_pValues.txt")
  
  write.table(pValues_v, file = fichier_pvalues,
              quote = F, row.names = T, col.names = F, sep = "\\t")
  
  print(paste("P-values saved in :", fichier_pvalues))

  # Add p-value significance stars to the plot using ggpubr
  p_final <- p + stat_compare_means(aes(group = ase), label = "p.signif")
  
  return(p_final)
}}





# Define a dynamic subtitle with current model parameters
modelSubtitle <- paste0("h: ", h_val, " | CRE: ", cre_val, "kb | REC: ", rec_val, "kb | sM2: ", sm2_val)
pdf_name <- paste0(outputPDF, "/SinglePlot_h", h_val , "_rec", rec_val, "kb_cre", cre_val, "kb_s", sm2_val, ".pdf")
pdf(file = pdf_name, width = 8, height = 6)

#plotting time
print(plot_stat_vs_sm3(df_timeFix, "Average Time to Fixation in CRE", "Generations", modelSubtitle))
print(plot_stat_vs_sm3(df_fixedM3, "Fixed Mutations in CRE", "Number of mutations", modelSubtitle))
print(plot_stat_vs_sm3(df_segM3, "Segregating Mutations in CRE", "Number of mutations", modelSubtitle))
print(plot_stat_vs_sm3(df_piCRE, expression("Nucleotide Diversity ("*pi*") in CRE"), expression(pi), modelSubtitle))

dev.off()
print(paste("PDF generated :", pdf_name))



EOF
        """

        

        
